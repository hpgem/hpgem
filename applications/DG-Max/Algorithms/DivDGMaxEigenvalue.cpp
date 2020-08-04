/*
This file forms part of hpGEM. This package has been developed over a number of
years by various people at the University of Twente and a full list of
contributors can be found at http://hpgem.org/about-the-code/team

This code is distributed using BSD 3-Clause License. A copy of which can found
below.


Copyright (c) 2018, Univesity of Twenete
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "DivDGMaxEigenvalue.h"

#include <valarray>

#include "Utilities/GlobalMatrix.h"
#include "Utilities/GlobalVector.h"

#include "Utils/KPhaseShift.h"

#include "DGMaxLogger.h"

using namespace hpgem;

template <std::size_t DIM>
DivDGMaxEigenvalue<DIM>::DivDGMaxEigenvalue(
    Base::MeshManipulator<DIM>& mesh, std::size_t order,
    typename DivDGMaxDiscretization<DIM>::Stab stab)
    : mesh_(mesh), order_(order), stab_(stab) {}

template <std::size_t DIM>
std::unique_ptr<AbstractEigenvalueResult<DIM>> DivDGMaxEigenvalue<DIM>::solve(
    const EigenvalueProblem<DIM>& input) {
    // Sometimes the solver finds more eigenvalues & vectors than requested, so
    // reserve some extra space for them.
    std::size_t numberOfEigenvalues = input.getNumberOfEigenvalues();
    const PetscInt numberOfEigenVectors =
        std::max(2 * numberOfEigenvalues, numberOfEigenvalues + 10);
    const KSpacePath<DIM>& kpath = input.getPath();

    PetscErrorCode error;
    DGMaxLogger(INFO, "Starting assembly");
    discretization.initializeBasisFunctions(mesh_, order_);
    discretization.computeElementIntegrands(mesh_, false, nullptr, nullptr,
                                            nullptr);
    discretization.computeFaceIntegrals(mesh_, nullptr, stab_);

    Utilities::GlobalIndexing indexing(&mesh_);
    Utilities::GlobalPetscMatrix massMatrix(
        indexing, DivDGMaxDiscretization<DIM>::ELEMENT_MASS_MATRIX_ID, -1),
        stiffnessMatrix(
            indexing, DivDGMaxDiscretization<DIM>::ELEMENT_STIFFNESS_MATRIX_ID,
            DivDGMaxDiscretization<DIM>::FACE_STIFFNESS_MATRIX_ID);
    DGMaxLogger(INFO, "Mass and stiffness matrices assembled");
    Utilities::GlobalPetscVector globalVector(indexing, -1, -1);
    DGMaxLogger(INFO, "Dummy vector assembled");

    // Setup the boundary block shifting //
    ///////////////////////////////////////
    DGMax::KPhaseShifts<DIM> kphaseshifts;
    {
        DGMax::FaceMatrixKPhaseShiftBuilder<DIM> builder;
        builder.setMatrixExtractor([&](const Base::Face* face) {
            const Base::FaceMatrix& faceMatrix = face->getFaceMatrix(
                DivDGMaxDiscretization<DIM>::FACE_STIFFNESS_MATRIX_ID);
            LinearAlgebra::MiddleSizeMatrix block1, block2;
            block1 = faceMatrix.getElementMatrix(Base::Side::LEFT,
                                                 Base::Side::RIGHT);
            block2 = faceMatrix.getElementMatrix(Base::Side::RIGHT,
                                                 Base::Side::LEFT);

            return std::make_pair(block1, block2);
        });
        builder.setIndexing(&indexing);
        kphaseshifts = builder.build();
        DGMaxLogger(VERBOSE, "Initialized boundary shifting");
    }

    // Setting up eigenvalue solver //
    //////////////////////////////////

    DGMaxLogger(INFO, "Configuring Eigenvalue solver");
    EPS eigenSolver;
    error = EPSCreate(PETSC_COMM_WORLD, &eigenSolver);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    // SH 180221
    error = EPSSetOperators(eigenSolver, massMatrix, stiffnessMatrix);
    // MatSetOption(product, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);
    // MatSetOption(stiffnessMatrix, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    // SH 180212
    // error = EPSSetUp(eigenSolver_);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = EPSSetDimensions(eigenSolver, numberOfEigenvalues, PETSC_DECIDE,
                             PETSC_DECIDE);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    // As the stiffness and mass matrix are switched (i.e. M u = lambda S u),
    // the eigenvalues correspond to omega^{-2}, thus the lowest bands give the
    // largest eigenvalues.
    EPSSetWhichEigenpairs(eigenSolver, EPS_LARGEST_MAGNITUDE);

    EPSSetProblemType(eigenSolver, EPS_GNHEP);
    // everything that is set in the code, but before this line is overridden by
    // command-line options
    error = EPSSetFromOptions(eigenSolver);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    // everything that is set in the code, but after this line overrides the
    // comand-line options
    DGMaxLogger(INFO, "Eigenvalue solver configured");

    // Storage vectors //
    /////////////////////

    Vec* eigenVectors;
    eigenVectors =
        new Vec[numberOfEigenVectors];  // a few extra in case SLEPc finds more
                                        // than the requested amount of
                                        // eigenvalues
    error = VecDuplicateVecs(globalVector, numberOfEigenVectors, &eigenVectors);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    int converged = 0;

    Vec waveVec, waveVecConjugate;
    error = VecDuplicate(globalVector, &waveVec);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = VecSetUp(waveVec);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    DGMaxLogger(INFO, "Storage vectors assembled");

    LinearAlgebra::SmallVector<DIM> dk = kpath.dk(1);

    //    makeShiftMatrix(dk, stiffnessMatrix.getGlobalIndex(), waveVec);
    //    error = VecAssemblyBegin(waveVec);
    //    CHKERRABORT(PETSC_COMM_WORLD, error);
    //    error = VecAssemblyEnd(waveVec);
    //    CHKERRABORT(PETSC_COMM_WORLD, error);
    //    error = VecDuplicate(waveVec, &waveVecConjugate);
    //    CHKERRABORT(PETSC_COMM_WORLD, error);
    //    error = VecCopy(waveVec, waveVecConjugate);
    //    CHKERRABORT(PETSC_COMM_WORLD, error);
    //    error = VecConjugate(waveVecConjugate);
    //    CHKERRABORT(PETSC_COMM_WORLD, error);

    std::size_t outputId = 0;
    std::size_t maxSteps = kpath.totalNumberOfSteps();
    // For testing
    //    maxSteps = 21;

    std::vector<std::vector<PetscScalar>> eigenvalues(maxSteps);

    DGMaxLogger(INFO, "Starting k-vector walk");
    for (int i = 0; i < maxSteps; ++i) {
        DGMaxLogger(INFO, "Solving for k-vector %/%", i + 1, maxSteps);

        if (kpath.dkDidChange(i)) {
            dk = kpath.dk(i);
            //            //recompute the shifts
            //            makeShiftMatrix(k, stiffnessMatrix.getGlobalIndex(),
            //            waveVec); error = VecAssemblyBegin(waveVec);
            //            CHKERRABORT(PETSC_COMM_WORLD, error);
            //            error = VecAssemblyEnd(waveVec);
            //            CHKERRABORT(PETSC_COMM_WORLD, error);
            //            error = VecCopy(waveVec, waveVecConjugate);
            //            CHKERRABORT(PETSC_COMM_WORLD, error);
            //            error = VecConjugate(waveVecConjugate);
            //            CHKERRABORT(PETSC_COMM_WORLD, error);
        }
        // SH 180216 turned this off
        // error = EPSGetInvariantSubspace(eigenSolver_, eigenVectors); //Must
        // be put after EPSSolve? CHKERRABORT(PETSC_COMM_WORLD, error);
        // SH 180221
        // error = MatDiagonalScale(product, waveVec, waveVecConjugate);
        // error = MatDiagonalScale(massMatrix, waveVec, waveVecConjugate);
        // error = MatDiagonalScale(stiffnessMatrix, waveVec, waveVecConjugate);
        // error = MatDiagonalScale(massMatrix, waveVecConjugate, waveVec);
        // error = MatDiagonalScale(stiffnessMatrix, waveVecConjugate, waveVec);
        // CHKERRABORT(PETSC_COMM_WORLD, error);

        // To match the AssemblyBegin and AssemblyEnd of Mat, we need to call
        // them the same number of times on each node. The current, slightly
        // inefficient approach is to call it max(BF_i) times, with BF_i the
        // number of boundary faces on node i. This can be further optimized
        // taking into account that not all boundary faces need a shift, and by
        // shifting multiple boundary faces in one step. But that is for a later
        // optimization round.
        DGMaxLogger(VERBOSE, "Starting boundary k-shift");
        kphaseshifts.apply(kpath.k(i), stiffnessMatrix);

        error = EPSSetOperators(eigenSolver, massMatrix, stiffnessMatrix);
        CHKERRABORT(PETSC_COMM_WORLD, error);
        // SH 180212
        if (i > 1) {
            error = EPSSetInitialSpace(eigenSolver, converged, eigenVectors);
            CHKERRABORT(PETSC_COMM_WORLD, error);
        }
        DGMaxLogger(INFO, "Setting up solve for k-vector %/%", i + 1, maxSteps);
        error = EPSSetUp(eigenSolver);
        CHKERRABORT(PETSC_COMM_WORLD, error);

        // Jelmer: Search eigenvalue in neighborhood of this number:
        // EPSSetTarget(eigenSolver_,40.0);
        EPSSetWhichEigenpairs(eigenSolver, EPS_LARGEST_MAGNITUDE);
        // EPSSetFromOptions(eigenSolver_);
        // CHKERRABORT(PETSC_COMM_WORLD);

        DGMaxLogger(INFO, "Solving for k-vector %/%", i + 1, maxSteps);
        error = EPSSolve(eigenSolver);
        CHKERRABORT(PETSC_COMM_WORLD, error);

        extractEigenvalues(eigenSolver, eigenvalues[i]);
        error = EPSGetEigenvector(eigenSolver, 0, globalVector, nullptr);
        CHKERRABORT(PETSC_COMM_WORLD, error);
        // globalVector.writeTimeIntegrationVector(outputId);
        outputId++;
    }
    DGMaxLogger(INFO, "Finished k-vector walk");

    error = EPSDestroy(&eigenSolver);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    return std::make_unique<Result>(input, eigenvalues);
}

template <std::size_t DIM>
void DivDGMaxEigenvalue<DIM>::makeShiftMatrix(
    LinearAlgebra::SmallVector<DIM>& direction,
    const Utilities::GlobalIndexing& index, Vec& waveVecMatrix) {
    PetscErrorCode error;

    for (typename Base::MeshManipulator<DIM>::ElementIterator it =
             mesh_.elementColBegin();
         it != mesh_.elementColEnd(); ++it) {
        Geometry::PointPhysical<DIM> centerPhys;
        const Geometry::PointReference<DIM>& center =
            (*it)->getReferenceGeometry()->getCenter();
        centerPhys = (*it)->referenceToPhysical(center);
        PetscScalar shift = exp(
            std::complex<double>(0, direction * centerPhys.getCoordinates()));

        std::size_t basisOffset = index.getGlobalIndex(*it, 0);
        for (std::size_t j = 0; j < (*it)->getNumberOfBasisFunctions(0); ++j) {
            error = VecSetValue(waveVecMatrix, basisOffset + j, shift,
                                INSERT_VALUES);
            CHKERRABORT(PETSC_COMM_WORLD, error);
        }

        basisOffset = index.getGlobalIndex(*it, 1);
        for (std::size_t j = 0; j < (*it)->getNumberOfBasisFunctions(1); ++j) {
            error = VecSetValue(waveVecMatrix, basisOffset + j, shift,
                                INSERT_VALUES);
            CHKERRABORT(PETSC_COMM_WORLD, error);
        }
    }
}

template <std::size_t DIM>
void DivDGMaxEigenvalue<DIM>::extractEigenvalues(
    const EPS& solver, std::vector<PetscScalar>& result) const {
    int converged;
    PetscErrorCode err;
    PetscScalar eigenvalue, neededOnlyForRealPetsc;

    err = EPSGetConverged(solver, &converged);
    CHKERRABORT(PETSC_COMM_WORLD, err);

    // Retrieve all non zero eigenvalues from the solver.
    result.resize(converged);
    for (int i = 0; i < converged; ++i) {
        // Note, the last parameter is only used for a PETSc compiled using real
        // numbers, where we need two output parameters for a complex number.
        err = EPSGetEigenvalue(solver, i, &eigenvalue, &neededOnlyForRealPetsc);
        CHKERRABORT(PETSC_COMM_WORLD, err);

        result[i] = eigenvalue;
    }

    DGMaxLogger(INFO, "Number of eigenvalues:  %.", result.size());
    // Sort eigenvalues in descending order with respect to the real part of the
    // eigenvalue and using the imaginary part as tie breaker.
    // Note descending, as this gives ascending values for the frequency.
    std::sort(result.begin(), result.end(),
              [](const PetscScalar& a, const PetscScalar& b) {
                  if (a.real() != b.real()) {
                      return a.real() > b.real();
                  }
                  return a.imag() > b.imag();
              });
}

template <std::size_t DIM>
DivDGMaxEigenvalue<DIM>::Result::Result(
    EigenvalueProblem<DIM> problem,
    std::vector<std::vector<PetscScalar>> eigenvalues)
    : problem_(problem), eigenvalues_(eigenvalues) {}

template <std::size_t DIM>
const EigenvalueProblem<DIM>& DivDGMaxEigenvalue<DIM>::Result::originalProblem()
    const {
    return problem_;
}

template <std::size_t DIM>
const std::vector<double> DivDGMaxEigenvalue<DIM>::Result::frequencies(
    std::size_t point) const {
    logger.assert_always(
        point >= 0 && point < problem_.getPath().totalNumberOfSteps(),
        "Invalid point");
    std::vector<double> frequencies(eigenvalues_[point].size());
    for (std::size_t i = 0; i < frequencies.size(); ++i) {
        frequencies[i] = std::sqrt(1.0 / eigenvalues_[point][i].real());
    }
    return frequencies;
}

template class DivDGMaxEigenvalue<2>;
template class DivDGMaxEigenvalue<3>;
