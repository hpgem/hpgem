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

#include "DivDGMaxEigenValue.h"

#include "../DGMaxLogger.h"

#include "Utilities/GlobalMatrix.h"
#include "Utilities/GlobalVector.h"

#include <slepceps.h>

#include <valarray>

template <std::size_t DIM>
DivDGMaxEigenValue<DIM>::DivDGMaxEigenValue(Base::MeshManipulator<DIM>& mesh)
    : mesh_(mesh) {}

template <std::size_t DIM>
typename DivDGMaxEigenValue<DIM>::Result DivDGMaxEigenValue<DIM>::solve(
    EigenValueProblem<DIM> input,
    typename DivDGMaxDiscretization<DIM>::Stab stab, std::size_t order) {
    // Sometimes the solver finds more eigenvalues & vectors than requested, so
    // reserve some extra space for them.
    std::size_t numberOfEigenvalues = input.getNumberOfEigenvalues();
    const PetscInt numberOfEigenVectors =
        std::max(2 * numberOfEigenvalues, numberOfEigenvalues + 10);
    const KSpacePath<DIM>& kpath = input.getPath();

    PetscErrorCode error;
    DGMaxLogger(INFO, "Starting assembly");
    discretization.initializeBasisFunctions(mesh_, order);
    discretization.computeElementIntegrands(mesh_, false, nullptr, nullptr,
                                            nullptr);
    discretization.computeFaceIntegrals(mesh_, nullptr, stab);

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

    const std::vector<Base::Face*> boundaryFaces = findPeriodicBoundaryFaces();
    unsigned long maxBoundaryFaces = boundaryFaces.size();
    // We need to know the maximum number of boundary faces on any node
    MPI_Allreduce(MPI_IN_PLACE, &maxBoundaryFaces, 1, MPI_UNSIGNED_LONG,
                  MPI_MAX, PETSC_COMM_WORLD);
    DGMaxLogger(INFO, "Max local number of periodic boundary faces: %",
                maxBoundaryFaces);

    // Implicitly assume that each element has the same number of DOFs
    std::size_t dofsUPerElement =
                    (*(mesh_.elementColBegin()))->getNumberOfBasisFunctions(0),
                dofsPPerElement =
                    (*(mesh_.elementColBegin()))->getNumberOfBasisFunctions(1),
                totalDofsPerElement = dofsUPerElement + dofsPPerElement;
    // Storage for shifting
    std::valarray<PetscScalar> lrBlockValues(totalDofsPerElement *
                                             totalDofsPerElement);
    std::valarray<PetscScalar> rlBlockValues(totalDofsPerElement *
                                             totalDofsPerElement);

    // Prepare indexing arrays.
    std::valarray<PetscInt> leftNumbers(totalDofsPerElement);
    std::valarray<PetscInt> rightNumbers(totalDofsPerElement);
    std::slice uslice(0, dofsUPerElement, 1),
        pslice(dofsUPerElement, dofsPPerElement, 1);

    for (PetscInt i = 0; i < dofsUPerElement; ++i) {
        leftNumbers[i] = i;
        rightNumbers[i] = i;
    }
    for (PetscInt i = 0; i < dofsPPerElement; ++i) {
        leftNumbers[i + dofsUPerElement] = i;
        rightNumbers[i + dofsUPerElement] = i;
    }

    // Last ids used for offsetting
    PetscInt lastLeftUOffset = 0, lastLeftPOffset = 0, lastRightUOffset = 0,
             lastRightPOffset = 0;
    DGMaxLogger(VERBOSE, "Initialized boundary shifting");

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
        DGMaxLogger(VERBOSE, "Starting boundary k-shift for max % faces",
                    maxBoundaryFaces);
        auto faceIter = boundaryFaces.begin();
        for (unsigned long j = 0; j < maxBoundaryFaces; ++j) {
            double kshift = 0;
            Base::Face* face;
            while (faceIter != boundaryFaces.end() &&
                   std::abs(kshift) <= 1e-12) {
                LinearAlgebra::SmallVector<DIM> shift =
                    boundaryFaceShift(*faceIter);
                kshift = dk * shift;
                // Store before using.
                face = *faceIter;
                faceIter++;
            }
            if (std::abs(kshift) > 1e-12) {
                PetscInt leftUOffset =
                    massMatrix.getGlobalIndex().getGlobalIndex(
                        face->getPtrElementLeft(), 0);
                PetscInt leftPOffset =
                    massMatrix.getGlobalIndex().getGlobalIndex(
                        face->getPtrElementLeft(), 1);
                PetscInt rightUOffset =
                    massMatrix.getGlobalIndex().getGlobalIndex(
                        face->getPtrElementRight(), 0);
                PetscInt rightPOffset =
                    massMatrix.getGlobalIndex().getGlobalIndex(
                        face->getPtrElementRight(), 1);

                PetscInt leftId = face->getPtrElementLeft()->getID(),
                         rightId = face->getPtrElementRight()->getID();

                // Unfortunately slices only support adding with a valarray, not
                // scalars.
                leftNumbers[uslice] += std::valarray<PetscInt>(
                    (leftUOffset - lastLeftUOffset), dofsUPerElement);
                rightNumbers[uslice] += std::valarray<PetscInt>(
                    (rightUOffset - lastRightUOffset), dofsUPerElement);
                leftNumbers[pslice] += std::valarray<PetscInt>(
                    (leftPOffset - lastLeftPOffset), dofsPPerElement);
                rightNumbers[pslice] += std::valarray<PetscInt>(
                    (rightPOffset - lastRightPOffset), dofsPPerElement);
                lastLeftUOffset = leftUOffset;
                lastLeftPOffset = leftPOffset;
                lastRightUOffset = rightUOffset;
                lastRightPOffset = rightPOffset;

                bool ownLeft =
                    face->getPtrElementLeft()->isOwnedByCurrentProcessor();
                bool ownRight =
                    face->getPtrElementRight()->isOwnedByCurrentProcessor();

                // Note that we interleave the left and right code to first do
                // the gets and then the sets and finally the assembly on the
                // global matrix. When doing left (get+set) and then right
                // (get+set) we would need a second assembly step in between.
                if (ownLeft) {
                    error = MatGetValues(stiffnessMatrix, totalDofsPerElement,
                                         &leftNumbers[0], totalDofsPerElement,
                                         &rightNumbers[0], &lrBlockValues[0]);
                    CHKERRABORT(PETSC_COMM_WORLD, error);
                }
                if (ownRight) {
                    error = MatGetValues(stiffnessMatrix, totalDofsPerElement,
                                         &rightNumbers[0], totalDofsPerElement,
                                         &leftNumbers[0], &rlBlockValues[0]);
                    CHKERRABORT(PETSC_COMM_WORLD, error);
                }

                if (ownLeft) {
                    lrBlockValues *= exp(std::complex<double>(0, kshift));
                    error = MatSetValues(stiffnessMatrix, totalDofsPerElement,
                                         &leftNumbers[0], totalDofsPerElement,
                                         &rightNumbers[0], &lrBlockValues[0],
                                         INSERT_VALUES);
                    CHKERRABORT(PETSC_COMM_WORLD, error);
                }

                if (ownRight) {
                    rlBlockValues *= exp(std::complex<double>(0, -kshift));
                    error = MatSetValues(stiffnessMatrix, totalDofsPerElement,
                                         &rightNumbers[0], totalDofsPerElement,
                                         &leftNumbers[0], &rlBlockValues[0],
                                         INSERT_VALUES);
                    CHKERRABORT(PETSC_COMM_WORLD, error);
                }
            }
            // Reassemble.
            error = MatAssemblyBegin(stiffnessMatrix, MAT_FINAL_ASSEMBLY);
            CHKERRABORT(PETSC_COMM_WORLD, error);
            error = MatAssemblyEnd(stiffnessMatrix, MAT_FINAL_ASSEMBLY);
            CHKERRABORT(PETSC_COMM_WORLD, error);
        }

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

    Result result(input, eigenvalues);
    return result;
}

template <std::size_t DIM>
void DivDGMaxEigenValue<DIM>::makeShiftMatrix(
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
void DivDGMaxEigenValue<DIM>::extractEigenvalues(
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
std::vector<Base::Face*> DivDGMaxEigenValue<DIM>::findPeriodicBoundaryFaces()
    const {
    std::vector<Base::Face*> result;
    auto end = mesh_.faceColEnd();
    for (Base::TreeIterator<Base::Face*> it = mesh_.faceColBegin(); it != end;
         ++it) {
        // To check if the face is on a periodic boundary we compare the
        // coordinates of the center of the face according to the elements on
        // each side of the face. For an internal face both elements touch in
        // the mesh, so they should give the same coordinates. For a periodic
        // boundary face, they should differ as they are on different sides of
        // the mesh (for example, one on the top and the other on the bottom).
        // As this should be zero for internal faces and of the size of the mesh
        // for boundary faces, we can use a very sloppy bound.

        // FIXME: The parallel meshes from the preprocessor contain 'boundary
        //  faces' for the ghost cells. These cells are probably there because
        //  the other side of the face is on the inside of the neighbouring
        //  domain, which is an element that is not known to 'this' processor.
        //  Note, the same hack is used in DGMaxEigenvalue
        if ((*it)->isInternal() &&
            Base::L2Norm(boundaryFaceShift(*it)) > 1e-3) {
            result.emplace_back(*it);
        }
    }
    return result;
}

template <std::size_t DIM>
LinearAlgebra::SmallVector<DIM> DivDGMaxEigenValue<DIM>::boundaryFaceShift(
    const Base::Face* face) const {
    logger.assert_always(face->isInternal(), "Internal face boundary");
    const Geometry::PointReference<DIM - 1>& p =
        face->getReferenceGeometry()->getCenter();
    const Geometry::PointPhysical<DIM> pLeftPhys =
        face->getPtrElementLeft()->referenceToPhysical(
            face->mapRefFaceToRefElemL(p));
    const Geometry::PointPhysical<DIM> pRightPhys =
        face->getPtrElementRight()->referenceToPhysical(
            face->mapRefFaceToRefElemR(p));
    return pLeftPhys.getCoordinates() - pRightPhys.getCoordinates();
}

template <std::size_t DIM>
DivDGMaxEigenValue<DIM>::Result::Result(
    EigenValueProblem<DIM> problem,
    std::vector<std::vector<PetscScalar>> eigenvalues)
    : problem_(problem), eigenvalues_(eigenvalues) {}

template <std::size_t DIM>
const EigenValueProblem<DIM>& DivDGMaxEigenValue<DIM>::Result::originalProblem()
    const {
    return problem_;
}

template <std::size_t DIM>
const std::vector<double> DivDGMaxEigenValue<DIM>::Result::frequencies(
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

template class DivDGMaxEigenValue<2>;
template class DivDGMaxEigenValue<3>;
