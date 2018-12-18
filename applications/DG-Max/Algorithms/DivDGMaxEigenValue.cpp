/*
This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
http://hpgem.org/about-the-code/team

This code is distributed using BSD 3-Clause License. A copy of which can found below.


Copyright (c) 2018, Univesity of Twenete
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "DivDGMaxEigenValue.h"

#include "Utilities/GlobalMatrix.h"
#include "Utilities/GlobalVector.h"

#include <slepceps.h>

#include <valarray>

DivDGMaxEigenValue::DivDGMaxEigenValue(hpGemUIExtentions& base)
    : base_ (base)
{
}

void DivDGMaxEigenValue::solve(EigenValueProblem input, DivDGMaxDiscretization::Stab stab)
{
    const PetscInt NUMBER_OF_EIGEN_VALUES = 24, NUMBER_OF_EIGEN_VECTORS = 40;

    PetscErrorCode error;
    std::cout << "finding a bunch of eigenvalues" << std::endl;
    discretization.initializeBasisFunctions(*base_.getMesh(0));
    discretization.computeElementIntegrands(
            *base_.getMesh(0),
            false,
            nullptr, nullptr, nullptr
    );
    discretization.computeFaceIntegrals(
            *base_.getMesh(0),
            nullptr,
            stab
    );

    Utilities::GlobalPetscMatrix
        massMatrix(base_.getMesh(0), DivDGMaxDiscretization::ELEMENT_MASS_MATRIX_ID, -1),
        stiffnessMatrix(base_.getMesh(0), DivDGMaxDiscretization::ELEMENT_STIFFNESS_MATRIX_ID, DivDGMaxDiscretization::FACE_STIFFNESS_MATRIX_ID);
    std::cout << "GlobalPetscMatrix initialised" << std::endl;
    Utilities::GlobalPetscVector globalVector(base_.getMesh(0), -1, -1);
    std::cout << "globalVector initialised" << std::endl;
    globalVector.assemble();
    std::cout << "globalVector assembled" << std::endl;

    // Setup the boundary block shifting //
    ///////////////////////////////////////

    const std::vector<Base::Face*> boundaryFaces = findPeriodicBoundaryFaces();
    // Implicitly assume that each element has the same number of DOFs
    std::size_t dofsUPerElement = (*base_.getMesh(0)->elementColBegin())->getNumberOfBasisFunctions(0),
                dofsPPerElement = (*base_.getMesh(0)->elementColBegin())->getNumberOfBasisFunctions(1),
                totalDofsPerElement = dofsUPerElement + dofsPPerElement;
    // Storage for shifting
    std::valarray<PetscScalar> lrBlockValues (totalDofsPerElement * totalDofsPerElement);
    std::valarray<PetscScalar> rlBlockValues (totalDofsPerElement * totalDofsPerElement);


    // Prepare indexing arrays.
    std::valarray<PetscInt> leftNumbers (totalDofsPerElement);
    std::valarray<PetscInt> rightNumbers (totalDofsPerElement);
    std::slice uslice (0, dofsUPerElement, 1), pslice (dofsUPerElement, dofsPPerElement, 1);

    for (PetscInt i = 0; i < dofsUPerElement; ++i)
    {
        leftNumbers[i] = i;
        rightNumbers[i] = i;
    }
    for (PetscInt i = 0; i < dofsPPerElement; ++i)
    {
        leftNumbers[i + dofsUPerElement] = i;
        rightNumbers[i + dofsUPerElement] = i;
    }

    // Last ids used for offsetting
    PetscInt lastLeftUOffset = 0, lastLeftPOffset = 0, lastRightUOffset = 0, lastRightPOffset = 0;

    EPS eigenSolver;
    error = EPSCreate(PETSC_COMM_WORLD, &eigenSolver);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    // SH 180221
    error = EPSSetOperators(eigenSolver, massMatrix, stiffnessMatrix);
    //MatSetOption(product, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);
    //MatSetOption(stiffnessMatrix, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    // SH 180212
    //error = EPSSetUp(eigenSolver_);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = EPSSetDimensions(eigenSolver, NUMBER_OF_EIGEN_VALUES, PETSC_DECIDE, PETSC_DECIDE);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    // As the stiffness and mass matrix are switched (i.e. M u = lambda S u),
    // the eigenvalues correspond to omega^{-2}, thus the lowest bands give the
    // largest eigenvalues.
    EPSSetWhichEigenpairs(eigenSolver, EPS_LARGEST_MAGNITUDE);

    EPSSetProblemType(eigenSolver, EPS_GNHEP);
    //everything that is set in the code, but before this line is overridden by command-line options
    error = EPSSetFromOptions(eigenSolver);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    //everything that is set in the code, but after this line overrides the comand-line options
    // SH 180212
    //error = EPSSolve(eigenSolver_);
    std::cout << "EPSSolve for k=0" << std::endl;

    CHKERRABORT(PETSC_COMM_WORLD, error);

    Vec *eigenVectors;
    eigenVectors = new Vec[NUMBER_OF_EIGEN_VECTORS]; //a few extra in case SLEPc finds more than the requested amount of eigenvalues
    error = VecDuplicateVecs(globalVector, NUMBER_OF_EIGEN_VECTORS, &eigenVectors);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    int converged = 0;

    Vec waveVec, waveVecConjugate;
    error = VecDuplicate(globalVector, &waveVec);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = VecSetUp(waveVec);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    LinearAlgebra::SmallVector<DIM> k;
    k[0] = M_PI / 20.0;
    k[1] = 0;

    makeShiftMatrix(k, stiffnessMatrix.getGlobalIndex(), waveVec);
    error = VecAssemblyBegin(waveVec);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = VecAssemblyEnd(waveVec);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = VecDuplicate(waveVec, &waveVecConjugate);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = VecCopy(waveVec, waveVecConjugate);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = VecConjugate(waveVecConjugate);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    std::size_t outputId = 0;
    std::size_t maxSteps;
    if (DIM == 2)
    {
        maxSteps = 41;
    }
    else if (DIM == 3)
    {
        maxSteps = 61;
    }
    // For testing
//    maxSteps = 21;

    error = EPSSolve(eigenSolver);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    std::vector<std::vector<PetscScalar>> eigenvalues (maxSteps);
    // The discretization_ is not correct for k = 0.
    //extractEigenvalues(eigenSolver, eigenvalues[0]);


    for (int i = 1; i < maxSteps; ++i)
    {
        std::cout << i << "/60" << std::endl;

        if (i == 21)
        {
            //these are only increments, actually this make the wavevector move from pi,0,0 to pi,pi,0
            k[0] = 0;
            k[1] = M_PI / 20.;
            //recompute the shifts
            makeShiftMatrix(k, stiffnessMatrix.getGlobalIndex(), waveVec);
            error = VecAssemblyBegin(waveVec);
            CHKERRABORT(PETSC_COMM_WORLD, error);
            error = VecAssemblyEnd(waveVec);
            CHKERRABORT(PETSC_COMM_WORLD, error);
            error = VecCopy(waveVec, waveVecConjugate);
            CHKERRABORT(PETSC_COMM_WORLD, error);
            error = VecConjugate(waveVecConjugate);
            CHKERRABORT(PETSC_COMM_WORLD, error);
        }
        else if (i == 41)
        {
            //these are only increments, actually this make the wavevector move from pi,pi,0 to pi,pi,pi
            k[1] = 0;
            k[2] = M_PI / 20.;
            //recompute the shifts
            makeShiftMatrix(k, stiffnessMatrix.getGlobalIndex(), waveVec);
            error = VecAssemblyBegin(waveVec);
            CHKERRABORT(PETSC_COMM_WORLD, error);
            error = VecAssemblyEnd(waveVec);
            CHKERRABORT(PETSC_COMM_WORLD, error);
            error = VecCopy(waveVec, waveVecConjugate);
            CHKERRABORT(PETSC_COMM_WORLD, error);
            error = VecConjugate(waveVecConjugate);
            CHKERRABORT(PETSC_COMM_WORLD, error);
        }
        // SH 180216 turned this off
        //error = EPSGetInvariantSubspace(eigenSolver_, eigenVectors); //Must be put after EPSSolve?
        CHKERRABORT(PETSC_COMM_WORLD, error);
        // SH 180221
        // error = MatDiagonalScale(product, waveVec, waveVecConjugate);
        //error = MatDiagonalScale(massMatrix, waveVec, waveVecConjugate);
        //error = MatDiagonalScale(stiffnessMatrix, waveVec, waveVecConjugate);
        //error = MatDiagonalScale(massMatrix, waveVecConjugate, waveVec);
        //error = MatDiagonalScale(stiffnessMatrix, waveVecConjugate, waveVec);
        //CHKERRABORT(PETSC_COMM_WORLD, error);

        for (Base::Face* face : boundaryFaces)
        {
            LinearAlgebra::SmallVector<DIM> shift = boundaryFaceShift(face);
            double kshift = k * shift;
            if (std::abs(kshift) > 1e-12)
            {
                PetscInt leftUOffset = massMatrix.getGlobalIndex().getGlobalIndex(face->getPtrElementLeft(), 0);
                PetscInt leftPOffset = massMatrix.getGlobalIndex().getGlobalIndex(face->getPtrElementLeft(), 1);
                PetscInt rightUOffset = massMatrix.getGlobalIndex().getGlobalIndex(face->getPtrElementRight(), 0);
                PetscInt rightPOffset = massMatrix.getGlobalIndex().getGlobalIndex(face->getPtrElementRight(), 1);

                PetscInt leftId = face->getPtrElementLeft()->getID(),
                        rightId = face->getPtrElementRight()->getID();

                // Unfortunately slices only support adding with a valarray, not scalars.
                leftNumbers[uslice]  += std::valarray<PetscInt>((leftUOffset - lastLeftUOffset), dofsUPerElement);
                rightNumbers[uslice] += std::valarray<PetscInt>((rightUOffset - lastRightUOffset), dofsUPerElement);
                leftNumbers[pslice]  += std::valarray<PetscInt>((leftPOffset - lastLeftPOffset) , dofsPPerElement);
                rightNumbers[pslice] += std::valarray<PetscInt>((rightPOffset - lastRightPOffset), dofsPPerElement);
                lastLeftUOffset = leftUOffset;
                lastLeftPOffset = leftPOffset;
                lastRightUOffset = rightUOffset;
                lastRightPOffset = rightPOffset;

                bool ownLeft = stiffnessMatrix.getGlobalIndex().isLocallyOwned(face->getPtrElementLeft(), 0);
                bool ownRight = stiffnessMatrix.getGlobalIndex().isLocallyOwned(face->getPtrElementRight(), 0);

                // Note that we interleave the left and right code to first do
                // the gets and then the sets and finally the assembly on the
                // global matrix. When doing left (get+set) and then right
                // (get+set) we would need a second assembly step in between.
                if (ownLeft)
                {
                    error = MatGetValues(stiffnessMatrix, totalDofsPerElement, &leftNumbers[0], totalDofsPerElement,
                                         &rightNumbers[0], &lrBlockValues[0]);
                    CHKERRABORT(PETSC_COMM_WORLD, error);
                }
                if (ownRight)
                {
                    error = MatGetValues(stiffnessMatrix, totalDofsPerElement, &rightNumbers[0], totalDofsPerElement,
                                         &leftNumbers[0], &rlBlockValues[0]);
                    CHKERRABORT(PETSC_COMM_WORLD, error);
                }

                if (ownLeft)
                {
                    lrBlockValues *= exp(std::complex<double>(0, kshift));
                    error = MatSetValues(stiffnessMatrix, totalDofsPerElement, &leftNumbers[0], totalDofsPerElement,
                                         &rightNumbers[0], &lrBlockValues[0], INSERT_VALUES);
                    CHKERRABORT(PETSC_COMM_WORLD, error);
                }

                if (ownRight)
                {
                    rlBlockValues *= exp(std::complex<double>(0, -kshift));
                    error = MatSetValues(stiffnessMatrix, totalDofsPerElement, &rightNumbers[0], totalDofsPerElement,
                                         &leftNumbers[0], &rlBlockValues[0], INSERT_VALUES);
                    CHKERRABORT(PETSC_COMM_WORLD, error);
                }
                // Reassemble.
                error = MatAssemblyBegin(stiffnessMatrix, MAT_FINAL_ASSEMBLY);
                CHKERRABORT(PETSC_COMM_WORLD, error);
                error = MatAssemblyEnd(stiffnessMatrix, MAT_FINAL_ASSEMBLY);
                CHKERRABORT(PETSC_COMM_WORLD, error);
            }
        }

        error = EPSSetOperators(eigenSolver, massMatrix, stiffnessMatrix);
        CHKERRABORT(PETSC_COMM_WORLD, error);
        // SH 180212
        if (i > 1)
        {
            error = EPSSetInitialSpace(eigenSolver, converged, eigenVectors);
            CHKERRABORT(PETSC_COMM_WORLD, error);
        }
        std::cout << k << std::endl;
        error = EPSSetUp(eigenSolver);
        CHKERRABORT(PETSC_COMM_WORLD, error);

        //Jelmer: Search eigenvalue in neighborhood of this number:
        //EPSSetTarget(eigenSolver_,40.0);
        EPSSetWhichEigenpairs(eigenSolver,EPS_LARGEST_MAGNITUDE);
        //EPSSetFromOptions(eigenSolver_);
        //CHKERRABORT(PETSC_COMM_WORLD);

        error = EPSSolve(eigenSolver);
        CHKERRABORT(PETSC_COMM_WORLD, error);
        extractEigenvalues(eigenSolver, eigenvalues[i]);
        //Jelmer: To plot Eigenvectors:
        error = EPSGetEigenvector(eigenSolver, 0, globalVector, NULL);
        //globalVector.writeTimeIntegrationVector(outputId);
        outputId++;

    }

    error = EPSDestroy(&eigenSolver);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    if (Base::MPIContainer::Instance().getProcessorID() == 0)
    {
        // Shows eigenvalues:
        for (std::size_t kPoint = 0; kPoint < maxSteps; kPoint++)
        {
            std::cout << "k-point " << kPoint << std::endl;
            for (PetscScalar eigenvalue : eigenvalues[kPoint])
            {
                std::cout << eigenvalue << std::endl;
            }
        }
        // Printing table of real parts
        for (std::size_t kPoint = 0; kPoint < maxSteps; kPoint++)
        {
            std::cout << kPoint;
            for (PetscScalar eigenValue : eigenvalues[kPoint])
            {
                std::cout << "\t" << std::sqrt(double(1.0) / eigenValue.real());
            }
            std::cout << std::endl;
        }
    }
}

void DivDGMaxEigenValue::makeShiftMatrix(LinearAlgebra::SmallVector<DIM>& direction, const Utilities::GlobalIndexing& index, Vec& waveVecMatrix)
{
    PetscErrorCode error;

    for (Base::MeshManipulator<DIM>::ElementIterator it = base_.getMesh(0)->elementColBegin();
            it != base_.getMesh(0)->elementColEnd(); ++it)
    {
        PointPhysicalT centerPhys;
        const PointElementReferenceT& center = (*it)->getReferenceGeometry()->getCenter();
        centerPhys = (*it)->referenceToPhysical(center);
        PetscScalar shift = exp(std::complex<double>(0, direction * centerPhys.getCoordinates()));

        std::size_t basisOffset = index.getGlobalIndex(*it, 0);
        for (std::size_t j = 0; j < (*it)->getNumberOfBasisFunctions(0); ++j)
        {
            error = VecSetValue(waveVecMatrix, basisOffset + j, shift, INSERT_VALUES);
            CHKERRABORT(PETSC_COMM_WORLD, error);
        }

        basisOffset = index.getGlobalIndex(*it, 1);
        for (std::size_t j = 0; j < (*it)->getNumberOfBasisFunctions(1); ++j)
        {
            error = VecSetValue(waveVecMatrix, basisOffset + j, shift, INSERT_VALUES);
            CHKERRABORT(PETSC_COMM_WORLD, error);
        }
    }
}

void DivDGMaxEigenValue::extractEigenvalues(const EPS &solver, std::vector<PetscScalar> &result) const
{
    int converged;
    PetscErrorCode err;
    PetscScalar eigenvalue, neededOnlyForRealPetsc;

    err = EPSGetConverged(solver, &converged);
    CHKERRABORT(PETSC_COMM_WORLD, err);

    // Retrieve all non zero eigenvalues from the solver.
    result.resize(converged);
    for (int i = 0; i < converged; ++i)
    {
        // Note, the last parameter is only used for a PETSc compiled using real numbers,
        // where we need two output parameters for a complex number.
        err = EPSGetEigenvalue(solver, i, &eigenvalue, &neededOnlyForRealPetsc);
        CHKERRABORT(PETSC_COMM_WORLD, err);

        result[i] = eigenvalue;
    }

    logger(INFO, "Number of eigenvalues:  %.", result.size());
    // Sort eigenvalues in descending order with respect to the real part of the
    // eigenvalue and using the imaginary part as tie breaker.
    // Note descending, as this gives ascending values for the frequency.
    std::sort(result.begin(), result.end(), [](const PetscScalar& a, const PetscScalar& b) {
        if (a.real() != b.real()) {
            return a.real() > b.real();
        } else {
            return a.imag() > b.imag();
        }
    });
}

std::vector<Base::Face*> DivDGMaxEigenValue::findPeriodicBoundaryFaces() const
{
    std::vector<Base::Face*> result;
    for (Base::TreeIterator<Base::Face*> it = base_.getMesh(0)->faceColBegin(); it != base_.getMesh(0)->faceColEnd(); ++it)
    {
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
        if ((*it)->isInternal() && Base::L2Norm(boundaryFaceShift(*it)) > 1e-3)
        {
            result.emplace_back(*it);
        }
    }
    return result;
}

LinearAlgebra::SmallVector<DIM> DivDGMaxEigenValue::boundaryFaceShift(const Base::Face *face) const
{
    logger.assert_always(face->isInternal(), "Internal face boundary");
    const PointFaceReferenceT& p = face->getReferenceGeometry()->getCenter();
    const PointPhysicalT pLeftPhys = face->getPtrElementLeft()
            ->referenceToPhysical(face->mapRefFaceToRefElemL(p));
    const PointPhysicalT pRightPhys = face->getPtrElementRight()
            ->referenceToPhysical(face->mapRefFaceToRefElemR(p));
    return pLeftPhys.getCoordinates() - pRightPhys.getCoordinates();
}
