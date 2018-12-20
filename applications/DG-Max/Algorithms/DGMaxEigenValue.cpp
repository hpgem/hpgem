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

#include "DGMaxEigenValue.h"

#include <iostream>
#include <valarray>
#include <vector>

#include <petscmat.h>
#include <petscvec.h>
#include <slepceps.h>

#include "LinearAlgebra/SmallVector.h"
#include "Utilities/GlobalMatrix.h"
#include "Utilities/GlobalVector.h"

DGMaxEigenValue::DGMaxEigenValue(hpGemUIExtentions &base)
    : base_ (base)
{
    discretization_.initializeBasisFunctions(*(base_.getMesh(0)), base_.getConfigData());
}

EPS DGMaxEigenValue::createEigenSolver()
{
    EPS solver;
    PetscErrorCode err = EPSCreate(PETSC_COMM_WORLD, &solver);
    CHKERRABORT(PETSC_COMM_WORLD, err);

    err = EPSSetProblemType(solver, EPS_NHEP);
    CHKERRABORT(PETSC_COMM_WORLD, err);
    err = EPSSetWhichEigenpairs(solver, EPS_TARGET_REAL);
    CHKERRABORT(PETSC_COMM_WORLD, err);
    // Note, this is to find the bottom band, which should run from omega = 0 to pi (Gamma-X)
    // Thus eigen values 0 to pi^2
//    double target = 2 * M_PI * M_PI;
    double target = 40; // Original target used in the older codes.
    err = EPSSetTarget(solver, target);
    CHKERRABORT(PETSC_COMM_WORLD, err);

    // So far we have configured the the parameters of the eigenvalue solver in
    // code. This overrides these settings with the values that are in SLEPc's
    // options database (if there are any). This can be used for commandline
    // overrides of the standard values.
    err = EPSSetFromOptions(solver);
    CHKERRABORT(PETSC_COMM_WORLD, err);

    return solver;
}

void DGMaxEigenValue::destroyEigenSolver(EPS& eps)
{
    PetscErrorCode err = EPSDestroy(&eps);
    CHKERRABORT(PETSC_COMM_WORLD, err);
}

void DGMaxEigenValue::initializeMatrices(double stab)
{
    //TODO: Which of these InputFunctions are needed?
    discretization_.computeElementIntegrands(*(base_.getMesh(0)), true, nullptr, nullptr, nullptr);
    discretization_.computeFaceIntegrals(*(base_.getMesh(0)), nullptr, stab);
}

void DGMaxEigenValue::solve(const EigenValueProblem& input, double stab)
{

    PetscErrorCode error;
    EPS eigenSolver = createEigenSolver();

    std::cout << "finding a bunch of eigenvalues" << std::endl;
    int measureAmount = 0;

    //For IP-DG solving a general eigenproblem is slightly faster,
    //but the Brezzi formulation needs an inverse of each block in the mass matrix anyway
    //so there the general eigenproblem is just more work
//    base_.MHasToBeInverted_ = true;
//    base_.assembler->fillMatrices(&base_);
    initializeMatrices(stab);

    Utilities::GlobalPetscMatrix massMatrix(base_.getMesh(0), DGMaxDiscretization::MASS_MATRIX_ID, -1),
            stiffnessMatrix(base_.getMesh(0), DGMaxDiscretization::STIFFNESS_MATRIX_ID, DGMaxDiscretization::FACE_MATRIX_ID);
    std::cout << "GlobalPetscMatrix initialised" << std::endl;
    Utilities::GlobalPetscVector
            sampleGlobalVector(base_.getMesh(0), -1, -1);
    std::cout << "GlobalPetscVector initialised" << std::endl;
    sampleGlobalVector.assemble();
    std::cout << "sampleGlobalVector assembled" << std::endl;

    Mat product;
    error = MatMatMult(massMatrix, stiffnessMatrix, MAT_INITIAL_MATRIX, 1.0, &product);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    //Setup the EPS eigen value solver of SLEPC to find the eigenvalues of `product`.
    const PetscInt NUMBER_OF_EIGEN_VALUES = 24;
    error = EPSSetOperators(eigenSolver, product, NULL);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = EPSSetUp(eigenSolver);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    // TODO: It seems that it should use PETSC_DEFAULT instead of PETSC_DECIDE.
    error = EPSSetDimensions(eigenSolver, NUMBER_OF_EIGEN_VALUES, PETSC_DEFAULT, PETSC_DEFAULT);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    // Setup is done, solve for the eigenvalues.
    error = EPSSolve(eigenSolver);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    // Setup eigen vector storage
    const std::size_t NUMBER_OF_EIGEN_VECTORS = 40; //a few extra in case SLEPc finds more than the requested amount of eigenvalues
    Vec *eigenVectors = new Vec[NUMBER_OF_EIGEN_VECTORS];
    error = VecDuplicateVecs(sampleGlobalVector, NUMBER_OF_EIGEN_VECTORS, &eigenVectors);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    // Setup initial wave vector and its conjugate.
    Vec waveVec, waveVecConjugate;
    const int *rows, *columns;
    error = VecDuplicate(sampleGlobalVector, &waveVec);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = VecSetUp(waveVec);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    LinearAlgebra::SmallVector<DIM> k;
    k[0] = M_PI / 20.;
    k[1] = 0;
    if (DIM == 3)
    {
        k[2] = 0;
    }
    makeShiftMatrix(*(base_.getMesh(0)), massMatrix.getGlobalIndex(), k, waveVec);
    error = VecDuplicate(waveVec, &waveVecConjugate);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = VecCopy(waveVec, waveVecConjugate);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = VecConjugate(waveVecConjugate);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    // Setup the boundary block shifting //
    ///////////////////////////////////////

    // For each periodic boundary face the elements L and R disagree on the
    // position of the face. Therefore the factor exp(i k x)

    std::vector<Base::Face*> periodicBoundaryFaces = findPeriodicBoundaryFaces();
    // Assume that the number of basis functions is constant over all the elements.
    std::size_t degreesOfFreedomPerElement = (*base_.getMesh(0)->elementColBegin())->getNumberOfBasisFunctions();

    // Storage for shifting the boundary blocks
    std::valarray<PetscScalar> lrblockvalues(degreesOfFreedomPerElement * degreesOfFreedomPerElement);
    std::valarray<PetscScalar> rlblockvalues(degreesOfFreedomPerElement * degreesOfFreedomPerElement);

    std::valarray<PetscInt> leftNumbers (degreesOfFreedomPerElement);
    std::valarray<PetscInt> rightNumbers (degreesOfFreedomPerElement);
    for(PetscInt i = 0; i < degreesOfFreedomPerElement; ++i)
    {
        leftNumbers[i] = i;
        rightNumbers[i] = i;
    }
    // Last id used for offsetting leftNumber/rightNumber
    PetscInt lastLeftId = 0, lastRightId = 0;

    std::size_t maxStep = 0;
    if (DIM == 2)
    {
        maxStep = 41;
    } else if (DIM == 3)
    {
        maxStep = 61;
    }
    // For testing
    // maxStep = 21;

    std::vector<std::vector<PetscScalar>> eigenvalues (maxStep);

    extractEigenValues(eigenSolver, eigenvalues[0]);
    // extractEigenValues removes all zero eigen values. This is not correct for
    // k = 0, where 0 is a physically interesting eigenvalue.
    eigenvalues[0].insert(eigenvalues[0].begin(), 0);

    for (int i = 1; i < maxStep; ++i)
    {
        std::cout << "Computing eigenvalues for k-point " << i << std::endl;
        // On steps 21 and 41 change the direction. Note that this gives 20 steps in each of the three directions
        if (i == 21)
        {
            //these are only increments, actually this make the wavevector move from pi,0,0 to pi,pi,0
            k[0] = 0;
            k[1] = M_PI / 20.;

            //recompute the shifts
            makeShiftMatrix(*(base_.getMesh(0)), massMatrix.getGlobalIndex(), k, waveVec);
            error = VecCopy(waveVec, waveVecConjugate);
            CHKERRABORT(PETSC_COMM_WORLD, error);
            error = VecConjugate(waveVecConjugate);
            CHKERRABORT(PETSC_COMM_WORLD, error);
        }
        else if (i == 41)
        {
            k[1] = 0;
            k[2] = M_PI / 20.;
            makeShiftMatrix(*(base_.getMesh(0)), massMatrix.getGlobalIndex(), k, waveVec);
            error = VecCopy(waveVec, waveVecConjugate);
            CHKERRABORT(PETSC_COMM_WORLD, error);
            error = VecConjugate(waveVecConjugate);
            CHKERRABORT(PETSC_COMM_WORLD, error);
        }
        error = EPSGetInvariantSubspace(eigenSolver, eigenVectors);
        CHKERRABORT(PETSC_COMM_WORLD, error);
        error = MatDiagonalScale(product, waveVec, waveVecConjugate);
        CHKERRABORT(PETSC_COMM_WORLD, error);


        for (Base::Face* face : periodicBoundaryFaces)
        {
            LinearAlgebra::SmallVector<DIM> shift = boundaryFaceShift(face);
            double kshift = k * shift;
            PetscInt leftId = face->getPtrElementLeft()->getID(),
                     rightId = face->getPtrElementRight()->getID();
            //Skip entries with no noticeable shift
            if (std::abs(kshift) > 1e-12)
            {
                // Compute the global dof-indices for the degrees of freedoms on both
                // sides of the face. Additionally we have to undo the offset from the
                // previous boundary face, hence the lastLeft/RightId.
                leftNumbers += (leftId - lastLeftId) * degreesOfFreedomPerElement;
                rightNumbers += (rightId - lastRightId) * degreesOfFreedomPerElement;
                lastLeftId = leftId;
                lastRightId = rightId;
                // Retrieve, shift and reinsert the matrix elements.
                error = MatGetValues(product, degreesOfFreedomPerElement, &leftNumbers[0], degreesOfFreedomPerElement, &rightNumbers[0], &lrblockvalues[0]);
                CHKERRABORT(PETSC_COMM_WORLD, error);
                error = MatGetValues(product, degreesOfFreedomPerElement, &rightNumbers[0], degreesOfFreedomPerElement, &leftNumbers[0], &rlblockvalues[0]);
                CHKERRABORT(PETSC_COMM_WORLD, error);
                lrblockvalues *= exp(std::complex<double>(0,  kshift));
                rlblockvalues *= exp(std::complex<double>(0, -kshift));
                error = MatSetValues(product, degreesOfFreedomPerElement, &leftNumbers[0], degreesOfFreedomPerElement, &rightNumbers[0], &lrblockvalues[0], INSERT_VALUES);
                CHKERRABORT(PETSC_COMM_WORLD, error);
                error = MatSetValues(product, degreesOfFreedomPerElement, &rightNumbers[0], degreesOfFreedomPerElement, &leftNumbers[0], &rlblockvalues[0], INSERT_VALUES);
                CHKERRABORT(PETSC_COMM_WORLD, error);
                error = MatAssemblyBegin(product, MAT_FINAL_ASSEMBLY);
                CHKERRABORT(PETSC_COMM_WORLD, error);
                error = MatAssemblyEnd(product, MAT_FINAL_ASSEMBLY);
                CHKERRABORT(PETSC_COMM_WORLD, error);
            }
        }

        //outputs 'old' data
        if (i % 20 == 1)
        {
            sampleGlobalVector.writeTimeIntegrationVector(measureAmount);
            measureAmount++;
        }

        int converged;
        error = EPSGetConverged(eigenSolver, &converged);
        CHKERRABORT(PETSC_COMM_WORLD, error);
        error = EPSSetOperators(eigenSolver, product, NULL);
        CHKERRABORT(PETSC_COMM_WORLD, error);
        error = EPSSetInitialSpace(eigenSolver, converged, eigenVectors);
        CHKERRABORT(PETSC_COMM_WORLD, error);
        error = EPSSetUp(eigenSolver);
        CHKERRABORT(PETSC_COMM_WORLD, error);
        error = EPSSolve(eigenSolver);
        CHKERRABORT(PETSC_COMM_WORLD, error);

        extractEigenValues(eigenSolver, eigenvalues[i]);
    }
    // Cleanup
    // Assembly of the eigen value vector

    for (std::size_t kPoint = 0; kPoint < maxStep; kPoint++)
    {
        std::cout << "k-point " << kPoint << "/" << maxStep << std::endl;
        for (PetscScalar eigenValue : eigenvalues[kPoint])
        {
            std::cout << eigenValue << std::endl;
        }
    }

    std::cout << "Real parts for plotting band structure." << std::endl;
    // Table form for plotting
    for (std::size_t kPoint = 0; kPoint < maxStep; kPoint++)
    {
        std::cout << kPoint;
        for (PetscScalar eigenvalue : eigenvalues[kPoint])
        {
            std::cout << "\t" << std::sqrt(eigenvalue.real());
        }
        std::cout << std::endl;
    }


    error = EPSGetInvariantSubspace(eigenSolver, eigenVectors);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    sampleGlobalVector.writeTimeIntegrationVector(measureAmount);

    error = VecDestroy(&waveVec);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = VecDestroy(&waveVecConjugate);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    //always clean up after you are done
    error = MatDestroy(&product);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    destroyEigenSolver(eigenSolver);
}

void DGMaxEigenValue::extractEigenValues(const EPS &solver, std::vector<PetscScalar> &result)
{
    const double ZERO_TOLLERANCE = 1e-10;

    PetscInt converged; //number of converged eigenpairs
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
    // Sort eigen values in ascending order with respect to the real part of the
    // eigenvalue and using the imaginairy part as tie breaker.
    std::sort(result.begin(), result.end(), [](const PetscScalar& a, const PetscScalar& b) {
        if (a.real() != b.real()) {
            return a.real() < b.real();
        } else {
            return a.imag() < b.imag();
        }
    });
}

std::vector<Base::Face*> DGMaxEigenValue::findPeriodicBoundaryFaces() const
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

        // TODO: temporary fix for internal faces, see DivDGMaxEigenvalue
        if ((*it)->isInternal() && Base::L2Norm(boundaryFaceShift(*it)) > 1e-3)
        {
            result.emplace_back(*it);
        }
    }
    return result;
}

LinearAlgebra::SmallVector<DIM> DGMaxEigenValue::boundaryFaceShift(const Base::Face *face) const
{
    logger.assert_always(face->isInternal(), "Internal face boundary");
    const PointFaceReferenceT& p = face->getReferenceGeometry()->getCenter();
    const PointPhysicalT pLeftPhys = face->getPtrElementLeft()
            ->referenceToPhysical(face->mapRefFaceToRefElemL(p));
    const PointPhysicalT pRightPhys = face->getPtrElementRight()
            ->referenceToPhysical(face->mapRefFaceToRefElemR(p));
    return pLeftPhys.getCoordinates() - pRightPhys.getCoordinates();
}

void DGMaxEigenValue::makeShiftMatrix(const Base::MeshManipulator<DIM>& mesh, const Utilities::GlobalIndexing& indexing,
                                       const LinearAlgebra::SmallVector<DIM>& direction, Vec& waveVecMatrix) const
{
    PetscErrorCode err = 0;
    for (Base::MeshManipulator<DIM>::ConstElementIterator it = mesh.elementColBegin(); it != mesh.elementColEnd(); ++it)
    {
        // Note this implicitly assumes we only uses DGBasisFunctions
        const std::size_t basisOffset = indexing.getGlobalIndex(*it, 0);
        for (int j = 0; j < (*it)->getNrOfBasisFunctions(); ++j)
        {
            PointPhysicalT centerPhys;
            const Geometry::PointReference<DIM>& center = (*it)->getReferenceGeometry()->getCenter();
            centerPhys = (*it)->referenceToPhysical(center);
            //this extra accuracy is probably irrelevant and a lot of extra ugly to get it working

            double imPart = direction * centerPhys.getCoordinates();
            PetscScalar value = exp(std::complex<double>(0, imPart));
            // Note this seems inefficient to call this function for each value.
            err = VecSetValue(waveVecMatrix, basisOffset + j, value, INSERT_VALUES);
            CHKERRABORT(PETSC_COMM_WORLD, err);
        }
    }
    VecAssemblyBegin(waveVecMatrix);
    CHKERRABORT(PETSC_COMM_WORLD, err);
    VecAssemblyEnd(waveVecMatrix);
    CHKERRABORT(PETSC_COMM_WORLD, err);
}