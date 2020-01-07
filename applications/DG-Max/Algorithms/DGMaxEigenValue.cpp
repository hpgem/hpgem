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

#include "Base/MeshManipulator.h"
#include "LinearAlgebra/SmallVector.h"
#include "Utilities/GlobalMatrix.h"
#include "Utilities/GlobalVector.h"

template<std::size_t DIM>
DGMaxEigenValue<DIM>::DGMaxEigenValue(Base::MeshManipulator<DIM> &mesh, std::size_t order)
    : mesh_ (mesh)
{
    discretization_.initializeBasisFunctions(mesh_, order);
}

template<std::size_t DIM>
EPS DGMaxEigenValue<DIM>::createEigenSolver()
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

template<std::size_t DIM>
void DGMaxEigenValue<DIM>::destroyEigenSolver(EPS& eps)
{
    PetscErrorCode err = EPSDestroy(&eps);
    CHKERRABORT(PETSC_COMM_WORLD, err);
}

template<std::size_t DIM>
void DGMaxEigenValue<DIM>::initializeMatrices(double stab)
{
    //TODO: Which of these InputFunctions are needed?
    discretization_.computeElementIntegrands(mesh_, true, nullptr, nullptr, nullptr);
    discretization_.computeFaceIntegrals(mesh_, nullptr, stab);
}

template<std::size_t DIM>
typename DGMaxEigenValue<DIM>::Result DGMaxEigenValue<DIM>::solve(
        const EigenValueProblem<DIM>& input, double stab)
{
    // Sometimes the solver finds more eigenvalues & vectors than requested, so
    // reserve some extra space for them.
    std::size_t numberOfEigenvalues = input.getNumberOfEigenvalues();
    const PetscInt numberOfEigenVectors = std::max(2 * numberOfEigenvalues, numberOfEigenvalues + 10);
    const KSpacePath<DIM>& kpath = input.getPath();

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

    Utilities::GlobalPetscMatrix massMatrix(&mesh_, DGMaxDiscretization<DIM>::MASS_MATRIX_ID, -1),
            stiffnessMatrix(&mesh_, DGMaxDiscretization<DIM>::STIFFNESS_MATRIX_ID, DGMaxDiscretization<DIM>::FACE_MATRIX_ID);
    std::cout << "GlobalPetscMatrix initialised" << std::endl;
    Utilities::GlobalPetscVector
            sampleGlobalVector(&mesh_, -1, -1);
    std::cout << "GlobalPetscVector initialised" << std::endl;
    sampleGlobalVector.assemble();
    std::cout << "sampleGlobalVector assembled" << std::endl;

    Mat product;
    error = MatMatMult(massMatrix, stiffnessMatrix, MAT_INITIAL_MATRIX, 1.0, &product);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    //Setup the EPS eigen value solver of SLEPC to find the eigenvalues of `product`.
    error = EPSSetOperators(eigenSolver, product, NULL);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = EPSSetUp(eigenSolver);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    // TODO: It seems that it should use PETSC_DEFAULT instead of PETSC_DECIDE.
    error = EPSSetDimensions(eigenSolver, numberOfEigenvalues, PETSC_DEFAULT, PETSC_DEFAULT);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    // Setup is done, solve for the eigenvalues.
    error = EPSSolve(eigenSolver);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    // Setup eigen vector storage
    Vec *eigenVectors = new Vec[numberOfEigenVectors];
    error = VecDuplicateVecs(sampleGlobalVector, numberOfEigenVectors, &eigenVectors);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    // Setup initial wave vector and its conjugate.
    Vec waveVec, waveVecConjugate;
    const int *rows, *columns;
    error = VecDuplicate(sampleGlobalVector, &waveVec);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = VecSetUp(waveVec);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    LinearAlgebra::SmallVector<DIM> dk = kpath.dk(1);

    makeShiftMatrix(mesh_, massMatrix.getGlobalIndex(), dk, waveVec);
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
    unsigned long maxBoundaryFaces = periodicBoundaryFaces.size();
    // We need to know the maximum number of boundary faces on any node
    MPI_Allreduce(MPI_IN_PLACE, &maxBoundaryFaces, 1, MPI_UNSIGNED_LONG, MPI_MAX, PETSC_COMM_WORLD);


    // Assume that the number of basis functions is constant over all the elements.
    std::size_t degreesOfFreedomPerElement = (*mesh_.elementColBegin())->getNumberOfBasisFunctions();

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
    PetscInt lastLeftOffset = 0, lastRightOffset = 0;

    std::size_t maxStep = kpath.totalNumberOfSteps();
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
        if (kpath.dkDidChange(i))
        {
            dk = kpath.dk(i);
            //recompute the shifts
            makeShiftMatrix(mesh_, massMatrix.getGlobalIndex(), dk, waveVec);
            error = VecCopy(waveVec, waveVecConjugate);
            CHKERRABORT(PETSC_COMM_WORLD, error);
            error = VecConjugate(waveVecConjugate);
            CHKERRABORT(PETSC_COMM_WORLD, error);
        }
        error = EPSGetInvariantSubspace(eigenSolver, eigenVectors);
        CHKERRABORT(PETSC_COMM_WORLD, error);
        error = MatDiagonalScale(product, waveVec, waveVecConjugate);
        CHKERRABORT(PETSC_COMM_WORLD, error);


        // To match the AssemblyBegin and AssemblyEnd of Mat, we need to call
        // them the same number of times on each node. The current, slightly
        // inefficient approach is to call it max(BF_i) times, with BF_i the
        // number of boundary faces on node i. This can be further optimized
        // taking into account that not all boundary faces need a shift, and by
        // shifting multiple boundary faces in one step. But that is for a later
        // optimization round.
        auto faceIter = periodicBoundaryFaces.begin();
        for(unsigned long j = 0; j < maxBoundaryFaces; ++j)
        {
            double kshift = 0;
            Base::Face* face;
            while (faceIter != periodicBoundaryFaces.end() && std::abs(kshift) <= 1e-12)
            {
                LinearAlgebra::SmallVector<DIM> shift = boundaryFaceShift(*faceIter);
                kshift = dk * shift;
                // Store before using.
                face = *faceIter;
                faceIter++;
            }
            //Skip entries with no noticeable shift
            if (std::abs(kshift) > 1e-12)
            {
                // Compute the global dof-indices for the degrees of freedoms on both
                // sides of the face. Additionally we have to undo the offset from the
                // previous boundary face, hence the lastLeft/RightId.

                PetscInt leftOffset = massMatrix.getGlobalIndex().getGlobalIndex(face->getPtrElementLeft(), 0);
                PetscInt rightOffset = massMatrix.getGlobalIndex().getGlobalIndex(face->getPtrElementRight(), 0);

                leftNumbers += leftOffset - lastLeftOffset;
                rightNumbers += rightOffset - lastRightOffset;
                lastLeftOffset = leftOffset;
                lastRightOffset = rightOffset;
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
            }
            error = MatAssemblyBegin(product, MAT_FINAL_ASSEMBLY);
            CHKERRABORT(PETSC_COMM_WORLD, error);
            error = MatAssemblyEnd(product, MAT_FINAL_ASSEMBLY);
            CHKERRABORT(PETSC_COMM_WORLD, error);
        }

        //outputs 'old' data
        if (i % 20 == 1)
        {
            //sampleGlobalVector.writeTimeIntegrationVector(measureAmount);
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

    error = EPSGetInvariantSubspace(eigenSolver, eigenVectors);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    //sampleGlobalVector.writeTimeIntegrationVector(measureAmount);

    error = VecDestroy(&waveVec);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = VecDestroy(&waveVecConjugate);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    //always clean up after you are done
    error = MatDestroy(&product);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    destroyEigenSolver(eigenSolver);

    Result result (input, eigenvalues);
    return result;
}

template<std::size_t DIM>
void DGMaxEigenValue<DIM>::extractEigenValues(const EPS &solver, std::vector<PetscScalar> &result)
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

template<std::size_t DIM>
std::vector<Base::Face*> DGMaxEigenValue<DIM>::findPeriodicBoundaryFaces() const
{
    std::vector<Base::Face*> result;
    for (Base::TreeIterator<Base::Face*> it = mesh_.faceColBegin(); it != mesh_.faceColEnd(); ++it)
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

template<std::size_t DIM>
LinearAlgebra::SmallVector<DIM> DGMaxEigenValue<DIM>::boundaryFaceShift(const Base::Face *face) const
{
    logger.assert_always(face->isInternal(), "Internal face boundary");
    const Geometry::PointReference<DIM-1>& p = face->getReferenceGeometry()->getCenter();
    const Geometry::PointPhysical<DIM> pLeftPhys = face->getPtrElementLeft()
            ->referenceToPhysical(face->mapRefFaceToRefElemL(p));
    const Geometry::PointPhysical<DIM> pRightPhys = face->getPtrElementRight()
            ->referenceToPhysical(face->mapRefFaceToRefElemR(p));
    return pLeftPhys.getCoordinates() - pRightPhys.getCoordinates();
}

template<std::size_t DIM>
void DGMaxEigenValue<DIM>::makeShiftMatrix(const Base::MeshManipulator<DIM>& mesh, const Utilities::GlobalIndexing& indexing,
                                       const LinearAlgebra::SmallVector<DIM>& direction, Vec& waveVecMatrix) const
{
    PetscErrorCode err = 0;
    for (typename Base::MeshManipulator<DIM>::ConstElementIterator it = mesh.elementColBegin();
            it != mesh.elementColEnd(); ++it)
    {
        // Note this implicitly assumes we only uses DGBasisFunctions
        const std::size_t basisOffset = indexing.getGlobalIndex(*it, 0);
        for (int j = 0; j < (*it)->getNrOfBasisFunctions(); ++j)
        {
            Geometry::PointPhysical<DIM> centerPhys;
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

template<std::size_t DIM>
DGMaxEigenValue<DIM>::Result::Result(
        EigenValueProblem<DIM> problem, std::vector<std::vector<PetscScalar>> values)
    : problem_ (problem)
    , eigenvalues_ (values)
{
    logger.assert_always(problem.getPath().totalNumberOfSteps() == values.size(),
        "Eigenvalues are not provided for each k-point.");
}

template<std::size_t DIM>
const EigenValueProblem<DIM>& DGMaxEigenValue<DIM>::Result::originalProblem() const
{
    return problem_;
}

template<std::size_t DIM>
const std::vector<double> DGMaxEigenValue<DIM>::Result::frequencies(std::size_t point) const
{
    logger.assert_debug(point >= 0 && point < problem_.getPath().totalNumberOfSteps(),
        "Point number outside of valid range for the path");

    std::vector<double> result;
    result.reserve(eigenvalues_[point].size());
    for(PetscScalar eigenvalue : eigenvalues_[point])
    {
        result.emplace_back(std::sqrt(eigenvalue.real()));
    }
    return result;
}

template class DGMaxEigenValue<2>;
template class DGMaxEigenValue<3>;
