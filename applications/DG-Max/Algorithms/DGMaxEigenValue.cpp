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
#include <DGMaxLogger.h>

#include "Base/MeshManipulator.h"
#include "LinearAlgebra/SmallVector.h"
#include "Utilities/GlobalMatrix.h"
#include "Utilities/GlobalVector.h"

template<std::size_t DIM>
DGMaxEigenValue<DIM>::DGMaxEigenValue(Base::MeshManipulator<DIM> &mesh, std::size_t order)
    : mesh_ (mesh)
    , discretization_ (false)
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
class KShift;

struct KShiftStorage
{
    /// Global indices of the rows
    std::vector<PetscInt> rowIndex_;
    /// Global indices of the columns
    std::vector<PetscInt> colIndex_;
    /// The values to be inserted
    std::vector<std::complex<double>> values_;
};

template<std::size_t DIM>
struct KShift
{
private:
    /// First global index for the rows
    const std::size_t startRowIndex_;
    /// The number of subsequent rows
    const std::size_t lengthRowIndex_;
    /// First global index for the columns
    const std::size_t startColIndex_;
    /// The number of subsequent columns
    const std::size_t lengthColIndex_;
    /// Whether the hermitian counter part also needs to be shifted
    const bool hermitian_;
    /// Distance between the two sides for the shift
    const LinearAlgebra::SmallVector<DIM> dx_;
    /// Block that needs to be shifted and inserted
    const LinearAlgebra::MiddleSizeMatrix block1_;
    /// Block that needs to be shifted and inserted for the hermitian part.
    const LinearAlgebra::MiddleSizeMatrix block2_;
public:
    KShift (std::size_t startRowIndex, std::size_t lengthRowIndex,
            std::size_t startColIndex, std::size_t lengthColIndex,
            bool hermitian, const LinearAlgebra::SmallVector<DIM>& dx,
            const LinearAlgebra::MiddleSizeMatrix& block1, const LinearAlgebra::MiddleSizeMatrix& block2)
        : startRowIndex_ (startRowIndex), lengthRowIndex_ (lengthRowIndex)
        , startColIndex_ (startColIndex), lengthColIndex_ (lengthColIndex)
        , hermitian_ (hermitian), dx_ (dx)
        , block1_ (block1), block2_ (block2)
    {}

    /// Create a KShift for a Periodic boundary face
    ///
    /// \param face The face
    /// \param indexing The indexing used in the matrix
    /// \param dx The difference x_L - x_R, where x_L is the position of the face from the
    ///     left element and x_R that as seen from the right.
    /// \return The corresponding shift.
    static KShift<DIM> faceShift(const Base::Face* face,
            const Utilities::GlobalIndexing& indexing, LinearAlgebra::SmallVector<DIM> dx)
    {

        // Element that this processor owns
        const Base::Element* ownedElement = face->getPtrElementLeft();
        // Element that this processor might own (owned == hermitian_)
        const Base::Element* otherElement;
        Base::Side ownedElementSide, otherElementSide;
        bool hermitian;
        if (!ownedElement->isOwnedByCurrentProcessor())
        {
            otherElement = ownedElement;
            ownedElement = face->getPtrElementRight();
            ownedElementSide = Base::Side::RIGHT;
            otherElementSide = Base::Side::LEFT;
            hermitian = false;
        }
        else
        {
            otherElement = face->getPtrElementRight();
            ownedElementSide = Base::Side::LEFT;
            otherElementSide = Base::Side::RIGHT;
            hermitian = otherElement->isOwnedByCurrentProcessor();
        }
        {
            // Rows are scaled by e^(ikx) and columns by e^(-ikx) where x is the centre
            // of the element owning the row/column. As the matrices are inserted from
            // scratch we need to add this factor.
            Geometry::PointPhysical<DIM> centerPhys;
            const Geometry::PointReference<DIM>& center1 = ownedElement->getReferenceGeometry()->getCenter();
            centerPhys = ownedElement->referenceToPhysical(center1);
            // Owned element corresponds to the rows -> +x
            dx += centerPhys.getCoordinates();
            const Geometry::PointReference<DIM>& center2 = otherElement->getReferenceGeometry()->getCenter();
            centerPhys = otherElement->referenceToPhysical(center2);
            // The other element is for the columns -> -x
            dx -= centerPhys.getCoordinates();
        }

        const Base::FaceMatrix& faceMatrix = face->getFaceMatrix(DGMaxDiscretization<DIM>::FACE_MATRIX_ID);
        LinearAlgebra::MiddleSizeMatrix block1, block2;
        // Note ordering, we need the block where the rows are from the ownedElement. The
        // rows correspond to the test functions, hence ownedElementSide as first
        // argument.
        block1 = faceMatrix.getElementMatrix(ownedElementSide, otherElementSide);
        block2 = faceMatrix.getElementMatrix(otherElementSide, ownedElementSide);
        // The stiffness matrix is premultiplied by the (inverse) mass matrix. Thus
        // multiply by the mass matrix of the element corresponding to the rows.
        block1 = ownedElement->getElementMatrix(DGMaxDiscretization<DIM>::MASS_MATRIX_ID) * block1;
        block2 = otherElement->getElementMatrix(DGMaxDiscretization<DIM>::MASS_MATRIX_ID) * block2;
        KShift result (
                indexing.getGlobalIndex(ownedElement, 0), ownedElement->getNumberOfBasisFunctions(0),
                indexing.getGlobalIndex(otherElement, 0), otherElement->getNumberOfBasisFunctions(0),
                hermitian, dx,
                block1, block2
        );
        return result;
    }

    bool shiftNeeded(LinearAlgebra::SmallVector<DIM> k) const
    {
        return std::abs(dx_ * k) > 1e-12;
    }

    /// Shift the block in the given matrix
    void shift(LinearAlgebra::SmallVector<DIM> k, KShiftStorage& storage, Mat mat) const
    {
        if (shiftNeeded(k))
        {
            setupShift(storage);
            insertShiftedBlock(k*dx_, false, storage, mat);
            if (hermitian_)
            {
                insertShiftedBlock(k * dx_, true, storage, mat);
            }
        }
    }

private:
    /// Setup the storage for shifing
    void setupShift(KShiftStorage& storage) const
    {
        // Setup the indices
        if (storage.rowIndex_.size() < lengthRowIndex_)
        {
            storage.rowIndex_.resize(lengthRowIndex_);
        }
        if (storage.colIndex_.size() < lengthColIndex_)
        {
            storage.colIndex_.resize(lengthColIndex_);
        }
        for (PetscInt i = 0; i < lengthRowIndex_; ++i)
        {
            storage.rowIndex_[i] = i + startRowIndex_;
        }
        for (PetscInt i = 0; i < lengthColIndex_; ++i)
        {
            storage.colIndex_[i] = i + startColIndex_;
        }
        // Ensure that there are enough entries
        std::size_t entries = lengthRowIndex_ * lengthColIndex_;
        if (storage.values_.size() < entries)
        {
            storage.values_.resize(entries);
        }
    }

    /// Insert the shifted block into the matrix
    ///
    /// \param kshift The shift k*dx.
    /// \param hermitianPart Whether to insert the Hermitian transposed block (true) or the normal one.
    /// \param storage The storage space to use (prepared with setupShift(KShiftStorage&))
    /// \param mat The matrix to insert into
    void insertShiftedBlock(double kshift, bool hermitianPart, KShiftStorage& storage, Mat mat) const
    {
        const std::complex<double> phase = exp(std::complex<double>(0, hermitianPart ? -kshift : kshift));
        // Copy data to temp matrix
        const LinearAlgebra::MiddleSizeMatrix& block = hermitianPart ? block2_ : block1_;
        for (std::size_t i = 0; i < block.size(); ++i)
        {
            storage.values_[i] = phase * block[i];
        }

        PetscErrorCode err;
#ifdef HPGEM_ASSERTS
        // Check if the matrix is column oriented, as MiddleSizeMatrix stores it in
        // column order.
        PetscBool isRowOriented = PETSC_FALSE;
        //TODO: Does not seem to do anything
        err = MatGetOption(mat, MAT_ROW_ORIENTED, &isRowOriented);
        CHKERRABORT(PETSC_COMM_WORLD, err);
        logger.assert_debug(isRowOriented == PETSC_FALSE, "Requires column orientation");
#endif
        if (!hermitianPart)
        {
            err = MatSetValues(mat,
                    lengthRowIndex_, storage.rowIndex_.data(),
                    lengthColIndex_, storage.colIndex_.data(),
                    storage.values_.data(), INSERT_VALUES);
        }
        else
        {
            err = MatSetValues(mat,
                    lengthColIndex_, storage.colIndex_.data(),
                    lengthRowIndex_, storage.rowIndex_.data(),
                    storage.values_.data(), INSERT_VALUES);
        }
        CHKERRABORT(PETSC_COMM_WORLD, err);
    }
};


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

    int measureAmount = 0;

    //For IP-DG solving a general eigenproblem is slightly faster,
    //but the Brezzi formulation needs an inverse of each block in the mass matrix anyway
    //so there the general eigenproblem is just more work
//    base_.MHasToBeInverted_ = true;
//    base_.assembler->fillMatrices(&base_);
    initializeMatrices(stab);


    std::vector<std::size_t> efieldUnknowns ({0});
    Utilities::GlobalIndexing indexing (&mesh_, Utilities::GlobalIndexing::BLOCKED_PROCESSOR, &efieldUnknowns);
    Utilities::GlobalPetscMatrix massMatrix(indexing, DGMaxDiscretization<DIM>::MASS_MATRIX_ID, -1),
            stiffnessMatrix(indexing, DGMaxDiscretization<DIM>::STIFFNESS_MATRIX_ID, DGMaxDiscretization<DIM>::FACE_MATRIX_ID);
    DGMaxLogger(INFO, "GlobalPetscMatrix initialised");
    Utilities::GlobalPetscVector
            sampleGlobalVector(indexing, -1, -1);
    DGMaxLogger(INFO, "GlobalPetscVector initialised");
    sampleGlobalVector.assemble();
    DGMaxLogger(INFO, "sampleGlobalVector assembled");

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
    error = VecDuplicate(waveVec, &waveVecConjugate);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = VecCopy(waveVec, waveVecConjugate);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = VecConjugate(waveVecConjugate);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    // Setup the boundary block shifting //
    ///////////////////////////////////////

    const std::vector<KShift<DIM>> periodicShifts = findPeriodicShifts(indexing);

    LinearAlgebra::SmallVector<DIM> dk; // Step in k-space from previous solve
    std::size_t maxStep = kpath.totalNumberOfSteps();

    std::vector<std::vector<PetscScalar>> eigenvalues (maxStep);

    for (int i = 0; i < maxStep; ++i)
    {
        DGMaxLogger(INFO, "Computing eigenvalues for k-point %/%", i+1, maxStep);
        if (kpath.dkDidChange(i))
        {

            dk = kpath.dk(i);
            //recompute the shifts
            makeShiftMatrix(massMatrix.getGlobalIndex(), dk, waveVec);
            error = VecCopy(waveVec, waveVecConjugate);
            CHKERRABORT(PETSC_COMM_WORLD, error);
            error = VecConjugate(waveVecConjugate);
            CHKERRABORT(PETSC_COMM_WORLD, error);
        }
        int converged;
        if (i > 0)
        {
            // Extract solutions from previous iteration, needs to be done before
            // changing the matrices for the new k-vector.
            error = EPSGetInvariantSubspace(eigenSolver, eigenVectors);
            CHKERRABORT(PETSC_COMM_WORLD, error);
            error = EPSGetConverged(eigenSolver, &converged);
            CHKERRABORT(PETSC_COMM_WORLD, error);
        }
        error = MatDiagonalScale(product, waveVec, waveVecConjugate);
        CHKERRABORT(PETSC_COMM_WORLD, error);


        // Change to column orientation of MiddleSizeMatrix
        error = MatSetOption(product, MAT_ROW_ORIENTED, PETSC_FALSE);
        CHKERRABORT(PETSC_COMM_WORLD, error);
        KShiftStorage storage;
        LinearAlgebra::SmallVector<DIM> k = kpath.k(i);
        for (const KShift<DIM>& toShift : periodicShifts)
        {
            toShift.shift(k, storage, product);
        }
        // Go back to the default row orientation
        error = MatSetOption(product, MAT_ROW_ORIENTED, PETSC_TRUE);
        CHKERRABORT(PETSC_COMM_WORLD, error);

        // Assembly after inserts from face-phase-shifts
        error = MatAssemblyBegin(product, MAT_FINAL_ASSEMBLY);
        CHKERRABORT(PETSC_COMM_WORLD, error);
        error = MatAssemblyEnd(product, MAT_FINAL_ASSEMBLY);
        CHKERRABORT(PETSC_COMM_WORLD, error);

        //outputs 'old' data
        if (i % 20 == 1)
        {
            //sampleGlobalVector.writeTimeIntegrationVector(measureAmount);
            measureAmount++;
        }

        error = EPSSetOperators(eigenSolver, product, NULL);
        CHKERRABORT(PETSC_COMM_WORLD, error);
        if (i > 0)
        {
            // Use solution of previous time as starting point for the next one.
            error = EPSSetInitialSpace(eigenSolver, converged, eigenVectors);
            CHKERRABORT(PETSC_COMM_WORLD, error);
        }
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
std::vector<KShift<DIM>> DGMaxEigenValue<DIM>::findPeriodicShifts(const Utilities::GlobalIndexing& indexing) const
{
    std::vector<KShift<DIM>> result;
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

        LinearAlgebra::SmallVector<DIM> dx = boundaryFaceShift(*it);
        // TODO: temporary fix for internal faces, see DivDGMaxEigenvalue
        if ((*it)->isInternal() && dx.l2Norm() > 1e-3)
        {
            result.emplace_back(KShift<DIM>::faceShift(*it, indexing, dx));
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
void DGMaxEigenValue<DIM>::makeShiftMatrix(const Utilities::GlobalIndexing& indexing,
                                       const LinearAlgebra::SmallVector<DIM>& direction, Vec& waveVecMatrix) const
{
    PetscErrorCode err = 0;
    for (typename Base::MeshManipulator<DIM>::ConstElementIterator it = mesh_.elementColBegin();
            it != mesh_.elementColEnd(); ++it)
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
