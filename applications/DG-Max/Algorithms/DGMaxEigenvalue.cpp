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

#include "DGMaxEigenvalue.h"

#include <chrono>  // For timing
#include <iostream>
#include <utility>
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
#include "Utils/KPhaseShift.h"

template <std::size_t DIM>
DGMaxEigenvalue<DIM>::DGMaxEigenvalue(Base::MeshManipulator<DIM>& mesh,
                                      std::size_t order, SolverConfig config)
    : mesh_(mesh),
      order_(order),
      config_(config),
      discretization_(config.useProjector_) {
    discretization_.initializeBasisFunctions(mesh_, order);
}

template <std::size_t DIM>
void DGMaxEigenvalue<DIM>::initializeMatrices() {
    auto massMatrixHandling = config_.useHermitian_
                                  ? DGMaxDiscretizationBase::ORTHOGONALIZE
                                  : DGMaxDiscretizationBase::INVERT;
    discretization_.computeElementIntegrands(mesh_, massMatrixHandling, nullptr,
                                             nullptr, nullptr);
    discretization_.computeFaceIntegrals(mesh_, massMatrixHandling, nullptr,
                                         config_.stab_);
}

// SolverWorkspace //
/////////////////////

/// Workspace area for the solver
struct SolverWorkspace {
    SolverWorkspace(DGMaxEigenvalueBase::SolverConfig config)
        : config_(config),
          mesh_(nullptr),
          fieldIndex_(nullptr),
          projectorIndex_(nullptr),
          stiffnessMatrix_(fieldIndex_,
                           DGMaxDiscretizationBase::STIFFNESS_MATRIX_ID,
                           DGMaxDiscretizationBase::FACE_MATRIX_ID),
          massMatrix_(fieldIndex_, DGMaxDiscretizationBase::MASS_MATRIX_ID, -1),
          projectorMatrix_(projectorIndex_, fieldIndex_,
                           DGMaxDiscretizationBase::PROJECTOR_MATRIX_ID, -1),
          tempFieldVector_(fieldIndex_, -1, -1),
          tempProjectorVector_(projectorIndex_, -1, -1),
          setupHasBeenRun_(false),
          targetFrequency_(1) {}

    /// Initialize the solver, should only be called once.
    void init(Base::MeshManipulatorBase* mesh, std::size_t numberOfEigenvalues,
              std::size_t numberOfEigenvectors);

    /// Update the blocks in the matrices that are affected by the k-shifted
    /// boundary conditions.
    template <std::size_t DIM>
    void shift(const DGMax::KPhaseShifts<DIM>& stiffnessMatrixShifts,
               const DGMax::KPhaseShifts<DIM>& projectorShifts,
               const LinearAlgebra::SmallVector<DIM>& k);
    /// Update the vectors used to adjust for the shifted basis functions.
    template <std::size_t DIM>
    void updateShiftVectors(const LinearAlgebra::SmallVector<DIM>& dk);

    /// Setup the solver for a solve, needs to be called after all the matrices
    /// are shifted.
    void setupSolver();
    /// Extract the eigenvectors after solving.
    void extractEigenVectors();

    /// Project a vector to remove the null-space
    void project(Vec vec);

    /// Cleanup
    void cleanup();

    /// Helper function for getting the Stiffness matrix to use as basis in the
    /// eigenvalue problem.
    Mat getActualStiffnessMatrix() {
        return config_.useHermitian_ ? stiffnessMatrix_ : product_;
    }

    DGMaxEigenvalueBase::SolverConfig config_;

    Base::MeshManipulatorBase* mesh_;
    Utilities::GlobalIndexing fieldIndex_;
    Utilities::GlobalIndexing projectorIndex_;

    // Matrices //
    //////////////
    Utilities::GlobalPetscMatrix stiffnessMatrix_;
    // Inverted
    Utilities::GlobalPetscMatrix massMatrix_;
    Utilities::GlobalPetscMatrix projectorMatrix_;
    // Temporary storage vectors
    Utilities::GlobalPetscVector tempFieldVector_;
    Vec tempFieldVector2_;
    Utilities::GlobalPetscVector tempProjectorVector_;

    // Product matrix, massMatrix * stiffnessMatrix
    Mat product_;
    /// Shell matrix P*product_, where P is the projection operator
    Mat shell_;
    /// Stiffness matrix used in the projection operator
    Mat projectionStiffness_;

    // Solver
    EPS solver_;
    /// Solver for the projection stiffness matrix
    KSP projectionSolver_;

    // Vectors corresponding to shifted basis functions
    Vec waveVec_, waveVecConjugate_;
    // Eigenvector storage
    PetscInt convergedEigenValues_;
    Vec* eigenVectors_;
    PetscInt numberOfEigenVectors_;

    // Setup has been run at least once, to allow reusing matrices
    bool setupHasBeenRun_;

    double targetFrequency_;

   private:
    void initMatrices();
    /// Initialize the Shell matrix
    void initShell();
    void initShiftVectors();
    void initSolver(std::size_t numberOfEigenvalues);
    void initEigenvectorStorage(std::size_t numberOfEigenvectors);

    void shellMultiply(Vec in, Vec out);
    static PetscErrorCode staticShellMultiply(Mat mat, Vec in, Vec out);
};

void SolverWorkspace::init(Base::MeshManipulatorBase* mesh,
                           std::size_t numberOfEigenvalues,
                           std::size_t numberOfEigenvectors) {
    mesh_ = mesh;
    initMatrices();
    DGMaxLogger(INFO, "Matrices assembled");
    initShell();
    initShiftVectors();
    initSolver(numberOfEigenvalues);
    initEigenvectorStorage(numberOfEigenvectors);
    DGMaxLogger(INFO, "Solver workspace init completed");
}

void SolverWorkspace::initMatrices() {
    DGMaxLogger(INFO, "DGMaxEigenvalue workspace init start");
    std::vector<std::size_t> fieldUnknowns({0});
    fieldIndex_.reset(mesh_,
                      Utilities::GlobalIndexing::Layout::BLOCKED_PROCESSOR,
                      &fieldUnknowns);

    // Reinit matrices after indices have be updated
    // This also assembles them from the local matrices.
    stiffnessMatrix_.reinit();
    massMatrix_.reinit();
    tempFieldVector_.reinit();

    if (config_.useProjector_) {
        std::vector<std::size_t> projectorUnknowns({1});
        projectorIndex_.reset(
            mesh_, Utilities::GlobalIndexing::Layout::BLOCKED_PROCESSOR,
            &projectorUnknowns);
        projectorMatrix_.reinit();
        tempProjectorVector_.reinit();
    }
    // Initialize the product matrix
    if (!config_.useHermitian_) {
        PetscErrorCode error;
        error = MatMatMult(massMatrix_, stiffnessMatrix_, MAT_INITIAL_MATRIX,
                           1.0, &product_);
        CHKERRABORT(PETSC_COMM_WORLD, error);
    }
    DGMaxLogger(INFO, "DGMaxEigenvalue workspace init finished");
}

void SolverWorkspace::initShell() {
    PetscErrorCode error;
    PetscInt rows = fieldIndex_.getNumberOfLocalBasisFunctions();
    error = MatCreateShell(PETSC_COMM_WORLD, rows, rows, PETSC_DETERMINE,
                           PETSC_DETERMINE, this, &shell_);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = MatShellSetOperation(shell_, MATOP_MULT,
                                 (void (*)(void))staticShellMultiply);
    if (config_.useProjector_) {
        CHKERRABORT(PETSC_COMM_WORLD, error);
        error = KSPCreate(PETSC_COMM_WORLD, &projectionSolver_);
        CHKERRABORT(PETSC_COMM_WORLD, error);
        error = KSPSetType(projectionSolver_, KSPPREONLY);
        CHKERRABORT(PETSC_COMM_WORLD, error);
        PC pc;
        error = KSPGetPC(projectionSolver_, &pc);
        CHKERRABORT(PETSC_COMM_WORLD, error);
        error = PCSetType(pc, PCLU);
        CHKERRABORT(PETSC_COMM_WORLD, error);
    }
}

PetscErrorCode SolverWorkspace::staticShellMultiply(Mat mat, Vec in, Vec out) {
    PetscErrorCode error;
    // TODO Check if the context is correct
    SolverWorkspace* workspace;
    error = MatShellGetContext(mat, &workspace);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    workspace->shellMultiply(in, out);
    return 0;
}

void SolverWorkspace::shellMultiply(Vec in, Vec out) {
    PetscErrorCode error;
    error = MatMult(getActualStiffnessMatrix(), in, out);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    if (config_.useProjector_) {
        project(out);
    }
}

void SolverWorkspace::project(Vec vec) {
    logger.assert_always(config_.useProjector_, "Projecting without projector");
    // Projection P of a vector u, this is
    // P u = u - M^{-1} * B^H * C^{-1} * B * u
    // where
    //   B is the projectionMatrix,
    //   C is projectionStiffness and
    //   M is the inverse mass matrix (only needed for the non Hermitian case)
    // The inverse of M is stored in massMatrix, while C^{-1} is done through
    // the KSP object 'projectorSolver_'
    PetscErrorCode error;

    // t = B*u
    error = MatMult(projectorMatrix_, vec, tempProjectorVector_);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    //    // Diagnostics
    //    PetscReal projectionNorm;
    //    VecNorm(tempProjectorVector_, NORM_2, &projectionNorm);

    // t =  C^{-1}t
    error =
        KSPSolve(projectionSolver_, tempProjectorVector_, tempProjectorVector_);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    // Multiply by -1, as the result of M^{-1} * B^H * C^{-1} * B * u needs to
    // be subtracted from u. By multiplying it by -1, we can later uses
    // MatMultAdd
    error = VecScale(tempProjectorVector_, -1.0);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    // The multiply and add version of the Hermitian case does not work on older
    // version of petsc (the non Hermitian was not tested, but probably also
    // fails). So the 1-function commented out version is replaced by a two step
    // version, first multiplying with a matrix and then adding to the vector.
    if (config_.useHermitian_) {
        // No need to multiply with the mass matrix in the Hermitian case
        // u = u + B^H (-t)
        // error = MatMultHermitianTransposeAdd(projectorMatrix_,
        //                                      tempProjectorVector_, vec, vec);
        error = MatMultHermitianTranspose(
            projectorMatrix_, tempProjectorVector_, tempFieldVector2_);
        CHKERRABORT(PETSC_COMM_WORLD, error);
    } else {
        // v = B^H (-t)
        error = MatMultHermitianTranspose(
            projectorMatrix_, tempProjectorVector_, tempFieldVector_);
        CHKERRABORT(PETSC_COMM_WORLD, error);
        // Finally compute u + M^{-1}v
        // error = MatMultAdd(massMatrix_, tempFieldVector_, vec, vec);
        error = MatMult(massMatrix_, tempFieldVector_, tempFieldVector2_);
        CHKERRABORT(PETSC_COMM_WORLD, error);
    }
    error = VecAXPY(vec, 1.0, tempFieldVector2_);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    //    // Diagnostics
    //    PetscReal correctionNorm, newProjectionNorm;
    //    VecNorm(tempFieldVector_, NORM_2, &correctionNorm);
    //    MatMult(projectorMatrix_, vec, tempProjectorVector_);
    //    VecNorm(tempProjectorVector_, NORM_2, &newProjectionNorm);
    //
    //    std::cout << "After projection " << newProjectionNorm << std::endl;
    //    std::cout
    //            << "Projection norm " << projectionNorm
    //            << ", correction norm " << correctionNorm
    //            << ", new projection norm " << newProjectionNorm
    //            << std::endl;
}

void SolverWorkspace::initShiftVectors() {
    // Setup initial wave vector and its conjugate.
    PetscErrorCode error;
    error = VecDuplicate(tempFieldVector_, &waveVec_);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = VecSetUp(waveVec_);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = VecDuplicate(waveVec_, &waveVecConjugate_);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = VecCopy(waveVec_, waveVecConjugate_);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = VecConjugate(waveVecConjugate_);
    CHKERRABORT(PETSC_COMM_WORLD, error);
}

PetscErrorCode compareEigen(PetscScalar ar, PetscScalar ai, PetscScalar br,
                            PetscScalar bi, PetscInt* res, void* ctx) {

    SolverWorkspace* ctxx = (SolverWorkspace*)ctx;

    // Custom eigenvalue comparison looking for those whose log is closest to a
    // certain target. Negative eigenvalues should not occur, but are sorted as
    // larger (further away from the target) than any positive number.

    // Target frequency, factor 2 as the eigenvalues are frequency squared
    const double target = 2 * std::log(ctxx->targetFrequency_);

    // Documentation is unclear on whether ai and bi are zero.
    double res1 = std::norm(ar) + std::norm(ai);
    double res2 = std::norm(br) + std::norm(bi);
    if (res1 <= 0 && res2 > 0) {
        (*res) = 1;  // Res 2 is positive and thus preferable
    } else if (res1 > 0 && res2 <= 0) {
        (*res) = -1;  // Res1 is positive and thus preferable
    } else if (res1 <= 0 && res2 <= 0) {
        // Both are negative, sort them according to standard order
        if (res1 < res2) {
            (*res) = -1;
        } else if (res2 < res1) {
            (*res) = 1;
        } else {
            (*res) = 0;
        }
    } else {
        // Both positive, order by |log(lambda) - target|
        double lres1 = std::abs(std::log(res1) - target);
        double lres2 = std::abs(std::log(res2) - target);
        if (lres1 < lres2) {
            (*res) = -1;
        } else if (lres2 < lres1) {
            (*res) = 1;
        } else {
            (*res) = 0;
        }
    }
    return 0;
}

void SolverWorkspace::initSolver(std::size_t numberOfEigenvalues) {
    PetscErrorCode err = EPSCreate(PETSC_COMM_WORLD, &solver_);
    CHKERRABORT(PETSC_COMM_WORLD, err);

    err =
        EPSSetProblemType(solver_, config_.useHermitian_ ? EPS_HEP : EPS_NHEP);
    CHKERRABORT(PETSC_COMM_WORLD, err);
    //    err = EPSSetWhichEigenpairs(solver_, EPS_SMALLEST_REAL);
    //    CHKERRABORT(PETSC_COMM_WORLD, err);
    err = EPSSetWhichEigenpairs(solver_, EPS_TARGET_REAL);
    err = EPSSetEigenvalueComparison(solver_, compareEigen, this);
    err = EPSSetTarget(solver_, targetFrequency_ * targetFrequency_);
    CHKERRABORT(PETSC_COMM_WORLD, err);

    err = EPSSetDimensions(solver_, numberOfEigenvalues, PETSC_DEFAULT,
                           PETSC_DEFAULT);
    CHKERRABORT(PETSC_COMM_WORLD, err);

    // So far we have configured the the parameters of the eigenvalue solver in
    // code. This overrides these settings with the values that are in SLEPc's
    // options database (if there are any). This can be used for commandline
    // overrides of the standard values.
    err = EPSSetFromOptions(solver_);
    CHKERRABORT(PETSC_COMM_WORLD, err);
}

void SolverWorkspace::initEigenvectorStorage(std::size_t numberOfEigenvectors) {
    convergedEigenValues_ = 0;
    numberOfEigenVectors_ = numberOfEigenvectors;
    PetscErrorCode error = VecDuplicateVecs(
        tempFieldVector_, numberOfEigenvectors, &eigenVectors_);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = VecDuplicate(tempFieldVector_, &tempFieldVector2_);
    CHKERRABORT(PETSC_COMM_WORLD, error);
}

void SolverWorkspace::setupSolver() {
    DGMaxLogger(INFO, "Setting up solver");
    PetscErrorCode error;
    // Setup the EPS eigen value solver of SLEPC to find the eigenvalues of
    // `product`.
    error = EPSSetOperators(solver_, shell_, NULL);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    if (config_.useProjector_) {
        Mat projectionH;
        error = MatHermitianTranspose(projectorMatrix_, MAT_INITIAL_MATRIX,
                                      &projectionH);
        CHKERRABORT(PETSC_COMM_WORLD, error);

        if (setupHasBeenRun_) {
            // Note MAT_REUSE_MATRIX has been used in the following matrix
            // creation but caused segfaults. Probably from slighly different
            // sparsity patterns.
            error = MatDestroy(&projectionStiffness_);
            CHKERRABORT(PETSC_COMM_WORLD, error);
        }
        if (config_.useHermitian_) {
            // No mass matrix needed.
            error =
                MatMatMult(projectorMatrix_, projectionH, MAT_INITIAL_MATRIX,
                           PETSC_DEFAULT, &projectionStiffness_);
            CHKERRABORT(PETSC_COMM_WORLD, error);
        } else {
            error = MatMatMatMult(projectorMatrix_, massMatrix_, projectionH,
                                  MAT_INITIAL_MATRIX, PETSC_DEFAULT,
                                  &projectionStiffness_);
            CHKERRABORT(PETSC_COMM_WORLD, error);
        }
        error = KSPSetOperators(projectionSolver_, projectionStiffness_,
                                projectionStiffness_);
        CHKERRABORT(PETSC_COMM_WORLD, error);
        error = KSPSetUp(projectionSolver_);
        CHKERRABORT(PETSC_COMM_WORLD, error);
        error = MatDestroy(&projectionH);
        CHKERRABORT(PETSC_COMM_WORLD, error);
        DGMaxLogger(INFO, "Projection solver setup completed");
    }
    // Mark as setup for future calls
    setupHasBeenRun_ = true;
}

void SolverWorkspace::extractEigenVectors() {
    PetscErrorCode error;
    error = EPSGetConverged(solver_, &convergedEigenValues_);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    if (convergedEigenValues_ > numberOfEigenVectors_) {
        // Reallocate
        DGMaxLogger(INFO,
                    "Reallocating eigenvector storage to store % instead of % "
                    "eigenvectors",
                    convergedEigenValues_, numberOfEigenVectors_);
        error = VecDestroyVecs(numberOfEigenVectors_, &eigenVectors_);
        CHKERRABORT(PETSC_COMM_WORLD, error);
        error = VecDuplicateVecs(tempFieldVector_, convergedEigenValues_,
                                 &eigenVectors_);
        CHKERRABORT(PETSC_COMM_WORLD, error);
        numberOfEigenVectors_ = convergedEigenValues_;
    }
    error = EPSGetInvariantSubspace(solver_, eigenVectors_);
    CHKERRABORT(PETSC_COMM_WORLD, error);
}

template <std::size_t DIM>
void SolverWorkspace::shift(
    const DGMax::KPhaseShifts<DIM>& stiffnessMatrixShifts,
    const DGMax::KPhaseShifts<DIM>& projectorShifts,
    const LinearAlgebra::SmallVector<DIM>& k) {
    // Maybe move to KShift as static function
    stiffnessMatrixShifts.apply(k, getActualStiffnessMatrix());
    if (config_.useProjector_) {
        projectorShifts.apply(k, projectorMatrix_);
    }
}

template <std::size_t DIM>
void SolverWorkspace::updateShiftVectors(
    const LinearAlgebra::SmallVector<DIM>& dk) {
    if (!config_.usesShifts()) {
        return;
    }

    auto* mesh = fieldIndex_.getMesh();
    PetscErrorCode err = 0;
    auto end = mesh->elementColEnd();
    for (auto it = mesh->elementColBegin(); it != end; ++it) {
        // Note this implicitly assumes we only uses DGBasisFunctions
        const std::size_t basisOffset = fieldIndex_.getGlobalIndex(*it, 0);
        for (int j = 0; j < (*it)->getNrOfBasisFunctions(0); ++j) {
            Geometry::PointPhysical<DIM> centerPhys;
            const Geometry::PointReference<DIM>& center =
                (*it)->getReferenceGeometry()->getCenter();
            centerPhys = (*it)->referenceToPhysical(center);
            // this extra accuracy is probably irrelevant and a lot of extra
            // ugly to get it working

            double imPart =
                dk * centerPhys.getCoordinates() * config_.shiftFactor_;
            PetscScalar value = exp(std::complex<double>(0, imPart));
            // Note this seems inefficient to call this function for each value.
            err = VecSetValue(waveVec_, basisOffset + j, value, INSERT_VALUES);
            CHKERRABORT(PETSC_COMM_WORLD, err);
        }
    }
    VecAssemblyBegin(waveVec_);
    CHKERRABORT(PETSC_COMM_WORLD, err);
    VecAssemblyEnd(waveVec_);
    CHKERRABORT(PETSC_COMM_WORLD, err);
    // Update the conjugate wave vector
    err = VecCopy(waveVec_, waveVecConjugate_);
    CHKERRABORT(PETSC_COMM_WORLD, err);
    err = VecConjugate(waveVecConjugate_);
    CHKERRABORT(PETSC_COMM_WORLD, err);
}

void SolverWorkspace::cleanup() {
    // TODO: Can this be moved to a destructor
    PetscErrorCode error;
    error = VecDestroyVecs(numberOfEigenVectors_, &eigenVectors_);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = VecDestroy(&tempFieldVector2_);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    error = EPSDestroy(&solver_);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    if (config_.useProjector_) {
        error = KSPDestroy(&projectionSolver_);
        CHKERRABORT(PETSC_COMM_WORLD, error);
        error = MatDestroy(&projectionStiffness_);
        CHKERRABORT(PETSC_COMM_WORLD, error);
    }

    error = VecDestroy(&waveVec_);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = VecDestroy(&waveVecConjugate_);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = MatDestroy(&shell_);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    // always clean up after you are done
    if (!config_.useHermitian_) {
        error = MatDestroy(&product_);
        CHKERRABORT(PETSC_COMM_WORLD, error);
    }
}

///

template <std::size_t DIM>
std::unique_ptr<AbstractEigenvalueResult<DIM>> DGMaxEigenvalue<DIM>::solve(
    const EigenvalueProblem<DIM>& input) {
    const std::size_t numberOfEigenvalues = input.getNumberOfEigenvalues();
    const KSpacePath<DIM>& kpath = input.getPath();

    PetscErrorCode error;
    initializeMatrices();

    SolverWorkspace workspace(config_);
    // Leave a bit room for extra converged eigenvectors
    workspace.init(&mesh_, numberOfEigenvalues,
                   std::max(2 * numberOfEigenvalues, numberOfEigenvalues + 10));

    // Setup the boundary block shifting //
    ///////////////////////////////////////

    DGMax::FaceMatrixKPhaseShiftBuilder<DIM> builder;
    builder.setMatrixExtractor([&](const Base::Face* face) {
        const Base::FaceMatrix& faceMatrix =
            face->getFaceMatrix(DGMaxDiscretizationBase::FACE_MATRIX_ID);
        LinearAlgebra::MiddleSizeMatrix block1, block2;
        block1 =
            faceMatrix.getElementMatrix(Base::Side::LEFT, Base::Side::RIGHT);
        block2 =
            faceMatrix.getElementMatrix(Base::Side::RIGHT, Base::Side::LEFT);

        if (!config_.useHermitian_) {
            // In the non Hermitian version we assemble the matrices S and
            // M^{-1} but solve the system with the product matrix M^{-1}S. So
            // when inserting the blocks for a k-phase shift we also need the
            // blocks for M^{-1}S. The block diagonal structure of M^{-1} makes
            // this easy, as we need to rescale the rows S by the block on the
            // same row in M^{-1}.
            block1 = face->getPtrElementLeft()->getElementMatrix(
                         DGMaxDiscretizationBase::MASS_MATRIX_ID) *
                     block1;
            block2 = face->getPtrElementRight()->getElementMatrix(
                         DGMaxDiscretizationBase::MASS_MATRIX_ID) *
                     block2;
        }

        return std::make_pair(block1, block2);
    });

    if (config_.usesShifts()) {
        builder.setExtraShift([&](const Base::Face* face) {
            LinearAlgebra::SmallVector<DIM> dx;
            // Rows are scaled by e^(ikx) and columns by e^(-ikx) where x is the
            // centre of the element owning the row/column. As the matrices are
            // inserted from scratch we need to add this factor. The factor for
            // this function is constructed considering a row from the left
            // element with a column for the right element.
            Geometry::PointPhysical<DIM> centerPhys;
            const Geometry::PointReference<DIM>& center1 =
                face->getPtrElementLeft()->getReferenceGeometry()->getCenter();
            centerPhys =
                face->getPtrElementLeft()->referenceToPhysical(center1);
            // Left element corresponds to the rows -> +x
            dx += centerPhys.getCoordinates() * config_.shiftFactor_;
            const Geometry::PointReference<DIM>& center2 =
                face->getPtrElementRight()->getReferenceGeometry()->getCenter();
            centerPhys =
                face->getPtrElementRight()->referenceToPhysical(center2);
            // The Right element is for the columns -> -x
            dx -= centerPhys.getCoordinates() * config_.shiftFactor_;

            return dx;
        });
    }

    const DGMax::KPhaseShifts<DIM> periodicShifts =
        builder.build(workspace.fieldIndex_);

    DGMax::KPhaseShifts<DIM> projectorShifts;
    if (config_.useProjector_) {
        DGMax::CGDGMatrixKPhaseShiftBuilder<DIM> projectorBuilder;

        projectorBuilder.setMatrixExtractor([&](const Base::Element* element) {
            return element->getElementMatrix(
                DGMaxDiscretizationBase::PROJECTOR_MATRIX_ID);
        });

        if (config_.usesShifts()) {
            // TODO: Think about where the -sign for this shift should go
            projectorBuilder.setExtraShift([&](const Base::Element* element) {
                const Geometry::PointReference<DIM>& center =
                    element->getReferenceGeometry()->getCenter();
                // Rescaling of the columns, as the shifts for the field basis
                // functions are e^{-ikx}. Where the convention is used that the
                // trial basis functions use the - sign.
                return element->referenceToPhysical(center).getCoordinates() *
                       config_.shiftFactor_;
            });
        }
        projectorShifts = projectorBuilder.build(workspace.projectorIndex_,
                                                 workspace.fieldIndex_);
    }

    LinearAlgebra::SmallVector<DIM> dk;  // Step in k-space from previous solve
    std::size_t maxStep = kpath.totalNumberOfSteps();

    std::vector<std::vector<PetscScalar>> eigenvalues(maxStep);

    for (int i = 0; i < maxStep; ++i) {
        DGMaxLogger(INFO, "Computing eigenvalues for k-point %/%", i + 1,
                    maxStep);
        if (kpath.dkDidChange(i)) {
            workspace.updateShiftVectors(kpath.dk(i));
        }
        if (i > 0) {
            workspace.extractEigenVectors();
        }
        if (config_.usesShifts()) {
            error = MatDiagonalScale(workspace.getActualStiffnessMatrix(),
                                     workspace.waveVec_,
                                     workspace.waveVecConjugate_);
            CHKERRABORT(PETSC_COMM_WORLD, error);
            // TODO: Check
            if (config_.useProjector_) {
                error = MatDiagonalScale(workspace.projectorMatrix_, nullptr,
                                         workspace.waveVecConjugate_);
            }
            CHKERRABORT(PETSC_COMM_WORLD, error);
        }

        workspace.shift(periodicShifts, projectorShifts, kpath.k(i));

        workspace.setupSolver();

        PetscInt usableInitialVectors;
        if (i == 0) {
            DGMaxLogger(INFO, "Generating initial vector");
            error = VecSetRandom(workspace.eigenVectors_[0], nullptr);
            CHKERRABORT(PETSC_COMM_WORLD, error);
            usableInitialVectors = 1;
        } else {
            DGMaxLogger(INFO, "Combining previous eigen vectors");
            for (PetscInt j = 1; j < workspace.convergedEigenValues_; ++j) {
                // Some eigenvalue solvers only uses a single starting vector.
                // Mix the eigenvalue spaces from the previous k-point in the
                // hope that these are rich in the eigenvectors for the next
                // space.
                error = VecAYPX(workspace.eigenVectors_[0], 1,
                                workspace.eigenVectors_[j]);
                CHKERRABORT(PETSC_COMM_WORLD, error);
            }
            // Add all previous eigenvectors to the eigenvalue solver, even
            // if they are not all used.
            usableInitialVectors = workspace.convergedEigenValues_;
        }
        if (config_.useProjector_) {
            for (std::size_t j = 0; j < usableInitialVectors; ++j) {
                workspace.project(workspace.eigenVectors_[j]);
            }
            DGMaxLogger(INFO, "Projected initial vector");
        }

        // Use solution of previous time as starting point for the next one.
        error = EPSSetInitialSpace(workspace.solver_, usableInitialVectors,
                                   workspace.eigenVectors_);
        CHKERRABORT(PETSC_COMM_WORLD, error);

        DGMaxLogger(INFO, "Solving eigenvalue problem");

        error = EPSSetUp(workspace.solver_);
        CHKERRABORT(PETSC_COMM_WORLD, error);
        DGMaxLogger(INFO, "Solver setup completed");

        auto start = std::chrono::high_resolution_clock::now();

        error = EPSSolve(workspace.solver_);
        CHKERRABORT(PETSC_COMM_WORLD, error);

        // Some basic statistics
        std::chrono::duration<double> time =
            std::chrono::high_resolution_clock::now() - start;
        PetscInt numEigenvalues, iterations;
        error = EPSGetConverged(workspace.solver_, &numEigenvalues);
        CHKERRABORT(PETSC_COMM_WORLD, error);
        error = EPSGetIterationNumber(workspace.solver_, &iterations);
        CHKERRABORT(PETSC_COMM_WORLD, error);

        DGMaxLogger(INFO,
                    "Eigenvalue solver stopped after % iterations with % "
                    "eigenvalues in %s",
                    iterations, numEigenvalues, time.count());

        // Actual result processing
        extractEigenValues(workspace.solver_, eigenvalues[i]);
    }

    workspace.extractEigenVectors();
    // Diagnostics
    // Useful for debugging, but not in a real situation.
    //    if (config_.useProjector_) {
    //        std::cout << "Test projection on results" << std::endl;
    //        for (PetscInt i = 0; i < workspace.convergedEigenValues_; ++i) {
    //            // Diagnostics on the projection operator. Theoretically we
    //            have
    //            // Bu == 0 <=> lambda != 0. For each of the eigenpairs
    //            (u,lambda) we
    //            // test whether Bu == 0 <=> lambda != 0 holds approximately.
    //            PetscScalar eigenValue;
    //            error = EPSGetEigenpair(workspace.solver_, i, &eigenValue,
    //            nullptr,
    //                                    workspace.eigenVectors_[i], nullptr);
    //            CHKERRABORT(PETSC_COMM_WORLD, error);
    //            std::cout << "Eigenvalue " << i << ":\t" <<
    //            std::abs(eigenValue)
    //                      << "\t" << eigenValue << std::endl;
    //            bool isPracticallyZero = std::abs(eigenValue) < 1e-5;
    //            PetscReal normB;
    //            MatMult(workspace.projectorMatrix_,
    //            workspace.eigenVectors_[i],
    //                    workspace.tempProjectorVector_);
    //            VecNorm(workspace.tempProjectorVector_, NORM_2, &normB);
    //            std::cout << "\t Bu=" << normB << " expected "
    //                      << (isPracticallyZero ? "!=" : "==") << " 0" <<
    //                      std::endl;
    //        }
    //    }

    workspace.cleanup();

    return std::make_unique<Result>(input, eigenvalues);
}

template <std::size_t DIM>
void DGMaxEigenvalue<DIM>::extractEigenValues(
    const EPS& solver, std::vector<PetscScalar>& result) {
    const double ZERO_TOLLERANCE = 1e-10;

    PetscInt converged;  // number of converged eigenpairs
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

    logger(INFO, "Number of eigenvalues:  %.", result.size());
    // Sort eigen values in ascending order with respect to the real part of the
    // eigenvalue and using the imaginairy part as tie breaker.
    std::sort(result.begin(), result.end(),
              [](const PetscScalar& a, const PetscScalar& b) {
                  if (a.real() != b.real()) {
                      return a.real() < b.real();
                  }
                  return a.imag() < b.imag();
              });
}

template <std::size_t DIM>
DGMaxEigenvalue<DIM>::Result::Result(
    EigenvalueProblem<DIM> problem,
    std::vector<std::vector<PetscScalar>> values)
    : problem_(problem), eigenvalues_(values) {
    logger.assert_always(
        problem.getPath().totalNumberOfSteps() == values.size(),
        "Eigenvalues are not provided for each k-point.");
}

template <std::size_t DIM>
const EigenvalueProblem<DIM>& DGMaxEigenvalue<DIM>::Result::originalProblem()
    const {
    return problem_;
}

template <std::size_t DIM>
const std::vector<double> DGMaxEigenvalue<DIM>::Result::frequencies(
    std::size_t point) const {
    logger.assert_debug(
        point >= 0 && point < problem_.getPath().totalNumberOfSteps(),
        "Point number outside of valid range for the path");

    std::vector<double> result;
    result.reserve(eigenvalues_[point].size());
    for (PetscScalar eigenvalue : eigenvalues_[point]) {
        result.emplace_back(std::sqrt(std::abs(eigenvalue.real())));
    }
    return result;
}

template class DGMaxEigenvalue<2>;
template class DGMaxEigenvalue<3>;
