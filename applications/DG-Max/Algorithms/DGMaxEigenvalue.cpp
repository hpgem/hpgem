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
#include <ElementInfos.h>

#include "Base/MeshManipulator.h"
#include "LinearAlgebra/SmallVector.h"
#include "Utilities/Eigenpairs.h"
#include "Utilities/GlobalMatrix.h"
#include "Utilities/GlobalVector.h"
#include "Utils/CGDGMatrixKPhaseShiftBuilder.h"
#include "Utils/FaceKPhaseShiftBuilder.h"

using namespace hpgem;

template <std::size_t DIM>
DGMaxEigenvalue<DIM>::DGMaxEigenvalue(Base::MeshManipulator<DIM>& mesh,
                                      std::size_t order, SolverConfig config)
    : mesh_(mesh),
      order_(order),
      config_(config),
      discretization_(config.useProjector_ != DGMaxEigenvalueBase::NONE) {
    discretization_.initializeBasisFunctions(mesh_, order);
}

template <std::size_t DIM>
void DGMaxEigenvalue<DIM>::initializeMatrices() {
    auto massMatrixHandling = config_.useHermitian_
                                  ? DGMaxDiscretizationBase::ORTHOGONALIZE
                                  : DGMaxDiscretizationBase::INVERT;
    // No element vectors
    std::map<std::size_t, typename DGMaxDiscretization<DIM>::InputFunction>
        elementVectors;
    discretization_.computeElementIntegrands(mesh_, massMatrixHandling,
                                             elementVectors);
    // No face vectors
    std::map<std::size_t, typename DGMaxDiscretization<DIM>::FaceInputFunction>
        faceVectors;
    discretization_.computeFaceIntegrals(mesh_, massMatrixHandling, faceVectors,
                                         config_.stab_);
}

// SolverWorkspace //
/////////////////////

template <std::size_t DIM>
class DGMaxEigenvalue<DIM>::SolverWorkspace {
   public:
    /**
     * Initialize the workspace
     * @param config  The configuration of the method
     * @param mesh The mesh to solve on
     * @param targetNumberOfEigenvalues Estimate of the target number of
     * eigenvalues (to allocate space).
     */
    SolverWorkspace(DGMaxEigenvalueBase::SolverConfig config,
                    Base::MeshManipulatorBase* mesh,
                    std::size_t targetNumberOfEigenvalues);

    ~SolverWorkspace();

    /// Update the blocks in the matrices that are affected by the k-shifted
    /// boundary conditions.
    void updateKPoint(const LinearAlgebra::SmallVector<DIM>& newK);
    /// Extract the eigenvectors after solving.
    void extractEigenVectors();

    void solve(std::size_t targetNumberOfEigenvalues);

    const LinearAlgebra::SmallVector<DIM>& getKPoint() const {
        return currentK_;
    }

    const DGMaxEigenvalueBase::SolverConfig& getConfig() const {
        return config_;
    }

    const Utilities::Eigenpairs& getEigenpairs() const { return eigenpairs_; }

    LinearAlgebra::MiddleSizeMatrix computeOverlapIntegrals();

    /**
     * Take an eigenvector and write it to the element local timeintegration
     * vector.
     * @param eigenvectorId The index of the eigenvector
     * @param timeIntegrationVectorId The index of the time integration to which
     * to write it.
     */
    void writeAsTimeIntegrationVector(std::size_t eigenvectorId,
                                      std::size_t timeIntegrationVectorId);

   private:
    /// Helper function for getting the Stiffness matrix to use as basis in the
    /// eigenvalue problem.
    Mat getActualStiffnessMatrix() {
        return config_.useHermitian_ ? stiffnessMatrix_ : product_;
    }

    DGMaxEigenvalueBase::SolverConfig config_;

    Base::MeshManipulatorBase* mesh_;
    Utilities::GlobalIndexing fieldIndex_;

    // Matrices //
    //////////////
    Utilities::GlobalPetscMatrix stiffnessMatrix_;
    // Inverted
    Utilities::GlobalPetscMatrix massMatrix_;
    // Temporary storage vectors
    Utilities::GlobalPetscVector tempFieldVector_;
    Vec tempFieldVector2_;

    // Product matrix, massMatrix * stiffnessMatrix
    Mat product_;
    /// Shell matrix P*product_, where P is the projection operator
    Mat shell_;

    // Solver
    EPS solver_;

    // Phase offset shifts
    std::unique_ptr<typename DGMaxEigenvalue<DIM>::ShiftWorkspace> shifts;
    // Projector
    std::unique_ptr<typename DGMaxEigenvalue<DIM>::ProjectorWorkspace>
        projector;
    double targetFrequency_;

    LinearAlgebra::SmallVector<DIM> currentK_;
    DGMax::KPhaseShifts<DIM> stiffnessMatrixShifts_;

    // Eigenvector storage
    Utilities::Eigenpairs eigenpairs_;
    Utilities::Eigenpairs previousEigenpairs_;

    void initStiffnessMatrixShifts();
    void initMatrices();
    /// Initialize the Shell matrix
    void initStiffnessShellMatrix();
    void initSolver();
    void initEigenvectorStorage(std::size_t targetNumberOfEigenvalues);

    void shellMultiply(Vec in, Vec out);
    static PetscErrorCode staticShellMultiply(Mat mat, Vec in, Vec out);

    friend class DGMaxEigenvalue<DIM>::ProjectorWorkspace;
    friend class DGMaxEigenvalue<DIM>::ShiftWorkspace;
};

template <std::size_t DIM>
class DGMaxEigenvalue<DIM>::ShiftWorkspace {
   public:
    ShiftWorkspace(Vec example, DGMaxEigenvalueBase::SolverConfig& config);
    // Probably need to remove some default constructors
    ~ShiftWorkspace();

    void updateShiftVectors(const LinearAlgebra::SmallVector<DIM>& dk,
                            const Utilities::GlobalIndexing& fieldIndex);

    /// Apply the the phase shift (as configured) in twosided fashion
    void applyToStiffnessMatrix(Mat mat);
    /// Apply the phase shift (as configured) in single sided fashion.
    void applyToProjector(Mat mat);

   private:
    // Vectors corresponding to shifted basis functions
    LinearAlgebra::SmallVector<DIM> currentDK_;
    Vec waveVec_, waveVecConjugate_;
    DGMaxEigenvalueBase::SolverConfig& config_;
};

/// Extra workspace for the Projector that removes the zero eigenvalue subspace.
template <std::size_t DIM>
class DGMaxEigenvalue<DIM>::ProjectorWorkspace {
   public:
    explicit ProjectorWorkspace(
        DGMaxEigenvalue<DIM>::SolverWorkspace& workspace);
    // TODO Remove constructors
    ~ProjectorWorkspace();

    /// Project a vector removing the subspace corresponding to the zero
    /// eigenvalues.
    void project(Vec vec);

    /// Update the matrices for a new k-point.
    void updateKPoint(const LinearAlgebra::SmallVector<DIM>& k);

    Utilities::GlobalIndexing projectorIndex_;
    Utilities::GlobalPetscMatrix projectorMatrix_;

   private:
    void initKPhaseShifts();

    DGMaxEigenvalue<DIM>::SolverWorkspace& workspace_;

    Utilities::GlobalPetscVector tempProjectorVector_;
    /// Stiffness matrix used in the projection operator
    Mat projectionStiffness_;
    /// Solver for the projection stiffness matrix
    KSP projectionSolver_;
    DGMax::KPhaseShifts<DIM> phaseShifts_;
};

template <std::size_t DIM>
class DGMaxEigenvalue<DIM>::Result final
    : public AbstractEigenvalueResult<DIM> {

   public:
    Result(SolverWorkspace& workspace, const Base::MeshManipulator<DIM>* mesh,
           const DGMaxDiscretization<DIM>& discretization)
        : workspace_(workspace), mesh_(mesh), discretization_(discretization){};

    std::vector<double> getFrequencies() final {
        std::vector<double> frequencies(workspace_.getEigenpairs().size());
        for (std::size_t i = 0; i < frequencies.size(); ++i) {
            frequencies[i] = std::sqrt(std::abs(
                PetscRealPart(workspace_.getEigenpairs().getEigenvalue(i))));
        }
        return frequencies;
    }

    const LinearAlgebra::SmallVector<DIM>& getKPoint() const final {
        return workspace_.getKPoint();
    }

    const Base::MeshManipulator<DIM>* getMesh() const final { return mesh_; }

    void writeField(std::size_t eigenvalue,
                    Output::VTKSpecificTimeWriter<DIM>& writer) final;

    LinearAlgebra::MiddleSizeMatrix computeFieldOverlap() const {
        return workspace_.computeOverlapIntegrals();
    }

   private:
    SolverWorkspace& workspace_;
    const Base::MeshManipulator<DIM>* mesh_;
    const DGMaxDiscretization<DIM>& discretization_;
};

template <std::size_t DIM>
void DGMaxEigenvalue<DIM>::Result::writeField(
    std::size_t eigenvalue, Output::VTKSpecificTimeWriter<DIM>& writer) {
    const std::size_t VECTOR_ID = 0;
    workspace_.writeAsTimeIntegrationVector(eigenvalue, VECTOR_ID);
    // When using the Hermitian system we applied a rescaling of the
    // solution coefficients to use y = L^H x (LL^H = M is the Cholesky
    // decomposition of the mass matrix and x the actual coefficients). Undo
    // this transformation to correctly compute the fields.
    if (workspace_.getConfig().useHermitian_) {
        for (Base::Element* element : mesh_->getElementsList()) {
            // Note: this all happens inplace.
            LinearAlgebra::MiddleSizeVector& coeffs =
                element->getTimeIntegrationVector(VECTOR_ID);
            element->getElementMatrix(DGMaxDiscretizationBase::MASS_MATRIX_ID)
                .solveLowerTriangular(
                    coeffs, LinearAlgebra::Side::OP_LEFT,
                    LinearAlgebra::Transpose::HERMITIAN_TRANSPOSE);
        }
    }
    // Write the actual fields
    writer.write(
        [&](const Base::Element* element,
            const Geometry::PointReference<DIM>& pref, std::size_t) {
            const LinearAlgebra::MiddleSizeVector& coefficients =
                element->getTimeIntegrationVector(VECTOR_ID);
            auto fields =
                discretization_.computeFields(element, pref, coefficients);
            return std::sqrt(fields.realEField.l2NormSquared() +
                             fields.imagEField.l2NormSquared());
        },
        "Emag");
    writer.write(
        [&](const Base::Element* element,
            const Geometry::PointReference<DIM>& pref, std::size_t) {
            const LinearAlgebra::MiddleSizeVector& coefficients =
                element->getTimeIntegrationVector(VECTOR_ID);
            auto fields =
                discretization_.computeFields(element, pref, coefficients);
            return fields.realEField;
        },
        "Ereal");
    writer.write(
        [&](const Base::Element* element,
            const Geometry::PointReference<DIM>& pref, std::size_t) {
            const LinearAlgebra::MiddleSizeVector& coefficients =
                element->getTimeIntegrationVector(VECTOR_ID);
            auto fields =
                discretization_.computeFields(element, pref, coefficients);
            return fields.imagEField;
        },
        "Eimag");
    // Also write epsilon for post processing
    writer.write(
        [&](const Base::Element* element, const Geometry::PointReference<DIM>&,
            std::size_t) {
            auto* userData = element->getUserData();
            const ElementInfos* elementInfo =
                dynamic_cast<ElementInfos*>(userData);
            if (elementInfo != nullptr) {
                return elementInfo->epsilon_;
            } else {
                return -1.0;  // Clearly invalid value
            }
        },
        "epsilon");
}

///

template <std::size_t DIM>
void DGMaxEigenvalue<DIM>::solve(AbstractEigenvalueSolverDriver<DIM>& driver) {
    std::size_t numberOfEigenvalues = driver.getTargetNumberOfEigenvalues();

    initializeMatrices();

    SolverWorkspace workspace(config_, &mesh_, numberOfEigenvalues);

    // Setup the boundary block shifting //
    ///////////////////////////////////////

    LinearAlgebra::SmallVector<DIM> dk;  // Step in k-space from previous solve
    std::size_t expectedNumberOfSteps = driver.getNumberOfKPoints();

    std::vector<PetscScalar> eigenvalues(numberOfEigenvalues);

    for (std::size_t solve = 0; !driver.stop(); driver.nextKPoint(), ++solve) {
        DGMaxLogger(INFO, "Computing eigenvalues for k-point %/%", solve + 1,
                    expectedNumberOfSteps);
        const LinearAlgebra::SmallVector<DIM>& currentK =
            driver.getCurrentKPoint();
        workspace.updateKPoint(currentK);

        workspace.solve(driver.getTargetNumberOfEigenvalues());
        // Actual result processing
        DGMaxEigenvalue<DIM>::Result result(workspace, &mesh_, discretization_);
        driver.handleResult(result);
    }
}

// SolverWorkspace //
/////////////////////

template <std::size_t DIM>
DGMaxEigenvalue<DIM>::SolverWorkspace::SolverWorkspace(
    DGMaxEigenvalueBase::SolverConfig config, Base::MeshManipulatorBase* mesh,
    std::size_t targetNumberOfEigenvalues)
    : config_(config),
      mesh_(mesh),
      fieldIndex_(nullptr),  // Initialized in initMatrices()
      stiffnessMatrix_(fieldIndex_,
                       DGMaxDiscretizationBase::STIFFNESS_MATRIX_ID,
                       DGMaxDiscretizationBase::FACE_MATRIX_ID),
      massMatrix_(fieldIndex_, DGMaxDiscretizationBase::MASS_MATRIX_ID, -1),
      tempFieldVector_(fieldIndex_, -1, -1),
      targetFrequency_(1) {

    initMatrices();
    DGMaxLogger(INFO, "Matrices assembled");
    initStiffnessShellMatrix();
    if (config_.usesShifts()) {
        shifts = std::make_unique<DGMaxEigenvalue<DIM>::ShiftWorkspace>(
            tempFieldVector_, config_);
    }
    if (config_.useProjector_ != DGMaxEigenvalueBase::NONE) {
        projector =
            std::make_unique<DGMaxEigenvalue<DIM>::ProjectorWorkspace>(*this);
    }
    initSolver();
    initEigenvectorStorage(targetNumberOfEigenvalues);
    initStiffnessMatrixShifts();
    DGMaxLogger(INFO, "Solver workspace init completed");
}

template <std::size_t DIM>
DGMaxEigenvalue<DIM>::SolverWorkspace::~SolverWorkspace() {

    // Force cleanup of the projector & shifts workspace (if available)
    projector = nullptr;
    shifts = nullptr;

    PetscErrorCode error;
    error = VecDestroy(&tempFieldVector2_);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    error = EPSDestroy(&solver_);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    error = MatDestroy(&shell_);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    // always clean up after you are done
    if (!config_.useHermitian_) {
        error = MatDestroy(&product_);
        CHKERRABORT(PETSC_COMM_WORLD, error);
    }
}

template <std::size_t DIM>
void DGMaxEigenvalue<DIM>::SolverWorkspace::initMatrices() {
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

    // Initialize the product matrix
    if (!config_.useHermitian_) {
        PetscErrorCode error;
        error = MatMatMult(massMatrix_, stiffnessMatrix_, MAT_INITIAL_MATRIX,
                           1.0, &product_);
        CHKERRABORT(PETSC_COMM_WORLD, error);
    }
    DGMaxLogger(INFO, "DGMaxEigenvalue workspace init finished");
}

template <std::size_t DIM>
void DGMaxEigenvalue<DIM>::SolverWorkspace::initStiffnessShellMatrix() {
    PetscErrorCode error;
    PetscInt rows = fieldIndex_.getNumberOfLocalBasisFunctions();
    error = MatCreateShell(PETSC_COMM_WORLD, rows, rows, PETSC_DETERMINE,
                           PETSC_DETERMINE, this, &shell_);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = MatShellSetOperation(shell_, MATOP_MULT,
                                 (void (*)(void))staticShellMultiply);
    CHKERRABORT(PETSC_COMM_WORLD, error);
}

template <std::size_t DIM>
PetscErrorCode DGMaxEigenvalue<DIM>::SolverWorkspace::staticShellMultiply(
    Mat mat, Vec in, Vec out) {
    PetscErrorCode error;
    // TODO Check if the context is correct
    DGMaxEigenvalue<DIM>::SolverWorkspace* workspace;
    error = MatShellGetContext(mat, &workspace);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    workspace->shellMultiply(in, out);
    return 0;
}

template <std::size_t DIM>
void DGMaxEigenvalue<DIM>::SolverWorkspace::shellMultiply(Vec in, Vec out) {
    PetscErrorCode error;
    error = MatMult(getActualStiffnessMatrix(), in, out);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    if (config_.useProjector_ == DGMaxEigenvalueBase::ALL) {
        projector->project(out);
    }
}

PetscErrorCode compareEigen(PetscScalar ar, PetscScalar ai, PetscScalar br,
                            PetscScalar bi, PetscInt* res, void* ctx) {

    // Custom eigenvalue comparison looking for those whose log is closest to a
    // certain target. Negative eigenvalues should not occur, but are sorted as
    // larger (further away from the target) than any positive number.

    // Target frequency, factor 2 as the eigenvalues are frequency squared
    const double target = 2 * std::log(*(double*)ctx);

    // Documentation is unclear on whether ai and bi are zero.
    double res1 = std::abs(ar) + std::abs(ai);
    double res2 = std::abs(br) + std::abs(bi);
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

template <std::size_t DIM>
void DGMaxEigenvalue<DIM>::SolverWorkspace::initSolver() {
    PetscErrorCode err = EPSCreate(PETSC_COMM_WORLD, &solver_);
    CHKERRABORT(PETSC_COMM_WORLD, err);

    err =
        EPSSetProblemType(solver_, config_.useHermitian_ ? EPS_HEP : EPS_NHEP);
    CHKERRABORT(PETSC_COMM_WORLD, err);
    err = EPSSetWhichEigenpairs(solver_, EPS_WHICH_USER);
    CHKERRABORT(PETSC_COMM_WORLD, err);
    err = EPSSetEigenvalueComparison(solver_, compareEigen,
                                     &(this->targetFrequency_));
    CHKERRABORT(PETSC_COMM_WORLD, err);
    err = EPSSetTarget(solver_, targetFrequency_ * targetFrequency_);
    CHKERRABORT(PETSC_COMM_WORLD, err);

    // So far we have configured the the parameters of the eigenvalue solver in
    // code. This overrides these settings with the values that are in SLEPc's
    // options database (if there are any). This can be used for commandline
    // overrides of the standard values.
    err = EPSSetFromOptions(solver_);
    CHKERRABORT(PETSC_COMM_WORLD, err);
    err = EPSSetOperators(solver_, shell_, nullptr);
    CHKERRABORT(PETSC_COMM_WORLD, err);
}

template <std::size_t DIM>
void DGMaxEigenvalue<DIM>::SolverWorkspace::initEigenvectorStorage(
    std::size_t targetNumberOfEigenvalues) {
    PetscErrorCode error;
    error = VecDuplicate(tempFieldVector_, &tempFieldVector2_);
    CHKERRABORT(PETSC_COMM_WORLD, error);
}

template <std::size_t DIM>
void DGMaxEigenvalue<DIM>::SolverWorkspace::initStiffnessMatrixShifts() {
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

    stiffnessMatrixShifts_ = builder.build(fieldIndex_);
}

template <std::size_t DIM>
void DGMaxEigenvalue<DIM>::SolverWorkspace::extractEigenVectors() {
    std::swap(eigenpairs_, previousEigenpairs_);
    eigenpairs_.loadEigenpairs(solver_, tempFieldVector_);
    // Reorder
    std::vector<std::size_t> ordering(eigenpairs_.size());
    std::iota(ordering.begin(), ordering.end(), 0);
    std::sort(ordering.begin(), ordering.end(),
              [&](const std::size_t& i1, const std::size_t& i2) {
                  PetscScalar e1 = eigenpairs_.getEigenvalue(ordering[i1]);
                  PetscScalar e2 = eigenpairs_.getEigenvalue(ordering[i2]);
                  if (PetscRealPart(e1) != PetscRealPart(e2)) {
                      return PetscRealPart(e1) < PetscRealPart(e2);
                  } else {
                      return PetscImaginaryPart(e1) < PetscImaginaryPart(e2);
                  }
              });
    eigenpairs_.reorder(ordering);
}

template <std::size_t DIM>
void DGMaxEigenvalue<DIM>::SolverWorkspace::updateKPoint(
    const LinearAlgebra::SmallVector<DIM>& newK) {
    LinearAlgebra::SmallVector<DIM> dk = newK - currentK_;
    if (config_.usesShifts()) {
        shifts->updateShiftVectors(dk, fieldIndex_);
        // TODO: Check with NON HERMITIAN
        shifts->applyToStiffnessMatrix(getActualStiffnessMatrix());
        if (config_.useProjector_ != DGMaxEigenvalueBase::NONE) {
            shifts->applyToProjector(projector->projectorMatrix_);
        }
    }

    // TODO: Check with NON HERMITIAN
    stiffnessMatrixShifts_.apply(newK, getActualStiffnessMatrix());

    if (config_.useProjector_ != DGMaxEigenvalueBase::NONE) {
        projector->updateKPoint(newK);
    }

    // Note, as we use a shell matrix for EPS, we don't have to call
    // EPSSetOperators, as the matrix used has not changed.

    currentK_ = newK;
}

template <std::size_t DIM>
void DGMaxEigenvalue<DIM>::SolverWorkspace::solve(
    std::size_t targetNumberOfEigenvalues) {
    PetscInt usableInitialVectors;
    PetscErrorCode error;

    // Setup search space //
    ////////////////////////

    if (eigenpairs_.size() == 0) {
        // Generate fresh starting vector
        eigenpairs_.reserve(1, tempFieldVector_);
        DGMaxLogger(INFO, "Generating initial vector");
        error = VecSetRandom(eigenpairs_.getRawEigenvectors()[0], nullptr);
        CHKERRABORT(PETSC_COMM_WORLD, error);
    } else {
        DGMaxLogger(INFO, "Combining previous eigen vectors");
        for (PetscInt j = 1; j < eigenpairs_.size(); ++j) {
            // Some eigenvalue solvers only uses a single starting vector.
            // Mix the eigenvalue spaces from the previous k-point in the
            // hope that these are rich in the eigenvectors for the next
            // space.
            error = VecAYPX(eigenpairs_.getRawEigenvectors()[0], 1,
                            eigenpairs_.getRawEigenvectors()[j]);
            CHKERRABORT(PETSC_COMM_WORLD, error);
        }
        // Add all previous eigenvectors to the eigenvalue solver, even
        // if they are not all used.
    }
    if (config_.useProjector_ != DGMaxEigenvalueBase::NONE) {
        for (std::size_t j = 0; j < usableInitialVectors; ++j) {
            projector->project(eigenpairs_.getRawEigenvectors()[j]);
        }
        DGMaxLogger(INFO, "Projected initial vector");
    }

    // Use solution of previous time as starting point for the next one.
    error = EPSSetInitialSpace(solver_, eigenpairs_.size(),
                               eigenpairs_.getRawEigenvectors());
    CHKERRABORT(PETSC_COMM_WORLD, error);

    // Final Options //
    ///////////////////

    // Set the actual target number of eigenvalues for this solve
    error = EPSSetDimensions(solver_, targetNumberOfEigenvalues, PETSC_DECIDE,
                             PETSC_DECIDE);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    // Final overrides from the command line just before solving.
    error = EPSSetFromOptions(solver_);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    // Actual Solve //
    //////////////////

    DGMaxLogger(INFO, "Solving eigenvalue problem");

    error = EPSSetUp(solver_);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    DGMaxLogger(INFO, "Solver setup completed");

    auto start = std::chrono::high_resolution_clock::now();

    error = EPSSolve(solver_);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    // Some basic statistics
    std::chrono::duration<double> time =
        std::chrono::high_resolution_clock::now() - start;
    PetscInt numEigenvalues, iterations;
    error = EPSGetConverged(solver_, &numEigenvalues);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = EPSGetIterationNumber(solver_, &iterations);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    DGMaxLogger(INFO,
                "Eigenvalue solver stopped after % iterations with % "
                "eigenvalues in %s",
                iterations, numEigenvalues, time.count());

    // Post processing //
    /////////////////////
    extractEigenVectors();
}

template <std::size_t DIM>
void DGMaxEigenvalue<DIM>::SolverWorkspace::writeAsTimeIntegrationVector(
    std::size_t eigenvectorId, std::size_t timeIntegrationVectorId) {
    logger.assert_debug(eigenvectorId < eigenpairs_.size(),
                        "Eigenvalue % is more than converged %", eigenvectorId,
                        eigenpairs_.size());
    // Distribute the solution coefficients to the local element vectors
    PetscErrorCode err;
    err = VecCopy(eigenpairs_.getEigenvector(eigenvectorId), tempFieldVector_);
    CHKERRABORT(PETSC_COMM_WORLD, err);
    tempFieldVector_.writeTimeIntegrationVector(timeIntegrationVectorId);
}

template <std::size_t DIM>
LinearAlgebra::MiddleSizeMatrix
    DGMaxEigenvalue<DIM>::SolverWorkspace::computeOverlapIntegrals() {
    LinearAlgebra::MiddleSizeMatrix result(eigenpairs_.size(),
                                           previousEigenpairs_.size());
    PetscErrorCode err;
    for (std::size_t i = 0; i < eigenpairs_.size(); ++i) {
        err = MatMult(massMatrix_, eigenpairs_.getEigenvector(i),
                      tempFieldVector_);
        CHKERRABORT(PETSC_COMM_WORLD, err);
        for (std::size_t j = 0; j < previousEigenpairs_.size(); ++j) {
            // Note VecMDot might be a bit faster, but
            // a) this is not needed
            // b) this is complicated by the possible reordering of
            // eigenvalues.
            err = VecDot(tempFieldVector_,
                         previousEigenpairs_.getEigenvector(j), &result(i, j));
            CHKERRABORT(PETSC_COMM_WORLD, err);
        }
    }
    return result;
}

// ShiftWorkspace //
////////////////////

template <std::size_t DIM>
DGMaxEigenvalue<DIM>::ShiftWorkspace::ShiftWorkspace(
    Vec example, DGMaxEigenvalueBase::SolverConfig& config)
    : config_(config), waveVec_(nullptr), waveVecConjugate_(nullptr) {

    PetscErrorCode error;
    error = VecDuplicate(example, &waveVec_);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = VecDuplicate(example, &waveVecConjugate_);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    // No dk yet, so set it to exp(i*0) = 1
    error = VecSet(waveVec_, 1.0);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = VecSet(waveVecConjugate_, 1.0);
    CHKERRABORT(PETSC_COMM_WORLD, error);
}

template <std::size_t DIM>
DGMaxEigenvalue<DIM>::ShiftWorkspace::~ShiftWorkspace() {
    PetscErrorCode error;
    error = VecDestroy(&waveVec_);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = VecDestroy(&waveVecConjugate_);
    CHKERRABORT(PETSC_COMM_WORLD, error);
}

template <std::size_t DIM>
void DGMaxEigenvalue<DIM>::ShiftWorkspace::applyToStiffnessMatrix(Mat mat) {
    PetscErrorCode error;
    error = MatDiagonalScale(mat, waveVec_, waveVecConjugate_);
    CHKERRABORT(PETSC_COMM_WORLD, error);
}

template <std::size_t DIM>
void DGMaxEigenvalue<DIM>::ShiftWorkspace::applyToProjector(Mat mat) {
    PetscErrorCode error;
    error = MatDiagonalScale(mat, nullptr, waveVecConjugate_);
    CHKERRABORT(PETSC_COMM_WORLD, error);
}

template <std::size_t DIM>
void DGMaxEigenvalue<DIM>::ShiftWorkspace::updateShiftVectors(
    const LinearAlgebra::SmallVector<DIM>& dk,
    const Utilities::GlobalIndexing& fieldIndex) {

    if ((dk - currentDK_).l2Norm() < 1e-12) {
        return;
    }
    currentDK_ = dk;

    auto* mesh = fieldIndex.getMesh();
    PetscErrorCode err = 0;
    auto end = mesh->elementColEnd();
    for (auto it = mesh->elementColBegin(); it != end; ++it) {
        // Note this implicitly assumes we only uses DGBasisFunctions
        const std::size_t basisOffset = fieldIndex.getGlobalIndex(*it, 0);
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

// ProjectorWorkspace //
////////////////////////

template <std::size_t DIM>
DGMaxEigenvalue<DIM>::ProjectorWorkspace::ProjectorWorkspace(
    DGMaxEigenvalue<DIM>::SolverWorkspace& workspace)
    : workspace_(workspace),
      projectorIndex_(nullptr),
      projectorMatrix_(projectorIndex_, workspace.fieldIndex_,
                       DGMaxDiscretizationBase::PROJECTOR_MATRIX_ID),
      tempProjectorVector_(projectorIndex_, -1, -1),
      projectionStiffness_(nullptr),
      projectionSolver_(nullptr) {

    // Setup vectors
    std::vector<std::size_t> projectorUnknowns({1});
    projectorIndex_.reset(workspace.mesh_,
                          Utilities::GlobalIndexing::Layout::BLOCKED_PROCESSOR,
                          &projectorUnknowns);
    projectorMatrix_.reinit();
    tempProjectorVector_.reinit();

    PetscErrorCode error;
    error = MatCreate(PETSC_COMM_WORLD, &projectionStiffness_);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    // Create the required solver
    error = KSPCreate(PETSC_COMM_WORLD, &projectionSolver_);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = KSPSetType(projectionSolver_, KSPPREONLY);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    PC pc;
    error = KSPGetPC(projectionSolver_, &pc);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = PCSetType(pc, PCLU);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    initKPhaseShifts();
}

template <std::size_t DIM>
DGMaxEigenvalue<DIM>::ProjectorWorkspace::~ProjectorWorkspace() {
    PetscErrorCode error;
    error = KSPDestroy(&projectionSolver_);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = MatDestroy(&projectionStiffness_);
    CHKERRABORT(PETSC_COMM_WORLD, error);
}

template <std::size_t DIM>
void DGMaxEigenvalue<DIM>::ProjectorWorkspace::project(Vec vec) {
    logger.assert_always(
        workspace_.config_.useProjector_ != DGMaxEigenvalueBase::NONE,
        "Projecting without projector");
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
    if (workspace_.config_.useHermitian_) {
        // No need to multiply with the mass matrix in the Hermitian case
        // u = u + B^H (-t)
        // error = MatMultHermitianTransposeAdd(projectorMatrix_,
        //                                      tempProjectorVector_, vec, vec);
        error =
            MatMultHermitianTranspose(projectorMatrix_, tempProjectorVector_,
                                      workspace_.tempFieldVector2_);
        CHKERRABORT(PETSC_COMM_WORLD, error);
    } else {
        // v = B^H (-t)
        error =
            MatMultHermitianTranspose(projectorMatrix_, tempProjectorVector_,
                                      workspace_.tempFieldVector_);
        CHKERRABORT(PETSC_COMM_WORLD, error);
        // Finally compute u + M^{-1}v
        // error = MatMultAdd(massMatrix_, tempFieldVector_, vec, vec);
        error = MatMult(workspace_.massMatrix_, workspace_.tempFieldVector_,
                        workspace_.tempFieldVector2_);
        CHKERRABORT(PETSC_COMM_WORLD, error);
    }
    error = VecAXPY(vec, 1.0, workspace_.tempFieldVector2_);
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

template <std::size_t DIM>
void DGMaxEigenvalue<DIM>::ProjectorWorkspace::updateKPoint(
    const LinearAlgebra::SmallVector<DIM>& k) {
    // Update the matrix
    phaseShifts_.apply(k, projectorMatrix_);

    // Update the KSP & the inner projectionMatrix.
    PetscErrorCode error;
    Mat projectionH;
    error = MatHermitianTranspose(projectorMatrix_, MAT_INITIAL_MATRIX,
                                  &projectionH);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    // Note MAT_REUSE_MATRIX has been used in the following matrix
    // creation but caused segfaults. Probably from slightly different
    // sparsity patterns. Hence we just destroy the previous matrix.
    error = MatDestroy(&projectionStiffness_);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    if (workspace_.config_.useHermitian_) {
        // No mass matrix needed.
        error = MatMatMult(projectorMatrix_, projectionH, MAT_INITIAL_MATRIX,
                           PETSC_DEFAULT, &projectionStiffness_);
        CHKERRABORT(PETSC_COMM_WORLD, error);
    } else {
        error = MatMatMatMult(projectorMatrix_, workspace_.massMatrix_,
                              projectionH, MAT_INITIAL_MATRIX, PETSC_DEFAULT,
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

template <std::size_t DIM>
void DGMaxEigenvalue<DIM>::ProjectorWorkspace::initKPhaseShifts() {
    DGMax::CGDGMatrixKPhaseShiftBuilder<DIM> projectorBuilder;

    projectorBuilder.setMatrixExtractor([&](const Base::Element* element) {
        return element->getElementMatrix(
            DGMaxDiscretizationBase::PROJECTOR_MATRIX_ID);
    });

    if (workspace_.config_.usesShifts()) {
        projectorBuilder.setExtraShift([&](const Base::Element* element) {
            const Geometry::PointReference<DIM>& center =
                element->getReferenceGeometry()->getCenter();
            // Rescaling of the columns, as the shifts for the field basis
            // functions are e^{-ikx}. Where the convention is used that the
            // trial basis functions use the - sign.
            return element->referenceToPhysical(center).getCoordinates() *
                   workspace_.config_.shiftFactor_;
        });
    }
    phaseShifts_ =
        projectorBuilder.build(projectorIndex_, workspace_.fieldIndex_);
}

template class DGMaxEigenvalue<2>;
template class DGMaxEigenvalue<3>;
