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

#include "Utilities/Eigenpairs.h"
#include "Utilities/GlobalMatrix.h"
#include "Utilities/GlobalVector.h"

#include "Utils/FaceKPhaseShiftBuilder.h"
#include "ElementInfos.h"

#include "DGMaxLogger.h"

using namespace hpgem;

/// Internal storage for the algorithm
template <std::size_t DIM>
class DivDGMaxEigenvalue<DIM>::SolverWorkspace {
   public:
    explicit SolverWorkspace(Base::MeshManipulator<DIM>* mesh)
        : indexing_(nullptr),
          stiffnessMatrix_(
              indexing_, DivDGMaxDiscretizationBase::STIFFNESS_MATRIX_ID,
              DivDGMaxDiscretizationBase::FACE_STIFFNESS_MATRIX_ID),
          massMatrix_(indexing_, DivDGMaxDiscretizationBase::MASS_MATRIX_ID,
                      -1),
          tempVector_(indexing_, -1, -1),
          kderivativeMat_(nullptr),
          solver_(nullptr) {
        // Separate from initializer list to allow for more flexibility
        init(mesh);
    };

    ~SolverWorkspace() {
        PetscErrorCode err;
        err = EPSDestroy(&solver_);
        CHKERRABORT(PETSC_COMM_WORLD, err);
        MatDestroy(&kderivativeMat_);
    }

    void init(Base::MeshManipulator<DIM>* mesh) {
        DGMaxLogger(DEBUG, "Initializing global index");
        indexing_.reset(mesh, Utilities::GlobalIndexing::BLOCKED_PROCESSOR);
        DGMaxLogger(DEBUG, "Assembling matrices");
        stiffnessMatrix_.reinit();
        massMatrix_.reinit();
        DGMaxLogger(DEBUG, "Building temporary global vector");
        tempVector_.reinit();

        PetscErrorCode error = MatDuplicate(
            stiffnessMatrix_, MAT_DO_NOT_COPY_VALUES, &kderivativeMat_);
        CHKERRABORT(PETSC_COMM_WORLD, error);

        initKShifts();
        initSolver();
    }

    void solve(LinearAlgebra::SmallVector<DIM> k,
               std::size_t numberOfEigenvalues) {

        kphaseshifts_.apply(k, stiffnessMatrix_);
        kpoint_ = k;

        PetscErrorCode err;
        err = EPSSetOperators(solver_, stiffnessMatrix_, massMatrix_);
        CHKERRABORT(PETSC_COMM_WORLD, err);
        err = EPSSetDimensions(solver_, numberOfEigenvalues, PETSC_DECIDE,
                               PETSC_DECIDE);
        CHKERRABORT(PETSC_COMM_WORLD, err);

        EPSSetInitialSpace(solver_, eigenpairs_.size(),
                           eigenpairs_.getRawEigenvectors());

        // Allow for overrides
        err = EPSSetFromOptions(solver_);
        CHKERRABORT(PETSC_COMM_WORLD, err);

        DGMaxLogger(INFO, "Setting up solve");
        err = EPSSetUp(solver_);
        CHKERRABORT(PETSC_COMM_WORLD, err);

        DGMaxLogger(INFO, "Solving");
        err = EPSSolve(solver_);
        CHKERRABORT(PETSC_COMM_WORLD, err);

        // Post solve
        PetscInt iterations;
        err = EPSGetIterationNumber(solver_, &iterations);
        CHKERRABORT(PETSC_COMM_WORLD, err);
        std::swap(eigenpairs_, previousEigenpairs_);
        eigenpairs_.loadEigenpairs(solver_, tempVector_);
        DGMaxLogger(INFO, "Number of eigenvalues % converged in % iterations",
                    eigenpairs_.size(), iterations);
    }

    const LinearAlgebra::SmallVector<DIM>& getKPoint() const { return kpoint_; }

    const Utilities::Eigenpairs& getEigenpairs() const { return eigenpairs_; }

    /**
     * Write the (global) eigenvector to a (local) time integration vector, for
     * example for plotting the solution.
     * @param eigenvalue The index of the eigenvalue
     * @param timeIntegrationVector The index of the time integration vector
     */
    void writeEigenvectorAsTimeIntegrationVector(
        std::size_t eigenvalue, std::size_t timeIntegrationVector) {
        logger.assert_debug(eigenvalue < eigenpairs_.size(),
                            "Eigenvalue index too large");
        PetscErrorCode err;
        err = VecCopy(eigenpairs_.getEigenvector(eigenvalue), tempVector_);
        CHKERRABORT(PETSC_COMM_WORLD, err);
        tempVector_.writeTimeIntegrationVector(timeIntegrationVector);
    }

    LinearAlgebra::MiddleSizeMatrix computeOverlapIntegrals() {
        LinearAlgebra::MiddleSizeMatrix result(eigenpairs_.size(),
                                               previousEigenpairs_.size());
        PetscErrorCode err;
        for (std::size_t i = 0; i < eigenpairs_.size(); ++i) {
            err = MatMult(massMatrix_, eigenpairs_.getEigenvector(i),
                          tempVector_);
            CHKERRABORT(PETSC_COMM_WORLD, err);
            for (std::size_t j = 0; j < previousEigenpairs_.size(); ++j) {
                // Note VecMDot might be a bit faster, but
                // a) this is not needed
                // b) this is complicated by the possible reordering of
                // eigenvalues.
                err = VecDot(tempVector_, previousEigenpairs_.getEigenvector(j),
                             &result(i, j));
                CHKERRABORT(PETSC_COMM_WORLD, err);
            }
        }
        return result;
    }

    std::array<LinearAlgebra::MiddleSizeMatrix, DIM>
        computeWaveVectorDerivatives() {
        PetscErrorCode err;
        std::array<LinearAlgebra::MiddleSizeMatrix, DIM> result;
        std::size_t numEigenvectors = eigenpairs_.size();
        LinearAlgebra::SmallVector<DIM> dk;
        for (std::size_t kdir = 0; kdir < DIM; ++kdir) {
            result[kdir].resize(numEigenvectors, numEigenvectors);
            // Setup derivative matrix in the i-th direction
            // Assumes that the kderivativeMat_ is only used for these
            // derivatives
            dk.set(0.0);
            dk[kdir] = 1.0;
            err = MatZeroEntries(kderivativeMat_);
            CHKERRABORT(PETSC_COMM_WORLD, err);
            MPI_Barrier(PETSC_COMM_WORLD);
            kphaseshifts_.applyDeriv(kpoint_, dk, kderivativeMat_);
            for(NormType nt : {NORM_1, NORM_FROBENIUS, NORM_INFINITY}) {
                PetscReal normval;
                MatNorm(kderivativeMat_, nt, &normval);
                DGMaxLogger(INFO, "Norm %: %", NormTypes[nt], normval);
            }


            MPI_Barrier(PETSC_COMM_WORLD);
            for (std::size_t i = 0; i < numEigenvectors; ++i) {
                err = MatMult(kderivativeMat_, eigenpairs_.getEigenvector(i),
                        tempVector_);
                CHKERRABORT(PETSC_COMM_WORLD, err);
                for (std::size_t j = i; j < numEigenvectors; ++j) {
                    std::complex<double> derivative;
                    err = VecDot(tempVector_, eigenpairs_.getEigenvector(j),
                           &derivative);
                    CHKERRABORT(PETSC_COMM_WORLD, err);
                    result[kdir](i, j) = derivative;
                    result[kdir](j, i) = std::conj(derivative);

                    if (i <= 1 && j <= 1) {
                        std::cout << "Term " << i << j << ":" << derivative << std::endl;
                    }
                }
            }
        }
        return result;
    }

   private:
    void initKShifts() {
        DGMaxLogger(VERBOSE, "Initializing boundary shifting");
        DGMax::FaceMatrixKPhaseShiftBuilder<DIM> builder;
        builder.setMatrixExtractor([&](const Base::Face* face) {
            const Base::FaceMatrix& faceMatrix = face->getFaceMatrix(
                DivDGMaxDiscretizationBase::FACE_STIFFNESS_MATRIX_ID);
            LinearAlgebra::MiddleSizeMatrix block1, block2;
            block1 = faceMatrix.getElementMatrix(Base::Side::LEFT,
                                                 Base::Side::RIGHT);
            block2 = faceMatrix.getElementMatrix(Base::Side::RIGHT,
                                                 Base::Side::LEFT);

            return std::make_pair(block1, block2);
        });
        kphaseshifts_ = builder.build(indexing_);
    }

    void initSolver() {
        PetscErrorCode err;
        err = EPSCreate(PETSC_COMM_WORLD, &solver_);
        CHKERRABORT(PETSC_COMM_WORLD, err);
        // Setting operators, but these will be set again after k-shifting
        err = EPSSetOperators(solver_, stiffnessMatrix_, massMatrix_);
        CHKERRABORT(PETSC_COMM_WORLD, err);

        // The mass matrix is positive semidefinite, the stiffness matrix is
        // Hermitian indefinite (k-phase-shifts are also Hermitian for real k).
        err = EPSSetProblemType(solver_, EPS_GHEP);
        CHKERRABORT(PETSC_COMM_WORLD, err);

        // By default configure to use shift & invert targeted to 0
        err = EPSSetWhichEigenpairs(solver_, EPS_TARGET_REAL);
        CHKERRABORT(PETSC_COMM_WORLD, err);
        err = EPSSetTarget(solver_, 0.0);
        CHKERRABORT(PETSC_COMM_WORLD, err);
        {
            ST st;
            err = EPSGetST(solver_, &st);
            CHKERRABORT(PETSC_COMM_WORLD, err);
            err = STSetType(st, STSINVERT);
            CHKERRABORT(PETSC_COMM_WORLD, err);
        }

        DGMaxLogger(INFO, "Eigenvalue solver configured");
    }

    Utilities::GlobalIndexing indexing_;
    Utilities::GlobalPetscMatrix stiffnessMatrix_;
    Utilities::GlobalPetscMatrix massMatrix_;
    /// Vector used as intermediate storage and sample
    Utilities::GlobalPetscVector tempVector_;

    /// exp(ikx) shifts used in the stiffness matrix
    DGMax::KPhaseShifts<DIM> kphaseshifts_;

    LinearAlgebra::SmallVector<DIM> kpoint_;

    /// Matrix used in computing the wave vector derivatives
    Mat kderivativeMat_;

    /// Eigenvalue solver
    EPS solver_;

    Utilities::Eigenpairs eigenpairs_;
    Utilities::Eigenpairs previousEigenpairs_;
};

template <std::size_t DIM>
class DivDGMaxEigenvalue<DIM>::Result final
    : public AbstractEigenvalueResult<DIM> {
   public:
    Result(typename DivDGMaxEigenvalue<DIM>::SolverWorkspace& workspace,
           const Base::MeshManipulator<DIM>* mesh,
           const DivDGMaxDiscretization<DIM>& discretization)
        : workspace_(workspace), mesh_(mesh), discretization_(discretization) {
        const auto& eigenpairs = workspace.getEigenpairs();
        frequencies_.resize(eigenpairs.size());
        for (std::size_t i = 0; i < frequencies_.size(); ++i) {
            frequencies_[i] =
                std::sqrt(PetscRealPart(eigenpairs.getEigenvalue(i)));
        }
    };

    std::vector<double> getFrequencies() final { return frequencies_; }

    const LinearAlgebra::SmallVector<DIM>& getKPoint() const final {
        return workspace_.getKPoint();
    }

    void writeField(std::size_t eigenvalue,
                    Output::VTKSpecificTimeWriter<DIM>& writer) final;

    const Base::MeshManipulator<DIM>* getMesh() const final { return mesh_; }

    LinearAlgebra::MiddleSizeMatrix computeFieldOverlap() const final {
        return workspace_.computeOverlapIntegrals();
    }

    bool supportsWaveVectorDerivatives() const final { return true; }
    std::array<LinearAlgebra::MiddleSizeMatrix, DIM>
        computeWaveVectorDerivatives() const override {
        return workspace_.computeWaveVectorDerivatives();
    }

   private:
    typename DivDGMaxEigenvalue<DIM>::SolverWorkspace& workspace_;
    const Base::MeshManipulator<DIM>* mesh_;
    const DivDGMaxDiscretization<DIM>& discretization_;
    std::vector<double> frequencies_;
};

template <std::size_t DIM>
void DivDGMaxEigenvalue<DIM>::Result::writeField(
    std::size_t eigenvalue, Output::VTKSpecificTimeWriter<DIM>& writer) {
    const std::size_t VECTOR_ID = 0;
    // Write to a time integration vector so that we can access it on each
    // element.
    workspace_.writeEigenvectorAsTimeIntegrationVector(eigenvalue, VECTOR_ID);
    // write all field components
    discretization_.writeFields(writer, VECTOR_ID);
}

template <std::size_t DIM>
DivDGMaxEigenvalue<DIM>::DivDGMaxEigenvalue(
    Base::MeshManipulator<DIM>& mesh, std::size_t order,
    DivDGMaxDiscretizationBase::Stab stab)
    : mesh_(mesh), discretization(order_, stab), order_(order), stab_(stab) {}

template <std::size_t DIM>
void DivDGMaxEigenvalue<DIM>::solve(
    AbstractEigenvalueSolverDriver<DIM>& driver) {
    // Sometimes the solver finds more eigenvalues & vectors than requested, so
    // reserve some extra space for them.
    std::size_t numberOfEigenvalues = driver.getTargetNumberOfEigenvalues();
    const PetscInt numberOfEigenVectors =
        std::max(2 * numberOfEigenvalues, numberOfEigenvalues + 10);

    PetscErrorCode error;
    DGMaxLogger(INFO, "Starting assembly");

    // Using dispersive materials is not supported
    const double dispersionFrequency =
        std::numeric_limits<double>::signaling_NaN();
    discretization.initializeBasisFunctions(mesh_);
    discretization.computeElementIntegrals(mesh_, {}, dispersionFrequency);
    discretization.computeFaceIntegrals(mesh_, {}, dispersionFrequency);

    SolverWorkspace workspace(&mesh_);

    std::size_t outputId = 0;
    std::size_t expectedNumberOfSteps = driver.getNumberOfKPoints();

    DGMaxLogger(INFO, "Starting k-vector walk");
    for (std::size_t solve = 0; !driver.stop(); driver.nextKPoint(), ++solve) {
        LinearAlgebra::SmallVector<DIM> currentK = driver.getCurrentKPoint();
        numberOfEigenvalues = driver.getTargetNumberOfEigenvalues();
        DGMaxLogger(INFO, "Solving for k-vector %/%", solve + 1,
                    expectedNumberOfSteps);
        workspace.solve(currentK, numberOfEigenvalues);

        Result result(workspace, &mesh_, discretization);

        driver.handleResult(result);

        // globalVector.writeTimeIntegrationVector(outputId);
        outputId++;
    }
    DGMaxLogger(INFO, "Finished k-vector walk");
}

template class DivDGMaxEigenvalue<2>;
template class DivDGMaxEigenvalue<3>;
