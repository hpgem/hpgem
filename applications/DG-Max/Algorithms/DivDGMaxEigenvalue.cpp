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

#include "Utils/FaceKPhaseShiftBuilder.h"

#include "DGMaxLogger.h"

using namespace hpgem;

/// Internal storage for the algorithm
template <std::size_t DIM>
struct Workspace {
    explicit Workspace(Base::MeshManipulator<DIM>* mesh)
        : indexing_(nullptr),
          stiffnessMatrix_(
              indexing_,
              DivDGMaxDiscretization<DIM>::ELEMENT_STIFFNESS_MATRIX_ID,
              DivDGMaxDiscretization<DIM>::FACE_STIFFNESS_MATRIX_ID),
          massMatrix_(indexing_,
                      DivDGMaxDiscretization<DIM>::ELEMENT_MASS_MATRIX_ID, -1),
          tempVector_(indexing_, -1, -1),
          solver_(nullptr),
          numberOfConvergedEigenpairs(0),
          numberOfEigenvectors_(0),
          eigenvectors_(nullptr) {
        // Separate from initializer list to allow for more flexibility
        init(mesh);
    };

    ~Workspace() {
        PetscErrorCode err;

        VecDestroyVecs(numberOfEigenvectors_, &eigenvectors_);
        err = EPSDestroy(&solver_);
        CHKERRABORT(PETSC_COMM_WORLD, err);
    }

    void init(Base::MeshManipulator<DIM>* mesh) {
        DGMaxLogger(DEBUG, "Initializing global index");
        indexing_.reset(mesh, Utilities::GlobalIndexing::BLOCKED_PROCESSOR);
        DGMaxLogger(DEBUG, "Assembling matrices");
        stiffnessMatrix_.reinit();
        massMatrix_.reinit();
        DGMaxLogger(DEBUG, "Building temporary global vector");
        tempVector_.reinit();

        initKShifts();
        initSolver();
    }

    void solve(LinearAlgebra::SmallVector<DIM> k,
               std::size_t numberOfEigenvalues) {

        kphaseshifts_.apply(k, stiffnessMatrix_);

        PetscErrorCode err;
        err = EPSSetOperators(solver_, massMatrix_, stiffnessMatrix_);
        CHKERRABORT(PETSC_COMM_WORLD, err);
        err = EPSSetDimensions(solver_, numberOfEigenvalues, PETSC_DECIDE,
                               PETSC_DECIDE);
        CHKERRABORT(PETSC_COMM_WORLD, err);

        EPSSetInitialSpace(solver_, numberOfConvergedEigenpairs, eigenvectors_);

        // Allow for overrides
        err = EPSSetFromOptions(solver_);
        CHKERRABORT(PETSC_COMM_WORLD, err);

        DGMaxLogger(INFO, "Setting up solve");
        err = EPSSetUp(solver_);
        CHKERRABORT(PETSC_COMM_WORLD, err);
        EPSSetWhichEigenpairs(solver_, EPS_LARGEST_REAL);
        DGMaxLogger(INFO, "Solving");

        err = EPSSolve(solver_);
        CHKERRABORT(PETSC_COMM_WORLD, err);

        // Post solve
        err = EPSGetConverged(solver_, &numberOfConvergedEigenpairs);
        CHKERRABORT(PETSC_COMM_WORLD, err);
        DGMaxLogger(INFO, "Number of eigenvalues %", numberOfEigenvalues);
        if (numberOfConvergedEigenpairs > numberOfEigenvectors_) {
            err = VecDestroyVecs(numberOfEigenvectors_, &eigenvectors_);
            CHKERRABORT(PETSC_COMM_WORLD, err);
            err = VecDuplicateVecs(tempVector_, numberOfConvergedEigenpairs,
                                   &eigenvectors_);
            CHKERRABORT(PETSC_COMM_WORLD, err);
        }

        eigenvalues_.resize(numberOfConvergedEigenpairs);
        for (PetscInt i = 0; i < numberOfConvergedEigenpairs; ++i) {
            err = EPSGetEigenpair(solver_, i, &eigenvalues_[i], nullptr,
                                  eigenvectors_[i], nullptr);
            CHKERRABORT(PETSC_COMM_WORLD, err);
        }
    }

   private:
    void initKShifts() {
        DGMaxLogger(VERBOSE, "Initializing boundary shifting");
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
        kphaseshifts_ = builder.build(indexing_);
    }

    void initSolver() {
        PetscErrorCode err;
        err = EPSCreate(PETSC_COMM_WORLD, &solver_);
        CHKERRABORT(PETSC_COMM_WORLD, err);
        // Setting operators, but these will be set again after k-shifting
        // Note the order Mx = (1/omega^2) Sx and NOT Sx = omega^2 M x and doing
        // a spectral transformation.
        // It may be possible to remove this relic with the current version of
        // SLEPc.
        err = EPSSetOperators(solver_, massMatrix_, stiffnessMatrix_);
        CHKERRABORT(PETSC_COMM_WORLD, err);

        err = EPSSetWhichEigenpairs(solver_, EPS_LARGEST_MAGNITUDE);
        CHKERRABORT(PETSC_COMM_WORLD, err);
        DGMaxLogger(INFO, "Eigenvalue solver configured");
    }

    Utilities::GlobalIndexing indexing_;
    Utilities::GlobalPetscMatrix stiffnessMatrix_;
    Utilities::GlobalPetscMatrix massMatrix_;
    /// Vector used as intermediate storage and sample
    Utilities::GlobalPetscVector tempVector_;

    /// exp(ikx) shifts used in the stiffness matrix
    DGMax::KPhaseShifts<DIM> kphaseshifts_;

    /// Eigenvalue solver
    EPS solver_;

    /// Number of converged eigenpairs
    PetscInt numberOfConvergedEigenpairs;
    /// Number of vectors available in eigenvectors_, this may be larger than
    /// the number of converged eigenvalues.
    std::size_t numberOfEigenvectors_;
    /// Storage for eigenvectors
    Vec* eigenvectors_;
    /// Storage for the most recent eigenvalues
    std::vector<PetscScalar> eigenvalues_;
};

template <std::size_t DIM>
DivDGMaxEigenvalue<DIM>::DivDGMaxEigenvalue(
    Base::MeshManipulator<DIM>& mesh, std::size_t order,
    typename DivDGMaxDiscretization<DIM>::Stab stab)
    : mesh_(mesh), order_(order), stab_(stab) {}

template <std::size_t DIM>
struct DivDGMaxResult : AbstractEigenvalueResult<DIM> {
    std::vector<double> getFrequencies() final { return frequencies_; }

    const LinearAlgebra::SmallVector<DIM>& getKPoint() const final {
        return kpoint_;
    }

    std::vector<double> frequencies_;
    LinearAlgebra::SmallVector<DIM> kpoint_;
};

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
    discretization.initializeBasisFunctions(mesh_, order_);
    discretization.computeElementIntegrands(mesh_, false, nullptr, nullptr,
                                            nullptr);
    discretization.computeFaceIntegrals(mesh_, nullptr, stab_);

    Workspace<DIM> workspace(&mesh_);

    std::size_t outputId = 0;
    std::size_t expectedNumberOfSteps = driver.getNumberOfKPoints();
    DivDGMaxResult<DIM> result;

    DGMaxLogger(INFO, "Starting k-vector walk");
    for (std::size_t solve = 0; !driver.stop(); driver.nextKPoint(), ++solve) {
        LinearAlgebra::SmallVector<DIM> currentK = driver.getCurrentKPoint();
        numberOfEigenvalues = driver.getTargetNumberOfEigenvalues();
        DGMaxLogger(INFO, "Solving for k-vector %/%", solve + 1,
                    expectedNumberOfSteps);
        workspace.solve(currentK, numberOfEigenvalues);

        result.kpoint_ = currentK;
        result.frequencies_.resize(workspace.numberOfConvergedEigenpairs);
        for (std::size_t i = 0; i < workspace.numberOfConvergedEigenpairs;
             ++i) {
            result.frequencies_[i] =
                std::sqrt(1 / workspace.eigenvalues_[i].real());
        }

        driver.handleResult(result);

        // globalVector.writeTimeIntegrationVector(outputId);
        outputId++;
    }
    DGMaxLogger(INFO, "Finished k-vector walk");
}

template class DivDGMaxEigenvalue<2>;
template class DivDGMaxEigenvalue<3>;
