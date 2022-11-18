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

#include "HarmonicSolver.h"

#include <petscksp.h>

#include <complex>
#include "Output/TecplotDiscontinuousSolutionWriter.h"
#include "Utilities/GlobalMatrix.h"
#include "Utilities/GlobalVector.h"

using namespace hpgem;

namespace DGMax {

template <std::size_t DIM>
class HarmonicSolver<DIM>::Result : public AbstractHarmonicResult<DIM> {
   public:
    Result(const HarmonicProblem<DIM>& problem,
           HarmonicSolver<DIM>::Workspace& workspace)
        : problem_(&problem), workspace_(&workspace){};

    const HarmonicProblem<DIM>& solvedProblem() final { return *problem_; }

    void writeVTK(Output::VTKSpecificTimeWriter<DIM>& output) final {
        workspace_->writeVTK(output);
    }
    double computeL2Error(const ExactHarmonicProblem<DIM>& solution) final {
        return workspace_->computeL2Error(solution);
    }
    Base::MeshManipulator<DIM>& getMesh() final {
        return workspace_->getMesh();
    }

    LinearAlgebra::SmallVector<4> computeEnergyFlux(
        Base::Face& face, hpgem::Base::Side side, double wavenumber,
        const FieldPattern<DIM>* background) final {
        return workspace_->computeEnergyFlux(face, side, wavenumber,
                                             background);
    }

    virtual LinearAlgebra::SmallVectorC<DIM> computeField(
        const Base::Element* element,
        const Geometry::PointReference<DIM>& p) final {
        return workspace_->computeField(element, p);
    }

    virtual LinearAlgebra::SmallVectorC<DIM> computeFieldCurl(
        const Base::Element* element,
        const Geometry::PointReference<DIM>& p) final {
        return workspace_->computeFieldCurl(element, p);
    }

    double computeFieldL2Integral(Base::Face& face, Base::Side side) final {
        return workspace_->computeFieldL2Integral(face, side);
    }

   private:
    const HarmonicProblem<DIM>* problem_;
    HarmonicSolver<DIM>::Workspace* workspace_;
};

template <std::size_t DIM>
class HarmonicSolver<DIM>::Workspace {
   public:
    Workspace(AbstractDiscretization<DIM>& discretization,
              Base::MeshManipulator<DIM>& mesh);
    ~Workspace();

    /**
     * (re) compute the element and face integrands
     */
    void computeIntegrals(AbstractHarmonicSolverDriver<DIM>& driver);

    /**
     * Assemble the solverMatrix out of the individual matrices
     * @param waveNumber The wave number
     */
    void assembleSolverMatrix(double waveNumber);

    /**
     * Runs the solver:
     *  - Assumes that the solverMatrix and loadVector has been set up
     *  - Result is written to resultVector_
     *  - resultVector_ is distributed
     */
    void solve();

    // Exposition to the solution
    Base::MeshManipulator<DIM>& getMesh() { return *mesh_; }

    LinearAlgebra::SmallVectorC<DIM> computeField(
        const Base::Element* element,
        const Geometry::PointReference<DIM>& p) const {
        return discretization_->computeField(
            element, p, element->getTimeIntegrationVector(VECTOR_ID));
    }

    LinearAlgebra::SmallVectorC<DIM> computeFieldCurl(
        const Base::Element* element,
        const Geometry::PointReference<DIM>& p) const {
        return discretization_->computeCurlField(
            element, p, element->getTimeIntegrationVector(VECTOR_ID));
    }

    double computeL2Error(const ExactHarmonicProblem<DIM>& solution) {
        return discretization_->computeL2Error(
            *mesh_, VECTOR_ID,
            [&solution](const Base::Element&,
                        const Geometry::PointPhysical<DIM>& p) {
                return solution.exactSolution(p);
            });
    }

    void writeVTK(Output::VTKSpecificTimeWriter<DIM>& output) {
        discretization_->writeFields(output, VECTOR_ID);
    }

    LinearAlgebra::SmallVector<4> computeEnergyFlux(
        Base::Face& face, hpgem::Base::Side side, double wavenumber,
        const FieldPattern<DIM>* background) {
        return discretization_->computeEnergyFluxes(face, side, wavenumber,
                                                    VECTOR_ID, background);
    }

    double computeFieldL2Integral(Base::Face& face, Base::Side side) {
        return discretization_->computeFieldL2Integral(face, side, VECTOR_ID);
    }

   private:
    static constexpr const std::size_t VECTOR_ID = 0;

    void configureSolver();
    bool checkDispersion();

    AbstractDiscretization<DIM>* discretization_;
    Base::MeshManipulator<DIM>* mesh_;
    /**
     * Whether any material in the mesh is frequency dependent
     *
     * Note for distributed computations this will include the for all
     * processors.
     */
    bool dispersion_;

    Utilities::GlobalIndexing indexing_;
    Utilities::GlobalPetscMatrix massMatrix_;
    Utilities::GlobalPetscMatrix stiffnessMatrix_;
    Utilities::GlobalPetscMatrix stiffnessImpedanceMatrix_;
    Utilities::GlobalPetscVector resultVector_;
    Utilities::GlobalPetscVector loadVector_;

    Mat solverMatrix_;
    KSP solver_;
};

template <std::size_t DIM>
HarmonicSolver<DIM>::Workspace::Workspace(
    AbstractDiscretization<DIM>& discretization,
    Base::MeshManipulator<DIM>& mesh)
    : discretization_(&discretization),
      mesh_(&mesh),
      dispersion_(checkDispersion()),
      indexing_(nullptr),
      massMatrix_(indexing_, AbstractDiscretizationBase::MASS_MATRIX_ID, -1),
      stiffnessMatrix_(indexing_,
                       AbstractDiscretizationBase::STIFFNESS_MATRIX_ID,
                       AbstractDiscretizationBase::FACE_STIFFNESS_MATRIX_ID),
      stiffnessImpedanceMatrix_(
          indexing_, -1, AbstractDiscretizationBase::FACE_IMPEDANCE_MATRIX_ID),
      resultVector_(indexing_, -1, -1),
      loadVector_(indexing_, AbstractDiscretizationBase::ELEMENT_VECTOR_ID,
                  AbstractDiscretizationBase::FACE_VECTOR_ID),
      solverMatrix_(nullptr),
      solver_(nullptr) {
    // Initialize the basis functions and indices after creating the matrices,
    // so that the matrices are not yet assembled.
    discretization_->initializeBasisFunctions(*mesh_);
    indexing_.reset(mesh_, Utilities::GlobalIndexing::BLOCKED_PROCESSOR);
    resultVector_.reinit();

    PetscErrorCode err;
    err = KSPCreate(PETSC_COMM_WORLD, &solver_);
    CHKERRABORT(PETSC_COMM_WORLD, err);

    configureSolver();
}

template <std::size_t DIM>
HarmonicSolver<DIM>::Workspace::~Workspace() {
    if (solverMatrix_ != nullptr) {
        MatDestroy(&solverMatrix_);
    }
    KSPDestroy(&solver_);
}

template <std::size_t DIM>
void HarmonicSolver<DIM>::Workspace::configureSolver() {
    PetscErrorCode error;

    error = KSPSetTolerances(solver_, 1e-8, PETSC_DEFAULT, PETSC_DEFAULT,
                             PETSC_DEFAULT);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    error = KSPSetFromOptions(solver_);
    CHKERRABORT(PETSC_COMM_WORLD, error);
}

template <std::size_t DIM>
bool HarmonicSolver<DIM>::Workspace::checkDispersion() {
    char has_dispersion = 0;
    // NOTE: If a local dispersion is desired, this may have to change to use
    // IteratorType::GLOBAL
    for (const Base::Element* element : mesh_->getElementsList()) {
        if (ElementInfos::get(*element).isDispersive()) {
            has_dispersion = 1;
            break;
        }
    }

#ifdef HPGEM_USE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &has_dispersion, 1, MPI_CHAR, MPI_MAX,
                  PETSC_COMM_WORLD);
#endif
    return has_dispersion;
}

template <std::size_t DIM>
void HarmonicSolver<DIM>::Workspace::computeIntegrals(
    AbstractHarmonicSolverDriver<DIM>& driver) {

    logger.assert_always(mesh_ != nullptr,
                         "Computing integrals without a mesh");

    const HarmonicProblem<DIM>& problem = driver.currentProblem();

    // Note: There are some possible improvements by splitting the dispersion
    // property into local (known elements of the mesh) and global (all the
    // elements).
    bool fullRecompute =
        driver.hasChanged(HarmonicProblemChanges::BOUNDARY_CONDITION_TYPE) ||
        (dispersion_ && driver.hasChanged(HarmonicProblemChanges::OMEGA));

    if (fullRecompute ||
        driver.hasChanged(HarmonicProblemChanges::CURRENT_SOURCE)) {
        std::map<std::size_t,
                 typename AbstractDiscretization<DIM>::InputFunction>
            elementVectors;
        elementVectors[AbstractDiscretizationBase::ELEMENT_VECTOR_ID] =
            [&problem](const Base::Element& element,
                       const Geometry::PointPhysical<DIM>& p) {
                return problem.sourceTerm(element, p);
            };
        discretization_->computeElementIntegrals(
            *mesh_, elementVectors,
            fullRecompute
                ? AbstractDiscretizationBase::LocalIntegrals::ALL
                : AbstractDiscretizationBase::LocalIntegrals::ONLY_VECTORS);
    }
    if (fullRecompute ||
        driver.hasChanged(HarmonicProblemChanges::BOUNDARY_CONDITION_VALUE)) {
        std::map<std::size_t,
                 typename AbstractDiscretization<DIM>::FaceInputFunction>
            faceVectors;
        faceVectors[AbstractDiscretizationBase::FACE_VECTOR_ID] =
            [&problem](Base::PhysicalFace<DIM>& pface) {
                return problem.boundaryCondition(pface);
            };
        discretization_->computeFaceIntegrals(
            *mesh_, faceVectors,
            [&problem](const Base::Face& face) {
                return problem.getBoundaryConditionType(face);
            },
            fullRecompute
                ? AbstractDiscretizationBase::LocalIntegrals::ALL
                : AbstractDiscretizationBase::LocalIntegrals::ONLY_VECTORS);
    }
    if (fullRecompute) {
        DGMaxLogger(INFO, "Assembling global matrices");
        massMatrix_.reinit();
        stiffnessMatrix_.reinit();
        stiffnessImpedanceMatrix_.reinit();

        PetscErrorCode error;
        if (solverMatrix_ == nullptr) {
            error =
                MatDuplicate(stiffnessMatrix_, MAT_COPY_VALUES, &solverMatrix_);
        } else {
            error = MatCopy(stiffnessMatrix_, solverMatrix_,
                            DIFFERENT_NONZERO_PATTERN);
        }
        CHKERRABORT(PETSC_COMM_WORLD, error);
    }
    DGMaxLogger(INFO, "Assembling load vector");
    loadVector_.reinit();
}

template <std::size_t DIM>
void HarmonicSolver<DIM>::Workspace::assembleSolverMatrix(double waveNumber) {
    PetscErrorCode error;
    error = MatCopy(stiffnessMatrix_, solverMatrix_, DIFFERENT_NONZERO_PATTERN);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    // Usually this is a subset, but sometimes the rounding is slightly
    // different and is not a subset. So using SUBSET_NONZERO would crash.

    MatStructure structure;
#if PETSC_VERSION_GE(3, 15, 0)
    structure = UNKNOWN_NONZERO_PATTERN;
#else
    structure = DIFFERENT_NONZERO_PATTERN;
#endif

    error = MatAXPY(solverMatrix_, waveNumber, stiffnessImpedanceMatrix_,
                    structure);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = MatAXPY(solverMatrix_, -waveNumber * waveNumber, massMatrix_,
                    structure);
    CHKERRABORT(PETSC_COMM_WORLD, error);
}

template <std::size_t DIM>
void HarmonicSolver<DIM>::Workspace::solve() {
    PetscErrorCode error;
    DGMaxLogger(INFO, "Setting up solver");
    error = KSPSetOperators(solver_, solverMatrix_, solverMatrix_);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = KSPSetUp(solver_);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    DGMaxLogger(INFO, "Solving harmonic problem");
    error = KSPSolve(solver_, loadVector_, resultVector_);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    {
        // Convergence diagnostics
        PetscInt niters;
        KSPGetIterationNumber(solver_, &niters);
        KSPConvergedReason converged;
        KSPGetConvergedReason(solver_, &converged);
        const char* convergedReason;

#if PETSC_VERSION_GE(3, 15, 0)
        KSPGetConvergedReasonString(solver_, &convergedReason);
#else
        convergedReason = KSPConvergedReasons[converged];
#endif

        if (converged > 0) {
            // Successful
            DGMaxLogger(INFO,
                        "Successfully converged in % iterations with reason %",
                        niters, convergedReason);
        } else {
            DGMaxLogger(WARN,
                        "Failed to converge in % iterations with reason %",
                        niters, convergedReason);
        }
    }
    resultVector_.writeTimeIntegrationVector(VECTOR_ID);
}

template <std::size_t DIM>
void HarmonicSolver<DIM>::solve(
    Base::MeshManipulator<DIM>& mesh,
    DGMax::AbstractHarmonicSolverDriver<DIM>& driver) {
    Workspace workspace(*discretization_, mesh);

    int index = -1;

    while (!driver.stop()) {
        index++;
        std::size_t expectedCount = driver.getExpectedNumberOfProblems();
        DGMaxLogger(INFO, "Starting problem %/%", index, expectedCount);
        driver.nextProblem();
        const HarmonicProblem<DIM>& problem = driver.currentProblem();

        if (driver.hasChanged(HarmonicProblemChanges::OMEGA)) {
            setDispersionWavenumber(problem.omega());
        }

        workspace.computeIntegrals(driver);

        if (driver.hasChanged(HarmonicProblemChanges::OMEGA) ||
            driver.hasChanged(
                HarmonicProblemChanges::BOUNDARY_CONDITION_TYPE)) {
            workspace.assembleSolverMatrix(problem.omega());
        }
        workspace.solve();

        Result result(problem, workspace);
        driver.handleResult(result);
    }
}

template class HarmonicSolver<2>;
template class HarmonicSolver<3>;

}  // namespace DGMax