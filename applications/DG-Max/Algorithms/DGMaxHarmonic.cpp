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

#include "DGMaxHarmonic.h"

#include <petscksp.h>

#include <complex>
#include "Output/TecplotDiscontinuousSolutionWriter.h"
#include "Utilities/GlobalMatrix.h"
#include "Utilities/GlobalVector.h"

using namespace hpgem;

namespace DGMax {

template <std::size_t DIM>
DGMaxHarmonic<DIM>::DGMaxHarmonic(Base::MeshManipulator<DIM>& mesh, double stab,
                                  std::size_t order)
    : mesh_(mesh), discretization(order, stab) {
    discretization.initializeBasisFunctions(mesh_);
    discretization.setMatrixHandling(DGMaxDiscretizationBase::NORMAL);
}

template <std::size_t DIM>
class DGMaxHarmonic<DIM>::Result : public AbstractHarmonicResult<DIM> {
   public:
    Result(const HarmonicProblem<DIM>& problem, DGMaxHarmonic<DIM>& solver)
        : problem_(&problem), solver_(&solver){};

    const HarmonicProblem<DIM>& solvedProblem() final { return *problem_; }

    void writeVTK(Output::VTKSpecificTimeWriter<DIM>& output) final {
        solver_->writeVTK(output);
    }
    double computeL2Error(const ExactHarmonicProblem<DIM>& solution) final {
        return solver_->computeL2Error(solution);
    }
    Base::MeshManipulator<DIM>& getMesh() final { return solver_->mesh_; }

   private:
    const HarmonicProblem<DIM>* problem_;
    DGMaxHarmonic<DIM>* solver_;
};

template <std::size_t DIM>
class DGMaxHarmonic<DIM>::Workspace {
   public:
    static constexpr const std::size_t VECTOR_ID = 0;

    Workspace(DGMaxDiscretization<DIM>& discretization,
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

   public:
    void configureSolver();

    DGMaxDiscretization<DIM>* discretization_;
    Base::MeshManipulator<DIM>* mesh_;

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
DGMaxHarmonic<DIM>::Workspace::Workspace(
    DGMaxDiscretization<DIM>& discretization, Base::MeshManipulator<DIM>& mesh)
    : discretization_(&discretization),
      mesh_(&mesh),
      indexing_(nullptr),
      massMatrix_(indexing_, DGMaxDiscretizationBase::MASS_MATRIX_ID, -1),
      stiffnessMatrix_(indexing_, DGMaxDiscretizationBase::STIFFNESS_MATRIX_ID,
                       DGMaxDiscretizationBase::FACE_STIFFNESS_MATRIX_ID),
      stiffnessImpedanceMatrix_(
          indexing_, -1, DGMaxDiscretizationBase::FACE_IMPEDANCE_MATRIX_ID),
      resultVector_(indexing_, -1, -1),
      loadVector_(indexing_, DGMaxDiscretizationBase::ELEMENT_VECTOR_ID,
                  DGMaxDiscretizationBase::FACE_VECTOR_ID),
      solverMatrix_(nullptr),
      solver_(nullptr) {
    indexing_.reset(mesh_, Utilities::GlobalIndexing::BLOCKED_GLOBAL);
    resultVector_.reinit();

    PetscErrorCode err;
    err = KSPCreate(PETSC_COMM_WORLD, &solver_);
    CHKERRABORT(PETSC_COMM_WORLD, err);

    configureSolver();
}

template <std::size_t DIM>
DGMaxHarmonic<DIM>::Workspace::~Workspace() {
    if (solverMatrix_ != nullptr) {
        MatDestroy(&solverMatrix_);
    }
    KSPDestroy(&solver_);
}

template <std::size_t DIM>
void DGMaxHarmonic<DIM>::Workspace::configureSolver() {
    PetscErrorCode error;

    error = KSPSetTolerances(solver_, 1e-8, PETSC_DEFAULT, PETSC_DEFAULT,
                             PETSC_DEFAULT);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    error = KSPSetFromOptions(solver_);
    CHKERRABORT(PETSC_COMM_WORLD, error);
}

template <std::size_t DIM>
void DGMaxHarmonic<DIM>::Workspace::computeIntegrals(
    AbstractHarmonicSolverDriver<DIM>& driver) {

    logger.assert_always(mesh_ != nullptr,
                         "Computing integrals without a mesh");

    const HarmonicProblem<DIM>& problem = driver.currentProblem();

    bool bctChanged =
        driver.hasChanged(HarmonicProblemChanges::BOUNDARY_CONDITION_TYPE);

    if (bctChanged ||
        driver.hasChanged(HarmonicProblemChanges::CURRENT_SOURCE)) {
        std::map<std::size_t, typename DGMaxDiscretization<DIM>::InputFunction>
            elementVectors;
        elementVectors[DGMaxDiscretization<DIM>::ELEMENT_VECTOR_ID] =
            [&problem](const Geometry::PointPhysical<DIM>& p) {
                return problem.sourceTerm(p);
            };
        discretization_->computeElementIntegrals(
            *mesh_, elementVectors,
            bctChanged ? DGMaxDiscretizationBase::LocalIntegrals::ALL
                       : DGMaxDiscretizationBase::LocalIntegrals::ONLY_VECTORS);
    }
    if (bctChanged ||
        driver.hasChanged(HarmonicProblemChanges::BOUNDARY_CONDITION_VALUE)) {
        std::map<std::size_t,
                 typename DGMaxDiscretization<DIM>::FaceInputFunction>
            faceVectors;
        faceVectors[DGMaxDiscretization<DIM>::FACE_VECTOR_ID] =
            [&problem](Base::PhysicalFace<DIM>& pface) {
                return problem.boundaryCondition(pface);
            };
        discretization_->computeFaceIntegrals(
            *mesh_, faceVectors,
            bctChanged ? DGMaxDiscretizationBase::LocalIntegrals::ALL
                       : DGMaxDiscretizationBase::LocalIntegrals::ONLY_VECTORS);
    }
    if (bctChanged) {
        DGMaxLogger(INFO, "Assembling global matrices vector");
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
void DGMaxHarmonic<DIM>::Workspace::assembleSolverMatrix(double waveNumber) {
    PetscErrorCode error;
    error = MatCopy(stiffnessMatrix_, solverMatrix_, DIFFERENT_NONZERO_PATTERN);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    // Usually this is a subset, but sometimes the rounding is slightly
    // different and is not a subset. So using SUBSET_NONZERO would crash.
    error = MatAXPY(solverMatrix_, waveNumber, stiffnessImpedanceMatrix_,
                    UNKNOWN_NONZERO_PATTERN);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = MatAXPY(solverMatrix_, -waveNumber * waveNumber, massMatrix_,
                    SUBSET_NONZERO_PATTERN);
    CHKERRABORT(PETSC_COMM_WORLD, error);
}

template <std::size_t DIM>
void DGMaxHarmonic<DIM>::Workspace::solve() {
    PetscErrorCode error;
    DGMaxLogger(INFO, "Setting up solver");
    error = KSPSetOperators(solver_, solverMatrix_, solverMatrix_);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = KSPSetUp(solver_);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    DGMaxLogger(INFO, "Solving harmonic problem");
    error = KSPSolve(solver_, loadVector_, resultVector_);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    DGMaxLogger(INFO, "Solved harmonic problem");
    resultVector_.writeTimeIntegrationVector(VECTOR_ID);
}

template <std::size_t DIM>
void DGMaxHarmonic<DIM>::solve(
    DGMax::AbstractHarmonicSolverDriver<DIM>& driver) {
    Workspace workspace(discretization, mesh_);

    while (!driver.stop()) {
        driver.nextProblem();
        const HarmonicProblem<DIM>& problem = driver.currentProblem();

        discretization.setBoundaryIndicator([&problem](const Base::Face& face) {
            return problem.getBoundaryConditionType(face);
        });

        workspace.computeIntegrals(driver);

        if (driver.hasChanged(HarmonicProblemChanges::OMEGA) ||
            driver.hasChanged(
                HarmonicProblemChanges::BOUNDARY_CONDITION_TYPE)) {
            workspace.assembleSolverMatrix(problem.omega());
        }
        workspace.solve();

        Result result(problem, *this);
        driver.handleResult(result);
    }
}

template <std::size_t DIM>
std::map<typename DGMaxDiscretization<DIM>::NormType, double>
    DGMaxHarmonic<DIM>::computeError(
        const std::set<typename DGMaxDiscretization<DIM>::NormType>& norms,
        const typename DGMaxDiscretization<DIM>::InputFunction& exactSolution,
        const typename DGMaxDiscretization<DIM>::InputFunction&
            exactSolutionCurl) {
    // Note, this only works by grace of distributing the solution as
    // timeIntegrationVector
    return discretization.computeError(mesh_,
                                       0,  // The abused time vector
                                       exactSolution, exactSolutionCurl, norms);
}

template <std::size_t DIM>
std::map<typename DGMaxDiscretization<DIM>::NormType, double>
    DGMaxHarmonic<DIM>::computeError(
        const std::set<typename DGMaxDiscretization<DIM>::NormType>& norms,
        const ExactHarmonicProblem<DIM>& problem) {
    return computeError(norms,
                        std::bind(&ExactHarmonicProblem<DIM>::exactSolution,
                                  std::ref(problem), std::placeholders::_1),
                        std::bind(&ExactHarmonicProblem<DIM>::exactSolutionCurl,
                                  std::ref(problem), std::placeholders::_1));
}

template <std::size_t DIM>
void DGMaxHarmonic<DIM>::writeTec(std::string fileName) const {
    std::ofstream fileWriter;
    fileWriter.open(fileName);
    Output::TecplotDiscontinuousSolutionWriter<DIM> tecplotWriter(
        fileWriter, "The electric field", "012", "E0,E1,E2,H0,H1,H2");

    logger.log(Log::WARN, "The magnetic field output is incorrectly computed.");

    std::string zoneName("field");  // Dummy
    tecplotWriter.write(
        &mesh_, zoneName, false,
        [&](const Base::Element* element,
            const Geometry::PointReference<DIM>& point, std::ostream& stream) {
            const LinearAlgebra::MiddleSizeVector coefficients =
                element->getTimeIntegrationVector(0);
            LinearAlgebra::SmallVector<DIM> electricField =
                discretization.computeField(element, point, coefficients);
            // TODO: Note that we computeCurlField already converts to real
            // numbers, so we can not compute H = i/(omega mu) curl E
            LinearAlgebra::SmallVector<DIM> curlField =
                discretization.computeCurlField(element, point, coefficients);
            stream << electricField[0] << " " << electricField[1] << " "
                   << electricField[2] << " " << curlField[0] << " "
                   << curlField[1] << " " << curlField[2] << std::endl;
        });

    fileWriter.close();
}

template <std::size_t DIM>
void DGMaxHarmonic<DIM>::writeVTK(
    Output::VTKSpecificTimeWriter<DIM>& output) const {
    discretization.writeFields(output, 0);
}

template class DGMaxHarmonic<2>;
template class DGMaxHarmonic<3>;

}  // namespace DGMax