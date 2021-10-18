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

#include "DivDGMaxHarmonic.h"

#include "DGMaxLogger.h"
#include <Output/TecplotDiscontinuousSolutionWriter.h>
#include "Utilities/GlobalMatrix.h"
#include "Utilities/GlobalVector.h"

#include <iostream>
#include <petscksp.h>

using namespace hpgem;

template <std::size_t DIM>
DivDGMaxHarmonic<DIM>::DivDGMaxHarmonic(Base::MeshManipulator<DIM>& mesh,
                                        DivDGMaxDiscretizationBase::Stab stab,
                                        std::size_t order)
    : mesh_(mesh), stab_(stab) {
    discretization_.initializeBasisFunctions(mesh_, order);
}

template <std::size_t DIM>
void DivDGMaxHarmonic<DIM>::solve(const HarmonicProblem<DIM>& input) {
    PetscErrorCode error;
    using Discretization = DivDGMaxDiscretization<DIM>;

    discretization_.setBoundaryIndicator(
        std::bind(&HarmonicProblem<DIM>::getBoundaryConditionType, &input,
                  std::placeholders::_1));
    {
        std::map<std::size_t, typename Discretization::InputFunction>
            elementVecs = {{Discretization::ELEMENT_SOURCE_VECTOR_ID,
                            std::bind(&HarmonicProblem<DIM>::sourceTerm,
                                      std::ref(input), std::placeholders::_1)}};
        discretization_.computeElementIntegrands(mesh_, elementVecs);
    }
    {
        std::map<std::size_t, typename Discretization::FaceInputFunction>
            faceVecs = {{Discretization::FACE_BOUNDARY_VECTOR_ID,
                         std::bind(&HarmonicProblem<DIM>::boundaryCondition,
                                   std::ref(input), std::placeholders::_1)}};
        discretization_.computeFaceIntegrals(mesh_, faceVecs, stab_);
    }

    Utilities::GlobalIndexing indexing(&mesh_);
    Utilities::GlobalPetscMatrix massMatrix(
        indexing, DivDGMaxDiscretizationBase::ELEMENT_MASS_MATRIX_ID, -1),
        stiffnessMatrix(indexing,
                        DivDGMaxDiscretizationBase::ELEMENT_STIFFNESS_MATRIX_ID,
                        DivDGMaxDiscretizationBase::FACE_STIFFNESS_MATRIX_ID),
        impedanceMatrix(
            indexing, -1,
            DivDGMaxDiscretizationBase::FACE_STIFFNESS_IMPEDANCE_MATRIX_ID);
    Utilities::GlobalPetscVector rhs(
        indexing, DivDGMaxDiscretizationBase::ELEMENT_SOURCE_VECTOR_ID,
        DivDGMaxDiscretizationBase::FACE_BOUNDARY_VECTOR_ID),
        result(indexing, -1, -1);

    rhs.assemble();

    std::complex<double> omega2 = input.omega();
    omega2 *= omega2;
    error = MatAXPY(stiffnessMatrix, -1.0 * omega2, massMatrix,
                    SUBSET_NONZERO_PATTERN);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = MatAXPY(stiffnessMatrix, input.omega(), impedanceMatrix,
                    SUBSET_NONZERO_PATTERN);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    KSP solver;
    error = KSPCreate(PETSC_COMM_WORLD, &solver);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    error = KSPSetType(solver, "minres");
    CHKERRABORT(PETSC_COMM_WORLD, error);

    error = KSPSetTolerances(solver, 1e-10, PETSC_DEFAULT, PETSC_DEFAULT,
                             PETSC_DEFAULT);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = KSPSetOperators(solver, stiffnessMatrix, stiffnessMatrix);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    error = KSPSetFromOptions(solver);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    error = KSPSetUp(solver);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = KSPSolve(solver, rhs, result);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    {
        // Convergence diagnostics
        PetscInt niters;
        KSPGetIterationNumber(solver, &niters);
        KSPConvergedReason converged;
        KSPGetConvergedReason(solver, &converged);
        const char* convergedReason;

#if PETSC_VERSION_GE(3, 15, 0)
        KSPGetConvergedReasonString(solver, &convergedReason);
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

    // TODO: This is a bit dirty to store it in the time integration vector.
    result.writeTimeIntegrationVector(0);

    // Cleanup
    error = KSPDestroy(&solver);
    CHKERRABORT(PETSC_COMM_WORLD, error);
}

template <std::size_t DIM>
void DivDGMaxHarmonic<DIM>::writeVTK(
    Output::VTKSpecificTimeWriter<DIM>& output) const {

    using Fields = typename DivDGMaxDiscretization<DIM>::Fields;

    // 4 fields to output, Ereal, Eimag, preal, pcomplex
    output.write(
        [this](Base::Element* element,
               const Geometry::PointReference<DIM>& point, std::size_t) {
            LinearAlgebra::MiddleSizeVector coefficients =
                element->getTimeIntegrationVector(0);
            Fields fields =
                discretization_.computeFields(element, point, coefficients);

            return fields.electricField.real();
        },
        "Ereal");
    output.write(
        [this](Base::Element* element,
               const Geometry::PointReference<DIM>& point, std::size_t) {
            LinearAlgebra::MiddleSizeVector coefficients =
                element->getTimeIntegrationVector(0);
            Fields fields =
                discretization_.computeFields(element, point, coefficients);

            return fields.electricField.imag();
        },
        "Eimag");
    output.write([this](Base::Element* element,
                        const Geometry::PointReference<DIM>& point,
                        std::size_t) {
        LinearAlgebra::MiddleSizeVector coefficients =
            element->getTimeIntegrationVector(0);
        Fields fields =
            discretization_.computeFields(element, point, coefficients);
        return fields.electricField.l2Norm();
    }, "Emag");

    output.write(
        [this](Base::Element* element,
               const Geometry::PointReference<DIM>& point, std::size_t) {
            LinearAlgebra::MiddleSizeVector coefficients =
                element->getTimeIntegrationVector(0);
            Fields fields =
                discretization_.computeFields(element, point, coefficients);

            return std::real(fields.potential);
        },
        "preal");

    output.write(
        [this](Base::Element* element,
               const Geometry::PointReference<DIM>& point, std::size_t) {
            LinearAlgebra::MiddleSizeVector coefficients =
                element->getTimeIntegrationVector(0);
            Fields fields =
                discretization_.computeFields(element, point, coefficients);

            return std::imag(fields.potential);
        },
        "pimag");
}

template <std::size_t DIM>
void DivDGMaxHarmonic<DIM>::writeTec(std::string fileName) const {
    std::ofstream fileWriter;
    fileWriter.open(fileName);
    Output::TecplotDiscontinuousSolutionWriter<DIM> tecplotWriter(
        fileWriter, "Electric field", "012", "E0,E1,E2");

    std::string zoneName("field");  // Dummy
    tecplotWriter.write(
        &mesh_, zoneName, false,
        [&](const Base::Element* element,
            const Geometry::PointReference<DIM>& point, std::ostream& stream) {
            const LinearAlgebra::MiddleSizeVector coefficients =
                element->getTimeIntegrationVector(0);
            LinearAlgebra::SmallVector<DIM> electricField =
                discretization_.computeField(element, point, coefficients);
            stream << electricField[0] << " " << electricField[1] << " "
                   << electricField[2] << std::endl;
        });

    fileWriter.close();
}

template <std::size_t DIM>
double DivDGMaxHarmonic<DIM>::computeL2Error(
    const typename DivDGMaxDiscretization<DIM>::InputFunction& exactSolution) {
    return discretization_.computeL2Error(mesh_, 0, exactSolution);
}

template <std::size_t DIM>
double DivDGMaxHarmonic<DIM>::computeL2Error(
    const ExactHarmonicProblem<DIM>& problem) {
    return computeL2Error(std::bind(&ExactHarmonicProblem<DIM>::exactSolution,
                                    std::ref(problem), std::placeholders::_1));
}

template class DivDGMaxHarmonic<2>;
template class DivDGMaxHarmonic<3>;
