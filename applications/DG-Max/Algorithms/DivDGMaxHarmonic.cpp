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

#include <Output/TecplotDiscontinuousSolutionWriter.h>
#include "Utilities/GlobalMatrix.h"
#include "Utilities/GlobalVector.h"

#include <iostream>
#include <petscksp.h>

using namespace hpgem;


template <std::size_t DIM>
DivDGMaxHarmonic<DIM>::DivDGMaxHarmonic(Base::MeshManipulator<DIM>& mesh)
    : mesh_(mesh) {}

template <std::size_t DIM>
void DivDGMaxHarmonic<DIM>::solve(
    const HarmonicProblem<DIM>& input,
    typename DivDGMaxDiscretization<DIM>::Stab stab, std::size_t order) {
    PetscErrorCode error;

    discretization_.initializeBasisFunctions(mesh_, order);
    discretization_.computeElementIntegrands(
        mesh_, false,
        std::bind(&HarmonicProblem<DIM>::sourceTerm, std::ref(input),
                  std::placeholders::_1, std::placeholders::_2),
        nullptr, nullptr);
    discretization_.computeFaceIntegrals(
        mesh_,
        std::bind(&HarmonicProblem<DIM>::boundaryCondition, std::ref(input),
                  std::placeholders::_1, std::placeholders::_2,
                  std::placeholders::_3),
        stab);

    Utilities::GlobalIndexing indexing(&mesh_);
    Utilities::GlobalPetscMatrix massMatrix(
        indexing, DivDGMaxDiscretization<DIM>::ELEMENT_MASS_MATRIX_ID, -1),
        stiffnessMatrix(
            indexing, DivDGMaxDiscretization<DIM>::ELEMENT_STIFFNESS_MATRIX_ID,
            DivDGMaxDiscretization<DIM>::FACE_STIFFNESS_MATRIX_ID);
    Utilities::GlobalPetscVector rhs(
        indexing, DivDGMaxDiscretization<DIM>::ELEMENT_SOURCE_VECTOR_ID,
        DivDGMaxDiscretization<DIM>::FACE_BOUNDARY_VECTOR_ID),
        result(indexing, -1, -1);

    rhs.assemble();

    std::complex<double> omega2 = input.omega();
    omega2 *= omega2;
    error = MatAXPY(stiffnessMatrix, -1.0 * omega2, massMatrix,
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

    // TODO: This is a bit dirty to store it in the time integration vector.
    result.writeTimeIntegrationVector(0);

    // Cleanup
    error = KSPDestroy(&solver);
    CHKERRABORT(PETSC_COMM_WORLD, error);
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
    const typename DivDGMaxDiscretization<DIM>::InputFunction& exactSolution)
    const {
    return discretization_.computeL2Error(mesh_, 0, exactSolution);
}

template <std::size_t DIM>
double DivDGMaxHarmonic<DIM>::computeL2Error(
    const ExactHarmonicProblem<DIM>& problem) const {
    return computeL2Error(std::bind(&ExactHarmonicProblem<DIM>::exactSolution,
                                    std::ref(problem), std::placeholders::_1,
                                    std::placeholders::_2));
}

template class DivDGMaxHarmonic<2>;
template class DivDGMaxHarmonic<3>;
