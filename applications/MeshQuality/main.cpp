/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2014, University of Twente
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.

 3. Neither the name of the copyright holder nor the names of its contributors
 may be used to endorse or promote products derived from this software without
 specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <chrono>
#include "petsc.h"

#include "Base/CommandLineOptions.h"
#include "Base/H1ConformingTransformation.h"
#include "Base/MeshManipulator.h"
#include "Integration/ElementIntegral.h"
#include "Utilities/GlobalIndexing.h"
#include "Utilities/GlobalMatrix.h"
#include "Utilities/GlobalVector.h"

#include "Logger.h"

using namespace hpgem;

auto& name = Base::register_argument<std::string>(
    'n', "meshName", "name of the mesh file", true);
auto& polynomialOrder = Base::register_argument<std::size_t>(
    'p', "order", "polynomial order of the solution", true);

template <std::size_t>
void runWithDimension();

int main(int argc, char** argv) {
    Base::parse_options(argc, argv);

    auto startTime = std::chrono::system_clock::now();

    runWithDimension<2>();

    auto endTime = std::chrono::system_clock::now();
    std::chrono::duration<double> duration = endTime - startTime;
    std::cout << "Elapsed time " << duration.count() << std::endl;

    return 0;
}

template <std::size_t DIM>
void runTestsWithDGBasis(Base::MeshManipulator<DIM>& mesh);
template <std::size_t DIM>
void runTestsWithConformingBasis(Base::MeshManipulator<DIM>& mesh);

template <std::size_t DIM>
void runWithDimension() {
    std::size_t numberOfVectors = 1;
    Base::ConfigurationData configData(1, 0);

    Base::MeshManipulator<DIM> mesh(&configData, 1, 1, 1, 1);
    mesh.readMesh(name.getValue());

    runTestsWithDGBasis(mesh);
    runTestsWithConformingBasis(mesh);
}

template <std::size_t DIM>
void runTestsWithDGBasis(Base::MeshManipulator<DIM>& mesh) {
    const std::size_t ELEM_MAT_ID = 0;
    const std::size_t ELEM_VEC_ID = 0;

    // Linear patch test
    mesh.useDefaultDGBasisFunctions(1);
    Integration::ElementIntegral<DIM> integral;
    integral.setTransformation(
        std::shared_ptr<Base::CoordinateTransformation<DIM>>(
            new Base::H1ConformingTransformation<DIM>()));

    std::vector<double> phis;

    for (Base::Element* element : mesh.getElementsList()) {
        LinearAlgebra::MiddleSizeMatrix mat;
        // Compute mass matrix
        mat = integral.integrate(
            element, [&phis](Base::PhysicalElement<DIM>& pelemnt) {
                std::size_t n = pelemnt.getNumberOfBasisFunctions();
                LinearAlgebra::MiddleSizeMatrix mat(n, n);
                phis.resize(n);
                for (std::size_t i = 0; i < n; ++i) {
                    phis[i] = pelemnt.basisFunction(i);
                }
                for (std::size_t i = 0; i < n; ++i) {
                    for (std::size_t j = i; j < n; ++j) {
                        mat(i, j) = phis[i] * phis[j];
                        mat(j, i) = mat(i, j);
                    }
                }
                return 1.0;
            });
        element->setElementMatrix(mat, ELEM_MAT_ID);
        // Compute source vector
        LinearAlgebra::MiddleSizeVector vec;
        vec = integral.integrate(
            element, [](Base::PhysicalElement<DIM>& pelement) {
                std::size_t n = pelement.getNumberOfBasisFunctions();
                LinearAlgebra::MiddleSizeVector vec(n);
                // Position for the patch test
                // x + 2y + 3z etc.
                const Geometry::PointPhysical<DIM>& p =
                    pelement.getPointPhysical();
                double functionValue = 0;
                for (std::size_t i = 0; i < DIM; ++i) {
                    functionValue += (i + 1) * p[i];
                }
                for (std::size_t i = 0; i < n; ++i) {
                    vec[i] = functionValue * pelement.basisFunction(i);
                }
                return vec;
            });
        element->setElementVector(vec, ELEM_VEC_ID);
    }
    // Setup global system
    Utilities::GlobalIndexing indexing(&mesh);
    Utilities::GlobalPetscMatrix massMatrix (indexing, ELEM_MAT_ID, -1);
    Utilities::GlobalPetscVector loadVector (indexing, ELEM_VEC_ID, -1);
    Utilities::GlobalPetscVector resultVector (indexing, -1, -1);
    // Solve system
    KSP ksp;
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetType(ksp, KSPPREONLY);

    KSPSetOperators(ksp, massMatrix, massMatrix);
    KSPSolve(ksp, loadVector, resultVector);

    // Computing the error
    resultVector.writeTimeIntegrationVector(0);

    for (Base::Element* element : mesh.getElementsList()) {

    }
}
