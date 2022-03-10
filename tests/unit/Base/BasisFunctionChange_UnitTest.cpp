/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2021, University of Twente
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

#include "../catch.hpp"
#include <cstdlib>

#include "Base/MeshManipulator.h"

// Showing that this is actually more like a self test
// but there we can't use catch yet.
#include "../../self/TestMeshes.h"

#ifdef HPGEM_USE_MPI
#include "mpi.h"
#endif

// Test to demonstrate that changing the basis function will result in correct
// basis function count afterwards. Note that switching between different basis
// functions / order of them may still have residual effects (e.g. the default
// order of the quadrature).

void initMPI() {
#ifdef HPGEM_USE_MPI
    int initialized;
    MPI_Initialized(&initialized);
    if (!initialized) {
        MPI_Init(nullptr, nullptr);
    }
    std::atexit([]() { MPI_Finalize(); });
#endif
}

using namespace hpgem::Base;

TEST_CASE("Reset CG -> DG", "[MeshManipulator]") {
    initMPI();

    // 2 unknowns
    ConfigurationData configData(2);
    MeshManipulator<3> mesh(&configData);
    // We get a small but not the smallest 3D mesh. This ensures we have:
    //  - All three of faces, edges and nodes
    //  - Both boundary and internal faces, etc.
    std::vector<std::string> meshes = hpgem::getUnitCubeCubeMeshes(1, 2);
    mesh.readMesh(meshes[0]);

    // First assign CG, so that we have basis functions on each part
    mesh.useDefaultDGBasisFunctions(0);
    // 4-th order so that we have basis functions on all parts
    mesh.useDefaultConformingBasisFunctions(4, 1);

    // Reset and use DG everywhere with order 0, which means 1 basis function on
    // each element per unknown.
    mesh.useDefaultDGBasisFunctions(0);
    mesh.useDefaultDGBasisFunctions(0, 1);

    // Check these counts
    INFO("Check Element basis function count");
    for (const Element* element : mesh.getElementsList()) {
        CHECK(element->getLocalNumberOfBasisFunctions(0) == 1);
        CHECK(element->getLocalNumberOfBasisFunctions(1) == 1);
        REQUIRE(element->getTotalNumberOfBasisFunctions() == 2);
    }
    INFO("Check face basis function count");
    for (const Face* face : mesh.getFacesList()) {
        REQUIRE(face->getTotalLocalNumberOfBasisFunctions() == 0);
    }

    INFO("Check edge basis function count");
    for (const Edge* edge : mesh.getEdgesList()) {
        REQUIRE(edge->getTotalLocalNumberOfBasisFunctions() == 0);
    }

    INFO("Check node basis function count");
    for (const Node* node : mesh.getNodesList()) {
        REQUIRE(node->getTotalLocalNumberOfBasisFunctions() == 0);
    }
}
