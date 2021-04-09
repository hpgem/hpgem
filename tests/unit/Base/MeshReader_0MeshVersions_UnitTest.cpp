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

/**
 * These unit tests are design to test the mesh reader with meshes with
 * different mesh-versions. This is done to ensure:
 *  - That it can read the intermediate meshes from different versions
 *  - That it correctly reads zones (which depends on the version)
 */

#include "Base/MeshManipulator.h"
#include "hpgem-cmake.h"
#include "../catch.hpp"

using namespace hpgem;

Base::ConfigurationData dummyConfig(1);

template <std::size_t DIM>
using MeshPtr = std::unique_ptr<Base::MeshManipulator<DIM>>;

template <std::size_t DIM>
MeshPtr<DIM> readMesh(const std::string& filename) {

    // Initialize MPI if needed
    int initialized;
    MPI_Initialized(&initialized);
    if (!initialized) {
        MPI_Init(nullptr, nullptr);
    }

    // Actual mesh reading
    auto mesh = std::make_unique<Base::MeshManipulator<DIM>>(&dummyConfig);
    mesh->readMesh(getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/fixedMeshes/" +
                   filename);
    return mesh;
}

// Checks that there is only 1 zone with the expected name, and that all
// elements belong to it.
template <std::size_t DIM>
void testSingleZone(MeshPtr<DIM>& mesh, std::string expectedZoneName) {
    const std::vector<Base::Zone>& zones = mesh->getZones();
    INFO("Exactly 1 zone")
    REQUIRE(zones.size() == 1);
    INFO("Zone has id 0")
    REQUIRE(zones[0].getZoneId() == 0);
    INFO("Correctly named zone");
    REQUIRE(zones[0].getName() == expectedZoneName);

    INFO("Elements have a zone")
    for (const Base::Element* elem : mesh->getElementsList()) {
        REQUIRE(elem->getZone() == 0);
    }
}

// Default zone name when reading the Version 1 mesh format
const std::string DEFAULT_V1_ZONENAME = "Main";

TEST_CASE("mesh format 1: 1D", "[Mesh reader - fixed meshes]") {
    auto mesh = readMesh<1>("meshD1N2v1.hpgem");
    testSingleZone<1>(mesh, DEFAULT_V1_ZONENAME);
}

TEST_CASE("mesh format 2: 1D single zone", "[Mesh reader - fixed meshes]") {
    auto mesh = readMesh<1>("meshD1N2v2-singlezone.hpgem");
    testSingleZone<1>(mesh, "TestZone");
}

TEST_CASE("mesh format 1: 3D", "[Mesh reader - fixed meshes]") {
    auto mesh = readMesh<3>("meshD3N2v1.hpgem");
    testSingleZone<3>(mesh, DEFAULT_V1_ZONENAME);
}

TEST_CASE("mesh format 2: 3D single zone", "[Mesh reader - fixed meshes]") {
    auto mesh = readMesh<3>("meshD3N2v2.hpgem");
    testSingleZone<3>(mesh, "TestZone");
}

TEST_CASE("mesh format 2: 1D two zones", "[Mesh reader - fixed meshes]") {
    auto mesh = readMesh<1>("meshD1N2v2-twozones.hpgem");

    INFO("Expect 2 elements");
    REQUIRE(mesh->getNumberOfElements() == 2);

    const std::size_t numZones = 2;
    const std::vector<Base::Zone>& zones = mesh->getZones();
    INFO("Expect 2 zones");
    REQUIRE(zones.size() == numZones);

    // Element id's are shared between all elements, not just those in a single
    // mesh. The zones are named 'Element0' and 'Element1'.
    std::size_t effectiveElementId = 0;
    for (const Base::Element* element : mesh->getElementsList()) {
        std::stringstream expectedZoneName;
        // From naming the zones manually
        expectedZoneName << "Element" << effectiveElementId;
        effectiveElementId++;

        INFO("Valid zone pointer")
        std::size_t zoneId = element->getZone();
        // Using pointer subtraction to get an idea of what actually might go
        // wrong. As raw addresses are hardly readable
        REQUIRE(zoneId >= 0);
        REQUIRE(zoneId < numZones);

        std::string expectedZoneNameStr = expectedZoneName.str();
        INFO("Zone matches ID for element " << element->getID());
        REQUIRE(zones[zoneId].getName() == expectedZoneNameStr);
    }
}