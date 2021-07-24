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

#include <fstream>
#include "hpgem-cmake.h"
#include "Base/MeshFileInformation.h"

using hpgem::Base::MeshFileInformation;

void testSingleMesh(const std::string& meshFile,
                    const MeshFileInformation& reference) {
    MeshFileInformation read = MeshFileInformation::readInformation(
        hpgem::getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/fixedMeshes/" +
        meshFile);

    INFO("Version")
    CHECK(read.version == reference.version);
    INFO("Dimension")
    CHECK(read.dimension == reference.dimension);
    INFO("Entity counts")
    CHECK(read.entityCount == reference.entityCount);
    INFO("Partition counts")
    CHECK(read.partitionNodeCounts == reference.partitionNodeCounts);
    INFO("Zone names")
    CHECK(read.zoneNames == reference.zoneNames);
}

// The reference data for the following test have been read manually from the
// mesh files.

TEST_CASE("Mesh format 1: 1D", "[MeshFileInformation - fixed meshes]") {
    MeshFileInformation reference;
    reference.version = 1;
    reference.dimension = 1;
    reference.entityCount = {4, 3};
    reference.partitionNodeCounts = {4};
    reference.zoneNames = {MeshFileInformation::MESH_V1_ZONENAME};
    testSingleMesh("meshD1N2v1.hpgem", reference);
}

TEST_CASE("mesh format 2: 1D single zone", "[MeshFileInformation]") {
    MeshFileInformation reference;
    reference.version = 2;
    reference.dimension = 1;
    reference.entityCount = {3, 2};
    reference.partitionNodeCounts = {3};
    reference.zoneNames = {"TestZone"};
    testSingleMesh("meshD1N2v2-singlezone.hpgem", reference);
}

TEST_CASE("mesh format 1: 3D", "[MeshFileInformation]") {
    MeshFileInformation reference;
    reference.version = 1;
    reference.dimension = 3;
    reference.entityCount = {27, 54, 36, 8};
    reference.partitionNodeCounts = {27};
    reference.zoneNames = {MeshFileInformation::MESH_V1_ZONENAME};
    testSingleMesh("meshD3N2v1.hpgem", reference);
}

TEST_CASE("mesh format 2: 3D single zone", "[MeshFileInformation]") {
    MeshFileInformation reference;
    reference.version = 2;
    reference.dimension = 3;
    reference.entityCount = {27, 54, 36, 8};
    reference.partitionNodeCounts = {27};
    reference.zoneNames = {"TestZone"};
    testSingleMesh("meshD3N2v2.hpgem", reference);
}

TEST_CASE("mesh format 2: 1D two zones", "[MeshFileInformation]") {
    MeshFileInformation reference;
    reference.version = 2;
    reference.dimension = 1;
    reference.entityCount = {3, 2};
    reference.partitionNodeCounts = {3};
    reference.zoneNames = {"Element0", "Element1"};
    testSingleMesh("meshD1N2v2-twozones.hpgem", reference);
}

// readInformation(std::istream) Post conditions
////////////////////////////////////////////////

TEST_CASE("Read information stream position - v1", "[MeshFileInformation]") {
    std::ifstream stream;
    stream.open(hpgem::getCMAKE_hpGEM_SOURCE_DIR() +
                "/tests/files/fixedMeshes/" + "meshD3N2v1.hpgem");
    MeshFileInformation information;
    information.readInformation(stream);
    // Post condition: The stream is now at the point just past the zones
    // section in the file. The next section is the nodes which starts with
    // node description. The first line is the number of partitions for which it
    // is relevant, followed by their indices. In this case that is 1 0.
    std::string line;
    std::getline(stream, line);

    // Extra trailing space is apparently present in the files.
    CHECK(line == "1 0 ");
}

TEST_CASE("Read information stream position - v2", "[MeshFileInformation]") {
    std::ifstream stream;
    stream.open(hpgem::getCMAKE_hpGEM_SOURCE_DIR() +
                "/tests/files/fixedMeshes/" + "meshD3N2v2.hpgem");
    MeshFileInformation information;
    information.readInformation(stream);
    // Post condition: The stream is now at the point just past the zones
    // section in the file. The next section is the nodes which starts with a
    // section header in the V2 format.
    std::string line;
    std::getline(stream, line);
    CHECK(line == "nodes");
}
