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

#include "MergePlan.h"

#include "../catch.hpp"

using namespace Preprocessor;

/**
 * Check the number of entities that are removed from the mesh after a merge.
 *
 * This is primarily a check on the execution of the merge.
 *
 * @tparam dim The dimension of the mesh
 * @param mesh The mesh after merging but before removing unused entities
 * @param plan The plan that was used in the merging
 */
template <std::size_t dim>
void checkRemovedEntityCount(Mesh<dim>& mesh, MergePlan<dim>& plan) {
    std::array<std::size_t, dim> beforeCounts {};
    for (std::size_t i = 0; i < dim; ++i) {
        beforeCounts[i] = mesh.getNumberOfEntities(i);
    }
    mesh.removeUnusedEntities();
    INFO("Checking removed entity count");
    const auto& merges = plan.getMerges();
    for (std::size_t i = 0; i < dim; ++i) {
        // Compute the number of entities that should be merged away
        std::size_t removals = 0;
        for (const auto& merge : merges[i]) {
            // If N entities are merged, N-1 are removed
            removals += merge.size() - 1;
        }

        CHECK(beforeCounts[i] - removals == mesh.getNumberOfEntities(i));
    }
}

TEST_CASE("1D line", "[MergePlan]") {
    // Test case: Create a line
    // 0 -- 1 -- 2
    // and merge both ends

    Mesh<1> mesh;

    for (int i = 0; i < 3; ++i) {
        EntityGId nodeId = mesh.addNode();
        mesh.addNodeCoordinate(nodeId,
                               LinearAlgebra::SmallVector<1>({1.0 * i}));
    }
    for (std::size_t i = 0; i < 2; ++i) {
        mesh.addElement({CoordId(i), CoordId(i + 1)});
    }
    // Setup the merge plan
    std::map<CoordId, CoordId> coordinatePairing;
    coordinatePairing[CoordId(0)] = CoordId(2);

    MergePlan<1> plan =
        MergePlan<1>::computeMergePlan(&mesh, coordinatePairing);
    auto& merges = plan.getMerges();
    INFO("Expect one merge")
    REQUIRE(merges[0].size() == 1);
    INFO("Expect nodes 0 and 2 to be merged");
    std::vector<EntityGId> merge = merges[0][0];
    std::sort(merge.begin(), merge.end());
    std::vector<EntityGId> expectedMerge({EntityGId(0), EntityGId(2)});
    REQUIRE(merge == expectedMerge);

    // Execute the merge
    plan.executeMerge();
    INFO("Coordinates use the same node");
    REQUIRE(mesh.getNodeCoordinates()[0].nodeIndex ==
            mesh.getNodeCoordinates()[2].nodeIndex);
    EntityGId targetNode = mesh.getNodeCoordinates()[0].nodeIndex;
    EntityGId unusedNode = EntityGId(2 - targetNode.id);
    INFO("Only one node is used");
    CHECK(mesh.getNode(targetNode).getNumberOfElements() == 2);
    CHECK(mesh.getNode(unusedNode).getNumberOfElements() == 0);
    INFO("Coordinates reference are unchanged");
    // Abusing that we know the coordinate indices from the original element
    // construction.
    CHECK(mesh.getElement(EntityGId(0)).getCoordinateIndex(0) == CoordId(0));
    CHECK(mesh.getElement(EntityGId(1)).getCoordinateIndex(1) == CoordId(2));

    checkRemovedEntityCount(mesh, plan);
}

TEST_CASE("2D squares", "[MergePlan]") {
    // Test case: Create 2 squares:
    // 0 - 2 - 4
    // |   |   |
    // 1 - 3 - 5
    //
    // With elements 0 (left) and 1 (right)

    Mesh<2> mesh;
    for (int i = 0; i < 3; ++i) {
        EntityGId topNode = mesh.addNode();
        mesh.addNodeCoordinate(topNode,
                               LinearAlgebra::SmallVector<2>({1.0 * i, 1.0}));
        EntityGId bottomNode = mesh.addNode();
        mesh.addNodeCoordinate(bottomNode,
                               LinearAlgebra::SmallVector<2>({1.0 * i, 0.0}));
    }
    for (std::size_t i = 0; i < 2; ++i) {
        mesh.addElement({CoordId(1 + 2 * i), CoordId(3 + 2 * i),
                         CoordId(0 + 2 * i), CoordId(2 + 2 * i)});
    }
    // Merge plan:
    // Merge coordinates 0 - 4, and 1 - 5
    // Expect that the edge gets merged too
    std::map<CoordId, CoordId> coordMapping;
    coordMapping[CoordId(0)] = CoordId(4);
    coordMapping[CoordId(1)] = CoordId(5);

    MergePlan<2> plan = MergePlan<2>::computeMergePlan(&mesh, coordMapping);

    auto& merges = plan.getMerges();
    INFO("Node merges")
    REQUIRE(merges[0].size() == 2);

    INFO("Edge/face merges")
    REQUIRE(merges[1].size() == 1);
    std::vector<EntityGId> merge = merges[1][0];
    std::sort(merge.begin(), merge.end());
    std::vector<EntityGId> expectedMerge;
    // To get the global Ids of the edges that we expect to be merged we use
    // that the left edge of a square is number 1, while the right edge is
    // number 2. Additionally, we use that the second element is constructed
    // after the first one and thus has the higher edge number.
    expectedMerge.push_back(mesh.getElement(EntityGId(0))
                                .template getIncidenceListAsIndices<1>()[1]);
    expectedMerge.push_back(mesh.getElement(EntityGId(1))
                                .template getIncidenceListAsIndices<1>()[2]);
    REQUIRE(merge == expectedMerge);

    plan.executeMerge();

    checkRemovedEntityCount(mesh, plan);
}

TEST_CASE("2D triangles", "[MergePlan]") {
    // Check the following case:
    //
    // 0 -- 1 -- 2
    //   \  |  /
    //    \ | /
    //      3
    // Where the following nodes with merges:
    // 0 -> 2
    // 1 -> 1
    // Where we expect the top edges to merge, but not the diagonals
    Mesh<2> mesh;
    for (int i = 0; i < 3; ++i) {
        EntityGId nodeId = mesh.addNode();
        mesh.addNodeCoordinate(nodeId,
                               LinearAlgebra::SmallVector<2>({1.0 * i, 0.0}));
    }
    {
        EntityGId nodeId = mesh.addNode();
        mesh.addNodeCoordinate(nodeId,
                               LinearAlgebra::SmallVector<2>({1.0, -1.0}));
    }
    mesh.addElement({CoordId(0), CoordId(1), CoordId(3)});
    mesh.addElement({CoordId(1), CoordId(2), CoordId(3)});

    std::map<CoordId, CoordId> coordMapping;
    coordMapping[CoordId(0)] = CoordId(2);
    coordMapping[CoordId(1)] = CoordId(1);

    MergePlan<2> plan = MergePlan<2>::computeMergePlan(&mesh, coordMapping);

    // Expected sizes
    const auto& merges = plan.getMerges();
    INFO("Merge counts")
    CHECK(merges[0].size() == 1);    // Only nodes 0 & 2
    REQUIRE(merges[1].size() == 1);  // The two top edges being merged
    INFO("Merged edges");
    // Using knowledge about the local edge numbers (from the reference
    // geometry) compute the expected edges to be merged
    std::vector<EntityGId> mergedEdges = merges[1][0];
    std::sort(mergedEdges.begin(), mergedEdges.end());
    std::vector<EntityGId> expectedEdges;
    // Note these should be sorted by construction order
    expectedEdges.push_back(mesh.getElement(EntityGId(0))
                                .template getIncidenceListAsIndices<1>()[0]);
    expectedEdges.push_back(mesh.getElement(EntityGId(1))
                                .template getIncidenceListAsIndices<1>()[0]);
    REQUIRE(mergedEdges == expectedEdges);

    plan.executeMerge();
    checkRemovedEntityCount(mesh, plan);
}