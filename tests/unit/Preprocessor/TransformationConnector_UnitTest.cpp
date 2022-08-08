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
#include <TransformationConnector.h>
#include <MeshSource.h>
#include <MeshFactory.h>

class TestingMeshSource : public Preprocessor::MeshSource2 {
   public:
    TestingMeshSource(std::size_t dim) : dimension_(dim){};

    const std::vector<Coord>& getCoordinates() override { return coords_; }
    const std::vector<Element>& getElements() override { return elements_; }

    size_t getDimension() const override { return dimension_; }

    std::vector<Coord> coords_;
    std::vector<Element> elements_;

   private:
    std::size_t dimension_;
};

using SCoord = Preprocessor::MeshSource2::Coord;
using SElement = Preprocessor::MeshSource2::Element;

TEST_CASE("1D line", "[TransformationConnector_UnitTest]") {
    // Create simple line mesh with 3 equally spaced elements of length 1.0
    TestingMeshSource source(1);
    for (std::size_t i = 0; i < 4; ++i) {
        source.coords_.push_back(SCoord(i, {(double)i}));
    }
    for (std::size_t i = 0; i < 3; ++i) {
        SElement element;
        element.coordinateIds = {i, i + 1};
        element.zoneName = "Dummy";
        source.elements_.push_back(std::move(element));
    }

    Preprocessor::Mesh<1> mesh = Preprocessor::fromMeshSource<1>(source);
    {
        INFO("Preconditions");
        REQUIRE(mesh.getNumberOfElements() == 3);
        REQUIRE(mesh.getNumberOfNodes() == 4);
    }

    // Actual merging of the outer nodes
    Preprocessor::translationConnect(mesh, {3.0});
    mesh.removeUnusedEntities();

    // Check the result
    INFO("After merge");
    REQUIRE(mesh.getNumberOfElements() == 3);
    REQUIRE(mesh.getNumberOfNodes() == 3);
    REQUIRE(mesh.getNodeCoordinates().size() == 4);
    // Check that the right coordinates have been merged
    Preprocessor::EntityGId left(-1), right(-1);
    for (auto& coord : mesh.getNodeCoordinates()) {
        if (coord.coordinate[0] == 0.0) {
            left = coord.nodeIndex;
        } else if (coord.coordinate[0] == 3.0) {
            right = coord.nodeIndex;
        }
    }
    REQUIRE(left == right);
}

TEST_CASE("2D square", "[TransformationConnector_UnitTest]") {
    // Test case: A square of 2x2 squares
    const std::size_t NX = 2, NY = 2;
    TestingMeshSource source(2);
    for (std::size_t xi = 0; xi < NX + 1; ++xi) {
        for (std::size_t yi = 0; yi < NY + 1; ++yi) {
            source.coords_.push_back(
                SCoord((NX + 1) * yi + xi, {(double)xi, (double)yi}));
        }
    }
    for (std::size_t xi = 0; xi < NX; ++xi) {
        for (std::size_t yi = 0; yi < NY; ++yi) {
            std::size_t bottomLeft = xi + (NX + 1) * yi;
            SElement element;
            element.zoneName = "Dummy";
            // Four corners in hpgem ordering
            element.coordinateIds = {bottomLeft, bottomLeft + 1,
                                     bottomLeft + (NX + 1),
                                     bottomLeft + (NX + 1) + 1};
            source.elements_.push_back(std::move(element));
        }
    }
    Preprocessor::Mesh<2> mesh = Preprocessor::fromMeshSource<2>(source);
    {
        INFO("Preconditions");
        REQUIRE(mesh.getNumberOfElements() == (NX * NY));
        REQUIRE(mesh.getNumberOfEdges() == (NX * (NY + 1) + (NX + 1) * NY));
        REQUIRE(mesh.getNumberOfNodes() == (NX + 1) * (NY + 1));
    }

    // Merge edges at x=0 and x=NX
    Preprocessor::translationConnect(mesh, {NX, 0.0});
    mesh.removeUnusedEntities();
    {
        INFO("After 1 merge");
        REQUIRE(mesh.getNumberOfElements() == (NX * NY));
        REQUIRE(mesh.getNumberOfEdges() == (NX * (NY + 1) + NX * NY));
        REQUIRE(mesh.getNumberOfNodes() == NX * (NY + 1));
    }
    // Also merge y=0 and y = NY
    Preprocessor::translationConnect(mesh, {0.0, NY});
    mesh.removeUnusedEntities();
    {
        INFO("After 2 merges");
        REQUIRE(mesh.getNumberOfElements() == (NX * NY));
        REQUIRE(mesh.getNumberOfEdges() == (NX * NY + NX * NY));
        REQUIRE(mesh.getNumberOfNodes() == NX * NY);
    }
}