/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2020, University of Twente
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

#include "gmsh.h"

#include "../catch.hpp"

TEST_CASE("ReadingFile", "[ReadingFile]") {

    Preprocessor::GmshReader reader(std::string(TEST_DATA_FOLDER) +
                                    "/gmsh2.2.gmsh");
    INFO("Dimension");
    REQUIRE(reader.getDimension() == 2);

    std::vector<Preprocessor::MeshSource2::Coord> ref_Coords;
    ref_Coords.reserve(6);
    ref_Coords.push_back({0, {0.0, 0.0}});
    ref_Coords.push_back({1, {1.0, 0.0}});
    ref_Coords.push_back({2, {1.0, 1.0}});
    ref_Coords.push_back({3, {0.0, 1.0}});
    ref_Coords.push_back({4, {2.0, 0.0}});
    ref_Coords.push_back({5, {2.0, 1.0}});

    const auto& coords = reader.getCoordinates();
    INFO("Nodes");
    REQUIRE(coords.size() == ref_Coords.size());
    for (size_t i = 0; i < coords.size(); i++) {
        INFO("Coordinate:" + std::to_string(i));
        REQUIRE(coords[i].nodeId == ref_Coords[i].nodeId);
        for (size_t dim = 0; dim < 2; dim++) {
            INFO("Dim:" + std::to_string(dim));
            REQUIRE(coords[i].coordinate[dim] ==
                    Approx(ref_Coords[i].coordinate[dim]));
        }
    }

    // Check periodic connections (Note: Fortran offset)
    // 5 - 1
    // 6 - 4
    std::vector<std::pair<std::size_t, std::size_t>> mergedNodes = {
        {4, 0},
        {5, 3}
    };
    const auto& actualMerges = reader.getMerges();
    CHECK(actualMerges.size() == 1);
    for(const auto& p : mergedNodes) {
        const auto pf = actualMerges[0].find(p.first);
        REQUIRE(pf != actualMerges[0].end());
        REQUIRE(pf->second == p.second);
    }

    std::vector<Preprocessor::MeshSource2::Element> ref_Elements;
    ref_Elements.reserve(2);
    ref_Elements.push_back({{0, 1, 3, 2}, "2", 2, 0});
    ref_Elements.push_back({{1, 4, 2, 5}, "1", 2, 1});

    const auto& elements = reader.getElements();

    INFO("Elements");
    REQUIRE(elements.size() == ref_Elements.size());

    for (size_t i = 0; i < elements.size(); i++) {
        INFO("Element:" + std::to_string(i));
        INFO("id");
        REQUIRE(elements[i].id == ref_Elements[i].id);
        INFO("dimension");
        REQUIRE(elements[i].dimension == ref_Elements[i].dimension);
        INFO("zone name");
        REQUIRE(elements[i].zoneName == ref_Elements[i].zoneName);
        INFO("coordinates");
        REQUIRE(elements[i].coordinateIds == ref_Elements[i].coordinateIds);
    }
}

TEST_CASE("ReadingFile_NOPBC", "[ReadingFile]") {

    Preprocessor::GmshReader reader(std::string(TEST_DATA_FOLDER) +
                                    "/gmsh2.2_nopbc.gmsh");
    INFO("Dimension");
    REQUIRE(reader.getDimension() == 2);

    std::vector<Preprocessor::MeshSource2::Coord> ref_Coords;
    ref_Coords.reserve(6);
    ref_Coords.push_back({0, {0.0, 0.0}});
    ref_Coords.push_back({1, {1.0, 0.0}});
    ref_Coords.push_back({2, {1.0, 1.0}});
    ref_Coords.push_back({3, {0.0, 1.0}});
    ref_Coords.push_back({4, {2.0, 0.0}});
    ref_Coords.push_back({5, {2.0, 1.0}});

    const auto& coords = reader.getCoordinates();
    INFO("Nodes");
    REQUIRE(coords.size() == ref_Coords.size());
    for (size_t i = 0; i < coords.size(); i++) {
        INFO("Coordinate:" + std::to_string(i));
        REQUIRE(coords[i].nodeId == ref_Coords[i].nodeId);
        for (size_t dim = 0; dim < 2; dim++) {
            INFO("Dim:" + std::to_string(dim));
            REQUIRE(coords[i].coordinate[dim] ==
                    Approx(ref_Coords[i].coordinate[dim]));
        }
    }

    REQUIRE(reader.getMerges().empty());

    std::vector<Preprocessor::MeshSource2::Element> ref_Elements;
    ref_Elements.reserve(2);
    ref_Elements.push_back({{0, 1, 3, 2}, "2", 2, 0});
    ref_Elements.push_back({{1, 4, 2, 5}, "1", 2, 1});

    const auto& elements = reader.getElements();

    INFO("Elements");
    REQUIRE(elements.size() == ref_Elements.size());

    for (size_t i = 0; i < elements.size(); i++) {
        INFO("Element:" + std::to_string(i));
        INFO("id");
        REQUIRE(elements[i].id == ref_Elements[i].id);
        INFO("dimension");
        REQUIRE(elements[i].dimension == ref_Elements[i].dimension);
        INFO("zone name");
        REQUIRE(elements[i].zoneName == ref_Elements[i].zoneName);
        INFO("coordinates");
        REQUIRE(elements[i].coordinateIds == ref_Elements[i].coordinateIds);
    }
}

TEST_CASE("Line mesh with physical names", "[GmshReader]") {
    Preprocessor::GmshReader reader(std::string(TEST_DATA_FOLDER) +
                                    "/gmsh2.physical_line_mesh.msh");
    {
        INFO("Dimension")
        REQUIRE(reader.getDimension() == 1);
    }
    // Coordinates at positions x = {0, 1, 2}
    const auto& coords = reader.getCoordinates();
    REQUIRE(coords.size() == 3);
    for (std::size_t i = 0; i < 3; ++i) {
        INFO("Coordinate " + std::to_string(i));
        REQUIRE(coords[i].coordinate.size() == 1);
        REQUIRE(coords[i].coordinate[0] == Approx(i));
    }
    std::vector<std::string> zones = {"line", "line2"};
    const auto& elements = reader.getElements();
    REQUIRE(elements.size() == 2);
    for (std::size_t i = 0; i < 2; ++i) {
        CHECK(elements[i].zoneName == zones[i]);
    }
}
