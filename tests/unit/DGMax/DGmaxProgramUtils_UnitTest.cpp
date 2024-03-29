/*
This file forms part of hpGEM. This package has been developed over a number of
years by various people at the University of Twente and a full list of
contributors can be found at http://hpgem.org/about-the-code/team

This code is distributed using BSD 3-Clause License. A copy of which can found
below.


Copyright (c) 2022, University of Twente
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

#include "../catch.hpp"
#include "DGMaxProgramUtils.h"

using DGMax::stringSplit;

TEST_CASE("Empty input", "[StringSplit]") {
    auto result = stringSplit("", ',');
    REQUIRE(result.empty());
}

TEST_CASE("Single input", "[StringSplit]") {
    auto result = stringSplit("test", ',');
    std::vector<std::string> expected{"test"};
    REQUIRE(result == expected);
}

TEST_CASE("Splitting", "[StringSplit]") {
    auto result = stringSplit("foo,bar,baz", ',');
    std::vector<std::string> expected = {"foo", "bar", "baz"};
    REQUIRE(result == expected);
}

TEST_CASE("Trailing separator", "[StringSplit") {
    auto result = stringSplit("foo,bar,", ',');
    std::vector<std::string> expected = {"foo", "bar", ""};
    REQUIRE(result == expected);
}

TEST_CASE("Separator", "[StringSplit]") {
    auto result = stringSplit("foo:bar:,:", ':');
    std::vector<std::string> expected = {"foo", "bar", ",", ""};
    REQUIRE(result == expected);
}

TEST_CASE("2D test", "[parsePMLZoneDescription]") {
    DGMax::PMLZoneDescription<2> pml =
        DGMax::parsePMLZoneDescription<2>("TEST,+-,1e-1,1e-2");
    CHECK(pml.zoneName_ == "TEST");
    LinearAlgebra::SmallVector<2> expectedDirection{1, -1};
    LinearAlgebra::SmallVector<2> expectedAttenuation{1e-1, 1e-2};
}

TEST_CASE("3D test", "[parsePMLZoneDescription[") {
    DGMax::PMLZoneDescription<3> pml =
        DGMax::parsePMLZoneDescription<3>("Dummy,0+0,1,1e-4,1");
    CHECK(pml.zoneName_ == "Dummy");
    LinearAlgebra::SmallVector<3> expectedDirection{0, 1, 0};
    LinearAlgebra::SmallVector<3> expectedAttenuation{1, 1e-4, 1};
}