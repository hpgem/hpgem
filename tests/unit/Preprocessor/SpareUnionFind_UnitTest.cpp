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

#include "SparseUnionFind.h"

#include "../catch.hpp"

TEST_CASE("Empty Set", "[Sparse Union Find]") {
    Preprocessor::SparseUnionFind empty;
    INFO("Empty iterator")
    REQUIRE(empty.begin() == empty.end());

    auto i = GENERATE(0, 2, 3, 10, 1000000);
    REQUIRE(empty.findSet(i) == i);
}

TEST_CASE("Single set", "[Sparse Union Find]") {
    Preprocessor::SparseUnionFind set;
    set.unionSets(0, 2);
    INFO("Unioned are in the same set")
    REQUIRE(set.inSameSet(0, 2));

    INFO("Distinct from other entries")
    CHECK_FALSE(set.inSameSet(0, 1));
    REQUIRE_FALSE(set.inSameSet(1, 2));

    // Add third element using union
    set.unionSets(0, 1);
    INFO("After merge in same set")
    CHECK(set.inSameSet(1, 2));
    REQUIRE(set.inSameSet(0, 1));
}

TEST_CASE("Two sets", "[Sparse Union Find]") {
    Preprocessor::SparseUnionFind set;
    // Setup two sets {3,4} and {10,11}
    set.unionSets(3, 4);
    set.unionSets(10, 11);

    INFO("Check initial connection")
    CHECK(set.inSameSet(3, 4));
    CHECK(set.inSameSet(10, 11));
    CHECK_FALSE(set.inSameSet(3, 10));
    REQUIRE_FALSE(set.inSameSet(11, 4));

    // Merge the sets
    set.unionSets(10, 3);
    INFO("Check post merge")
    auto i = GENERATE(3, 4, 10, 11);
    auto j = GENERATE(3, 4, 10, 11);
    REQUIRE(set.inSameSet(i, j));
}
