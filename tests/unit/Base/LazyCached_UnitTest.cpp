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
#include "Base/LazyCached.h"

using hpgem::Base::LazyCached;
using hpgem::Base::LazyVectorCached;

TEST_CASE("Basic test", "[LazyCached]") {

    // Counter to mirror the internal state
    std::size_t index = 0;
    LazyCached<std::size_t> cache = LazyCached<std::size_t>(
        [&index](std::size_t& value) { index = ++value; });

    REQUIRE(cache.get() == 1);  // Does the first computation
    REQUIRE(index == 1);        // Should been done once
    // No reset, so still the same value
    REQUIRE(cache.get() == 1);

    cache.reset();
    REQUIRE(index == 1);        // Computation is only run when calling get()
    REQUIRE(cache.get() == 2);  // Now the computation should run
    REQUIRE(index == 2);        // Second computation did actually run

    // Check that the value remains the same
    REQUIRE(cache.get() == 2);  // Value did not change
    REQUIRE(cache.get() == 2);  // Value did not change

    // Last reset and check
    cache.reset();
    REQUIRE(cache.get() == 3);
}

TEST_CASE("Basic Vector test", "[LazyVectorCached]") {

    std::array<std::size_t, 2> indices = {0, 0};
    LazyVectorCached<std::size_t> cache = LazyVectorCached<std::size_t>(
        [&indices](std::size_t& value, std::size_t index) {
            REQUIRE(indices[index] == value);
            value++;
            indices[index] = value;
        });

    INFO("Initial resize")
    REQUIRE(cache.size() == 0);
    cache.reset(2);
    REQUIRE(cache.size() == 2);

    // Second entry
    INFO("After first reset")
    REQUIRE(indices[1] == 0);    // No computation
    REQUIRE(cache.get(1) == 1);  // First computation
    CHECK(indices[0] == 0);      // Different entry -> Untouched
    REQUIRE(cache[1] == 1);      // Computation has happened
    REQUIRE(cache.get(1) == 1);  // No reset -> Same value

    INFO("After second reset");
    cache.reset();
    REQUIRE(indices[1] == 1);    // No computation yet
    REQUIRE(cache.get(1) == 2);  // Second computation
    CHECK(indices[0] == 0);      // Different entry -> Untouched
    REQUIRE(cache[1] == 2);      // Computation has happened
    REQUIRE(cache.get(1) == 2);  // No reset -> Same value

    // Check only now compute the first entry -> Should have value 1
    REQUIRE(cache.get(0) == 1);
    CHECK(indices[0] == 1);
    REQUIRE(indices[1] == 2);
    REQUIRE(cache.get(0) == 1);
}

TEST_CASE("Resize reset test", "[LazyVectorCached]") {
    std::array<std::size_t, 3> indices = {0, 0, 0};
    LazyVectorCached<std::size_t> cache = LazyVectorCached<std::size_t>(
        [&indices](std::size_t& value, std::size_t index) {
            REQUIRE(indices[index] == value);
            value++;
            indices[index] = value;
        });

    INFO("Initial sizing");
    cache.reset(2);
    REQUIRE(cache.size() == 2);

    INFO("Compute initial values")
    REQUIRE(cache[0] == 1);
    REQUIRE(cache[1] == 1);

    INFO("Resize-reset 1")
    cache.reset(1);
    REQUIRE(cache.size() == 1);
    REQUIRE(cache[0] == 2);

    INFO("Resize-reset 3")
    cache.reset(3);
    REQUIRE(cache.size() == 3);

    CHECK(cache[0] == 3);
    // As the cache was resized from 1 to 3 elements, the starting values of
    // element 2 & 3 is 0. Hence, we should reset the expected index for the
    // second entry.
    indices[1] = 0;
    CHECK(cache[1] == 1);
    REQUIRE(cache[2] == 1);
}
