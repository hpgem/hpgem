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

#include "utils/TemplateArray.h"
#include "../catch.hpp"

// Dummy template to test with
template <std::size_t>
struct V {
    V() = default;
    explicit V(std::size_t value) : value(value){};
    std::size_t value;
};

TEST_CASE("Constructor & get", "[TemplateArray]") {
    // Purely testing that it does not generate an error
    Preprocessor::TemplateArray<0, V> empty;
    // Single item
    Preprocessor::TemplateArray<1, V> singleton{V<0>{1}};
    INFO("Get singleton item")
    CHECK(singleton.get<0>().value == 1);

    // Multiple items
    Preprocessor::TemplateArray<3, V> triple{V<2>{3}, V<1>{2}, V<0>{1}};
    INFO("Check multiple item content")
    CHECK(triple.get<0>().value == 1);
    CHECK(triple.get<1>().value == 2);
    CHECK(triple.get<2>().value == 3);
}

TEST_CASE("Getters", "[TemplateArray]") {
    Preprocessor::TemplateArray<2, V> sample{V<1>{1}, V<0>{1}};
    REQUIRE(sample.get<0>().value == 1);
    REQUIRE(sample.get<1>().value == 1);

    // Other accessors
    INFO("Check accessors")
    REQUIRE(sample[Preprocessor::tag<0>{}].value == 1);
    REQUIRE(sample[Preprocessor::itag<0>{}].value == 1);
    // Overwrite
    sample[Preprocessor::tag<0>{}] = V<0>{2};
    // Check
    INFO("Check post modification");
    REQUIRE(sample.get<0>().value == 2);
    REQUIRE(sample.get<1>().value == 1);
}