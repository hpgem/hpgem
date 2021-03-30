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

#include "ElementReorder.h"

#include "../catch.hpp"

TEST_CASE("Reorder", "Reorder") {

    Preprocessor::ElementReorder r;

    r.addElementType(2, "square", {0, 2, 3, 1});
   
    std::vector<size_t> index1{5, 8, 7, 4};
    std::vector<size_t> index1_ref1{5, 7, 4, 8};
    std::vector<size_t> index1_ref2 = index1;

    r.reorderToHpGem(2,index1);
    INFO("square -> hpgem");
    REQUIRE(index1==index1_ref1);
    r.reorderFromHpGem(2,index1);
     INFO("square <- hpgem");
    REQUIRE(index1==index1_ref2);

    r.addElementType(3, "cube", {2, 3, 4, 1, 5, 6, 7, 0});

    std::vector<size_t> index2{1, 2, 3, 4, 5, 6, 7, 8};
    std::vector<size_t> index2_ref1{3, 4, 5, 2, 6, 7, 8, 1};
    std::vector<size_t> index2_ref2 = index2;


    r.reorderToHpGem(3,index2);
    INFO("cube -> hpgem");
    REQUIRE(index2==index2_ref1);
    r.reorderFromHpGem(3,index2);
    INFO("cube <- hpgem");
    REQUIRE(index2==index2_ref2);


}