
/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2019, University of Twente
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

#include "Base/ElementBasisFunctions.h"
#include "FE/BasisFunctionsMonomials.h"
#include "Logger.h"

#include "../catch.hpp"

using namespace hpgem;
TEST_CASE("035ElementBasisFunctions_UnitTest",
          "[035ElementBasisFunctions_UnitTest]") {
    using namespace Base;

    {
        ElementBasisFunctions empty;  // Empty basis function set
        INFO("Expecting zero local basis functions");
        CHECK(0 == empty.getTotalLocalNumberOfBasisFunctions());
        empty.validatePositions();
    }

    // Test case for testing with basis functions
    std::shared_ptr<FE::BasisFunctionSet> set =
        std::make_shared<FE::BasisFunctionSet>(3);
    FE::assembleMonomialBasisFunctions3D(*set, 3);
    ElementBasisFunctions::CollectionOfBasisFunctionSets sets;
    sets.emplace_back(set);

    {
        // Test that when we clear all basis functions the result is a
        // completely empty set.
        const std::size_t UNKNOWNS = 3;
        ElementBasisFunctions almostEmpty(&sets, UNKNOWNS);
        for (std::size_t i = 0; i < UNKNOWNS; ++i) {
            almostEmpty.validatePositions();
            INFO("Expecting no local basis functions");
            CHECK(0 == almostEmpty.getNumberOfLocalBasisFunctions(i));
            INFO("Expecting no basis functions");
            CHECK(0 == almostEmpty.getNumberOfBasisFunctions(i));
        }
        INFO("No basis functions, expecting order 0");
        CHECK(0 == almostEmpty.getMaximumOrder());
        INFO("Expecting total of 0 local basis functions");
        CHECK(0 == almostEmpty.getTotalLocalNumberOfBasisFunctions());

        // Register the basis function as first basis function for second
        // unknown
        almostEmpty.registerBasisFunctionPosition(1, 0, 0);
        almostEmpty.validatePositions();
        INFO("Matching basis function count");
        CHECK(set->size() == almostEmpty.getNumberOfBasisFunctions(1));
        INFO("Expecting registered for local functions");
        CHECK(set->size() == almostEmpty.getNumberOfLocalBasisFunctions(1));
        INFO("Expecting no registered other basis functions");
        CHECK(set->size() == almostEmpty.getTotalLocalNumberOfBasisFunctions());
        // Register it on only the second position for the zeroth unknown
        // hence not as local basis functions
        almostEmpty.clearBasisFunctionPosition(1);
        almostEmpty.registerBasisFunctionPosition(0, 1, 0);
        INFO("Expecting no local basis functions");
        CHECK(0 == almostEmpty.getNumberOfLocalBasisFunctions(0));
        INFO("Expecting non local basis function");
        CHECK(set->size() == almostEmpty.getNumberOfBasisFunctions(0));
        for (std::size_t i = 1; i < UNKNOWNS; ++i) {
            INFO("Expecting no basis functions for other unknowns");
            CHECK(0 == almostEmpty.getNumberOfBasisFunctions(i));
        }
    }

    {
        // Adding secondary lower order basis functions
        std::shared_ptr<FE::BasisFunctionSet> secondSet =
            std::make_shared<FE::BasisFunctionSet>(1);
        FE::assembleMonomialBasisFunctions3D(*secondSet, 1);
        sets.emplace_back(secondSet);

        const std::size_t UNKNOWNS = 2;
        ElementBasisFunctions filled(&sets, UNKNOWNS);
        for (std::size_t i = 0; i < UNKNOWNS; ++i) {
            filled.validatePositions();
        }
        filled.registerBasisFunctionPosition(0, 0, 0);  // set for unknown 0
        filled.registerBasisFunctionPosition(1, 0,
                                             1);  // secondSet for unknown 1
    }
}