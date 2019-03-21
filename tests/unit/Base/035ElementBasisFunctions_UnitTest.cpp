
/*
 This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
 http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found below.


 Copyright (c) 2019, University of Twente
 All rights reserved.

 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "Base/ElementBasisFunctions.h"
#include "Base/AssembleBasisFunctionSet.h"
#include "Logger.h"

int main()
{
    using namespace Base;

    {
        ElementBasisFunctions empty; // Empty basis function set
        logger.assert_always(0 == empty.getTotalLocalNumberOfBasisFunctions(), "Expecting zero local basis functions");
        empty.validatePositions();
    }

    // Test case for testing with basis functions
    std::shared_ptr<BasisFunctionSet> set = std::make_shared<BasisFunctionSet> (3);
    AssembleBasisFunctionSet_3D_Ord3_A1(*set);
    ElementBasisFunctions::CollectionOfBasisFunctionSets sets;
    sets.emplace_back(set);

    {
        // Test that when we clear all basis functions the result is a completely empty set.
        const std::size_t UNKNOWNS = 3;
        ElementBasisFunctions almostEmpty(&sets, UNKNOWNS);
        for (std::size_t i = 0; i < UNKNOWNS; ++i)
        {
            almostEmpty.validatePositions();
            logger.assert_always(0 == almostEmpty.getNumberOfLocalBasisFunctions(i), "Expecting no local basis functions");
            logger.assert_always(0 == almostEmpty.getNumberOfBasisFunctions(i), "Expecting no basis functions");
        }
        logger.assert_always(0 == almostEmpty.getMaximumOrder(), "No basis functions, expecting order 0");
        logger.assert_always(0 == almostEmpty.getTotalLocalNumberOfBasisFunctions(),
                             "Expecting total of 0 local basis functions");

        // Register the basis function as first basis function for second unknown
        almostEmpty.registerBasisFunctionPosition(1, 0, 0);
        almostEmpty.validatePositions();
        logger.assert_always(set->size() == almostEmpty.getNumberOfBasisFunctions(1),
                "Matching basis function count");
        logger.assert_always(set->size() == almostEmpty.getNumberOfLocalBasisFunctions(1),
                "Expecting registered for local functions");
        logger.assert_always(set->size() == almostEmpty.getTotalLocalNumberOfBasisFunctions(),
                "Expecting no registered other basis functions");
        // Register it on only the second position for the zeroth unknown
        // hence not as local basis functions
        almostEmpty.clearBasisFunctionPosition(1);
        almostEmpty.registerBasisFunctionPosition(0, 1, 0);
        logger.assert_always(0 == almostEmpty.getNumberOfLocalBasisFunctions(0), "Expecting no local basis functions");
        logger.assert_always(set->size() == almostEmpty.getNumberOfBasisFunctions(0),
                "Expecting non local basis function");
        for(std::size_t i = 1; i < UNKNOWNS; ++i)
        {
            logger.assert_always(0 == almostEmpty.getNumberOfBasisFunctions(i),
                    "Expecting no basis functions for other unknowns");
        }
    }

    {
        // Adding secondary lower order basis functions
        std::shared_ptr<BasisFunctionSet> secondSet = std::make_shared<BasisFunctionSet> (1);
        AssembleBasisFunctionSet_3D_Ord1_A1(*secondSet);
        sets.emplace_back(secondSet);

        const std::size_t UNKNOWNS = 2;
        ElementBasisFunctions filled(&sets, UNKNOWNS);
        for(std::size_t i = 0; i < UNKNOWNS; ++i)
        {
            filled.validatePositions();
        }
        filled.registerBasisFunctionPosition(0, 0, 0); // set for unknown 0
        filled.registerBasisFunctionPosition(1, 0, 1); // secondSet for unknown 1



    }

    return 0;
}