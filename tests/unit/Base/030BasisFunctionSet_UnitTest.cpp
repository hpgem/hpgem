/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2014, University of Twente
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

// naming convention: <Digit><ClassName>_UnitTest.cpp where <Digit> is a number
// that will make sure the unit tests are ordered such that the first failing
// unit test indicate the culprit class and other 'unit' tests may assume
// correct execution of all prior unit tests
#include "FE/BasisFunctionsMonomials.h"

#include "FE/BasisFunctionSet.h"
#include "Geometry/PointReference.h"
#include "LinearAlgebra/MiddleSizeVector.h"
#include "FE/BaseBasisFunction.h"

#include "Logger.h"

#include "../catch.hpp"

using namespace hpgem;
TEST_CASE("030BasisFunctionSet_UnitTest", "[030BasisFunctionSet_UnitTest]") {

    FE::BasisFunctionSet all1DbasisFunctions(5);

    FE::assembleMonomialBasisFunctions1D(all1DbasisFunctions, 5);
    Geometry::PointReference<1> point1D;
    for (std::size_t i = 0; i < all1DbasisFunctions.size(); ++i) {
        const FE::BaseBasisFunction* test = all1DbasisFunctions[i];
        for (point1D[0] = -1.5; point1D[0] < 1.51; point1D[0] += 0.1) {
            INFO("eval");
            CHECK((test->eval((point1D)) ==
                   all1DbasisFunctions.eval(i, (point1D))));
            INFO("derivative");
            CHECK((test->evalDeriv0((point1D)) ==
                   all1DbasisFunctions.evalDeriv(i, 0, (point1D))));
        }
    }

    FE::BasisFunctionSet all2DbasisFunctions(5);
    FE::assembleMonomialBasisFunctions2D(all2DbasisFunctions, 5);
    Geometry::PointReference<2> point2D;
    for (std::size_t i = 0; i < all2DbasisFunctions.size(); ++i) {
        const FE::BaseBasisFunction* test = all2DbasisFunctions[i];
        for (point2D[0] = -1.5; point2D[0] < 1.51; point2D[0] += 0.2) {
            for (point2D[1] = -1.5; point2D[1] < 1.51; point2D[1] += 0.2) {
                INFO("eval");
                CHECK((test->eval((point2D)) ==
                       all2DbasisFunctions.eval(i, (point2D))));
                INFO("derivative");
                CHECK((test->evalDeriv0((point2D)) ==
                       all2DbasisFunctions.evalDeriv(i, 0, (point2D))));
                INFO("derivative");
                CHECK((test->evalDeriv1((point2D)) ==
                       all2DbasisFunctions.evalDeriv(i, 1, (point2D))));
            }
        }
    }

    FE::BasisFunctionSet all3DbasisFunctions(5);
    FE::assembleMonomialBasisFunctions3D(all3DbasisFunctions, 5);
    Geometry::PointReference<3> point3D;
    for (std::size_t i = 0; i < all3DbasisFunctions.size(); ++i) {
        const FE::BaseBasisFunction* test = all3DbasisFunctions[i];
        for (point3D[0] = -1.5; point3D[0] < 1.51; point3D[0] += 0.6) {
            for (point3D[1] = -1.5; point3D[1] < 1.51; point3D[1] += 0.7) {
                for (point3D[2] = -1.5; point3D[2] < 1.51; point3D[2] += 0.8) {
                    INFO("eval");
                    CHECK((test->eval((point3D)) ==
                           all3DbasisFunctions.eval(i, (point3D))));
                    INFO("derivative");
                    CHECK((test->evalDeriv0((point3D)) ==
                           all3DbasisFunctions.evalDeriv(i, 0, (point3D))));
                    INFO("derivative");
                    CHECK((test->evalDeriv1((point3D)) ==
                           all3DbasisFunctions.evalDeriv(i, 1, (point3D))));
                    INFO("derivative");
                    CHECK((test->evalDeriv2((point3D)) ==
                           all3DbasisFunctions.evalDeriv(i, 2, (point3D))));
                }
            }
        }
    }
}
