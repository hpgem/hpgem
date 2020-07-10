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
#include "Integration/QuadratureRules/GaussQuadratureRulesForSquare.h"
#include "Logger.h"
#include <typeinfo>

#include "Utilities/BasisFunctions2DH1ConformingSquare.h"
#include "Geometry/ReferenceSquare.h"
#include "Base/BasisFunctionSet.h"
#include "Geometry/PointReference.h"
#include "LinearAlgebra/MiddleSizeVector.h"
#include <cmath>

#define CATCH_CONFIG_MAIN
#include "../catch.hpp"

using namespace hpgem;
void testRule(QuadratureRules::GaussQuadratureRule& test,
              std::size_t expectedOrder) {
    std::cout << test.getName() << std::endl;
    INFO("dimension");
    CHECK((test.dimension() == 2));
    INFO("order");
    CHECK((test.order() >= expectedOrder));
    const Geometry::ReferenceGeometry& refGeo = *test.forReferenceGeometry();
    INFO("forReferenceGeometry");
    CHECK((typeid(refGeo) == typeid(Geometry::ReferenceSquare)));

    Base::BasisFunctionSet* functions =
        Utilities::createDGBasisFunctionSet2DH1Square(expectedOrder);

    for (std::size_t i = 0; i < functions->size(); ++i) {
        double integrated = 0;
        for (std::size_t j = 0; j < test.getNumberOfPoints(); ++j) {
            const Geometry::PointReference<2>& point = test.getPoint(j);
            integrated += test.weight(j) * functions->eval(i, point);
        }
        if (i < 4) {
            INFO("integration");
            CHECK((std::abs(integrated - 1) < 1e-12));
        } else if (i < 8) {
            INFO("integration");
            CHECK((std::abs(integrated + std::sqrt(2. / 3.)) < 1e-12));
        } else if (i == 8) {
            INFO("integration");
            CHECK((std::abs(integrated - 2. / 3.) < 1e-12));
        } else {
            INFO("integration");
            CHECK((std::abs(integrated) < 1e-12));
        }
    }

    delete functions;
}

TEST_CASE("030GaussQuadratureRulesForSquare_UnitTest",
          "[030GaussQuadratureRulesForSquare_UnitTest]") {

    testRule(QuadratureRules::Cn2_1_1::Instance(), 1);
    testRule(QuadratureRules::Cn2_3_4::Instance(), 3);
    testRule(QuadratureRules::Cn2_5_9::Instance(), 5);
    testRule(QuadratureRules::C2_7_4::Instance(), 7);
    testRule(QuadratureRules::C2_9_5::Instance(), 9);
    testRule(QuadratureRules::C2_11_6::Instance(), 11);
}
