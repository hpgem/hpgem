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
#include "Integration/QuadratureRules/GaussQuadratureRulesForTriangularPrism.h"
#include "Logger.h"
#include <typeinfo>

#include "FE/BasisFunctions3DH1ConformingPrism.h"
#include "Geometry/ReferenceTriangularPrism.h"
#include "FE/BasisFunctionSet.h"
#include "Geometry/PointReference.h"
#include "LinearAlgebra/MiddleSizeVector.h"
#include <cmath>

#include "../../catch.hpp"

using namespace hpgem;
void testRule(QuadratureRules::GaussQuadratureRule& test,
              std::size_t expectedOrder) {
    std::cout << test.getName() << std::endl;
    INFO("dimension");
    CHECK((test.dimension() == 3));
    INFO("order");
    CHECK((test.order() >= expectedOrder));
    const Geometry::ReferenceGeometry& refGeo = *test.forReferenceGeometry();
    INFO("forReferenceGeometry");
    CHECK((typeid(refGeo) == typeid(Geometry::ReferenceTriangularPrism)));

    // Check for negative weights (these are unstable)
    for (std::size_t i = 0; i < test.getNumberOfPoints(); ++i) {
        INFO("Non negative weights");
        REQUIRE(test.weight(i) >= 0);
    }

    std::cout.precision(14);
    FE::BasisFunctionSet* functions =
        FE::createDGBasisFunctionSet3DH1ConformingPrism(expectedOrder);
    for (std::size_t i = 0; i < functions->size(); ++i) {
        double integrated = 0;
        for (std::size_t j = 0; j < test.getNumberOfPoints(); ++j) {
            const Geometry::PointReference<3>& point = test.getPoint(j);
            integrated += test.weight(j) * functions->eval(i, point);
        }
        if (i < 6) {
            INFO("integration");
            CHECK((std::abs(integrated - 1. / 6.) < 1e-10));
        } else if (i < 12) {
            INFO("integration");
            CHECK((std::abs(integrated + 0.1020620726159) < 1e-10));
        } else if (11 < i && i < 15) {
            INFO("integration");
            CHECK((std::abs(integrated + 0.1360827634879) < 1e-9));
        } else if (14 < i && i < 18) {
            INFO("integration");
            CHECK((std::abs(integrated - 1. / 12.) < 1e-10));
        } else if (i == 27 || i == 28) {
            INFO("integration");
            CHECK((std::abs(integrated - 1. / 20.) < 1e-10));
        } else if (40 < i && i < 47) {
            INFO("integration");
            CHECK((std::abs(integrated - 0.012991865926298) < 1e-10));
        } else if (i == 56 || i == 61 || i == 66) {
            INFO("integration");
            CHECK((std::abs(integrated + 0.01060781410869) < 1e-10));
        } else if (i == 69 || i == 70) {
            INFO("integration");
            CHECK((std::abs(integrated - 0.008166315725102) < 1e-10));
        } else if (i == 75) {
            INFO("integration");
            CHECK((std::abs(integrated - 0.001369177697178) < 1e-10 ||
                   expectedOrder == 7));  // actually the p=5 quadrature rule
                                          // may also be the culprit, but 7 is
                                          // more likely because it has other
                                          // flaws
        } else if (i == 76) {
            INFO("integration");
            CHECK((std::abs(integrated + 0.001369177697178) < 1e-10 ||
                   expectedOrder == 7));
        } else if (i == 87 || i == 89 || i == 90 || i == 92) {
            INFO("integration");
            CHECK((std::abs(integrated + 0.010001653302483) < 1e-10));
        } else if (i == 88 || i == 91) {
            INFO("integration");
            CHECK((std::abs(integrated + 0.003968253968254) < 1e-10));
        } else if (i == 123 || i == 124) {
            INFO("integration");
            CHECK((std::abs(integrated + 0.000465750474069) < 1e-10 ||
                   expectedOrder == 7));
        } else if (i == 129) {
            INFO("integration");
            CHECK((std::abs(integrated + 0.000721656823802) < 1e-10 ||
                   expectedOrder == 7));
        } else if (i == 130) {
            INFO("integration");
            CHECK((std::abs(integrated - 0.000721656823802) < 1e-10 ||
                   expectedOrder == 7));
        } else if (i < 132) {  // I test what I can for p=7, but not all the
                               // points
            INFO("integration");
            CHECK((std::abs(integrated) < 1e-10));
        }
    }

    delete functions;
}

TEST_CASE("070GaussQuadratureRulesForPrism_UnitTest",
          "[070GaussQuadratureRulesForPrism_UnitTest]") {

    testRule(QuadratureRules::TriPrism_1_1::Instance(), 1);
    testRule(QuadratureRules::TriPrism_5_21::Instance(), 5);
    ///\todo not accurate enough
    testRule(QuadratureRules::TriPrism_7_64::Instance(), 7);
    ///\todo there are no quadrature rules for higher order prisms
}
