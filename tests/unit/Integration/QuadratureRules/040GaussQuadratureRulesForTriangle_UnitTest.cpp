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
#include "Integration/QuadratureRules/GaussQuadratureRulesForTriangle.h"
#include "Logger.h"
#include <typeinfo>

#include "Utilities/BasisFunctions2DH1ConformingTriangle.h"
#include "Geometry/ReferenceTriangle.h"
#include "Base/BasisFunctionSet.h"
#include "Geometry/PointReference.h"
#include "LinearAlgebra/MiddleSizeVector.h"
#include <cmath>

#include "../../catch.hpp"

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
    CHECK((typeid(refGeo) == typeid(Geometry::ReferenceTriangle)));
    std::cout.precision(14);
    Base::BasisFunctionSet* functions =
        Utilities::createDGBasisFunctionSet2DH1Triangle(expectedOrder);
    for (std::size_t i = 0; i < functions->size(); ++i) {
        double integrated = 0;
        for (std::size_t j = 0; j < test.getNumberOfPoints(); ++j) {
            const Geometry::PointReference<2>& point = test.getPoint(j);
            integrated += test.weight(j) * functions->eval(i, point);
        }
        logger(INFO, "basisfunction %: integral equal to %", i, integrated);
        if (i < 3) {
            INFO("integration");
            CHECK((std::abs(integrated - 1. / 6.) < 1e-12));
        } else if (i < 6) {
            INFO("integration");
            CHECK((std::abs(integrated + 0.102062072616) < 1e-12));
        } else if (i == 6) {
            INFO("integration");
            CHECK((std::abs(integrated - 1. / 20.) < 1e-12));
        } else if (11 < i && i < 15) {
            INFO("integration");
            CHECK((std::abs(integrated - 0.0129918659263) < 1e-12));
        } else if (i == 15 || i == 17) {
            INFO("integration");
            CHECK((std::abs(integrated + 0.010001653302) < 1e-12));
        } else if (i == 16) {
            INFO("integration");
            CHECK((std::abs(integrated + 0.00396825396825) < 1e-12));
        } else if (i == 22) {
            INFO("integration");
            CHECK((std::abs(integrated - 0.001467281692237) < 1e-12));
        } else if (i == 23) {
            INFO("integration");
            CHECK((std::abs(integrated + 0.001467281692237) < 1e-12));
        } else if (24 < i && i < 28) {
            INFO("integration");
            CHECK((std::abs(integrated + 0.000814308291636) < 1e-12));
        } else if (i == 28 || i == 32) {
            INFO("integration");
            CHECK((std::abs(integrated - 0.001994639807826) < 1e-12));
        } else if (i == 29 || i == 31) {
            INFO("integration");
            CHECK((std::abs(integrated - 0.001663741054687) < 1e-12));
        } else if (i == 30) {
            INFO("integration");
            CHECK((std::abs(integrated - 0.0025173611111111) < 1e-12));
        } else if (i == 37) {
            INFO("integration");
            CHECK((std::abs(integrated + 0.0010300275676522) < 1e-12));
        } else if (i == 38) {
            INFO("integration");
            CHECK((std::abs(integrated - 0.0000984282481784) < 1e-12));
        } else if (i == 39) {
            INFO("integration");
            CHECK((std::abs(integrated + 0.0000984282481784) < 1e-12));
        } else if (i == 40) {
            INFO("integration");
            CHECK((std::abs(integrated - 0.0010300275676522) < 1e-12));
        } else if (41 < i && i < 45) {
            INFO("integration");
            CHECK((std::abs(integrated - 0.0001528243743039) < 1e-12));
        } else if (i == 45 || i == 51) {
            INFO("integration");
            CHECK((std::abs(integrated + 0.0003743417373047) < 1e-12));
        } else if (i == 46 || i == 50) {
            INFO("integration");
            CHECK((std::abs(integrated + 0.0045990061560235) < 1e-12));
        } else if (i == 47 || i == 49) {
            INFO("integration");
            CHECK((std::abs(integrated + 0.0005242977910035) < 1e-12));
        } else if (i == 48) {
            INFO("integration");
            CHECK((std::abs(integrated + 0.001020698051948) < 1e-12));
        } else if (i == 56) {
            INFO("integration");
            CHECK((std::abs(integrated - 0.0004027275873253) < 1e-12));
        } else if (i == 57) {
            INFO("integration");
            CHECK((std::abs(integrated + 0.001833588494787) < 1e-12));
        } else if (i == 58) {
            INFO("integration");
            CHECK((std::abs(integrated - 0.0002669631696782) < 1e-12));
        } else if (i == 59) {
            INFO("integration");
            CHECK((std::abs(integrated + 0.0002669631696782) < 1e-12));
        } else if (i == 60) {
            INFO("integration");
            CHECK((std::abs(integrated - 0.001833588494787) < 1e-12));
        } else if (i == 61) {
            INFO("integration");
            CHECK((std::abs(integrated + 0.0004027275873253) < 1e-12));
        } else if (62 < i && i < 66) {
            INFO("integration");
            CHECK((std::abs(integrated + 0.0000445921151835) < 1e-12));
        } else if (i == 66 || i == 74) {
            INFO("integration");
            CHECK((std::abs(integrated - 0.00010922792875103) < 1e-12));
        } else if (i == 67 || i == 73) {
            INFO("integration");
            CHECK((std::abs(integrated - 6.6692306548705e-05) < 1e-12));
        } else if (i == 68 || i == 72) {
            INFO("integration");
            CHECK((std::abs(integrated + 5.5593265874347e-05) < 1e-12));
        } else if (i == 69 || i == 71) {
            INFO("integration");
            CHECK((std::abs(integrated - 0.0018449672212449) < 1e-12));
        } else if (i == 70) {
            INFO("integration");
            CHECK((std::abs(integrated - 0.00021701388888889) < 1e-12));
        } else {
            INFO("integration");
            CHECK((std::abs(integrated) < 1e-12));
        }
    }

    delete functions;
}

TEST_CASE("040GaussQuadratureRulesForTriangle_UnitTest",
          "[040GaussQuadratureRulesForTriangle_UnitTest]") {

    testRule(QuadratureRules::Tn2_1_1::Instance(), 1);
    testRule(QuadratureRules::Tn2_2_3::Instance(), 2);
    testRule(QuadratureRules::Tn2_3_4::Instance(), 3);
    testRule(QuadratureRules::Tn2_4_6::Instance(), 4);
    testRule(QuadratureRules::T2_5_7::Instance(), 5);
    testRule(QuadratureRules::T2_6_12::Instance(), 6);
    testRule(QuadratureRules::T2_7_13::Instance(), 7);
    testRule(QuadratureRules::T2_8_16::Instance(), 8);
    testRule(QuadratureRules::T2_9_19::Instance(), 9);
    testRule(QuadratureRules::T2_10_25::Instance(), 10);
    testRule(QuadratureRules::T2_11_28::Instance(), 11);
}
