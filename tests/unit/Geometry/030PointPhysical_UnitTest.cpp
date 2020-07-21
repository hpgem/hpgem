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
#include "Geometry/PointPhysical.h"
#include <iostream>
#include "Logger.h"
#include <cmath>


#include "../catch.hpp"

using namespace hpgem;
using Geometry::PointPhysical;

TEST_CASE("030PointPhysical_UnitTest", "[030PointPhysical_UnitTest]") {

    double *coord0 = nullptr;
    double coord1[] = {1.1};
    double coord2[] = {1.2, 2.2};
    double coord3[] = {1.3, 2.3, 3.3};
    double coord4[] = {1.4, 2.4, 3.4, 4.4};

    LinearAlgebra::SmallVector<0> vec0(coord0);
    LinearAlgebra::SmallVector<1> vec1(coord1);
    LinearAlgebra::SmallVector<2> vec2(coord2);
    LinearAlgebra::SmallVector<3> vec3(coord3);
    LinearAlgebra::SmallVector<4> vec4(coord4);

    Geometry::Point<0> pb0;
    Geometry::Point<1> pb1;
    Geometry::Point<2> pb2;
    Geometry::Point<3> pb3;
    Geometry::Point<4> pb4;

    // testing constructors up to DIM=4
    Geometry::PointPhysical<0> p0, pp0(pb0), pv0(vec0);
    Geometry::PointPhysical<1> p1, pp1(pb1), pv1(vec1);
    Geometry::PointPhysical<2> p2, pp2(pb2), pv2(vec2);
    Geometry::PointPhysical<3> p3, pp3(pb3), pv3(vec3);
    Geometry::PointPhysical<4> p4, pp4(pb4), pv4(vec4);

    INFO("1D default constructor");
    CHECK((p1[0] == 0.));
    for (std::size_t i = 0; i < 2; ++i) {
        INFO("2D default constructor");
        CHECK((p2[i] == 0.));
    }
    for (std::size_t i = 0; i < 3; ++i) {
        INFO("3D default constructor");
        CHECK((p3[i] == 0.));
    }
    for (std::size_t i = 0; i < 4; ++i) {
        INFO("4D default constructor");
        CHECK((p4[i] == 0.));
    }

    INFO("1D copy constructor");
    CHECK((pp1[0] == 0.));
    for (std::size_t i = 0; i < 2; ++i) {
        INFO("2D copy constructor");
        CHECK((pp2[i] == 0.));
    }
    for (std::size_t i = 0; i < 3; ++i) {
        INFO("3D copy constructor");
        CHECK((pp3[i] == 0.));
    }
    for (std::size_t i = 0; i < 4; ++i) {
        INFO("4D copy constructor");
        CHECK((pp4[i] == 0.));
    }

    INFO("1D from NumericalVector constructor");
    CHECK((std::abs(pv1[0] - 1.1) < 1e-12));
    for (std::size_t i = 0; i < 2; ++i) {
        INFO("2D from NumericalVector constructor");
        CHECK((std::abs(pv2[i] - 1.2 - i) < 1e-12));
    }
    for (std::size_t i = 0; i < 3; ++i) {
        INFO("3D from NumericalVector constructor");
        CHECK((std::abs(pv3[i] - 1.3 - i) < 1e-12));
    }
    for (std::size_t i = 0; i < 4; ++i) {
        INFO("4D from NumericalVector constructor");
        CHECK((std::abs(pv4[i] - 1.4 - i) < 1e-12));
    }

    // testing operators

    const PointPhysical<0> pr0 = pp0 = pv0;
    const PointPhysical<1> pr1 = pp1 = pv1;
    const PointPhysical<2> pr2 = pp2 = pv2;
    const PointPhysical<3> pr3 = pp3 = pv3;
    const PointPhysical<4> pr4 = pp4 = pv4;

    pv0 * 6.;
    INFO("1D multiplication");
    CHECK((std::abs((pv1 * 5.)[0] - 5.5) < 1e-12));
    for (std::size_t i = 0; i < 2; ++i) {
        INFO("2D multiplication");
        CHECK((std::abs((pv2 * 4.)[i] - 4.8 - 4 * i) < 1e-12));
    }
    for (std::size_t i = 0; i < 3; ++i) {
        INFO("3D multiplication");
        CHECK((std::abs((pv3 * 3.)[i] - 3.9 - 3 * i) < 1e-12));
    }
    for (std::size_t i = 0; i < 4; ++i) {
        INFO("4D multiplication");
        CHECK((std::abs((pv4 * 2.)[i] - 2.8 - 2 * i) < 1e-12));
    }

    INFO("0D multiplication");
    CHECK(((pr0 * 0.) == p0));
    INFO("1D multiplication");
    CHECK(((pr1 * 0.) == p1));
    INFO("2D multiplication");
    CHECK(((pr2 * 0.) == p2));
    INFO("3D multiplication");
    CHECK(((pr3 * 0.) == p3));
    INFO("4D multiplication");
    CHECK(((pr4 * 0.) == p4));

    pp0 + pv0;
    INFO("1D addition");
    CHECK((std::abs((pp1 + pv1)[0] - 2.2) < 1e-12));
    for (std::size_t i = 0; i < 2; ++i) {
        INFO("2D addition");
        CHECK((std::abs((pp2 + pv2)[i] - 2.4 - 2 * i) < 1e-12));
    }
    for (std::size_t i = 0; i < 3; ++i) {
        INFO("3D addition");
        CHECK((std::abs((pp3 + pv3)[i] - 2.6 - 2 * i) < 1e-12));
    }
    for (std::size_t i = 0; i < 4; ++i) {
        INFO("4D addition");
        CHECK((std::abs((pp4 + pv4)[i] - 2.8 - 2 * i) < 1e-12));
    }

    pr0 + pv0;
    INFO("1D addition");
    CHECK((std::abs((pr1 + pv1)[0] - 2.2) < 1e-12));
    for (std::size_t i = 0; i < 2; ++i) {
        INFO("2D addition");
        CHECK((std::abs((pr2 + pv2)[i] - 2.4 - 2 * i) < 1e-12));
    }
    for (std::size_t i = 0; i < 3; ++i) {
        INFO("3D addition");
        CHECK((std::abs((pr3 + pv3)[i] - 2.6 - 2 * i) < 1e-12));
    }
    for (std::size_t i = 0; i < 4; ++i) {
        INFO("4D addition");
        CHECK((std::abs((pr4 + pv4)[i] - 2.8 - 2 * i) < 1e-12));
    }

    pp0 - pv0;
    INFO("1D subtraction");
    CHECK((std::abs((pp1 - pv1)[0]) < 1e-12));
    for (std::size_t i = 0; i < 2; ++i) {
        INFO("2D subtraction");
        CHECK((std::abs((pp2 - pv2)[i]) < 1e-12));
    }
    for (std::size_t i = 0; i < 3; ++i) {
        INFO("3D subtraction");
        CHECK((std::abs((pp3 - pv3)[i]) < 1e-12));
    }
    for (std::size_t i = 0; i < 4; ++i) {
        INFO("4D subtraction");
        CHECK((std::abs((pp4 - pv4)[i]) < 1e-12));
    }

    pr0 - pv0;
    INFO("1D subtraction");
    CHECK((std::abs((pr1 - pv1)[0]) < 1e-12));
    for (std::size_t i = 0; i < 2; ++i) {
        INFO("2D subtraction");
        CHECK((std::abs((pr2 - pv2)[i]) < 1e-12));
    }
    for (std::size_t i = 0; i < 3; ++i) {
        INFO("3D subtraction");
        CHECK((std::abs((pr3 - pv3)[i]) < 1e-12));
    }
    for (std::size_t i = 0; i < 4; ++i) {
        INFO("4D subtraction");
        CHECK((std::abs((pr4 - pv4)[i]) < 1e-12));
    }

    // and point already works so done
}
