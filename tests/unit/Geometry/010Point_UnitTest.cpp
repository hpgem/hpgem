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
#include "Geometry/Point.h"
#include "Logger.h"
#include <iostream>
#include <cmath>

#define CATCH_CONFIG_MAIN
#include "../catch.hpp"

using namespace hpgem;
using Geometry::Point;

TEST_CASE("010Point_UnitTest", "[010Point_UnitTest]") {
    // when using zero-dimensional points, take care that the 0-size double*
    // array is a language extension
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

    // testing constructors up to DIM=4
    Point<0> p0, pc0(coord0), pp0(p0), pv0(vec0);
    Point<1> p1, pc1(coord1), pp1(p1), pv1(vec1);
    Point<2> p2, pc2(coord2), pp2(p2), pv2(vec2);
    Point<3> p3, pc3(coord3), pp3(p3), pv3(vec3);
    Point<4> p4, pc4(coord4), pp4(p4), pv4(vec4);

    // testing operator[]

    INFO("1D default constructor or operator[] of Point");
    CHECK((p1[0] == 0.));
    for (std::size_t i = 0; i < 2; ++i) {
        INFO("2D default constructor or operator[] of Point");
        CHECK((p2[i] == 0.));
    }
    for (std::size_t i = 0; i < 3; ++i) {
        INFO("3D default constructor or operator[] of Point");
        CHECK((p3[i] == 0.));
    }
    for (std::size_t i = 0; i < 4; ++i) {
        INFO("4D default constructor or operator[] of Point");
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

    INFO("1D from array constructor or operator[] of Point");
    CHECK((std::abs(pc1[0] - 1.1) < 1e-12));
    for (std::size_t i = 0; i < 2; ++i) {
        INFO("2D from array constructor or operator[] of Point");
        CHECK((std::abs(pc2[i] - 1.2 - i) < 1e-12));
    }
    for (std::size_t i = 0; i < 3; ++i) {
        INFO("3D from array constructor or operator[] of Point");
        CHECK((std::abs(pc3[i] - 1.3 - i) < 1e-12));
    }
    for (std::size_t i = 0; i < 4; ++i) {
        INFO("4D from array constructor or operator[] of Point");
        CHECK((std::abs(pc4[i] - 1.4 - i) < 1e-12));
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

    // testing setCoordinates and setCoordinate

    p1.setCoordinates(vec1);
    INFO("1D setCoordinates");
    CHECK((std::abs(p1[0] - 1.1) < 1e-12));
    p2.setCoordinates(vec2);
    for (std::size_t i = 0; i < 2; ++i) {
        INFO("2D setCoordinates");
        CHECK((std::abs(p2[i] - 1.2 - i) < 1e-12));
    }
    p3.setCoordinates(vec3);
    for (std::size_t i = 0; i < 3; ++i) {
        INFO("3D setCoordinates");
        CHECK((std::abs(p3[i] - 1.3 - i) < 1e-12));
    }
    p4.setCoordinates(vec4);
    for (std::size_t i = 0; i < 4; ++i) {
        INFO("4D setCoordinates");
        CHECK((std::abs(p4[i] - 1.4 - i) < 1e-12));
    }

    p1.setCoordinate(0, 0.9);
    INFO("1D setCoordinate");
    CHECK((std::abs(p1[0] - 0.9) < 1e-12));
    for (std::size_t i = 0; i < 2; ++i) {
        p2.setCoordinate(i, 0.8 + double(i));
        for (std::size_t j = 0; j <= i; ++j) {
            INFO("2D setCoordinate");
            CHECK((std::abs(p2[j] - 0.8 - j) < 1e-12));
        }
        for (std::size_t j = i + 1; j < 2; ++j) {
            INFO("2D setCoordinate");
            CHECK((std::abs(p2[j] - 1.2 - j) < 1e-12));
        }
    }
    for (std::size_t i = 0; i < 3; ++i) {
        p3.setCoordinate(i, 0.7 + i);
        for (std::size_t j = 0; j <= i; ++j) {
            INFO("3D setCoordinate");
            CHECK((std::abs(p3[j] - 0.7 - j) < 1e-12));
        }
        for (std::size_t j = i + 1; j < 3; ++j) {
            INFO("3D setCoordinate");
            CHECK((std::abs(p3[j] - 1.3 - j) < 1e-12));
        }
    }
    for (std::size_t i = 0; i < 4; ++i) {
        p4.setCoordinate(i, 0.6 + i);
        for (std::size_t j = 0; j <= i; ++j) {
            INFO("4D setCoordinate");
            CHECK((std::abs(p4[j] - 0.6 - j) < 1e-12));
        }
        for (std::size_t j = i + 1; j < 4; ++j) {
            INFO("4D setCoordinate");
            CHECK((std::abs(p4[j] - 1.4 - j) < 1e-12));
        }
    }

    // testing operators

    const Point<0> pr0 = pc0 = p0;
    const Point<1> pr1 = pc1 = p1;
    INFO("1D assignment operator");
    CHECK((std::abs(pc1[0] - 0.9) < 1e-12));
    INFO("1D assignment operator");
    CHECK((std::abs(pr1[0] - 0.9) < 1e-12));
    const Point<2> pr2 = pc2 = p2;
    for (std::size_t i = 0; i < 2; ++i) {
        INFO("2D assignment operator");
        CHECK((std::abs(pc2[i] - 0.8 - i) < 1e-12));
        INFO("2D assignment operator");
        CHECK((std::abs(pr2[i] - 0.8 - i) < 1e-12));
    }
    const Point<3> pr3 = pc3 = p3;
    for (std::size_t i = 0; i < 3; ++i) {
        INFO("3D assignment operator");
        CHECK((std::abs(pc3[i] - 0.7 - i) < 1e-12));
        INFO("3D assignment operator");
        CHECK((std::abs(pr3[i] - 0.7 - i) < 1e-12));
    }
    const Point<4> pr4 = pc4 = p4;
    for (std::size_t i = 0; i < 4; ++i) {
        INFO("4D assignment operator");
        CHECK((std::abs(pc4[i] - 0.6 - i) < 1e-12));
        INFO("4D assignment operator");
        CHECK((std::abs(pr4[i] - 0.6 - i) < 1e-12));
    }

    INFO("0D equality operator");
    CHECK((pr0 == pc0 && pc0 == pr0 && pc0 == p0));
    INFO("1D equality operator");
    CHECK((pr1 == pc1 && pc1 == pr1 && pc1 == p1 &&
           !(pr1 == pv1 || pv1 == pr1 || p1 == pv1)));
    INFO("2D equality operator");
    CHECK((pr2 == pc2 && pc2 == pr2 && pc2 == p2 &&
           !(pr2 == pv2 || pv2 == pr2 || p2 == pv2)));
    INFO("3D equality operator");
    CHECK((pr3 == pc3 && pc3 == pr3 && pc3 == p3 &&
           !(pr3 == pv3 || pv3 == pr3 || p3 == pv3)));
    INFO("4D equality operator");
    CHECK((pr4 == pc4 && pc4 == pr4 && pc4 == p4 &&
           !(pr4 == pv4 || pv4 == pr4 || p4 == pv4)));

    pc0 += p0;
    pc1 += p1;
    INFO("1D increment operator");
    CHECK((std::abs(pc1[0] - 1.8) < 1e-12));
    pc2 += p2;
    for (std::size_t i = 0; i < 2; ++i) {
        INFO("2D increment operator");
        CHECK((std::abs(pc2[i] - 1.6 - 2 * i) < 1e-12));
    }
    pc3 += p3;
    for (std::size_t i = 0; i < 3; ++i) {
        INFO("3D increment operator");
        CHECK((std::abs(pc3[i] - 1.4 - 2 * i) < 1e-12));
    }
    pc4 += p4;
    for (std::size_t i = 0; i < 4; ++i) {
        INFO("4D increment operator");
        CHECK((std::abs(pc4[i] - 1.2 - 2 * i) < 1e-12));
    }

    pc0 -= pv0;
    pc1 -= pv1;
    INFO("1D decrement operator");
    CHECK((std::abs(pc1[0] - 0.7) < 1e-12));
    pc2 -= pv2;
    for (std::size_t i = 0; i < 2; ++i) {
        INFO("2D decrement operator");
        CHECK((std::abs(pc2[i] - 0.4 - i) < 1e-12));
    }
    pc3 -= pv3;
    for (std::size_t i = 0; i < 3; ++i) {
        INFO("3D decrement operator");
        CHECK((std::abs(pc3[i] - 0.1 - i) < 1e-12));
    }
    pc4 -= pv4;
    for (std::size_t i = 0; i < 4; ++i) {
        INFO("4D decrement operator or negative numbers");
        CHECK((std::abs(pc4[i] + 0.2 - i) < 1e-12));
    }

    pc0 *= 2.;
    pc1 *= 3.;
    INFO("1D multiply operator");
    CHECK((std::abs(pc1[0] - 2.1) < 1e-12));
    pc2 *= 4.;
    for (std::size_t i = 0; i < 2; ++i) {
        INFO("2D multiply operator");
        CHECK((std::abs(pc2[i] - 1.6 - 4 * i) < 1e-12));
    }
    pc3 *= 5.;
    for (std::size_t i = 0; i < 3; ++i) {
        INFO("3D multiply operator");
        CHECK((std::abs(pc3[i] - 0.5 - 5 * i) < 1e-12));
    }
    pc4 *= 6.;
    for (std::size_t i = 0; i < 4; ++i) {
        INFO("4D multiply operator");
        CHECK((std::abs(pc4[i] + 1.2 - 6 * i) < 1e-12));
    }

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
    CHECK(((pr0 * 0.) == pp0));
    INFO("1D multiplication");
    CHECK(((pr1 * 0.) == pp1));
    INFO("2D multiplication");
    CHECK(((pr2 * 0.) == pp2));
    INFO("3D multiplication");
    CHECK(((pr3 * 0.) == pp3));
    INFO("4D multiplication");
    CHECK(((pr4 * 0.) == pp4));

    pc0 + pv0;
    INFO("1D addition");
    CHECK((std::abs((pc1 + pv1)[0] - 3.2) < 1e-12));
    for (std::size_t i = 0; i < 2; ++i) {
        INFO("2D addition");
        CHECK((std::abs((pc2 + pv2)[i] - 2.8 - 5 * i) < 1e-12));
    }
    for (std::size_t i = 0; i < 3; ++i) {
        INFO("3D addition");
        CHECK((std::abs((pc3 + pv3)[i] - 1.8 - 6 * i) < 1e-12));
    }
    for (std::size_t i = 0; i < 4; ++i) {
        INFO("4D addition");
        CHECK((std::abs((pc4 + pv4)[i] - 0.2 - 7 * i) < 1e-12));
    }

    pr0 + pv0;
    INFO("1D addition");
    CHECK((std::abs((pr1 + pv1)[0] - 2.) < 1e-12));
    for (std::size_t i = 0; i < 2; ++i) {
        INFO("2D addition");
        CHECK((std::abs((pr2 + pv2)[i] - 2. - 2 * i) < 1e-12));
    }
    for (std::size_t i = 0; i < 3; ++i) {
        INFO("3D addition");
        CHECK((std::abs((pr3 + pv3)[i] - 2. - 2 * i) < 1e-12));
    }
    for (std::size_t i = 0; i < 4; ++i) {
        INFO("4D addition");
        CHECK((std::abs((pr4 + pv4)[i] - 2. - 2 * i) < 1e-12));
    }

    pc0 - pv0;
    INFO("1D subtraction");
    CHECK((std::abs((pc1 - pv1)[0] - 1.) < 1e-12));
    for (std::size_t i = 0; i < 2; ++i) {
        INFO("2D subtraction");
        CHECK((std::abs((pc2 - pv2)[i] - 0.4 - 3 * i) < 1e-12));
    }
    for (std::size_t i = 0; i < 3; ++i) {
        INFO("3D subtraction");
        CHECK((std::abs((pc3 - pv3)[i] + 0.8 - 4 * i) < 1e-12));
    }
    for (std::size_t i = 0; i < 4; ++i) {
        INFO("4D subtraction");
        CHECK((std::abs((pc4 - pv4)[i] + 2.6 - 5 * i) < 1e-12));
    }

    pr0 - pv0;
    INFO("1D subtraction");
    CHECK((std::abs((pr1 - pv1)[0] + 0.2) < 1e-12));
    for (std::size_t i = 0; i < 2; ++i) {
        INFO("2D subtraction");
        CHECK((std::abs((pr2 - pv2)[i] + 0.4) < 1e-12));
    }
    for (std::size_t i = 0; i < 3; ++i) {
        INFO("3D subtraction");
        CHECK((std::abs((pr3 - pv3)[i] + 0.6) < 1e-12));
    }
    for (std::size_t i = 0; i < 4; ++i) {
        INFO("4D subtraction");
        CHECK((std::abs((pr4 - pv4)[i] + 0.8) < 1e-12));
    }

    // testing size

    INFO("size of a 0D point");
    CHECK((p0.size() == 0 && pp0.size() == 0 && pr0.size() == 0 &&
           pc0.size() == 0 && pv0.size() == 0));
    INFO("size of a 1D point");
    CHECK((p1.size() == 1 && pp1.size() == 1 && pr1.size() == 1 &&
           pc1.size() == 1 && pv1.size() == 1));
    INFO("size of a 2D point");
    CHECK((p2.size() == 2 && pp2.size() == 2 && pr2.size() == 2 &&
           pc2.size() == 2 && pv2.size() == 2));
    INFO("size of a 3D point");
    CHECK((p3.size() == 3 && pp3.size() == 3 && pr3.size() == 3 &&
           pc3.size() == 3 && pv3.size() == 3));
    INFO("size of a 4D point");
    CHECK((p4.size() == 4 && pp4.size() == 4 && pr4.size() == 4 &&
           pc4.size() == 4 && pv4.size() == 4));

    // testing getCoordinate and getCoordinates

    INFO("1D getCoordinate");
    CHECK((p1.getCoordinate(0) == p1[0]));
    for (std::size_t i = 0; i < 2; ++i) {
        INFO("2D getCoordinate");
        CHECK((p2.getCoordinate(i) == p2[i]));
    }
    for (std::size_t i = 0; i < 3; ++i) {
        INFO("3D getCoordinate");
        CHECK((p3.getCoordinate(i) == p3[i]));
    }
    for (std::size_t i = 0; i < 4; ++i) {
        INFO("4D getCoordinate");
        CHECK((p4.getCoordinate(i) == p4[i]));
    }

    INFO("0D getCoordinates");
    CHECK((pv0.getCoordinates() == vec0));
    INFO("1D getCoordinates");
    CHECK((pv1.getCoordinates() == vec1));
    INFO("2D getCoordinates");
    CHECK((pv2.getCoordinates() == vec2));
    INFO("3D getCoordinates");
    CHECK((pv3.getCoordinates() == vec3));
    INFO("4D getCoordinates");
    CHECK((pv4.getCoordinates() == vec4));

    // testing friends

    -pc0;
    INFO("1D unary -");
    CHECK((std::abs((-pc1)[0] + 2.1) < 1e-12));
    for (std::size_t i = 0; i < 2; ++i) {
        INFO("2D unary -");
        CHECK((std::abs((-pc2)[i] + 1.6 + 4 * i) < 1e-12));
    }
    for (std::size_t i = 0; i < 3; ++i) {
        INFO("3D unary -");
        CHECK((std::abs((-pc3)[i] + 0.5 + 5 * i) < 1e-12));
    }
    for (std::size_t i = 0; i < 4; ++i) {
        INFO("4D unary -");
        CHECK((std::abs((-pc4)[i] - 1.2 + 6 * i) < 1e-12));
    }

    6. * pv0;
    INFO("1D left multiplication");
    CHECK((std::abs((5. * pv1)[0] - 5.5) < 1e-12));
    for (std::size_t i = 0; i < 2; ++i) {
        INFO("2D left multiplication");
        CHECK((std::abs((4. * pv2)[i] - 4.8 - 4 * i) < 1e-12));
    }
    for (std::size_t i = 0; i < 3; ++i) {
        INFO("3D left multiplication");
        CHECK((std::abs((3. * pv3)[i] - 3.9 - 3 * i) < 1e-12));
    }
    for (std::size_t i = 0; i < 4; ++i) {
        INFO("4D left multiplication");
        CHECK((std::abs((2. * pv4)[i] - 2.8 - 2 * i) < 1e-12));
    }

    std::cout << p0 << p1 << p2 << p3 << p4 << pr0 << pr1 << pr2 << pr3 << pr4
              << std::endl;
}
