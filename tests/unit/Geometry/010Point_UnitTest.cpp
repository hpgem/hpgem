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
// unit test indicate the culprit class and other 'unit' tests may assume correct
// execution of all prior unit tests
#include "Geometry/Point.h"
#include "Logger.h"
#include <iostream>
#include <cmath>
using Geometry::Point;

int main() {
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

    logger.assert_always((p1[0] == 0.),
                         "1D default constructor or operator[] of Point");
    for (std::size_t i = 0; i < 2; ++i) {
        logger.assert_always((p2[i] == 0.),
                             "2D default constructor or operator[] of Point");
    }
    for (std::size_t i = 0; i < 3; ++i) {
        logger.assert_always((p3[i] == 0.),
                             "3D default constructor or operator[] of Point");
    }
    for (std::size_t i = 0; i < 4; ++i) {
        logger.assert_always((p4[i] == 0.),
                             "4D default constructor or operator[] of Point");
    }

    logger.assert_always((pp1[0] == 0.), "1D copy constructor");
    for (std::size_t i = 0; i < 2; ++i) {
        logger.assert_always((pp2[i] == 0.), "2D copy constructor");
    }
    for (std::size_t i = 0; i < 3; ++i) {
        logger.assert_always((pp3[i] == 0.), "3D copy constructor");
    }
    for (std::size_t i = 0; i < 4; ++i) {
        logger.assert_always((pp4[i] == 0.), "4D copy constructor");
    }

    logger.assert_always((std::abs(pc1[0] - 1.1) < 1e-12),
                         "1D from array constructor or operator[] of Point");
    for (std::size_t i = 0; i < 2; ++i) {
        logger.assert_always(
            (std::abs(pc2[i] - 1.2 - i) < 1e-12),
            "2D from array constructor or operator[] of Point");
    }
    for (std::size_t i = 0; i < 3; ++i) {
        logger.assert_always(
            (std::abs(pc3[i] - 1.3 - i) < 1e-12),
            "3D from array constructor or operator[] of Point");
    }
    for (std::size_t i = 0; i < 4; ++i) {
        logger.assert_always(
            (std::abs(pc4[i] - 1.4 - i) < 1e-12),
            "4D from array constructor or operator[] of Point");
    }

    logger.assert_always((std::abs(pv1[0] - 1.1) < 1e-12),
                         "1D from NumericalVector constructor");
    for (std::size_t i = 0; i < 2; ++i) {
        logger.assert_always((std::abs(pv2[i] - 1.2 - i) < 1e-12),
                             "2D from NumericalVector constructor");
    }
    for (std::size_t i = 0; i < 3; ++i) {
        logger.assert_always((std::abs(pv3[i] - 1.3 - i) < 1e-12),
                             "3D from NumericalVector constructor");
    }
    for (std::size_t i = 0; i < 4; ++i) {
        logger.assert_always((std::abs(pv4[i] - 1.4 - i) < 1e-12),
                             "4D from NumericalVector constructor");
    }

    // testing setCoordinates and setCoordinate

    p1.setCoordinates(vec1);
    logger.assert_always((std::abs(p1[0] - 1.1) < 1e-12), "1D setCoordinates");
    p2.setCoordinates(vec2);
    for (std::size_t i = 0; i < 2; ++i) {
        logger.assert_always((std::abs(p2[i] - 1.2 - i) < 1e-12),
                             "2D setCoordinates");
    }
    p3.setCoordinates(vec3);
    for (std::size_t i = 0; i < 3; ++i) {
        logger.assert_always((std::abs(p3[i] - 1.3 - i) < 1e-12),
                             "3D setCoordinates");
    }
    p4.setCoordinates(vec4);
    for (std::size_t i = 0; i < 4; ++i) {
        logger.assert_always((std::abs(p4[i] - 1.4 - i) < 1e-12),
                             "4D setCoordinates");
    }

    p1.setCoordinate(0, 0.9);
    logger.assert_always((std::abs(p1[0] - 0.9) < 1e-12), "1D setCoordinate");
    for (std::size_t i = 0; i < 2; ++i) {
        p2.setCoordinate(i, 0.8 + double(i));
        for (std::size_t j = 0; j <= i; ++j) {
            logger.assert_always((std::abs(p2[j] - 0.8 - j) < 1e-12),
                                 "2D setCoordinate");
        }
        for (std::size_t j = i + 1; j < 2; ++j) {
            logger.assert_always((std::abs(p2[j] - 1.2 - j) < 1e-12),
                                 "2D setCoordinate");
        }
    }
    for (std::size_t i = 0; i < 3; ++i) {
        p3.setCoordinate(i, 0.7 + i);
        for (std::size_t j = 0; j <= i; ++j) {
            logger.assert_always((std::abs(p3[j] - 0.7 - j) < 1e-12),
                                 "3D setCoordinate");
        }
        for (std::size_t j = i + 1; j < 3; ++j) {
            logger.assert_always((std::abs(p3[j] - 1.3 - j) < 1e-12),
                                 "3D setCoordinate");
        }
    }
    for (std::size_t i = 0; i < 4; ++i) {
        p4.setCoordinate(i, 0.6 + i);
        for (std::size_t j = 0; j <= i; ++j) {
            logger.assert_always((std::abs(p4[j] - 0.6 - j) < 1e-12),
                                 "4D setCoordinate");
        }
        for (std::size_t j = i + 1; j < 4; ++j) {
            logger.assert_always((std::abs(p4[j] - 1.4 - j) < 1e-12),
                                 "4D setCoordinate");
        }
    }

    // testing operators

    const Point<0> pr0 = pc0 = p0;
    const Point<1> pr1 = pc1 = p1;
    logger.assert_always((std::abs(pc1[0] - 0.9) < 1e-12),
                         "1D assignment operator");
    logger.assert_always((std::abs(pr1[0] - 0.9) < 1e-12),
                         "1D assignment operator");
    const Point<2> pr2 = pc2 = p2;
    for (std::size_t i = 0; i < 2; ++i) {
        logger.assert_always((std::abs(pc2[i] - 0.8 - i) < 1e-12),
                             "2D assignment operator");
        logger.assert_always((std::abs(pr2[i] - 0.8 - i) < 1e-12),
                             "2D assignment operator");
    }
    const Point<3> pr3 = pc3 = p3;
    for (std::size_t i = 0; i < 3; ++i) {
        logger.assert_always((std::abs(pc3[i] - 0.7 - i) < 1e-12),
                             "3D assignment operator");
        logger.assert_always((std::abs(pr3[i] - 0.7 - i) < 1e-12),
                             "3D assignment operator");
    }
    const Point<4> pr4 = pc4 = p4;
    for (std::size_t i = 0; i < 4; ++i) {
        logger.assert_always((std::abs(pc4[i] - 0.6 - i) < 1e-12),
                             "4D assignment operator");
        logger.assert_always((std::abs(pr4[i] - 0.6 - i) < 1e-12),
                             "4D assignment operator");
    }

    logger.assert_always((pr0 == pc0 && pc0 == pr0 && pc0 == p0),
                         "0D equality operator");
    logger.assert_always((pr1 == pc1 && pc1 == pr1 && pc1 == p1 &&
                          !(pr1 == pv1 || pv1 == pr1 || p1 == pv1)),
                         "1D equality operator");
    logger.assert_always((pr2 == pc2 && pc2 == pr2 && pc2 == p2 &&
                          !(pr2 == pv2 || pv2 == pr2 || p2 == pv2)),
                         "2D equality operator");
    logger.assert_always((pr3 == pc3 && pc3 == pr3 && pc3 == p3 &&
                          !(pr3 == pv3 || pv3 == pr3 || p3 == pv3)),
                         "3D equality operator");
    logger.assert_always((pr4 == pc4 && pc4 == pr4 && pc4 == p4 &&
                          !(pr4 == pv4 || pv4 == pr4 || p4 == pv4)),
                         "4D equality operator");

    pc0 += p0;
    pc1 += p1;
    logger.assert_always((std::abs(pc1[0] - 1.8) < 1e-12),
                         "1D increment operator");
    pc2 += p2;
    for (std::size_t i = 0; i < 2; ++i) {
        logger.assert_always((std::abs(pc2[i] - 1.6 - 2 * i) < 1e-12),
                             "2D increment operator");
    }
    pc3 += p3;
    for (std::size_t i = 0; i < 3; ++i) {
        logger.assert_always((std::abs(pc3[i] - 1.4 - 2 * i) < 1e-12),
                             "3D increment operator");
    }
    pc4 += p4;
    for (std::size_t i = 0; i < 4; ++i) {
        logger.assert_always((std::abs(pc4[i] - 1.2 - 2 * i) < 1e-12),
                             "4D increment operator");
    }

    pc0 -= pv0;
    pc1 -= pv1;
    logger.assert_always((std::abs(pc1[0] - 0.7) < 1e-12),
                         "1D decrement operator");
    pc2 -= pv2;
    for (std::size_t i = 0; i < 2; ++i) {
        logger.assert_always((std::abs(pc2[i] - 0.4 - i) < 1e-12),
                             "2D decrement operator");
    }
    pc3 -= pv3;
    for (std::size_t i = 0; i < 3; ++i) {
        logger.assert_always((std::abs(pc3[i] - 0.1 - i) < 1e-12),
                             "3D decrement operator");
    }
    pc4 -= pv4;
    for (std::size_t i = 0; i < 4; ++i) {
        logger.assert_always((std::abs(pc4[i] + 0.2 - i) < 1e-12),
                             "4D decrement operator or negative numbers");
    }

    pc0 *= 2.;
    pc1 *= 3.;
    logger.assert_always((std::abs(pc1[0] - 2.1) < 1e-12),
                         "1D multiply operator");
    pc2 *= 4.;
    for (std::size_t i = 0; i < 2; ++i) {
        logger.assert_always((std::abs(pc2[i] - 1.6 - 4 * i) < 1e-12),
                             "2D multiply operator");
    }
    pc3 *= 5.;
    for (std::size_t i = 0; i < 3; ++i) {
        logger.assert_always((std::abs(pc3[i] - 0.5 - 5 * i) < 1e-12),
                             "3D multiply operator");
    }
    pc4 *= 6.;
    for (std::size_t i = 0; i < 4; ++i) {
        logger.assert_always((std::abs(pc4[i] + 1.2 - 6 * i) < 1e-12),
                             "4D multiply operator");
    }

    pv0 * 6.;
    logger.assert_always((std::abs((pv1 * 5.)[0] - 5.5) < 1e-12),
                         "1D multiplication");
    for (std::size_t i = 0; i < 2; ++i) {
        logger.assert_always((std::abs((pv2 * 4.)[i] - 4.8 - 4 * i) < 1e-12),
                             "2D multiplication");
    }
    for (std::size_t i = 0; i < 3; ++i) {
        logger.assert_always((std::abs((pv3 * 3.)[i] - 3.9 - 3 * i) < 1e-12),
                             "3D multiplication");
    }
    for (std::size_t i = 0; i < 4; ++i) {
        logger.assert_always((std::abs((pv4 * 2.)[i] - 2.8 - 2 * i) < 1e-12),
                             "4D multiplication");
    }

    logger.assert_always(((pr0 * 0.) == pp0), "0D multiplication");
    logger.assert_always(((pr1 * 0.) == pp1), "1D multiplication");
    logger.assert_always(((pr2 * 0.) == pp2), "2D multiplication");
    logger.assert_always(((pr3 * 0.) == pp3), "3D multiplication");
    logger.assert_always(((pr4 * 0.) == pp4), "4D multiplication");

    pc0 + pv0;
    logger.assert_always((std::abs((pc1 + pv1)[0] - 3.2) < 1e-12),
                         "1D addition");
    for (std::size_t i = 0; i < 2; ++i) {
        logger.assert_always((std::abs((pc2 + pv2)[i] - 2.8 - 5 * i) < 1e-12),
                             "2D addition");
    }
    for (std::size_t i = 0; i < 3; ++i) {
        logger.assert_always((std::abs((pc3 + pv3)[i] - 1.8 - 6 * i) < 1e-12),
                             "3D addition");
    }
    for (std::size_t i = 0; i < 4; ++i) {
        logger.assert_always((std::abs((pc4 + pv4)[i] - 0.2 - 7 * i) < 1e-12),
                             "4D addition");
    }

    pr0 + pv0;
    logger.assert_always((std::abs((pr1 + pv1)[0] - 2.) < 1e-12),
                         "1D addition");
    for (std::size_t i = 0; i < 2; ++i) {
        logger.assert_always((std::abs((pr2 + pv2)[i] - 2. - 2 * i) < 1e-12),
                             "2D addition");
    }
    for (std::size_t i = 0; i < 3; ++i) {
        logger.assert_always((std::abs((pr3 + pv3)[i] - 2. - 2 * i) < 1e-12),
                             "3D addition");
    }
    for (std::size_t i = 0; i < 4; ++i) {
        logger.assert_always((std::abs((pr4 + pv4)[i] - 2. - 2 * i) < 1e-12),
                             "4D addition");
    }

    pc0 - pv0;
    logger.assert_always((std::abs((pc1 - pv1)[0] - 1.) < 1e-12),
                         "1D subtraction");
    for (std::size_t i = 0; i < 2; ++i) {
        logger.assert_always((std::abs((pc2 - pv2)[i] - 0.4 - 3 * i) < 1e-12),
                             "2D subtraction");
    }
    for (std::size_t i = 0; i < 3; ++i) {
        logger.assert_always((std::abs((pc3 - pv3)[i] + 0.8 - 4 * i) < 1e-12),
                             "3D subtraction");
    }
    for (std::size_t i = 0; i < 4; ++i) {
        logger.assert_always((std::abs((pc4 - pv4)[i] + 2.6 - 5 * i) < 1e-12),
                             "4D subtraction");
    }

    pr0 - pv0;
    logger.assert_always((std::abs((pr1 - pv1)[0] + 0.2) < 1e-12),
                         "1D subtraction");
    for (std::size_t i = 0; i < 2; ++i) {
        logger.assert_always((std::abs((pr2 - pv2)[i] + 0.4) < 1e-12),
                             "2D subtraction");
    }
    for (std::size_t i = 0; i < 3; ++i) {
        logger.assert_always((std::abs((pr3 - pv3)[i] + 0.6) < 1e-12),
                             "3D subtraction");
    }
    for (std::size_t i = 0; i < 4; ++i) {
        logger.assert_always((std::abs((pr4 - pv4)[i] + 0.8) < 1e-12),
                             "4D subtraction");
    }

    // testing size

    logger.assert_always(
        (p0.size() == 0 && pp0.size() == 0 && pr0.size() == 0 &&
         pc0.size() == 0 && pv0.size() == 0),
        "size of a 0D point");
    logger.assert_always(
        (p1.size() == 1 && pp1.size() == 1 && pr1.size() == 1 &&
         pc1.size() == 1 && pv1.size() == 1),
        "size of a 1D point");
    logger.assert_always(
        (p2.size() == 2 && pp2.size() == 2 && pr2.size() == 2 &&
         pc2.size() == 2 && pv2.size() == 2),
        "size of a 2D point");
    logger.assert_always(
        (p3.size() == 3 && pp3.size() == 3 && pr3.size() == 3 &&
         pc3.size() == 3 && pv3.size() == 3),
        "size of a 3D point");
    logger.assert_always(
        (p4.size() == 4 && pp4.size() == 4 && pr4.size() == 4 &&
         pc4.size() == 4 && pv4.size() == 4),
        "size of a 4D point");

    // testing getCoordinate and getCoordinates

    logger.assert_always((p1.getCoordinate(0) == p1[0]), "1D getCoordinate");
    for (std::size_t i = 0; i < 2; ++i) {
        logger.assert_always((p2.getCoordinate(i) == p2[i]),
                             "2D getCoordinate");
    }
    for (std::size_t i = 0; i < 3; ++i) {
        logger.assert_always((p3.getCoordinate(i) == p3[i]),
                             "3D getCoordinate");
    }
    for (std::size_t i = 0; i < 4; ++i) {
        logger.assert_always((p4.getCoordinate(i) == p4[i]),
                             "4D getCoordinate");
    }

    logger.assert_always((pv0.getCoordinates() == vec0), "0D getCoordinates");
    logger.assert_always((pv1.getCoordinates() == vec1), "1D getCoordinates");
    logger.assert_always((pv2.getCoordinates() == vec2), "2D getCoordinates");
    logger.assert_always((pv3.getCoordinates() == vec3), "3D getCoordinates");
    logger.assert_always((pv4.getCoordinates() == vec4), "4D getCoordinates");

    // testing friends

    -pc0;
    logger.assert_always((std::abs((-pc1)[0] + 2.1) < 1e-12), "1D unary -");
    for (std::size_t i = 0; i < 2; ++i) {
        logger.assert_always((std::abs((-pc2)[i] + 1.6 + 4 * i) < 1e-12),
                             "2D unary -");
    }
    for (std::size_t i = 0; i < 3; ++i) {
        logger.assert_always((std::abs((-pc3)[i] + 0.5 + 5 * i) < 1e-12),
                             "3D unary -");
    }
    for (std::size_t i = 0; i < 4; ++i) {
        logger.assert_always((std::abs((-pc4)[i] - 1.2 + 6 * i) < 1e-12),
                             "4D unary -");
    }

    6. * pv0;
    logger.assert_always((std::abs((5. * pv1)[0] - 5.5) < 1e-12),
                         "1D left multiplication");
    for (std::size_t i = 0; i < 2; ++i) {
        logger.assert_always((std::abs((4. * pv2)[i] - 4.8 - 4 * i) < 1e-12),
                             "2D left multiplication");
    }
    for (std::size_t i = 0; i < 3; ++i) {
        logger.assert_always((std::abs((3. * pv3)[i] - 3.9 - 3 * i) < 1e-12),
                             "3D left multiplication");
    }
    for (std::size_t i = 0; i < 4; ++i) {
        logger.assert_always((std::abs((2. * pv4)[i] - 2.8 - 2 * i) < 1e-12),
                             "4D left multiplication");
    }

    std::cout << p0 << p1 << p2 << p3 << p4 << pr0 << pr1 << pr2 << pr3 << pr4
              << std::endl;

    return 0;
}
