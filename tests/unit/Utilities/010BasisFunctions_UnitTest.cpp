/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2014, University of Twente
 All rights reserved.

 Redistribution and use INFO($2);
CHECK($1); source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer.

 2. Redistributions INFO($2);
CHECK($1); binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer INFO($2);
CHECK($1); the documentation
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
// correct execution of all prior unit tests see also
// Base/???BasisFunction_UnitTest.cpp
// - running the same series of checks on different basisfunctions
#include "FE/BasisFunctions1DH1ConformingLine.h"
#include "FE/BasisFunctions2DH1ConformingSquare.h"
#include "FE/BasisFunctions2DH1ConformingTriangle.h"
#include "FE/BasisFunctions3DH1ConformingCube.h"
#include "FE/BasisFunctions3DH1ConformingTetrahedron.h"
#include "FE/BasisFunctions3DH1ConformingPrism.h"
#include "FE/BasisFunctions3DH1ConformingPyramid.h"
#include "FE/BasisFunctions3DNedelec.h"
#include "FE/BasisFunctions3DAinsworthCoyle.h"
#include "Logger.h"

#include "FE/BasisFunctionSet.h"
#include "Geometry/PointReference.h"
#include "LinearAlgebra/MiddleSizeVector.h"
#include "Geometry/ReferenceGeometry.h"
#include "Geometry/ReferenceLine.h"
#include <cmath>

#include "../catch.hpp"

using namespace hpgem;

TEST_CASE("BasisFunctions", "[BasisFunctions]") {

    // 1D
    FE::BasisFunctionSet* all1DbasisFunctions =
        FE::createDGBasisFunctionSet1DH1Line(5);
    Geometry::PointReference<1> point1D;
    LinearAlgebra::SmallVector<1> ret;
    for (auto test : *all1DbasisFunctions) {
        for (point1D[0] = -1.5; point1D[0] < 1.51; point1D[0] += 0.5) {
            point1D[0] += -1.e-8;
            double x0 = test->eval((point1D));
            point1D[0] += 2.e-8;
            double x1 = test->eval((point1D));

            point1D[0] += -1e-8;
            ret = test->evalDeriv((point1D));
            INFO("gradient");
            CHECK((std::abs(ret[0] - 5.e7 * (x1 - x0)) < 1e-5));
            INFO("derivative");
            CHECK((std::abs(test->evalDeriv0((point1D)) - 5.e7 * (x1 - x0)) <
                   1e-5));
        }
    }

    delete all1DbasisFunctions;

    // 2D
    FE::BasisFunctionSet* all2DbasisFunctions =
        FE::createDGBasisFunctionSet2DH1Square(5);
    Geometry::PointReference<2> point2D;
    LinearAlgebra::SmallVector<2> ret2;
    for (auto test : *all2DbasisFunctions) {
        for (point2D[0] = -1.5; point2D[0] < 1.51; point2D[0] += 0.8) {
            for (point2D[1] = -1.5; point2D[1] < 1.51; point2D[1] += 0.9) {
                point2D[0] += -1.e-8;
                double x0 = test->eval((point2D));
                point2D[0] += 2.e-8;
                double x1 = test->eval((point2D));

                point2D[0] += -1e-8;
                ret2 = test->evalDeriv((point2D));
                INFO("gradient");
                CHECK((std::abs(ret2[0] - 5.e7 * (x1 - x0)) < 1e-5));
                INFO("derivative");
                CHECK((std::abs(test->evalDeriv0((point2D)) -
                                5.e7 * (x1 - x0)) < 1e-5));

                point2D[1] += -1.e-8;
                x0 = test->eval((point2D));
                point2D[1] += 2.e-8;
                x1 = test->eval((point2D));

                point2D[1] += -1e-8;
                INFO("gradient");
                CHECK((std::abs(ret2[1] - 5.e7 * (x1 - x0)) < 1e-5));
                INFO("derivative");
                CHECK((std::abs(test->evalDeriv1((point2D)) -
                                5.e7 * (x1 - x0)) < 1e-5));
            }
        }
    }

    delete all2DbasisFunctions;

    all2DbasisFunctions = FE::createDGBasisFunctionSet2DH1Triangle(5);
    for (auto test : *all2DbasisFunctions) {
        for (point2D[0] = -1.5; point2D[0] < 1.51; point2D[0] += 0.8) {
            for (point2D[1] = -1.5; point2D[1] < 1.51; point2D[1] += 0.9) {
                point2D[0] += -1.e-8;
                double x0 = test->eval((point2D));
                point2D[0] += 2.e-8;
                double x1 = test->eval((point2D));

                point2D[0] += -1e-8;
                ret2 =
                    test->evalDeriv((point2D));  // exact to within absolute OR
                                                 // relative tolerace of 1e-5
                double derivative = test->evalDeriv0((point2D));
                INFO("gradient");
                CHECK((
                    std::abs(ret2[0] - 5.e7 * (x1 - x0)) < 1e-5 ||
                    (ret2.l2Norm() > 1 && std::abs(ret2[0] - 5.e7 * (x1 - x0)) <
                                              1e-5 * ret2.l2Norm())));
                INFO("derivative");
                CHECK((std::abs(derivative - 5.e7 * (x1 - x0)) < 1e-5 ||
                       (std::abs(derivative) > 1 &&
                        std::abs(derivative - 5.e7 * (x1 - x0)) <
                            1e-5 * std::abs(derivative))));

                point2D[1] += -1.e-8;
                x0 = test->eval((point2D));
                point2D[1] += 2.e-8;
                x1 = test->eval((point2D));

                point2D[1] += -1e-8;
                derivative = test->evalDeriv1((point2D));
                INFO("gradient");
                CHECK((
                    std::abs(ret2[1] - 5.e7 * (x1 - x0)) < 1e-5 ||
                    (ret2.l2Norm() > 1 && std::abs(ret2[1] - 5.e7 * (x1 - x0)) <
                                              1e-5 * ret2.l2Norm())));
                INFO("derivative");
                CHECK((std::abs(derivative - 5.e7 * (x1 - x0)) < 1e-5 ||
                       (std::abs(derivative) > 1 &&
                        std::abs(derivative - 5.e7 * (x1 - x0)) <
                            1e-5 * std::abs(derivative))));
            }
        }
    }

    delete all2DbasisFunctions;

    // 3D

    FE::BasisFunctionSet* all3DbasisFunctions =
        FE::createDGBasisFunctionSet3DH1Cube(5);
    Geometry::PointReference<3> point3D;
    LinearAlgebra::SmallVector<3> ret3;
    for (auto test : *all3DbasisFunctions) {
        for (point3D[0] = -1.5; point3D[0] < 1.51; point3D[0] += 0.6) {
            for (point3D[1] = -1.5; point3D[1] < 1.51; point3D[1] += 0.6) {
                for (point3D[2] = -1.5; point3D[2] < 1.51; point3D[2] += 1.2) {
                    point3D[0] += -1.e-8;
                    double x0 = test->eval((point3D));
                    point3D[0] += 2.e-8;
                    double x1 = test->eval((point3D));

                    point3D[0] += -1e-8;
                    ret3 = test->evalDeriv((point3D));
                    double derivative = test->evalDeriv0((point3D));
                    INFO("gradient");
                    CHECK((std::abs(ret3[0] - 5.e7 * (x1 - x0)) < 1e-5 ||
                           (ret3.l2Norm() > 1 &&
                            std::abs(ret3[0] - 5.e7 * (x1 - x0)) <
                                1e-5 * ret3.l2Norm())));
                    INFO("derivative");
                    CHECK((std::abs(derivative - 5.e7 * (x1 - x0)) < 1e-5 ||
                           (std::abs(derivative) > 1 &&
                            std::abs(derivative - 5.e7 * (x1 - x0)) <
                                1e-5 * std::abs(derivative))));

                    point3D[1] += -1.e-8;
                    x0 = test->eval((point3D));
                    point3D[1] += 2.e-8;
                    x1 = test->eval((point3D));

                    point3D[1] += -1e-8;
                    derivative = test->evalDeriv1((point3D));
                    INFO("gradient");
                    CHECK((std::abs(ret3[1] - 5.e7 * (x1 - x0)) < 1e-5 ||
                           (ret3.l2Norm() > 1 &&
                            std::abs(ret3[1] - 5.e7 * (x1 - x0)) <
                                1e-5 * ret3.l2Norm())));
                    INFO("derivative");
                    CHECK((std::abs(derivative - 5.e7 * (x1 - x0)) < 1e-5 ||
                           (std::abs(derivative) > 1 &&
                            std::abs(derivative - 5.e7 * (x1 - x0)) <
                                1e-5 * std::abs(derivative))));

                    point3D[2] += -1.e-8;
                    x0 = test->eval((point3D));
                    point3D[2] += 2.e-8;
                    x1 = test->eval((point3D));

                    point3D[2] += -1e-8;
                    derivative = test->evalDeriv2((point3D));
                    INFO("gradient");
                    CHECK((std::abs(ret3[2] - 5.e7 * (x1 - x0)) < 1e-5 ||
                           (ret3.l2Norm() > 1 &&
                            std::abs(ret3[2] - 5.e7 * (x1 - x0)) <
                                1e-5 * ret3.l2Norm())));
                    INFO("derivative");
                    CHECK((std::abs(derivative - 5.e7 * (x1 - x0)) < 1e-5 ||
                           (std::abs(derivative) > 1 &&
                            std::abs(derivative - 5.e7 * (x1 - x0)) <
                                1e-5 * std::abs(derivative))));
                }
            }
        }
    }

    delete all3DbasisFunctions;

    all3DbasisFunctions = FE::createDGBasisFunctionSet3DH1Tetrahedron(5);
    for (auto test : *all3DbasisFunctions) {
        for (point3D[0] = -1.5; point3D[0] < 1.51; point3D[0] += 0.6) {
            for (point3D[1] = -1.5; point3D[1] < 1.51; point3D[1] += 0.6) {
                for (point3D[2] = -1.5; point3D[2] < 1.51; point3D[2] += 1.2) {
                    point3D[0] += -1.e-8;
                    double x0 = test->eval((point3D));
                    point3D[0] += 2.e-8;
                    double x1 = test->eval((point3D));

                    point3D[0] += -1e-8;
                    ret3 = test->evalDeriv((point3D));
                    double derivative = test->evalDeriv0((point3D));
                    INFO("gradient");
                    CHECK((std::abs(ret3[0] - 5.e7 * (x1 - x0)) < 1e-5 ||
                           (ret3.l2Norm() > 1 &&
                            std::abs(ret3[0] - 5.e7 * (x1 - x0)) <
                                1e-5 * ret3.l2Norm())));
                    INFO("derivative");
                    CHECK((std::abs(derivative - 5.e7 * (x1 - x0)) < 1e-5 ||
                           (std::abs(derivative) > 1 &&
                            std::abs(derivative - 5.e7 * (x1 - x0)) <
                                1e-5 * std::abs(derivative))));

                    point3D[1] += -1.e-8;
                    x0 = test->eval((point3D));
                    point3D[1] += 2.e-8;
                    x1 = test->eval((point3D));

                    point3D[1] += -1e-8;
                    derivative = test->evalDeriv1((point3D));
                    INFO("gradient");
                    CHECK((std::abs(ret3[1] - 5.e7 * (x1 - x0)) < 1e-5 ||
                           (ret3.l2Norm() > 1 &&
                            std::abs(ret3[1] - 5.e7 * (x1 - x0)) <
                                1e-5 * ret3.l2Norm())));
                    INFO("derivative");
                    CHECK((std::abs(derivative - 5.e7 * (x1 - x0)) < 1e-5 ||
                           (std::abs(derivative) > 1 &&
                            std::abs(derivative - 5.e7 * (x1 - x0)) <
                                1e-5 * std::abs(derivative))));

                    point3D[2] += -1.e-8;
                    x0 = test->eval((point3D));
                    point3D[2] += 2.e-8;
                    x1 = test->eval((point3D));

                    point3D[2] += -1e-8;
                    derivative = test->evalDeriv2((point3D));
                    INFO("gradient");
                    CHECK((std::abs(ret3[2] - 5.e7 * (x1 - x0)) < 1e-5 ||
                           (ret3.l2Norm() > 1 &&
                            std::abs(ret3[2] - 5.e7 * (x1 - x0)) <
                                1e-5 * ret3.l2Norm())));
                    INFO("derivative");
                    CHECK((std::abs(derivative - 5.e7 * (x1 - x0)) < 1e-5 ||
                           (std::abs(derivative) > 1 &&
                            std::abs(derivative - 5.e7 * (x1 - x0)) <
                                1e-5 * std::abs(derivative))));
                }
            }
        }
    }

    delete all3DbasisFunctions;

    all3DbasisFunctions = FE::createDGBasisFunctionSet3DNedelec(5);
    for (auto test : *all3DbasisFunctions) {
        for (point3D[0] = -1.5; point3D[0] < 1.51; point3D[0] += 0.6) {
            for (point3D[1] = -1.5; point3D[1] < 1.51; point3D[1] += 0.6) {
                for (point3D[2] = -1.5; point3D[2] < 1.51; point3D[2] += 1.2) {
                    point3D[0] += -1.e-8;
                    LinearAlgebra::SmallVector<3> x0;
                    test->eval((point3D), x0);
                    point3D[0] += 2.e-8;
                    LinearAlgebra::SmallVector<3> x1;
                    test->eval((point3D), x1);
                    point3D[0] += -1e-8;

                    point3D[1] += -1.e-8;
                    LinearAlgebra::SmallVector<3> y0;
                    test->eval((point3D), y0);
                    point3D[1] += 2.e-8;
                    LinearAlgebra::SmallVector<3> y1;
                    test->eval((point3D), y1);
                    point3D[1] += -1e-8;

                    point3D[2] += -1.e-8;
                    LinearAlgebra::SmallVector<3> z0;
                    test->eval((point3D), z0);
                    point3D[2] += 2.e-8;
                    LinearAlgebra::SmallVector<3> z1;
                    test->eval((point3D), z1);
                    point3D[2] += -1e-8;

                    ret3 = test->evalCurl((point3D));
                    INFO("curl");
                    CHECK((std::abs(ret3[0] - 5.e7 * (y1[2] - y0[2]) +
                                    5.e7 * (z1[1] - z0[1])) +
                                   std::abs(ret3[1] + 5.e7 * (x1[2] - x0[2]) -
                                            5.e7 * (z1[0] - z0[0])) +
                                   std::abs(ret3[2] - 5.e7 * (x1[1] - x0[1]) +
                                            5.e7 * (y1[0] - y0[0])) <
                               1e-5 ||
                           (ret3.l2Norm() &&
                            std::abs(ret3[0] - 5.e7 * (y1[2] - y0[2]) +
                                     5.e7 * (z1[1] - z0[1])) +
                                    std::abs(ret3[1] + 5.e7 * (x1[2] - x0[2]) -
                                             5.e7 * (z1[0] - z0[0])) +
                                    std::abs(ret3[2] - 5.e7 * (x1[1] - x0[1]) +
                                             5.e7 * (y1[0] - y0[0])) <
                                1e-5 * ret3.l2Norm())));
                }
            }
        }
    }

    delete all3DbasisFunctions;

    all3DbasisFunctions = FE::createDGBasisFunctionSet3DAinsworthCoyle(5);
    for (auto test : *all3DbasisFunctions) {
        for (point3D[0] = -1.5; point3D[0] < 1.51; point3D[0] += 0.8) {
            for (point3D[1] = -1.5; point3D[1] < 1.51; point3D[1] += 1.6) {
                for (point3D[2] = -1.5; point3D[2] < 1.51; point3D[2] += 2.) {
                    point3D[0] += -1.e-8;
                    LinearAlgebra::SmallVector<3> x0;
                    test->eval((point3D), x0);
                    point3D[0] += 2.e-8;
                    LinearAlgebra::SmallVector<3> x1;
                    test->eval((point3D), x1);
                    point3D[0] += -1e-8;

                    point3D[1] += -1.e-8;
                    LinearAlgebra::SmallVector<3> y0;
                    test->eval((point3D), y0);
                    point3D[1] += 2.e-8;
                    LinearAlgebra::SmallVector<3> y1;
                    test->eval((point3D), y1);
                    point3D[1] += -1e-8;

                    point3D[2] += -1.e-8;
                    LinearAlgebra::SmallVector<3> z0;
                    test->eval((point3D), z0);
                    point3D[2] += 2.e-8;
                    LinearAlgebra::SmallVector<3> z1;
                    test->eval((point3D), z1);
                    point3D[2] += -1e-8;

                    ret3 = test->evalCurl((point3D));
                    INFO("derivative");
                    CHECK((std::abs(ret3[0] - 5.e7 * (y1[2] - y0[2]) +
                                    5.e7 * (z1[1] - z0[1])) +
                                   std::abs(ret3[1] + 5.e7 * (x1[2] - x0[2]) -
                                            5.e7 * (z1[0] - z0[0])) +
                                   std::abs(ret3[2] - 5.e7 * (x1[1] - x0[1]) +
                                            5.e7 * (y1[0] - y0[0])) <
                               1e-5 ||
                           (ret3.l2Norm() &&
                            std::abs(ret3[0] - 5.e7 * (y1[2] - y0[2]) +
                                     5.e7 * (z1[1] - z0[1])) +
                                    std::abs(ret3[1] + 5.e7 * (x1[2] - x0[2]) -
                                             5.e7 * (z1[0] - z0[0])) +
                                    std::abs(ret3[2] - 5.e7 * (x1[1] - x0[1]) +
                                             5.e7 * (y1[0] - y0[0])) <
                                1e-5 * ret3.l2Norm())));
                }
            }
        }
    }

    delete all3DbasisFunctions;

    all3DbasisFunctions = FE::createDGBasisFunctionSet3DH1ConformingPrism(5);
    for (auto test : *all3DbasisFunctions) {
        for (point3D[0] = -1.5; point3D[0] < 1.51; point3D[0] += 0.6) {
            for (point3D[1] = -1.5; point3D[1] < 1.51; point3D[1] += 0.6) {
                for (point3D[2] = -1.5; point3D[2] < 1.51; point3D[2] += 1.2) {
                    point3D[0] += -1.e-8;
                    double x0 = test->eval((point3D));
                    point3D[0] += 2.e-8;
                    double x1 = test->eval((point3D));

                    point3D[0] += -1e-8;
                    ret3 = test->evalDeriv((point3D));
                    double derivative = test->evalDeriv0((point3D));
                    INFO("gradient");
                    CHECK((std::abs(ret3[0] - 5.e7 * (x1 - x0)) < 1e-5 ||
                           (ret3.l2Norm() > 1 &&
                            std::abs(ret3[0] - 5.e7 * (x1 - x0)) <
                                1e-5 * ret3.l2Norm())));
                    INFO("derivative");
                    CHECK((std::abs(derivative - 5.e7 * (x1 - x0)) < 1e-5 ||
                           (std::abs(derivative) > 1 &&
                            std::abs(derivative - 5.e7 * (x1 - x0)) <
                                1e-5 * std::abs(derivative))));

                    point3D[1] += -1.e-8;
                    x0 = test->eval((point3D));
                    point3D[1] += 2.e-8;
                    x1 = test->eval((point3D));

                    point3D[1] += -1e-8;
                    derivative = test->evalDeriv1((point3D));
                    INFO("gradient");
                    CHECK((std::abs(ret3[1] - 5.e7 * (x1 - x0)) < 1e-5 ||
                           (ret3.l2Norm() > 1 &&
                            std::abs(ret3[1] - 5.e7 * (x1 - x0)) <
                                1e-5 * ret3.l2Norm())));
                    INFO("derivative");
                    CHECK((std::abs(derivative - 5.e7 * (x1 - x0)) < 1e-5 ||
                           (std::abs(derivative) > 1 &&
                            std::abs(derivative - 5.e7 * (x1 - x0)) <
                                1e-5 * std::abs(derivative))));

                    point3D[2] += -1.e-8;
                    x0 = test->eval((point3D));
                    point3D[2] += 2.e-8;
                    x1 = test->eval((point3D));

                    point3D[2] += -1e-8;
                    derivative = test->evalDeriv2((point3D));
                    INFO("gradient");
                    CHECK((std::abs(ret3[2] - 5.e7 * (x1 - x0)) < 1e-5 ||
                           (ret3.l2Norm() > 1 &&
                            std::abs(ret3[2] - 5.e7 * (x1 - x0)) <
                                1e-5 * ret3.l2Norm())));
                    INFO("derivative");
                    CHECK((std::abs(derivative - 5.e7 * (x1 - x0)) < 1e-5 ||
                           (std::abs(derivative) > 1 &&
                            std::abs(derivative - 5.e7 * (x1 - x0)) <
                                1e-5 * std::abs(derivative))));
                }
            }
        }
    }

    delete all3DbasisFunctions;

    all3DbasisFunctions = FE::createDGBasisFunctionSet3DH1ConformingPyramid(1);
    for (auto test : *all3DbasisFunctions) {
        for (point3D[0] = -1.5; point3D[0] < 1.51; point3D[0] += 0.6) {
            for (point3D[1] = -1.5; point3D[1] < 1.51; point3D[1] += 0.6) {
                for (point3D[2] = -1.5; point3D[2] < 1.51; point3D[2] += 1.2) {
                    point3D[0] += -1.e-8;
                    double x0 = test->eval((point3D));
                    point3D[0] += 2.e-8;
                    double x1 = test->eval((point3D));

                    point3D[0] += -1e-8;
                    ret3 = test->evalDeriv((point3D));
                    double derivative = test->evalDeriv0((point3D));
                    INFO("gradient");
                    CHECK((std::abs(ret3[0] - 5.e7 * (x1 - x0)) < 1e-5 ||
                           (ret3.l2Norm() > 1 &&
                            std::abs(ret3[0] - 5.e7 * (x1 - x0)) <
                                1e-5 * ret3.l2Norm())));
                    INFO("derivative");
                    CHECK((std::abs(derivative - 5.e7 * (x1 - x0)) < 1e-5 ||
                           (std::abs(derivative) > 1 &&
                            std::abs(derivative - 5.e7 * (x1 - x0)) <
                                1e-5 * std::abs(derivative))));

                    point3D[1] += -1.e-8;
                    x0 = test->eval((point3D));
                    point3D[1] += 2.e-8;
                    x1 = test->eval((point3D));

                    point3D[1] += -1e-8;
                    derivative = test->evalDeriv1((point3D));
                    INFO("gradient");
                    CHECK((std::abs(ret3[1] - 5.e7 * (x1 - x0)) < 1e-5 ||
                           (ret3.l2Norm() > 1 &&
                            std::abs(ret3[1] - 5.e7 * (x1 - x0)) <
                                1e-5 * ret3.l2Norm())));
                    INFO("derivative");
                    CHECK((std::abs(derivative - 5.e7 * (x1 - x0)) < 1e-5 ||
                           (std::abs(derivative) > 1 &&
                            std::abs(derivative - 5.e7 * (x1 - x0)) <
                                1e-5 * std::abs(derivative))));

                    point3D[2] += -1.e-8;
                    x0 = test->eval((point3D));
                    point3D[2] += 2.e-8;
                    x1 = test->eval((point3D));

                    point3D[2] += -1e-8;
                    derivative = test->evalDeriv2((point3D));
                    INFO("gradient");
                    CHECK((std::abs(ret3[2] - 5.e7 * (x1 - x0)) < 1e-5 ||
                           (ret3.l2Norm() > 1 &&
                            std::abs(ret3[2] - 5.e7 * (x1 - x0)) <
                                1e-5 * ret3.l2Norm())));
                    INFO("derivative");
                    CHECK((std::abs(derivative - 5.e7 * (x1 - x0)) < 1e-5 ||
                           (std::abs(derivative) > 1 &&
                            std::abs(derivative - 5.e7 * (x1 - x0)) <
                                1e-5 * std::abs(derivative))));
                }
            }
        }
    }

    delete all3DbasisFunctions;
}
