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

#include "Geometry/Mappings/RefinementMapsForPyramid.h"
#include "Geometry/PointReference.h"
#include "Geometry/ReferencePyramid.h"
#include <iostream>

#include "../catch.hpp"

using namespace hpgem;
TEST_CASE("240RefinementForPyramid_UnitTest",
          "[240RefinementForPyramid_UnitTest]") {
    Geometry::PointReference<3> refPoint;
    Geometry::PointReference<3> point, compare;
    LinearAlgebra::SmallMatrix<3, 3> jac;
    for (const Geometry::RefinementMapping* test :
         {Geometry::RefinementMapForPyramid0::instance()}) {
        std::cout << test->getName();
        for (std::size_t i = 0; i < test->getNumberOfSubElements(); ++i) {
            for (refPoint[0] = -1.51; refPoint[0] < 1.51; refPoint[0] += 0.6) {
                for (refPoint[1] = -1.51; refPoint[1] < 1.51;
                     refPoint[1] += 0.7) {
                    for (refPoint[2] = -1.51; refPoint[2] < 1.51;
                         refPoint[2] += 0.8) {
                        point = test->refinementTransform(i, (refPoint));

                        refPoint[0] += -1.e-8;
                        compare = test->refinementTransform(i, (refPoint));
                        refPoint[0] += 2.e-8;
                        point = test->refinementTransform(i, (refPoint));

                        refPoint[0] += -1e-8;
                        jac = test->getRefinementMappingMatrixL(i, (refPoint));
                        INFO("jacobian");
                        CHECK(
                            (std::abs(jac[0] - 5.e7 * (point[0] - compare[0])) <
                             1e-5));  // estimate is a bit rough, but should
                                      // work for most mappings
                        INFO("jacobian");
                        CHECK(
                            (std::abs(jac[1] - 5.e7 * (point[1] - compare[1])) <
                             1e-5));  // implementations are very strongly
                                      // recommended to be more accurate
                        INFO("jacobian");
                        CHECK((std::abs(jac[2] - 5.e7 * (point[2] -
                                                         compare[2])) < 1e-5));

                        refPoint[1] += -1.e-8;
                        compare = test->refinementTransform(i, (refPoint));
                        refPoint[1] += 2.e-8;
                        point = test->refinementTransform(i, (refPoint));

                        refPoint[1] += -1e-8;
                        jac = test->getRefinementMappingMatrixL(i, (refPoint));
                        INFO("jacobian");
                        CHECK((std::abs(jac[3] - 5.e7 * (point[0] -
                                                         compare[0])) < 1e-5));
                        INFO("jacobian");
                        CHECK((std::abs(jac[4] - 5.e7 * (point[1] -
                                                         compare[1])) < 1e-5));
                        INFO("jacobian");
                        CHECK((std::abs(jac[5] - 5.e7 * (point[2] -
                                                         compare[2])) < 1e-5));

                        refPoint[2] += -1.e-8;
                        compare = test->refinementTransform(i, (refPoint));
                        refPoint[2] += 2.e-8;
                        point = test->refinementTransform(i, (refPoint));

                        refPoint[2] += -1e-8;
                        jac = test->getRefinementMappingMatrixL(i, (refPoint));
                        INFO("jacobian");
                        CHECK((std::abs(jac[6] - 5.e7 * (point[0] -
                                                         compare[0])) < 1e-5));
                        INFO("jacobian");
                        CHECK((std::abs(jac[7] - 5.e7 * (point[1] -
                                                         compare[1])) < 1e-5));
                        INFO("jacobian");
                        CHECK((std::abs(jac[8] - 5.e7 * (point[2] -
                                                         compare[2])) < 1e-5));
                        jac *= test->getRefinementMappingMatrixR(i, (refPoint));
                        INFO("inverse of jacobian");
                        CHECK(std::abs(jac[0] - 1.) < 1e-12);
                        INFO("inverse of jacobian");
                        CHECK(std::abs(jac[1] - 0.) < 1e-12);
                        INFO("inverse of jacobian");
                        CHECK(std::abs(jac[2] - 0.) < 1e-12);
                        INFO("inverse of jacobian");
                        CHECK(std::abs(jac[3] - 0.) < 1e-12);
                        INFO("inverse of jacobian");
                        CHECK(std::abs(jac[4] - 1.) < 1e-12);
                        INFO("inverse of jacobian");
                        CHECK(std::abs(jac[5] - 0.) < 1e-12);
                        INFO("inverse of jacobian");
                        CHECK(std::abs(jac[6] - 0.) < 1e-12);
                        INFO("inverse of jacobian");
                        CHECK(std::abs(jac[7] - 0.) < 1e-12);
                        INFO("inverse of jacobian");
                        CHECK(std::abs(jac[8] - 1.) < 1e-12);
                    }
                }
            }
            for (std::size_t index : test->getSubElementLocalNodeIndices(i)) {
                INFO("local index out of bounds");
                CHECK(
                    index <
                    Geometry::ReferencePyramid::Instance().getNumberOfNodes() +
                        test->getNumberOfNewNodes());
            }
            for (const Geometry::RefinementMapping* faceMapping :
                 test->getCodim1RefinementMaps()) {
                INFO("face gets too much new nodes");
                CHECK(faceMapping->getNumberOfNewNodes() <=
                      test->getNumberOfNewNodes());
                INFO("face gets too much new elements");
                CHECK(faceMapping->getNumberOfSubElements() <=
                      test->getNumberOfSubElements());
            }
        }
    }
}
