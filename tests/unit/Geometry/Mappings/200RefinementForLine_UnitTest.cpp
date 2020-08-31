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

#include "Geometry/Mappings/RefinementMapsForLine.h"
#include "Geometry/PointReference.h"
#include "Geometry/ReferenceLine.h"
#include <iostream>

#include "../catch.hpp"

using namespace hpgem;
TEST_CASE("200RefinementForLine_UnitTest", "[200RefinementForLine_UnitTest]") {
    Geometry::PointReference<1> refPoint;
    Geometry::PointReference<1> point, compare;
    LinearAlgebra::SmallMatrix<1, 1> jac;
    for (const Geometry::RefinementMapping* test :
         std::initializer_list<const Geometry::RefinementMapping*>{
             Geometry::RefinementMapForLine0::instance(),
             Geometry::RefinementMapForLine1::instance()}) {
        std::cout << test->getName();
        for (std::size_t i = 0; i < test->getNumberOfSubElements(); ++i) {
            for (refPoint[0] = -2.8189; refPoint[0] < 3.141;
                 refPoint[0] += 0.1) {
                point = test->refinementTransform(i, (refPoint));
                refPoint[0] += -1.e-8;
                compare = test->refinementTransform(i, (refPoint));
                refPoint[0] += 2.e-8;
                point = test->refinementTransform(i, (refPoint));

                refPoint[0] += -1e-8;
                jac = test->getRefinementMappingMatrixL(i, (refPoint));
                INFO("jacobian");
                CHECK((std::abs(jac[0] - 5.e7 * (point[0] - compare[0])) <
                       1e-5));  // estimate is a bit rough, but should work
                                // for most mappings
                jac *= test->getRefinementMappingMatrixR(i, (refPoint));
                INFO("inverse of jacobian");
                CHECK(std::abs(jac[0] - 1.) < 1e-12);
            }
            for (std::size_t index : test->getSubElementLocalNodeIndices(i)) {
                INFO("local index out of bounds");
                CHECK(index <
                      Geometry::ReferenceLine::Instance().getNumberOfNodes() +
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
