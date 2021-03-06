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
#include "Geometry/Mappings/MappingToRefLineToTriangle.h"
#include "Logger.h"

#include "Geometry/ReferenceLine.h"
#include "Geometry/ReferenceTriangle.h"
#include "Geometry/PointReference.h"
#include "Geometry/Jacobian.h"
#include "LinearAlgebra/MiddleSizeVector.h"
#include <cmath>

#include "../catch.hpp"

using namespace hpgem;
TEST_CASE("120MappingToRefLineToTriangle_UnitTest",
          "[120MappingToRefLineToTriangle_UnitTest]") {

    Geometry::PointReference<1> refPoint;
    Geometry::PointReference<2> point, compare;

    Geometry::ReferenceTriangle& eGeom =
        Geometry::ReferenceTriangle::Instance();
    Geometry::ReferenceLine& fGeom = Geometry::ReferenceLine::Instance();

    Geometry::Jacobian<1, 2> jac;

    std::vector<std::size_t> nodesAfterTransformation(2);

    const Geometry::MappingReferenceToReference<1>* test =
        &Geometry::MappingToRefLineToTriangle0::Instance();
    nodesAfterTransformation[0] = 0;
    nodesAfterTransformation[1] = 1;

    for (refPoint[0] = -2.8189; refPoint[0] < 3.141; refPoint[0] += 0.1) {
        point = test->transform((refPoint));
        INFO("transform");
        CHECK((fGeom.isInternalPoint((refPoint)) ==
               eGeom.isInternalPoint((point))));

        refPoint[0] += -1.e-8;
        compare = test->transform((refPoint));
        refPoint[0] += 2.e-8;
        point = test->transform((refPoint));

        refPoint[0] += -1e-8;
        jac = test->calcJacobian((refPoint));
        INFO("jacobian");
        CHECK((std::abs(jac[0] - 5.e7 * (point[0] - compare[0])) <
               1e-5));  // estimate is a bit rough, but should work for most
                        // mappings
        INFO("jacobian");
        CHECK((std::abs(jac[1] - 5.e7 * (point[1] - compare[1])) <
               1e-5));  // implementations are very strongly recommended to be
                        // more accurate
    }

    for (std::size_t i = 0; i < fGeom.getNumberOfNodes(); ++i) {
        refPoint = fGeom.getReferenceNodeCoordinate(i);
        compare = eGeom.getReferenceNodeCoordinate(nodesAfterTransformation[i]);
        point = test->transform((refPoint));
        INFO("transform");
        CHECK((std::abs(point[0] - compare[0]) < 1e-12));
        INFO("transform");
        CHECK((std::abs(point[1] - compare[1]) < 1e-12));
    }

    INFO("getTargetDimension");
    CHECK((test->getTargetDimension() == 2));

    test = &Geometry::MappingToRefLineToTriangle1::Instance();
    nodesAfterTransformation[0] = 0;
    nodesAfterTransformation[1] = 2;

    for (refPoint[0] = -2.8189; refPoint[0] < 3.141; refPoint[0] += 0.1) {
        point = test->transform((refPoint));
        INFO("transform");
        CHECK((fGeom.isInternalPoint((refPoint)) ==
               eGeom.isInternalPoint((point))));

        refPoint[0] += -1.e-8;
        compare = test->transform((refPoint));
        refPoint[0] += 2.e-8;
        point = test->transform((refPoint));

        refPoint[0] += -1e-8;
        jac = test->calcJacobian((refPoint));
        INFO("jacobian");
        CHECK((std::abs(jac[0] - 5.e7 * (point[0] - compare[0])) <
               1e-5));  // estimate is a bit rough, but should work for most
                        // mappings
        INFO("jacobian");
        CHECK((std::abs(jac[1] - 5.e7 * (point[1] - compare[1])) <
               1e-5));  // implementations are very strongly recommended to be
                        // more accurate
    }

    for (std::size_t i = 0; i < fGeom.getNumberOfNodes(); ++i) {
        refPoint = fGeom.getReferenceNodeCoordinate(i);
        compare = eGeom.getReferenceNodeCoordinate(nodesAfterTransformation[i]);
        point = test->transform((refPoint));
        INFO("transform");
        CHECK((std::abs(point[0] - compare[0]) < 1e-12));
        INFO("transform");
        CHECK((std::abs(point[1] - compare[1]) < 1e-12));
    }

    INFO("getTargetDimension");
    CHECK((test->getTargetDimension() == 2));

    test = &Geometry::MappingToRefLineToTriangle2::Instance();
    nodesAfterTransformation[0] = 1;
    nodesAfterTransformation[1] = 2;

    for (refPoint[0] = -2.8189; refPoint[0] < 3.141; refPoint[0] += 0.1) {
        point = test->transform((refPoint));
        std::cout << refPoint << " " << point << std::endl;
        refPoint[0] += -1.e-8;
        compare = test->transform((refPoint));
        refPoint[0] += 2.e-8;
        point = test->transform((refPoint));

        refPoint[0] += -1e-8;
        jac = test->calcJacobian((refPoint));
        INFO("jacobian");
        CHECK((std::abs(jac[0] - 5.e7 * (point[0] - compare[0])) <
               1e-5));  // estimate is a bit rough, but should work for most
                        // mappings
        INFO("jacobian");
        CHECK((std::abs(jac[1] - 5.e7 * (point[1] - compare[1])) <
               1e-5));  // implementations are very strongly recommended to be
                        // more accurate
    }

    for (std::size_t i = 0; i < fGeom.getNumberOfNodes(); ++i) {
        refPoint = fGeom.getReferenceNodeCoordinate(i);
        compare = eGeom.getReferenceNodeCoordinate(nodesAfterTransformation[i]);
        point = test->transform((refPoint));
        INFO("transform");
        CHECK((std::abs(point[0] - compare[0]) < 1e-12));
        INFO("transform");
        CHECK((std::abs(point[1] - compare[1]) < 1e-12));
    }

    INFO("getTargetDimension");
    CHECK((test->getTargetDimension() == 2));
}
