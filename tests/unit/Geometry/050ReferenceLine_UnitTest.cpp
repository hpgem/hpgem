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
#include "Geometry/ReferenceLine.h"
#include "Geometry/ReferencePoint.h"
#include <iostream>

#include "Geometry/PointReference.h"
#include "Geometry/Mappings/MappingToRefLineToLine.h"
#include "Geometry/Mappings/MappingToRefPointToLine.h"
#include "Integration/QuadratureRules/GaussQuadratureRule.h"
#include "Logger.h"
#include <cmath>

#define CATCH_CONFIG_MAIN
#include "../catch.hpp"

using namespace hpgem;
using Geometry::ReferenceLine;

TEST_CASE("050ReferenceLine_UnitTest", "[050ReferenceLine_UnitTest]") {
    ReferenceLine& test = ReferenceLine::Instance();

    Geometry::PointReference<1> pTest;

    // testing basic functionality

    for (pTest[0] = -3.141; pTest[0] < -1.; pTest[0] += 0.1) {
        INFO("isInternalPoint");
        CHECK((!test.isInternalPoint(pTest)));
    }
    for (; pTest[0] < 1; pTest[0] += 0.1) {
        INFO("isInternalPoint");
        CHECK((test.isInternalPoint(pTest)));
    }
    for (; pTest[0] < 3.141; pTest[0] += 0.1) {
        INFO("isInternalPoint");
        CHECK((!test.isInternalPoint(pTest)));
    }

    pTest = test.getCenter();
    INFO("getCenter");
    CHECK((test.isInternalPoint(pTest) && std::abs(pTest[0]) < 1e-12));
    pTest = test.getReferenceNodeCoordinate(0);
    INFO("getNode 0");
    CHECK((std::abs(pTest[0] + 1) < 1e-12));
    pTest = test.getReferenceNodeCoordinate(1);
    INFO("getNode 1");
    CHECK((std::abs(pTest[0] - 1) < 1e-12));
    std::cout << test.getName();

    INFO("getLocalNodeIndex 0");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 0) == 0));
    INFO("getLocalNodeIndex 1");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 0) == 1));

    std::cout << test;

    // testing mappings and quadrature rules

    std::vector<std::size_t> base(2), transformed(2), faceIndices(1);
    for (std::size_t i = 0; i < 2; ++i) {
        base[i] = transformed[i] = i;
    }
    INFO("getCodim0MappingIndex&Ptr");
    CHECK((test.getCodim0MappingPtr(
               test.getCodim0MappingIndex(base, transformed)) ==
           &Geometry::MappingToRefLineToLine0::Instance()));
    INFO("getCodim0MappingIndex&Ptr");
    CHECK((test.getCodim0MappingPtr(base, transformed) ==
           &Geometry::MappingToRefLineToLine0::Instance()));
    transformed[0] = 1;
    transformed[1] = 0;
    INFO("getCodim0MappingIndex&Ptr");
    CHECK((test.getCodim0MappingPtr(
               test.getCodim0MappingIndex(base, transformed)) ==
           &Geometry::MappingToRefLineToLine1::Instance()));
    INFO("getCodim0MappingIndex&Ptr");
    CHECK((test.getCodim0MappingPtr(base, transformed) ==
           &Geometry::MappingToRefLineToLine1::Instance()));

    INFO("higher codimensional entities");
    CHECK((test.getNumberOfCodim1Entities() == 2 &&
           test.getNumberOfCodim2Entities() == 0) &&
          test.getNumberOfCodim3Entities() == 0);
    INFO("getCodim1ReferenceGeometry");
    CHECK((test.getCodim1ReferenceGeometry(0) ==
               &Geometry::ReferencePoint::Instance() &&
           test.getCodim1ReferenceGeometry(1) ==
               &Geometry::ReferencePoint::Instance()));
    INFO("getCodim1MappingPtr");
    CHECK((test.getCodim1MappingPtr(0) ==
           &Geometry::MappingToRefPointToLine0::Instance()));
    INFO("getCodim1MappingPtr");
    CHECK((test.getCodim1MappingPtr(1) ==
           &Geometry::MappingToRefPointToLine1::Instance()));
    faceIndices = test.getCodim1EntityLocalIndices(0);
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[0] == test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 0)));
    faceIndices = test.getCodim1EntityLocalIndices(1);
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[0] == test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 0)));

    INFO("quadrature rules");
    CHECK((test.getGaussQuadratureRule(3)->order() >= 3));
    INFO("quadrature rules");
    CHECK((test.getGaussQuadratureRule(5)->order() >= 5));
    INFO("quadrature rules");
    CHECK((test.getGaussQuadratureRule(7)->order() >= 7));
    INFO("quadrature rules");
    CHECK((test.getGaussQuadratureRule(9)->order() >= 9));
    INFO("quadrature rules");
    CHECK((test.getGaussQuadratureRule(11)->order() >= 11));

    // testing functionality of abstract parent classes

    INFO("number of nodes");
    CHECK((test.getNumberOfNodes() == 2));
    INFO("type of geometry");
    CHECK((test.getGeometryType() == Geometry::ReferenceGeometryType::LINE));

    ///\todo testing that the refinement maps behave exactly like the forwarded
    /// calls of this class
}
