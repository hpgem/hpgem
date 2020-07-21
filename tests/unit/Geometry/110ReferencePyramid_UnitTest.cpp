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
#include "Geometry/ReferencePyramid.h"
#include "Geometry/ReferenceSquare.h"
#include "Geometry/ReferenceLine.h"
#include "Geometry/ReferenceTriangle.h"
#include "Geometry/ReferencePoint.h"
#include <iostream>
#include "Logger.h"

#include "Geometry/PointReference.h"
#include "Geometry/Mappings/MappingToRefFaceToPyramid.h"
#include "Integration/QuadratureRules/GaussQuadratureRule.h"
#include "Integration/QuadratureRules/GaussQuadratureRule.h"
#include <cmath>

#define CATCH_CONFIG_MAIN
#include "../catch.hpp"

using namespace hpgem;
using Geometry::ReferencePyramid;

TEST_CASE("110ReferencePyramid_UnitTest", "[110ReferencePyramid_UnitTest]") {
    ReferencePyramid& test = ReferencePyramid::Instance();

    Geometry::PointReference<3> pTest;

    // testing basic functionality

    for (pTest[0] = -1.51; pTest[0] < -1.; pTest[0] += 0.2) {
        for (pTest[1] = -1.511; pTest[1] < 1.51; pTest[1] += 0.2) {
            for (pTest[2] = -1.512; pTest[2] < 1.51; pTest[2] += 0.2) {
                INFO("isInternalPoint");
                CHECK((!test.isInternalPoint((pTest))));
            }
        }
    }
    for (; pTest[0] < 1; pTest[0] += 0.2) {
        for (pTest[1] = -1.511; pTest[1] < -1.; pTest[1] += 0.2) {
            for (pTest[2] = -1.512; pTest[2] < 1.51; pTest[2] += 0.2) {
                INFO("isInternalPoint");
                CHECK((!test.isInternalPoint((pTest))));
            }
        }
        for (; pTest[1] < 1.; pTest[1] += 0.2) {
            for (pTest[2] = -1.512; pTest[2] < 0.; pTest[2] += 0.2) {
                INFO("isInternalPoint");
                CHECK((!test.isInternalPoint((pTest))));
            }
            for (; pTest[2] < 1 - pTest[0] && pTest[2] < 1 + pTest[0] &&
                   pTest[2] < 1 - pTest[1] && pTest[2] < 1 + pTest[1];
                 pTest[2] += 0.2) {
                INFO("isInternalPoint");
                CHECK((test.isInternalPoint((pTest))));
            }
            for (; pTest[2] < 1.512; pTest[2] += 0.2) {
                INFO("isInternalPoint");
                CHECK((!test.isInternalPoint((pTest))));
            }
        }
        for (; pTest[1] < 1.51; pTest[1] += 0.2) {
            for (pTest[2] = -1.512; pTest[2] < 1.51; pTest[2] += 0.2) {
                INFO("isInternalPoint");
                CHECK((!test.isInternalPoint((pTest))));
            }
        }
    }
    for (; pTest[0] < 1.51; pTest[0] += 0.2) {
        for (pTest[1] = -1.511; pTest[1] < 1.51; pTest[1] += 0.2) {
            for (pTest[2] = -1.512; pTest[2] < 1.51; pTest[2] += 0.2) {
                INFO("isInternalPoint");
                CHECK((!test.isInternalPoint((pTest))));
            }
        }
    }

    pTest = test.getCenter();  ///\BUG it is not very clear to me where the
                               /// center of a pyramid lies
    INFO("getCenter");
    CHECK(test.isInternalPoint((pTest)));
    CHECK(std::abs(pTest[0]) < 1e-12);
    CHECK(std::abs(pTest[1]) < 1e-12);
    CHECK(std::abs(pTest[2] - 1. / 4.) < 1e-12);
    pTest = test.getReferenceNodeCoordinate(0);
    INFO("getNode 0");
    CHECK(std::abs(pTest[0]) < 1e-12);
    CHECK(std::abs(pTest[1]) < 1e-12);
    CHECK(std::abs(pTest[2] - 1) < 1e-12);
    pTest = test.getReferenceNodeCoordinate(1);
    INFO("getNode 1");
    CHECK(std::abs(pTest[0] + 1) < 1e-12);
    CHECK(std::abs(pTest[1] + 1) < 1e-12);
    CHECK(std::abs(pTest[2]) < 1e-12);
    pTest = test.getReferenceNodeCoordinate(2);
    INFO("getNode 2");
    CHECK(std::abs(pTest[0] - 1) < 1e-12);
    CHECK(std::abs(pTest[1] + 1) < 1e-12);
    CHECK(std::abs(pTest[2]) < 1e-12);
    pTest = test.getReferenceNodeCoordinate(3);
    INFO("getNode 3");
    CHECK(std::abs(pTest[0] + 1) < 1e-12);
    CHECK(std::abs(pTest[1] - 1) < 1e-12);
    CHECK(std::abs(pTest[2]) < 1e-12);
    pTest = test.getReferenceNodeCoordinate(4);
    INFO("getNode 4");
    CHECK(std::abs(pTest[0] - 1) < 1e-12);
    CHECK(std::abs(pTest[1] - 1) < 1e-12);
    CHECK(std::abs(pTest[2]) < 1e-12);

    INFO("getLocalNodeIndex 0");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 0) ==
           3));  // specified IN THIS SPECIFIC ORDER
    INFO("getLocalNodeIndex 0");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 1) ==
           4));  // should at least verify
    INFO("getLocalNodeIndex 0");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 2) ==
           1));  // specified
                 // twice
    INFO("getLocalNodeIndex 0");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 3) ==
           2));  // ordering of the nodes is consistent
    INFO("getLocalNodeIndex 1");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 0) == 3));
    INFO("getLocalNodeIndex 1");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 1) == 1));
    INFO("getLocalNodeIndex 1");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 2) == 0));
    INFO("getLocalNodeIndex 2");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 0) == 2));
    INFO("getLocalNodeIndex 2");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 1) == 4));
    INFO("getLocalNodeIndex 2");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 2) == 0));
    INFO("getLocalNodeIndex 3");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 0) == 1));
    INFO("getLocalNodeIndex 3");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 1) == 2));
    INFO("getLocalNodeIndex 3");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 2) == 0));
    INFO("getLocalNodeIndex 4");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 0) == 4));
    INFO("getLocalNodeIndex 4");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 1) == 3));
    INFO("getLocalNodeIndex 4");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 2) == 0));
    std::cout << test;

    // testing mappings and quadrature rules

    std::vector<std::size_t> faceIndices(4);
    // codim0maps dont exist so they dont need to be found properly

    INFO("higher codimensional entities");
    CHECK(test.getNumberOfCodim1Entities() == 5);
    CHECK(test.getNumberOfCodim2Entities() == 8);
    CHECK(test.getNumberOfCodim3Entities() == 5);
    INFO("getCodim1ReferenceGeometry");
    CHECK(test.getCodim1ReferenceGeometry(0) ==
          &Geometry::ReferenceSquare::Instance());
    CHECK(test.getCodim1ReferenceGeometry(1) ==
          &Geometry::ReferenceTriangle::Instance());
    CHECK(test.getCodim1ReferenceGeometry(2) ==
          &Geometry::ReferenceTriangle::Instance());
    CHECK(test.getCodim1ReferenceGeometry(3) ==
          &Geometry::ReferenceTriangle::Instance());
    CHECK(test.getCodim1ReferenceGeometry(4) ==
          &Geometry::ReferenceTriangle::Instance());
    INFO("getCodim2ReferenceGeometry");
    CHECK(test.getCodim2ReferenceGeometry(0) ==
          &Geometry::ReferenceLine::Instance());
    CHECK(test.getCodim2ReferenceGeometry(1) ==
          &Geometry::ReferenceLine::Instance());
    CHECK(test.getCodim2ReferenceGeometry(2) ==
          &Geometry::ReferenceLine::Instance());
    CHECK(test.getCodim2ReferenceGeometry(3) ==
          &Geometry::ReferenceLine::Instance());
    CHECK(test.getCodim2ReferenceGeometry(4) ==
          &Geometry::ReferenceLine::Instance());
    CHECK(test.getCodim2ReferenceGeometry(5) ==
          &Geometry::ReferenceLine::Instance());
    CHECK(test.getCodim2ReferenceGeometry(6) ==
          &Geometry::ReferenceLine::Instance());
    CHECK(test.getCodim2ReferenceGeometry(7) ==
          &Geometry::ReferenceLine::Instance());
    INFO("getCodim1MappingPtr");
    CHECK((test.getCodim1MappingPtr(0) ==
           &Geometry::MappingToRefFaceToPyramid0::Instance()));
    INFO("getCodim1MappingPtr");
    CHECK((test.getCodim1MappingPtr(1) ==
           &Geometry::MappingToRefFaceToPyramid1::Instance()));
    INFO("getCodim1MappingPtr");
    CHECK((test.getCodim1MappingPtr(2) ==
           &Geometry::MappingToRefFaceToPyramid2::Instance()));
    INFO("getCodim1MappingPtr");
    CHECK((test.getCodim1MappingPtr(3) ==
           &Geometry::MappingToRefFaceToPyramid3::Instance()));
    INFO("getCodim1MappingPtr");
    CHECK((test.getCodim1MappingPtr(4) ==
           &Geometry::MappingToRefFaceToPyramid4::Instance()));
    faceIndices = test.getCodim1EntityLocalIndices(0);
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[0] == test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 0)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[1] == test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 1)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[2] == test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 2)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[3] == test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 3)));
    faceIndices.resize(3);
    faceIndices = test.getCodim1EntityLocalIndices(1);
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[0] == test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 0)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[1] == test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 1)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[2] == test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 2)));
    faceIndices = test.getCodim1EntityLocalIndices(2);
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[0] == test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 0)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[1] == test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 1)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[2] == test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 2)));
    faceIndices = test.getCodim1EntityLocalIndices(3);
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[0] == test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 0)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[1] == test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 1)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[2] == test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 2)));
    faceIndices = test.getCodim1EntityLocalIndices(4);
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[0] == test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 0)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[1] == test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 1)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[2] == test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 2)));
    faceIndices.resize(2);
    faceIndices = test.getCodim2EntityLocalIndices(0);
    INFO("getCodim2EntityLocalIndices");
    CHECK((faceIndices[0] == 0));
    INFO("getCodim2EntityLocalIndices");
    CHECK((faceIndices[1] == 1));
    faceIndices = test.getCodim2EntityLocalIndices(1);
    INFO("getCodim2EntityLocalIndices");
    CHECK((faceIndices[0] == 0));
    INFO("getCodim2EntityLocalIndices");
    CHECK((faceIndices[1] == 2));
    faceIndices = test.getCodim2EntityLocalIndices(2);
    INFO("getCodim2EntityLocalIndices");
    CHECK((faceIndices[0] == 0));
    INFO("getCodim2EntityLocalIndices");
    CHECK((faceIndices[1] == 3));
    faceIndices = test.getCodim2EntityLocalIndices(3);
    INFO("getCodim2EntityLocalIndices");
    CHECK((faceIndices[0] == 0));
    INFO("getCodim2EntityLocalIndices");
    CHECK((faceIndices[1] == 4));
    faceIndices = test.getCodim2EntityLocalIndices(4);
    INFO("getCodim2EntityLocalIndices");
    CHECK((faceIndices[0] == 1));
    INFO("getCodim2EntityLocalIndices");
    CHECK((faceIndices[1] == 2));
    faceIndices = test.getCodim2EntityLocalIndices(5);
    INFO("getCodim2EntityLocalIndices");
    CHECK((faceIndices[0] == 2));
    INFO("getCodim2EntityLocalIndices");
    CHECK((faceIndices[1] == 4));
    faceIndices = test.getCodim2EntityLocalIndices(6);
    INFO("getCodim2EntityLocalIndices");
    CHECK((faceIndices[0] == 4));
    INFO("getCodim2EntityLocalIndices");
    CHECK((faceIndices[1] == 3));
    faceIndices = test.getCodim2EntityLocalIndices(7);
    INFO("getCodim2EntityLocalIndices");
    CHECK((faceIndices[0] == 3));
    INFO("getCodim2EntityLocalIndices");
    CHECK((faceIndices[1] == 1));
    faceIndices.resize(1);
    faceIndices = test.getCodim3EntityLocalIndices(0);
    INFO("getCodim3EntityLocalIndices");
    CHECK((faceIndices[0] == 0));
    faceIndices = test.getCodim3EntityLocalIndices(1);
    INFO("getCodim3EntityLocalIndices");
    CHECK((faceIndices[0] == 1));
    faceIndices = test.getCodim3EntityLocalIndices(2);
    INFO("getCodim3EntityLocalIndices");
    CHECK((faceIndices[0] == 2));
    faceIndices = test.getCodim3EntityLocalIndices(3);
    INFO("getCodim3EntityLocalIndices");
    CHECK((faceIndices[0] == 3));
    faceIndices = test.getCodim3EntityLocalIndices(4);
    INFO("getCodim3EntityLocalIndices");
    CHECK((faceIndices[0] == 4));
    INFO("quadrature rules");
    CHECK((test.getGaussQuadratureRule(3)->order() >= 3));
    INFO("quadrature rules");
    CHECK((test.getGaussQuadratureRule(5)->order() >= 5));
    INFO("quadrature rules");
    CHECK((test.getGaussQuadratureRule(7)->order() >= 7));  ///\todo implement
                                                            /// more quadrature
                                                            /// rules
    // assert_debug(("quadrature
    // rules",test.getGaussQuadratureRule(9)->order()>=9));
    // assert_debug(("quadrature
    // rules",test.getGaussQuadratureRule(11)->order()>=11));

    // testing functionality of abstract parent classes

    INFO("number of nodes");
    CHECK((test.getNumberOfNodes() == 5));
    INFO("type of geometry");
    CHECK((test.getGeometryType() == Geometry::ReferenceGeometryType::PYRAMID));
    ///\todo testing that the refinement maps behave exactly like the forwarded
    /// calls of this class
}
