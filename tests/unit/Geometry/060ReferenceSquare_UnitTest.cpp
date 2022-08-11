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
#include "Geometry/ReferenceSquare.h"
#include "Geometry/ReferenceLine.h"
#include "Geometry/ReferencePoint.h"
#include "Geometry/Mappings/MappingToRefSquareToSquare.h"
#include "Geometry/Mappings/MappingToRefLineToSquare.h"
#include "ReferenceGeometryChecks.h"
#include <iostream>
#include "Logger.h"

#include "Geometry/PointReference.h"
#include "Integration/QuadratureRules/GaussQuadratureRule.h"
#include <cmath>

#include "../catch.hpp"

using namespace hpgem;
using Geometry::ReferenceSquare;

TEST_CASE("ReferenceSquare measure", "[060ReferenceSquare_UnitTest]") {
    testMeasure<2>(ReferenceSquare::Instance());
}

TEST_CASE("060ReferenceSquare_UnitTest", "[060ReferenceSquare_UnitTest]") {
    ReferenceSquare& test = ReferenceSquare::Instance();

    Geometry::PointReference<2> pTest;

    // testing basic functionality

    for (pTest[0] = -3.141; pTest[0] < -1; pTest[0] += 0.1) {
        for (pTest[1] = -3.1416; pTest[1] < 3.1416; pTest[1] += 0.1) {
            INFO("isInternalPoint");
            CHECK((!test.isInternalPoint(pTest)));
        }
    }
    for (; pTest[0] < 1; pTest[0] += 0.1) {
        for (pTest[1] = -3.1416; pTest[1] < -1; pTest[1] += 0.1) {
            INFO("isInternalPoint");
            CHECK((!test.isInternalPoint((pTest))));
        }
        for (; pTest[1] < 1.; pTest[1] += 0.1) {
            INFO("isInternalPoint");
            CHECK((test.isInternalPoint((pTest))));
        }
        for (; pTest[1] < 3.141; pTest[1] += 0.1) {
            INFO("isInternalPoint");
            CHECK((!test.isInternalPoint((pTest))));
        }
    }
    for (; pTest[0] < 3.141; pTest[0] += 0.1) {
        for (pTest[1] = -3.1416; pTest[1] < 3.1416; pTest[1] += 0.1) {
            INFO("isInternalPoint");
            CHECK((!test.isInternalPoint((pTest))));
        }
    }

    pTest = test.getCenter();
    INFO("getCenter");
    CHECK((test.isInternalPoint((pTest)) && std::abs(pTest[0]) < 1e-12 &&
           std::abs(pTest[1]) < 1e-12));
    pTest = test.getReferenceNodeCoordinate(0);
    INFO("getNode 0");
    CHECK((std::abs(pTest[0] + 1) < 1e-12 && std::abs(pTest[1] + 1) < 1e-12));
    pTest = test.getReferenceNodeCoordinate(1);
    INFO("getNode 1");
    CHECK((std::abs(pTest[0] - 1) < 1e-12 && std::abs(pTest[1] + 1) < 1e-12));
    pTest = test.getReferenceNodeCoordinate(2);
    INFO("getNode 2");
    CHECK((std::abs(pTest[0] + 1) < 1e-12 && std::abs(pTest[1] - 1) < 1e-12));
    pTest = test.getReferenceNodeCoordinate(3);
    INFO("getNode 3");
    CHECK((std::abs(pTest[0] - 1) < 1e-12 && std::abs(pTest[1] - 1) < 1e-12));
    std::cout << test.getName();

    INFO("getLocalNodeIndex 0");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 0) ==
           0));  // the nodes of the face must always be
                 // specified IN THIS SPECIFIC ORDER
    INFO("getLocalNodeIndex 0");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 1) ==
           1));  // this is needed because the outward pointing
                 // normal vector
    INFO("getLocalNodeIndex 1");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 0) ==
           0));  // will automatically point outward when
                 // compute it using this node ordering
    INFO("getLocalNodeIndex 1");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 1) ==
           2));  ///\BUG some orderings are wrong
    INFO("getLocalNodeIndex 2");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 0) == 1));
    INFO("getLocalNodeIndex 2");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 1) == 3));
    INFO("getLocalNodeIndex 3");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 0) == 2));
    INFO("getLocalNodeIndex 3");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 1) == 3));

    std::cout << test;

    // testing mappings and quadrature rules

    std::vector<std::size_t> base(4), transformed(4), faceIndices(2);
    for (std::size_t i = 0; i < 4; ++i) {  // doesnt test against reordering of
                                           // the nodes in the first vector
        base[i] = transformed[i] = i;
    }
    INFO("getCodim0MappingIndex&Ptr");
    CHECK((test.getCodim0MappingPtr(
               test.getCodim0MappingIndex(base, transformed)) ==
           &Geometry::MappingToRefSquareToSquare0::Instance()));
    INFO("getCodim0MappingIndex&Ptr");
    CHECK((test.getCodim0MappingPtr(base, transformed) ==
           &Geometry::MappingToRefSquareToSquare0::Instance()));
    transformed[0] = 1;
    transformed[1] = 3;
    transformed[2] = 0;
    transformed[3] = 2;
    INFO("getCodim0MappingIndex&Ptr");
    CHECK((test.getCodim0MappingPtr(
               test.getCodim0MappingIndex(base, transformed)) ==
           &Geometry::MappingToRefSquareToSquare1::Instance()));
    INFO("getCodim0MappingIndex&Ptr");
    CHECK((test.getCodim0MappingPtr(base, transformed) ==
           &Geometry::MappingToRefSquareToSquare1::Instance()));
    transformed[0] = 3;
    transformed[1] = 2;
    transformed[2] = 1;
    transformed[3] = 0;
    INFO("getCodim0MappingIndex&Ptr");
    CHECK((test.getCodim0MappingPtr(
               test.getCodim0MappingIndex(base, transformed)) ==
           &Geometry::MappingToRefSquareToSquare2::Instance()));
    INFO("getCodim0MappingIndex&Ptr");
    CHECK((test.getCodim0MappingPtr(base, transformed) ==
           &Geometry::MappingToRefSquareToSquare2::Instance()));
    transformed[0] = 2;
    transformed[1] = 0;
    transformed[2] = 3;
    transformed[3] = 1;
    INFO("getCodim0MappingIndex&Ptr");
    CHECK((test.getCodim0MappingPtr(
               test.getCodim0MappingIndex(base, transformed)) ==
           &Geometry::MappingToRefSquareToSquare3::Instance()));
    INFO("getCodim0MappingIndex&Ptr");
    CHECK((test.getCodim0MappingPtr(base, transformed) ==
           &Geometry::MappingToRefSquareToSquare3::Instance()));
    transformed[0] = 2;
    transformed[1] = 3;
    transformed[2] = 0;
    transformed[3] = 1;
    INFO("getCodim0MappingIndex&Ptr");
    CHECK((test.getCodim0MappingPtr(
               test.getCodim0MappingIndex(base, transformed)) ==
           &Geometry::MappingToRefSquareToSquare4::Instance()));
    INFO("getCodim0MappingIndex&Ptr");
    CHECK((test.getCodim0MappingPtr(base, transformed) ==
           &Geometry::MappingToRefSquareToSquare4::Instance()));
    transformed[0] = 1;
    transformed[1] = 0;
    transformed[2] = 3;
    transformed[3] = 2;
    INFO("getCodim0MappingIndex&Ptr");
    CHECK((test.getCodim0MappingPtr(
               test.getCodim0MappingIndex(base, transformed)) ==
           &Geometry::MappingToRefSquareToSquare5::Instance()));
    INFO("getCodim0MappingIndex&Ptr");
    CHECK((test.getCodim0MappingPtr(base, transformed) ==
           &Geometry::MappingToRefSquareToSquare5::Instance()));
    transformed[0] = 3;
    transformed[1] = 1;
    transformed[2] = 2;
    transformed[3] = 0;
    INFO("getCodim0MappingIndex&Ptr");
    CHECK((test.getCodim0MappingPtr(
               test.getCodim0MappingIndex(base, transformed)) ==
           &Geometry::MappingToRefSquareToSquare6::Instance()));
    INFO("getCodim0MappingIndex&Ptr");
    CHECK((test.getCodim0MappingPtr(base, transformed) ==
           &Geometry::MappingToRefSquareToSquare6::Instance()));
    transformed[0] = 0;
    transformed[1] = 2;
    transformed[2] = 1;
    transformed[3] = 3;
    INFO("getCodim0MappingIndex&Ptr");
    CHECK((test.getCodim0MappingPtr(
               test.getCodim0MappingIndex(base, transformed)) ==
           &Geometry::MappingToRefSquareToSquare7::Instance()));
    INFO("getCodim0MappingIndex&Ptr");
    CHECK((test.getCodim0MappingPtr(base, transformed) ==
           &Geometry::MappingToRefSquareToSquare7::Instance()));

    INFO("higher codimensional entities");
    CHECK(test.getNumberOfCodim1Entities() == 4);
    CHECK(test.getNumberOfCodim2Entities() == 4);
    CHECK(test.getNumberOfCodim3Entities() == 0);
    INFO("getCodim1ReferenceGeometry");
    CHECK(test.getCodim1ReferenceGeometry(0) ==
          &Geometry::ReferenceLine::Instance());
    CHECK(test.getCodim1ReferenceGeometry(1) ==
          &Geometry::ReferenceLine::Instance());
    CHECK(test.getCodim1ReferenceGeometry(2) ==
          &Geometry::ReferenceLine::Instance());
    CHECK(test.getCodim1ReferenceGeometry(3) ==
          &Geometry::ReferenceLine::Instance());
    INFO("getCodim2ReferenceGeometry");
    CHECK(test.getCodim2ReferenceGeometry(0) ==
          &Geometry::ReferencePoint::Instance());
    CHECK(test.getCodim2ReferenceGeometry(1) ==
          &Geometry::ReferencePoint::Instance());
    CHECK(test.getCodim2ReferenceGeometry(2) ==
          &Geometry::ReferencePoint::Instance());
    CHECK(test.getCodim2ReferenceGeometry(3) ==
          &Geometry::ReferencePoint::Instance());

    INFO("getCodim1MappingPtr");
    CHECK((test.getCodim1MappingPtr(0) ==
           &Geometry::MappingToRefLineToSquare0::Instance()));
    INFO("getCodim1MappingPtr");
    CHECK((test.getCodim1MappingPtr(1) ==
           &Geometry::MappingToRefLineToSquare1::Instance()));
    INFO("getCodim1MappingPtr");
    CHECK((test.getCodim1MappingPtr(2) ==
           &Geometry::MappingToRefLineToSquare2::Instance()));
    INFO("getCodim1MappingPtr");
    CHECK((test.getCodim1MappingPtr(3) ==
           &Geometry::MappingToRefLineToSquare3::Instance()));
    faceIndices = test.getCodim1EntityLocalIndices(0);
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[0] == test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 0)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[1] == test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 1)));
    faceIndices = test.getCodim1EntityLocalIndices(1);
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[0] == test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 0)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[1] == test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 1)));
    faceIndices = test.getCodim1EntityLocalIndices(2);
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[0] == test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 0)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[1] == test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 1)));
    faceIndices = test.getCodim1EntityLocalIndices(3);
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[0] == test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 0)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[1] == test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 1)));
    faceIndices.resize(1);
    faceIndices = test.getCodim2EntityLocalIndices(0);
    INFO("getCodim2EntityLocalIndices");
    CHECK((faceIndices[0] == 0));
    faceIndices = test.getCodim2EntityLocalIndices(1);
    INFO("getCodim2EntityLocalIndices");
    CHECK((faceIndices[0] == 1));
    faceIndices = test.getCodim2EntityLocalIndices(2);
    INFO("getCodim2EntityLocalIndices");
    CHECK((faceIndices[0] == 2));
    faceIndices = test.getCodim2EntityLocalIndices(3);
    INFO("getCodim2EntityLocalIndices");
    CHECK((faceIndices[0] == 3));

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
    CHECK((test.getNumberOfNodes() == 4));
    INFO("type of geometry");
    CHECK((test.getGeometryType() == Geometry::ReferenceGeometryType::SQUARE));

    ///\todo testing that the refinement maps behave exactly like the forwarded
    /// calls of this class
}
