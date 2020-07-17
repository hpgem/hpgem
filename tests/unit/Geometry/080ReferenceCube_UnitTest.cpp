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
#include "Geometry/ReferenceCube.h"
#include "Geometry/ReferenceSquare.h"
#include "Geometry/ReferenceLine.h"
#include "Geometry/ReferencePoint.h"
#include <iostream>
#include "Logger.h"

#include "Geometry/PointReference.h"
#include "Geometry/Mappings/MappingToRefSquareToCube.h"
#include "Geometry/Mappings/MappingToRefCubeToCube.h"
#include "Integration/QuadratureRules/GaussQuadratureRule.h"
#include <cmath>

#define CATCH_CONFIG_MAIN
#include "../catch.hpp"

using namespace hpgem;
using Geometry::ReferenceCube;

TEST_CASE("080ReferenceCube_UnitTest", "[080ReferenceCube_UnitTest]") {
    ReferenceCube& test = ReferenceCube::Instance();

    Geometry::PointReference<3> pTest;

    // testing basic functionality

    for (pTest[0] = -1.51; pTest[0] < -1.; pTest[0] += 0.2) {
        for (pTest[1] = -1.51; pTest[1] < 1.51; pTest[1] += 0.2) {
            for (pTest[2] = -1.51; pTest[2] < 1.51; pTest[2] += 0.2) {
                INFO("isInternalPoint");
                CHECK((!test.isInternalPoint((pTest))));
            }
        }
    }
    for (; pTest[0] < 1; pTest[0] += 0.2) {
        for (pTest[1] = -1.51; pTest[1] < -1.; pTest[1] += 0.2) {
            for (pTest[2] = -1.51; pTest[2] < 1.51; pTest[2] += 0.2) {
                INFO("isInternalPoint");
                CHECK((!test.isInternalPoint((pTest))));
            }
        }
        for (; pTest[1] < 1.; pTest[1] += 0.2) {
            for (pTest[2] = -1.51; pTest[2] < -1.; pTest[2] += 0.2) {
                INFO("isInternalPoint");
                CHECK((!test.isInternalPoint((pTest))));
            }
            for (; pTest[2] < 1.; pTest[2] += 0.2) {
                INFO("isInternalPoint");
                CHECK((test.isInternalPoint((pTest))));
            }
            for (; pTest[2] < 1.51; pTest[2] += 0.2) {
                INFO("isInternalPoint");
                CHECK((!test.isInternalPoint((pTest))));
            }
        }
        for (; pTest[1] < 1.51; pTest[1] += 0.2) {
            for (pTest[2] = -1.51; pTest[2] < 1.51; pTest[2] += 0.2) {
                INFO("isInternalPoint");
                CHECK((!test.isInternalPoint((pTest))));
            }
        }
    }
    for (; pTest[0] < 1.51; pTest[0] += 0.2) {
        for (pTest[1] = -1.51; pTest[1] < 3.1416; pTest[1] += 0.2) {
            for (pTest[2] = -1.51; pTest[2] < 3.1416; pTest[2] += 0.2) {
                INFO("isInternalPoint");
                CHECK((!test.isInternalPoint((pTest))));
            }
        }
    }

    pTest = test.getCenter();
    INFO("getCenter");
    CHECK(test.isInternalPoint(pTest));
    CHECK(std::abs(pTest[0]) < 1e-12);
    CHECK(std::abs(pTest[1]) < 1e-12);
    CHECK(std::abs(pTest[2]) < 1e-12);
    pTest = test.getReferenceNodeCoordinate(0);
    INFO("getNode 0");
    CHECK(std::abs(pTest[0] + 1) < 1e-12);
    CHECK(std::abs(pTest[1] + 1) < 1e-12);
    CHECK(std::abs(pTest[2] + 1) < 1e-12);
    pTest = test.getReferenceNodeCoordinate(1);
    INFO("getNode 1");
    CHECK(std::abs(pTest[0] - 1) < 1e-12);
    CHECK(std::abs(pTest[1] + 1) < 1e-12);
    CHECK(std::abs(pTest[2] + 1) < 1e-12);
    pTest = test.getReferenceNodeCoordinate(2);
    INFO("getNode 2");
    CHECK(std::abs(pTest[0] + 1) < 1e-12);
    CHECK(std::abs(pTest[1] - 1) < 1e-12);
    CHECK(std::abs(pTest[2] + 1) < 1e-12);
    pTest = test.getReferenceNodeCoordinate(3);
    INFO("getNode 3");
    CHECK(std::abs(pTest[0] - 1) < 1e-12);
    CHECK(std::abs(pTest[1] - 1) < 1e-12);
    CHECK(std::abs(pTest[2] + 1) < 1e-12);
    pTest = test.getReferenceNodeCoordinate(4);
    INFO("getNode 4");
    CHECK(std::abs(pTest[0] + 1) < 1e-12);
    CHECK(std::abs(pTest[1] + 1) < 1e-12);
    CHECK(std::abs(pTest[2] - 1) < 1e-12);
    pTest = test.getReferenceNodeCoordinate(5);
    INFO("getNode 5");
    CHECK(std::abs(pTest[0] - 1) < 1e-12);
    CHECK(std::abs(pTest[1] + 1) < 1e-12);
    CHECK(std::abs(pTest[2] - 1) < 1e-12);
    pTest = test.getReferenceNodeCoordinate(6);
    INFO("getNode 6");
    CHECK(std::abs(pTest[0] + 1) < 1e-12);
    CHECK(std::abs(pTest[1] - 1) < 1e-12);
    CHECK(std::abs(pTest[2] - 1) < 1e-12);
    pTest = test.getReferenceNodeCoordinate(7);
    INFO("getNode 7");
    CHECK(std::abs(pTest[0] - 1) < 1e-12);
    CHECK(std::abs(pTest[1] - 1) < 1e-12);
    CHECK(std::abs(pTest[2] - 1) < 1e-12);
    std::cout << test.getName();

    INFO("getLocalNodeIndex 0");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 0) ==
           0));  // specified IN THIS SPECIFIC ORDER
    INFO("getLocalNodeIndex 0");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 1) ==
           1));  // should at least verify
    INFO("getLocalNodeIndex 0");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 2) == 2));  // specified
                                                                       // twice
    INFO("getLocalNodeIndex 0");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 3) ==
           3));  // ordering of the nodes is consistent
    INFO("getLocalNodeIndex 1");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 0) == 0));
    INFO("getLocalNodeIndex 1");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 1) == 1));
    INFO("getLocalNodeIndex 1");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 2) == 4));
    INFO("getLocalNodeIndex 1");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 3) == 5));
    INFO("getLocalNodeIndex 2");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 0) == 0));
    INFO("getLocalNodeIndex 2");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 1) == 2));
    INFO("getLocalNodeIndex 2");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 2) == 4));
    INFO("getLocalNodeIndex 2");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 3) == 6));
    INFO("getLocalNodeIndex 3");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 0) == 1));
    INFO("getLocalNodeIndex 3");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 1) == 3));
    INFO("getLocalNodeIndex 3");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 2) == 5));
    INFO("getLocalNodeIndex 3");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 3) == 7));
    INFO("getLocalNodeIndex 4");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 0) == 2));
    INFO("getLocalNodeIndex 4");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 1) == 3));
    INFO("getLocalNodeIndex 4");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 2) == 6));
    INFO("getLocalNodeIndex 4");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 3) == 7));
    INFO("getLocalNodeIndex 5");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(5, 0) == 4));
    INFO("getLocalNodeIndex 5");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(5, 1) == 5));
    INFO("getLocalNodeIndex 5");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(5, 2) == 6));
    INFO("getLocalNodeIndex 5");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(5, 3) == 7));
    std::cout << test;

    // testing mappings and quadrature rules

    std::vector<std::size_t> base(8), transformed(8), faceIndices(4);
    for (std::size_t i = 0; i < 8; ++i) {  // doesnt test against reordering of
                                           // the nodes in the first vector
        base[i] = transformed[i] = i;
    }
    INFO("getCodim0MappingIndex&Ptr");
    CHECK((test.getCodim0MappingPtr(
               test.getCodim0MappingIndex(base, transformed)) ==
           &Geometry::MappingToRefCubeToCube0::Instance()));
    INFO("getCodim0MappingIndex&Ptr");
    CHECK((test.getCodim0MappingPtr(base, transformed) ==
           &Geometry::MappingToRefCubeToCube0::Instance()));
    transformed[0] = 1;
    transformed[1] = 3;
    transformed[2] = 0;
    transformed[3] = 2;
    transformed[4] = 5;
    transformed[5] = 7;
    transformed[6] = 4;
    transformed[7] = 6;
    INFO("getCodim0MappingIndex&Ptr");
    CHECK((test.getCodim0MappingPtr(
               test.getCodim0MappingIndex(base, transformed)) ==
           &Geometry::MappingToRefCubeToCube1::Instance()));
    INFO("getCodim0MappingIndex&Ptr");
    CHECK((test.getCodim0MappingPtr(base, transformed) ==
           &Geometry::MappingToRefCubeToCube1::Instance()));
    transformed[0] = 3;
    transformed[1] = 2;
    transformed[2] = 1;
    transformed[3] = 0;
    transformed[4] = 7;
    transformed[5] = 6;
    transformed[6] = 5;
    transformed[7] = 4;
    INFO("getCodim0MappingIndex&Ptr");
    CHECK((test.getCodim0MappingPtr(
               test.getCodim0MappingIndex(base, transformed)) ==
           &Geometry::MappingToRefCubeToCube2::Instance()));
    INFO("getCodim0MappingIndex&Ptr");
    CHECK((test.getCodim0MappingPtr(base, transformed) ==
           &Geometry::MappingToRefCubeToCube2::Instance()));
    transformed[0] = 2;
    transformed[1] = 0;
    transformed[2] = 3;
    transformed[3] = 1;
    transformed[4] = 6;
    transformed[5] = 4;
    transformed[6] = 7;
    transformed[7] = 5;
    INFO("getCodim0MappingIndex&Ptr");
    CHECK((test.getCodim0MappingPtr(
               test.getCodim0MappingIndex(base, transformed)) ==
           &Geometry::MappingToRefCubeToCube3::Instance()));
    INFO("getCodim0MappingIndex&Ptr");
    CHECK((test.getCodim0MappingPtr(base, transformed) ==
           &Geometry::MappingToRefCubeToCube3::Instance()));
    transformed[0] = 2;
    transformed[1] = 3;
    transformed[2] = 0;
    transformed[3] = 1;
    transformed[4] = 6;
    transformed[5] = 7;
    transformed[6] = 4;
    transformed[7] = 5;
    INFO("getCodim0MappingIndex&Ptr");
    CHECK((test.getCodim0MappingPtr(
               test.getCodim0MappingIndex(base, transformed)) ==
           &Geometry::MappingToRefCubeToCube4::Instance()));
    INFO("getCodim0MappingIndex&Ptr");
    CHECK((test.getCodim0MappingPtr(base, transformed) ==
           &Geometry::MappingToRefCubeToCube4::Instance()));
    transformed[0] = 1;
    transformed[1] = 0;
    transformed[2] = 3;
    transformed[3] = 2;
    transformed[4] = 5;
    transformed[5] = 4;
    transformed[6] = 7;
    transformed[7] = 6;
    INFO("getCodim0MappingIndex&Ptr");
    CHECK((test.getCodim0MappingPtr(
               test.getCodim0MappingIndex(base, transformed)) ==
           &Geometry::MappingToRefCubeToCube5::Instance()));
    INFO("getCodim0MappingIndex&Ptr");
    CHECK((test.getCodim0MappingPtr(base, transformed) ==
           &Geometry::MappingToRefCubeToCube5::Instance()));
    transformed[0] = 3;
    transformed[1] = 1;
    transformed[2] = 2;
    transformed[3] = 0;
    transformed[4] = 7;
    transformed[5] = 5;
    transformed[6] = 6;
    transformed[7] = 4;
    INFO("getCodim0MappingIndex&Ptr");
    CHECK((test.getCodim0MappingPtr(
               test.getCodim0MappingIndex(base, transformed)) ==
           &Geometry::MappingToRefCubeToCube6::Instance()));
    INFO("getCodim0MappingIndex&Ptr");
    CHECK((test.getCodim0MappingPtr(base, transformed) ==
           &Geometry::MappingToRefCubeToCube6::Instance()));
    transformed[0] = 0;
    transformed[1] = 2;
    transformed[2] = 1;
    transformed[3] = 3;
    transformed[4] = 4;
    transformed[5] = 6;
    transformed[6] = 5;
    transformed[7] = 7;
    INFO("getCodim0MappingIndex&Ptr");
    CHECK((test.getCodim0MappingPtr(
               test.getCodim0MappingIndex(base, transformed)) ==
           &Geometry::MappingToRefCubeToCube7::Instance()));
    INFO("getCodim0MappingIndex&Ptr");
    CHECK((test.getCodim0MappingPtr(base, transformed) ==
           &Geometry::MappingToRefCubeToCube7::Instance()));
    INFO("higher codimensional entities");
    CHECK(test.getNumberOfCodim1Entities() == 6);
CHECK(test.getNumberOfCodim2Entities() == 12);
CHECK(test.getNumberOfCodim3Entities() == 8);
INFO("getCodim1ReferenceGeometry");
CHECK(test.getCodim1ReferenceGeometry(0) ==
      &Geometry::ReferenceSquare::Instance());
CHECK(test.getCodim1ReferenceGeometry(1) ==
      &Geometry::ReferenceSquare::Instance());
CHECK(test.getCodim1ReferenceGeometry(2) ==
      &Geometry::ReferenceSquare::Instance());
CHECK(test.getCodim1ReferenceGeometry(3) ==
      &Geometry::ReferenceSquare::Instance());
CHECK(test.getCodim1ReferenceGeometry(4) ==
      &Geometry::ReferenceSquare::Instance());
CHECK(test.getCodim1ReferenceGeometry(5) ==
      &Geometry::ReferenceSquare::Instance());
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
CHECK(test.getCodim2ReferenceGeometry(8) ==
      &Geometry::ReferenceLine::Instance());
CHECK(test.getCodim2ReferenceGeometry(9) ==
      &Geometry::ReferenceLine::Instance());
CHECK(test.getCodim2ReferenceGeometry(10) ==
      &Geometry::ReferenceLine::Instance());
CHECK(test.getCodim2ReferenceGeometry(11) ==
      &Geometry::ReferenceLine::Instance());
INFO("getCodim1MappingPtr");
CHECK((test.getCodim1MappingPtr(0) ==
       &Geometry::MappingToRefSquareToCube0::Instance()));
INFO("getCodim1MappingPtr");
CHECK((test.getCodim1MappingPtr(1) ==
       &Geometry::MappingToRefSquareToCube1::Instance()));
INFO("getCodim1MappingPtr");
CHECK((test.getCodim1MappingPtr(2) ==
       &Geometry::MappingToRefSquareToCube2::Instance()));
INFO("getCodim1MappingPtr");
CHECK((test.getCodim1MappingPtr(3) ==
       &Geometry::MappingToRefSquareToCube3::Instance()));
INFO("getCodim1MappingPtr");
CHECK((test.getCodim1MappingPtr(4) ==
       &Geometry::MappingToRefSquareToCube4::Instance()));
INFO("getCodim1MappingPtr");
CHECK((test.getCodim1MappingPtr(5) ==
       &Geometry::MappingToRefSquareToCube5::Instance()));
faceIndices = test.getCodim1EntityLocalIndices(0);
INFO("getCodim1EntityLocalIndices");
CHECK((faceIndices[0] == test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 0)));
INFO("getCodim1EntityLocalIndices");
CHECK((faceIndices[1] == test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 1)));
INFO("getCodim1EntityLocalIndices");
CHECK((faceIndices[2] == test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 2)));
INFO("getCodim1EntityLocalIndices");
CHECK((faceIndices[3] == test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 3)));
faceIndices = test.getCodim1EntityLocalIndices(1);
INFO("getCodim1EntityLocalIndices");
CHECK((faceIndices[0] == test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 0)));
INFO("getCodim1EntityLocalIndices");
CHECK((faceIndices[1] == test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 1)));
INFO("getCodim1EntityLocalIndices");
CHECK((faceIndices[2] == test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 2)));
INFO("getCodim1EntityLocalIndices");
CHECK((faceIndices[3] == test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 3)));
faceIndices = test.getCodim1EntityLocalIndices(2);
INFO("getCodim1EntityLocalIndices");
CHECK((faceIndices[0] == test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 0)));
INFO("getCodim1EntityLocalIndices");
CHECK((faceIndices[1] == test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 1)));
INFO("getCodim1EntityLocalIndices");
CHECK((faceIndices[2] == test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 2)));
INFO("getCodim1EntityLocalIndices");
CHECK((faceIndices[3] == test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 3)));
faceIndices = test.getCodim1EntityLocalIndices(3);
INFO("getCodim1EntityLocalIndices");
CHECK((faceIndices[0] == test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 0)));
INFO("getCodim1EntityLocalIndices");
CHECK((faceIndices[1] == test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 1)));
INFO("getCodim1EntityLocalIndices");
CHECK((faceIndices[2] == test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 2)));
INFO("getCodim1EntityLocalIndices");
CHECK((faceIndices[3] == test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 3)));
faceIndices = test.getCodim1EntityLocalIndices(4);
INFO("getCodim1EntityLocalIndices");
CHECK((faceIndices[0] == test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 0)));
INFO("getCodim1EntityLocalIndices");
CHECK((faceIndices[1] == test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 1)));
INFO("getCodim1EntityLocalIndices");
CHECK((faceIndices[2] == test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 2)));
INFO("getCodim1EntityLocalIndices");
CHECK((faceIndices[3] == test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 3)));
faceIndices = test.getCodim1EntityLocalIndices(5);
INFO("getCodim1EntityLocalIndices");
CHECK((faceIndices[0] == test.getLocalNodeIndexFromFaceAndIndexOnFace(5, 0)));
INFO("getCodim1EntityLocalIndices");
CHECK((faceIndices[1] == test.getLocalNodeIndexFromFaceAndIndexOnFace(5, 1)));
INFO("getCodim1EntityLocalIndices");
CHECK((faceIndices[2] == test.getLocalNodeIndexFromFaceAndIndexOnFace(5, 2)));
INFO("getCodim1EntityLocalIndices");
CHECK((faceIndices[3] == test.getLocalNodeIndexFromFaceAndIndexOnFace(5, 3)));
faceIndices.resize(2);
faceIndices = test.getCodim2EntityLocalIndices(0);
INFO("getCodim2EntityLocalIndices");
CHECK((faceIndices[0] == 0));
INFO("getCodim2EntityLocalIndices");
CHECK((faceIndices[1] == 1));
faceIndices = test.getCodim2EntityLocalIndices(1);
INFO("getCodim2EntityLocalIndices");
CHECK((faceIndices[0] == 2));
INFO("getCodim2EntityLocalIndices");
CHECK((faceIndices[1] == 3));
faceIndices = test.getCodim2EntityLocalIndices(2);
INFO("getCodim2EntityLocalIndices");
CHECK((faceIndices[0] == 4));
INFO("getCodim2EntityLocalIndices");
CHECK((faceIndices[1] == 5));
faceIndices = test.getCodim2EntityLocalIndices(3);
INFO("getCodim2EntityLocalIndices");
CHECK((faceIndices[0] == 6));
INFO("getCodim2EntityLocalIndices");
CHECK((faceIndices[1] == 7));
faceIndices = test.getCodim2EntityLocalIndices(4);
INFO("getCodim2EntityLocalIndices");
CHECK((faceIndices[0] == 0));
INFO("getCodim2EntityLocalIndices");
CHECK((faceIndices[1] == 2));
faceIndices = test.getCodim2EntityLocalIndices(5);
INFO("getCodim2EntityLocalIndices");
CHECK((faceIndices[0] == 1));
INFO("getCodim2EntityLocalIndices");
CHECK((faceIndices[1] == 3));
faceIndices = test.getCodim2EntityLocalIndices(6);
INFO("getCodim2EntityLocalIndices");
CHECK((faceIndices[0] == 4));
INFO("getCodim2EntityLocalIndices");
CHECK((faceIndices[1] == 6));
faceIndices = test.getCodim2EntityLocalIndices(7);
INFO("getCodim2EntityLocalIndices");
CHECK((faceIndices[0] == 5));
INFO("getCodim2EntityLocalIndices");
CHECK((faceIndices[1] == 7));
faceIndices = test.getCodim2EntityLocalIndices(8);
INFO("getCodim2EntityLocalIndices");
CHECK((faceIndices[0] == 0));
INFO("getCodim2EntityLocalIndices");
CHECK((faceIndices[1] == 4));
faceIndices = test.getCodim2EntityLocalIndices(9);
INFO("getCodim2EntityLocalIndices");
CHECK((faceIndices[0] == 1));
INFO("getCodim2EntityLocalIndices");
CHECK((faceIndices[1] == 5));
faceIndices = test.getCodim2EntityLocalIndices(10);
INFO("getCodim2EntityLocalIndices");
CHECK((faceIndices[0] == 2));
INFO("getCodim2EntityLocalIndices");
CHECK((faceIndices[1] == 6));
faceIndices = test.getCodim2EntityLocalIndices(11);
INFO("getCodim2EntityLocalIndices");
CHECK((faceIndices[0] == 3));
INFO("getCodim2EntityLocalIndices");
CHECK((faceIndices[1] == 7));
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
faceIndices = test.getCodim3EntityLocalIndices(5);
INFO("getCodim3EntityLocalIndices");
CHECK((faceIndices[0] == 5));
faceIndices = test.getCodim3EntityLocalIndices(6);
INFO("getCodim3EntityLocalIndices");
CHECK((faceIndices[0] == 6));
faceIndices = test.getCodim3EntityLocalIndices(7);
INFO("getCodim3EntityLocalIndices");
CHECK((faceIndices[0] == 7));
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
CHECK((test.getNumberOfNodes() == 8));
INFO("type of geometry");
CHECK((test.getGeometryType() == Geometry::ReferenceGeometryType::CUBE));
///\todo testing that the refinement maps behave exactly like the forwarded
/// calls of this class
}
