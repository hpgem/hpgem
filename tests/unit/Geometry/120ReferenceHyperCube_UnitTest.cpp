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
#include "Geometry/ReferenceHypercube.h"
#include "Geometry/ReferenceCube.h"
#include "Geometry/ReferenceSquare.h"
#include "Geometry/ReferenceLine.h"
#include "Geometry/ReferencePoint.h"
#include <iostream>
#include "Logger.h"

#include "Geometry/PointReference.h"
#include "Geometry/Mappings/MappingToRefCubeToHypercube.h"
#include "Integration/QuadratureRules/GaussQuadratureRule.h"
#include <cmath>


#include "../catch.hpp"

using namespace hpgem;
using Geometry::ReferenceHypercube;

TEST_CASE("120ReferenceHyperCube_UnitTest",
          "[120ReferenceHyperCube_UnitTest]") {
    ReferenceHypercube& test = ReferenceHypercube::Instance();

    Geometry::PointReference<4> pTest;

    // testing basic functionality

    for (pTest[0] = -1.51; pTest[0] < -1.; pTest[0] += 0.4) {
        for (pTest[1] = -1.51; pTest[1] < 1.51; pTest[1] += 0.4) {
            for (pTest[2] = -1.51; pTest[2] < 1.51; pTest[2] += 0.4) {
                for (pTest[3] = -1.51; pTest[3] < 1.51; pTest[3] += 0.4) {
                    INFO("isInternalPoint");
                    CHECK((!test.isInternalPoint((pTest))));
                }
            }
        }
    }
    for (; pTest[0] < 1; pTest[0] += 0.4) {
        for (pTest[1] = -1.51; pTest[1] < -1.; pTest[1] += 0.4) {
            for (pTest[2] = -1.51; pTest[2] < 1.51; pTest[2] += 0.4) {
                for (pTest[3] = -1.51; pTest[3] < 1.51; pTest[3] += 0.4) {
                    INFO("isInternalPoint");
                    CHECK((!test.isInternalPoint((pTest))));
                }
            }
        }
        for (; pTest[1] < 1.; pTest[1] += 0.4) {
            for (pTest[2] = -1.51; pTest[2] < -1.; pTest[2] += 0.4) {
                for (pTest[3] = -1.51; pTest[3] < 1.51; pTest[3] += 0.4) {
                    INFO("isInternalPoint");
                    CHECK((!test.isInternalPoint((pTest))));
                }
            }
            for (; pTest[2] < 1.; pTest[2] += 0.4) {
                for (pTest[3] = -1.51; pTest[3] < -1.; pTest[3] += 0.4) {
                    INFO("isInternalPoint");
                    CHECK((!test.isInternalPoint((pTest))));
                }
                for (; pTest[3] < 1.; pTest[3] += 0.4) {
                    INFO("isInternalPoint");
                    CHECK((test.isInternalPoint((pTest))));
                }
                for (; pTest[3] < 1.51; pTest[3] += 0.4) {
                    INFO("isInternalPoint");
                    CHECK((!test.isInternalPoint((pTest))));
                }
            }
            for (; pTest[2] < 1.51; pTest[2] += 0.4) {
                for (pTest[3] = -1.51; pTest[3] < 1.51; pTest[3] += 0.4) {
                    INFO("isInternalPoint");
                    CHECK((!test.isInternalPoint((pTest))));
                }
            }
        }
        for (; pTest[1] < 1.51; pTest[1] += 0.4) {
            for (pTest[2] = -1.51; pTest[2] < 1.51; pTest[2] += 0.4) {
                for (pTest[3] = -1.51; pTest[3] < 1.51; pTest[3] += 0.4) {
                    INFO("isInternalPoint");
                    CHECK((!test.isInternalPoint((pTest))));
                }
            }
        }
    }
    for (; pTest[0] < 1.51; pTest[0] += 0.4) {
        for (pTest[1] = -1.51; pTest[1] < 1.51; pTest[1] += 0.4) {
            for (pTest[2] = -1.51; pTest[2] < 1.51; pTest[2] += 0.4) {
                for (pTest[3] = -1.51; pTest[3] < 1.51; pTest[3] += 0.4) {
                    INFO("isInternalPoint");
                    CHECK((!test.isInternalPoint((pTest))));
                }
            }
        }
    }

    pTest = test.getCenter();
    INFO("getCenter");
    CHECK(test.isInternalPoint((pTest)));
    CHECK(std::abs(pTest[0]) < 1e-12);
    CHECK(std::abs(pTest[1]) < 1e-12);
    CHECK(std::abs(pTest[2]) < 1e-12);
    CHECK(std::abs(pTest[3]) < 1e-12);
    pTest = test.getReferenceNodeCoordinate(0);
    INFO("getNode 0");
    CHECK(std::abs(pTest[0] + 1) < 1e-12);
    CHECK(std::abs(pTest[1] + 1) < 1e-12);
    CHECK(std::abs(pTest[2] + 1) < 1e-12);
    CHECK(std::abs(pTest[3] + 1) < 1e-12);
    pTest = test.getReferenceNodeCoordinate(1);
    INFO("getNode 1");
    CHECK(std::abs(pTest[0] - 1) < 1e-12);
    CHECK(std::abs(pTest[1] + 1) < 1e-12);
    CHECK(std::abs(pTest[2] + 1) < 1e-12);
    CHECK(std::abs(pTest[3] + 1) < 1e-12);
    pTest = test.getReferenceNodeCoordinate(2);
    INFO("getNode 2");
    CHECK(std::abs(pTest[0] + 1) < 1e-12);
    CHECK(std::abs(pTest[1] - 1) < 1e-12);
    CHECK(std::abs(pTest[2] + 1) < 1e-12);
    CHECK(std::abs(pTest[3] + 1) < 1e-12);
    pTest = test.getReferenceNodeCoordinate(3);
    INFO("getNode 3");
    CHECK(std::abs(pTest[0] - 1) < 1e-12);
    CHECK(std::abs(pTest[1] - 1) < 1e-12);
    CHECK(std::abs(pTest[2] + 1) < 1e-12);
    CHECK(std::abs(pTest[3] + 1) < 1e-12);
    pTest = test.getReferenceNodeCoordinate(4);
    INFO("getNode 4");
    CHECK(std::abs(pTest[0] + 1) < 1e-12);
    CHECK(std::abs(pTest[1] + 1) < 1e-12);
    CHECK(std::abs(pTest[2] - 1) < 1e-12);
    CHECK(std::abs(pTest[3] + 1) < 1e-12);
    pTest = test.getReferenceNodeCoordinate(5);
    INFO("getNode 5");
    CHECK(std::abs(pTest[0] - 1) < 1e-12);
    CHECK(std::abs(pTest[1] + 1) < 1e-12);
    CHECK(std::abs(pTest[2] - 1) < 1e-12);
    CHECK(std::abs(pTest[3] + 1) < 1e-12);
    pTest = test.getReferenceNodeCoordinate(6);
    INFO("getNode 6");
    CHECK(std::abs(pTest[0] + 1) < 1e-12);
    CHECK(std::abs(pTest[1] - 1) < 1e-12);
    CHECK(std::abs(pTest[2] - 1) < 1e-12);
    CHECK(std::abs(pTest[3] + 1) < 1e-12);
    pTest = test.getReferenceNodeCoordinate(7);
    INFO("getNode 7");
    CHECK(std::abs(pTest[0] - 1) < 1e-12);
    CHECK(std::abs(pTest[1] - 1) < 1e-12);
    CHECK(std::abs(pTest[2] - 1) < 1e-12);
    CHECK(std::abs(pTest[3] + 1) < 1e-12);
    pTest = test.getReferenceNodeCoordinate(8);
    INFO("getNode 8");
    CHECK(std::abs(pTest[0] + 1) < 1e-12);
    CHECK(std::abs(pTest[1] + 1) < 1e-12);
    CHECK(std::abs(pTest[2] + 1) < 1e-12);
    CHECK(std::abs(pTest[3] - 1) < 1e-12);
    pTest = test.getReferenceNodeCoordinate(9);
    INFO("getNode 9");
    CHECK(std::abs(pTest[0] - 1) < 1e-12);
    CHECK(std::abs(pTest[1] + 1) < 1e-12);
    CHECK(std::abs(pTest[2] + 1) < 1e-12);
    CHECK(std::abs(pTest[3] - 1) < 1e-12);
    pTest = test.getReferenceNodeCoordinate(10);
    INFO("getNode 10");
    CHECK(std::abs(pTest[0] + 1) < 1e-12);
    CHECK(std::abs(pTest[1] - 1) < 1e-12);
    CHECK(std::abs(pTest[2] + 1) < 1e-12);
    CHECK(std::abs(pTest[3] - 1) < 1e-12);
    pTest = test.getReferenceNodeCoordinate(11);
    INFO("getNode 11");
    CHECK(std::abs(pTest[0] - 1) < 1e-12);
    CHECK(std::abs(pTest[1] - 1) < 1e-12);
    CHECK(std::abs(pTest[2] + 1) < 1e-12);
    CHECK(std::abs(pTest[3] - 1) < 1e-12);
    pTest = test.getReferenceNodeCoordinate(12);
    INFO("getNode 12");
    CHECK(std::abs(pTest[0] + 1) < 1e-12);
    CHECK(std::abs(pTest[1] + 1) < 1e-12);
    CHECK(std::abs(pTest[2] - 1) < 1e-12);
    CHECK(std::abs(pTest[3] - 1) < 1e-12);
    pTest = test.getReferenceNodeCoordinate(13);
    INFO("getNode 13");
    CHECK(std::abs(pTest[0] - 1) < 1e-12);
    CHECK(std::abs(pTest[1] + 1) < 1e-12);
    CHECK(std::abs(pTest[2] - 1) < 1e-12);
    CHECK(std::abs(pTest[3] - 1) < 1e-12);
    pTest = test.getReferenceNodeCoordinate(14);
    INFO("getNode 14");
    CHECK(std::abs(pTest[0] + 1) < 1e-12);
    CHECK(std::abs(pTest[1] - 1) < 1e-12);
    CHECK(std::abs(pTest[2] - 1) < 1e-12);
    CHECK(std::abs(pTest[3] - 1) < 1e-12);
    pTest = test.getReferenceNodeCoordinate(15);
    INFO("getNode 15");
    CHECK(std::abs(pTest[0] - 1) < 1e-12);
    CHECK(std::abs(pTest[1] - 1) < 1e-12);
    CHECK(std::abs(pTest[2] - 1) < 1e-12);
    CHECK(std::abs(pTest[3] - 1) < 1e-12);
    std::cout << test.getName();

    INFO("getLocalNodeIndex 0");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 0) ==
           0));  // specified IN THIS SPECIFIC ORDER
    INFO("getLocalNodeIndex 0");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 1) ==
           1));  // should at least verify
    INFO("getLocalNodeIndex 0");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 2) ==
           2));  // specified
                 // twice
    INFO("getLocalNodeIndex 0");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 3) ==
           3));  // ordering of the nodes is consistent
    INFO("getLocalNodeIndex 0");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 4) == 4));
    INFO("getLocalNodeIndex 0");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 5) == 5));
    INFO("getLocalNodeIndex 0");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 6) == 6));
    INFO("getLocalNodeIndex 0");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 7) == 7));
    INFO("getLocalNodeIndex 1");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 0) == 0));
    INFO("getLocalNodeIndex 1");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 1) == 1));
    INFO("getLocalNodeIndex 1");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 2) == 2));
    INFO("getLocalNodeIndex 1");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 3) == 3));
    INFO("getLocalNodeIndex 1");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 4) == 8));
    INFO("getLocalNodeIndex 1");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 5) == 9));
    INFO("getLocalNodeIndex 1");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 6) == 10));
    INFO("getLocalNodeIndex 1");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 7) == 11));
    INFO("getLocalNodeIndex 2");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 0) == 0));
    INFO("getLocalNodeIndex 2");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 1) == 1));
    INFO("getLocalNodeIndex 2");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 2) == 4));
    INFO("getLocalNodeIndex 2");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 3) == 5));
    INFO("getLocalNodeIndex 2");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 4) == 8));
    INFO("getLocalNodeIndex 2");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 5) == 9));
    INFO("getLocalNodeIndex 2");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 6) == 12));
    INFO("getLocalNodeIndex 2");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 7) == 13));
    INFO("getLocalNodeIndex 3");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 0) == 0));
    INFO("getLocalNodeIndex 3");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 1) == 2));
    INFO("getLocalNodeIndex 3");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 2) == 4));
    INFO("getLocalNodeIndex 3");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 3) == 6));
    INFO("getLocalNodeIndex 3");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 4) == 8));
    INFO("getLocalNodeIndex 3");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 5) == 10));
    INFO("getLocalNodeIndex 3");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 6) == 12));
    INFO("getLocalNodeIndex 3");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 7) == 14));
    INFO("getLocalNodeIndex 4");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 0) == 1));
    INFO("getLocalNodeIndex 4");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 1) == 3));
    INFO("getLocalNodeIndex 4");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 2) == 5));
    INFO("getLocalNodeIndex 4");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 3) == 7));
    INFO("getLocalNodeIndex 4");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 4) == 9));
    INFO("getLocalNodeIndex 4");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 5) == 11));
    INFO("getLocalNodeIndex 4");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 6) == 13));
    INFO("getLocalNodeIndex 4");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 7) == 15));
    INFO("getLocalNodeIndex 5");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(5, 0) == 2));
    INFO("getLocalNodeIndex 5");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(5, 1) == 3));
    INFO("getLocalNodeIndex 5");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(5, 2) == 6));
    INFO("getLocalNodeIndex 5");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(5, 3) == 7));
    INFO("getLocalNodeIndex 5");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(5, 4) == 10));
    INFO("getLocalNodeIndex 5");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(5, 5) == 11));
    INFO("getLocalNodeIndex 5");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(5, 6) == 14));
    INFO("getLocalNodeIndex 5");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(5, 7) == 15));
    INFO("getLocalNodeIndex 6");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(6, 0) == 4));
    INFO("getLocalNodeIndex 6");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(6, 1) == 5));
    INFO("getLocalNodeIndex 6");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(6, 2) == 6));
    INFO("getLocalNodeIndex 6");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(6, 3) == 7));
    INFO("getLocalNodeIndex 6");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(6, 4) == 12));
    INFO("getLocalNodeIndex 6");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(6, 5) == 13));
    INFO("getLocalNodeIndex 6");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(6, 6) == 14));
    INFO("getLocalNodeIndex 6");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(6, 7) == 15));
    INFO("getLocalNodeIndex 7");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(7, 0) == 8));
    INFO("getLocalNodeIndex 7");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(7, 1) == 9));
    INFO("getLocalNodeIndex 7");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(7, 2) == 10));
    INFO("getLocalNodeIndex 7");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(7, 3) == 11));
    INFO("getLocalNodeIndex 7");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(7, 4) == 12));
    INFO("getLocalNodeIndex 7");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(7, 5) == 13));
    INFO("getLocalNodeIndex 7");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(7, 6) == 14));
    INFO("getLocalNodeIndex 7");
    CHECK((test.getLocalNodeIndexFromFaceAndIndexOnFace(7, 7) == 15));
    std::cout << test;

    // testing mappings and quadrature rules

    std::vector<std::size_t> faceIndices(8);
    // there is no 5D element so codim0mappings are not needed

    INFO("higher codimensional entities");
    CHECK(test.getNumberOfCodim1Entities() == 8);
    CHECK(test.getNumberOfCodim2Entities() == 24);
    CHECK(test.getNumberOfCodim3Entities() == 32);
    INFO("getCodim1ReferenceGeometry");
    CHECK(test.getCodim1ReferenceGeometry(0) ==
          &Geometry::ReferenceCube::Instance());
    CHECK(test.getCodim1ReferenceGeometry(1) ==
          &Geometry::ReferenceCube::Instance());
    CHECK(test.getCodim1ReferenceGeometry(2) ==
          &Geometry::ReferenceCube::Instance());
    CHECK(test.getCodim1ReferenceGeometry(3) ==
          &Geometry::ReferenceCube::Instance());
    CHECK(test.getCodim1ReferenceGeometry(4) ==
          &Geometry::ReferenceCube::Instance());
    CHECK(test.getCodim1ReferenceGeometry(5) ==
          &Geometry::ReferenceCube::Instance());
    CHECK(test.getCodim1ReferenceGeometry(6) ==
          &Geometry::ReferenceCube::Instance());
    CHECK(test.getCodim1ReferenceGeometry(7) ==
          &Geometry::ReferenceCube::Instance());
    INFO("getCodim2ReferenceGeometry");
    CHECK(test.getCodim2ReferenceGeometry(0) ==
          &Geometry::ReferenceSquare::Instance());
    CHECK(test.getCodim2ReferenceGeometry(1) ==
          &Geometry::ReferenceSquare::Instance());
    CHECK(test.getCodim2ReferenceGeometry(2) ==
          &Geometry::ReferenceSquare::Instance());
    CHECK(test.getCodim2ReferenceGeometry(3) ==
          &Geometry::ReferenceSquare::Instance());
    CHECK(test.getCodim2ReferenceGeometry(4) ==
          &Geometry::ReferenceSquare::Instance());
    CHECK(test.getCodim2ReferenceGeometry(5) ==
          &Geometry::ReferenceSquare::Instance());
    CHECK(test.getCodim2ReferenceGeometry(6) ==
          &Geometry::ReferenceSquare::Instance());
    CHECK(test.getCodim2ReferenceGeometry(7) ==
          &Geometry::ReferenceSquare::Instance());
    CHECK(test.getCodim2ReferenceGeometry(8) ==
          &Geometry::ReferenceSquare::Instance());
    CHECK(test.getCodim2ReferenceGeometry(9) ==
          &Geometry::ReferenceSquare::Instance());
    CHECK(test.getCodim2ReferenceGeometry(10) ==
          &Geometry::ReferenceSquare::Instance());
    CHECK(test.getCodim2ReferenceGeometry(11) ==
          &Geometry::ReferenceSquare::Instance());
    CHECK(test.getCodim2ReferenceGeometry(12) ==
          &Geometry::ReferenceSquare::Instance());
    CHECK(test.getCodim2ReferenceGeometry(13) ==
          &Geometry::ReferenceSquare::Instance());
    CHECK(test.getCodim2ReferenceGeometry(14) ==
          &Geometry::ReferenceSquare::Instance());
    CHECK(test.getCodim2ReferenceGeometry(15) ==
          &Geometry::ReferenceSquare::Instance());
    CHECK(test.getCodim2ReferenceGeometry(16) ==
          &Geometry::ReferenceSquare::Instance());
    CHECK(test.getCodim2ReferenceGeometry(17) ==
          &Geometry::ReferenceSquare::Instance());
    CHECK(test.getCodim2ReferenceGeometry(18) ==
          &Geometry::ReferenceSquare::Instance());
    CHECK(test.getCodim2ReferenceGeometry(19) ==
          &Geometry::ReferenceSquare::Instance());
    CHECK(test.getCodim2ReferenceGeometry(20) ==
          &Geometry::ReferenceSquare::Instance());
    CHECK(test.getCodim2ReferenceGeometry(21) ==
          &Geometry::ReferenceSquare::Instance());
    CHECK(test.getCodim2ReferenceGeometry(22) ==
          &Geometry::ReferenceSquare::Instance());
    CHECK(test.getCodim2ReferenceGeometry(23) ==
          &Geometry::ReferenceSquare::Instance());
    INFO("getCodim1MappingPtr");
    CHECK((test.getCodim1MappingPtr(0) ==
           &Geometry::MappingToRefCubeToHypercube0::Instance()));
    INFO("getCodim1MappingPtr");
    CHECK((test.getCodim1MappingPtr(1) ==
           &Geometry::MappingToRefCubeToHypercube1::Instance()));
    INFO("getCodim1MappingPtr");
    CHECK((test.getCodim1MappingPtr(2) ==
           &Geometry::MappingToRefCubeToHypercube2::Instance()));
    INFO("getCodim1MappingPtr");
    CHECK((test.getCodim1MappingPtr(3) ==
           &Geometry::MappingToRefCubeToHypercube3::Instance()));
    INFO("getCodim1MappingPtr");
    CHECK((test.getCodim1MappingPtr(4) ==
           &Geometry::MappingToRefCubeToHypercube4::Instance()));
    INFO("getCodim1MappingPtr");
    CHECK((test.getCodim1MappingPtr(5) ==
           &Geometry::MappingToRefCubeToHypercube5::Instance()));
    INFO("getCodim1MappingPtr");
    CHECK((test.getCodim1MappingPtr(6) ==
           &Geometry::MappingToRefCubeToHypercube6::Instance()));
    INFO("getCodim1MappingPtr");
    CHECK((test.getCodim1MappingPtr(7) ==
           &Geometry::MappingToRefCubeToHypercube7::Instance()));
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
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[4] == test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 4)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[5] == test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 5)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[6] == test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 6)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[7] == test.getLocalNodeIndexFromFaceAndIndexOnFace(0, 7)));
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
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[3] == test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 3)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[4] == test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 4)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[5] == test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 5)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[6] == test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 6)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[7] == test.getLocalNodeIndexFromFaceAndIndexOnFace(1, 7)));
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
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[3] == test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 3)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[4] == test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 4)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[5] == test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 5)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[6] == test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 6)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[7] == test.getLocalNodeIndexFromFaceAndIndexOnFace(2, 7)));
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
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[3] == test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 3)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[4] == test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 4)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[5] == test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 5)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[6] == test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 6)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[7] == test.getLocalNodeIndexFromFaceAndIndexOnFace(3, 7)));
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
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[3] == test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 3)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[4] == test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 4)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[5] == test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 5)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[6] == test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 6)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[7] == test.getLocalNodeIndexFromFaceAndIndexOnFace(4, 7)));
    faceIndices = test.getCodim1EntityLocalIndices(5);
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[0] == test.getLocalNodeIndexFromFaceAndIndexOnFace(5, 0)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[1] == test.getLocalNodeIndexFromFaceAndIndexOnFace(5, 1)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[2] == test.getLocalNodeIndexFromFaceAndIndexOnFace(5, 2)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[3] == test.getLocalNodeIndexFromFaceAndIndexOnFace(5, 3)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[4] == test.getLocalNodeIndexFromFaceAndIndexOnFace(5, 4)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[5] == test.getLocalNodeIndexFromFaceAndIndexOnFace(5, 5)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[6] == test.getLocalNodeIndexFromFaceAndIndexOnFace(5, 6)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[7] == test.getLocalNodeIndexFromFaceAndIndexOnFace(5, 7)));
    faceIndices = test.getCodim1EntityLocalIndices(6);
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[0] == test.getLocalNodeIndexFromFaceAndIndexOnFace(6, 0)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[1] == test.getLocalNodeIndexFromFaceAndIndexOnFace(6, 1)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[2] == test.getLocalNodeIndexFromFaceAndIndexOnFace(6, 2)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[3] == test.getLocalNodeIndexFromFaceAndIndexOnFace(6, 3)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[4] == test.getLocalNodeIndexFromFaceAndIndexOnFace(6, 4)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[5] == test.getLocalNodeIndexFromFaceAndIndexOnFace(6, 5)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[6] == test.getLocalNodeIndexFromFaceAndIndexOnFace(6, 6)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[7] == test.getLocalNodeIndexFromFaceAndIndexOnFace(6, 7)));
    faceIndices = test.getCodim1EntityLocalIndices(7);
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[0] == test.getLocalNodeIndexFromFaceAndIndexOnFace(7, 0)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[1] == test.getLocalNodeIndexFromFaceAndIndexOnFace(7, 1)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[2] == test.getLocalNodeIndexFromFaceAndIndexOnFace(7, 2)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[3] == test.getLocalNodeIndexFromFaceAndIndexOnFace(7, 3)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[4] == test.getLocalNodeIndexFromFaceAndIndexOnFace(7, 4)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[5] == test.getLocalNodeIndexFromFaceAndIndexOnFace(7, 5)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[6] == test.getLocalNodeIndexFromFaceAndIndexOnFace(7, 6)));
    INFO("getCodim1EntityLocalIndices");
    CHECK(
        (faceIndices[7] == test.getLocalNodeIndexFromFaceAndIndexOnFace(7, 7)));
    // other codimensions are not implemented

    INFO("quadrature rules");
    CHECK((test.getGaussQuadratureRule(3)->order() >=
           3));  ///\todo implement more quadrature rules
    // assert_debug(("quadrature
    // rules",test.getGaussQuadratureRule(5)->order()>=5));
    // assert_debug(("quadrature
    // rules",test.getGaussQuadratureRule(7)->order()>=7));
    // assert_debug(("quadrature
    // rules",test.getGaussQuadratureRule(9)->order()>=9));
    // assert_debug(("quadrature
    // rules",test.getGaussQuadratureRule(11)->order()>=11));

    // testing functionality of abstract parent classes

    INFO("number of nodes");
    CHECK((test.getNumberOfNodes() == 16));
    INFO("type of geometry");
    CHECK(
        (test.getGeometryType() == Geometry::ReferenceGeometryType::HYPERCUBE));
    ///\todo testing that the refinement maps behave exactly like the forwarded
    /// calls of this class
}
