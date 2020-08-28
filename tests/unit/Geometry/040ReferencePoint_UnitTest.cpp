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
#include "Geometry/ReferencePoint.h"
#include <iostream>
#include "Logger.h"

#include "Geometry/PointReference.h"
#include "Geometry/Mappings/MappingToRefPointToPoint.h"
#include "Integration/QuadratureRules/GaussQuadratureRule.h"

#include <cmath>

#include "../catch.hpp"

using namespace hpgem;
using Geometry::ReferencePoint;

TEST_CASE("040ReferencePoint_UnitTest", "[040ReferencePoint_UnitTest]") {
    ReferencePoint& test = ReferencePoint::Instance();

    Geometry::PointReference<0> pTest = {};

    // testing basic functionality

    INFO("isInternalPoint");
    CHECK((test.isInternalPoint(pTest)));
    const Geometry::PointReference<0>& pTest1 = test.getCenter();
    const Geometry::PointReference<0>& pTest2 =
        test.getReferenceNodeCoordinate(0);
    std::cout << test.getName();

    // getLocalNodeIndex should always break since dimension -1 entities dont
    // exist

    std::vector<std::size_t> base(1), transformed(1);

    // testing mappings and quadrature rules

    INFO("getCodim0MappingIndex&Ptr");
    CHECK((test.getCodim0MappingPtr(
               test.getCodim0MappingIndex(base, transformed)) ==
           &Geometry::MappingToRefPointToPoint::Instance()));
    INFO("getCodim0MappingIndex&Ptr");
    CHECK((test.getCodim0MappingPtr(base, transformed) ==
           &Geometry::MappingToRefPointToPoint::Instance()));
    INFO("higher codimensional entities");
    CHECK(test.getNumberOfCodim1Entities() == 0);
    CHECK(test.getNumberOfCodim2Entities() == 0);
    CHECK(test.getNumberOfCodim3Entities() == 0);

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
    CHECK((test.getNumberOfNodes() == 1));
    INFO("type of geometry");
    CHECK((test.getGeometryType() == Geometry::ReferenceGeometryType::POINT));
    ///\todo testing that the refinement maps behave exactly like the forwarded
    /// calls of this class
}
