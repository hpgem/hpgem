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
// correct execution of all prior unit tests A lot of the functionality of
// Element, Face and Edge only makes sense if they are embedded in a larger mesh
// that is completely set up Setting up a testing case to validate this
// functionality is as hard as generating a mesh. As such this functionality
// will be tested as part of the mesh generation self tests

#include "Base/PhysicalElement.h"
#include "Base/Face.h"
#include "Logger.h"

#include "FE/BasisFunctionsMonomials.h"
#include "Integration/QuadratureRules/GaussQuadratureRulesForSquare.h"

#include "Geometry/PointPhysical.h"
#include "FE/BasisFunctionSet.h"
#include "Base/Element.h"
#include "Geometry/PointReference.h"
#include "FE/BaseBasisFunction.h"

#include "../catch.hpp"

using namespace hpgem;
TEST_CASE("050Face_UnitTest", "[050Face_UnitTest]") {

    std::vector<std::size_t> pointIndexes;
    std::vector<Geometry::PointPhysical<3>> nodes;

    Geometry::PointPhysical<3> point;

    // Cube
    pointIndexes.push_back(4);
    pointIndexes.push_back(7);
    pointIndexes.push_back(10);
    pointIndexes.push_back(11);
    pointIndexes.push_back(12);
    pointIndexes.push_back(13);
    pointIndexes.push_back(14);
    pointIndexes.push_back(15);

    for (double i = 0.; i < 1 - 1e-10; i += 0.1) {
        point[0] = 1. + i;
        point[1] = 2. + i;
        point[2] = 3. + i;
        nodes.push_back(point);
    }

    point[0] = 3.5;
    point[1] = 4.6;
    point[2] = 5.4;
    nodes.push_back(point);
    point[0] = 6.7;
    point[1] = 2.8;
    point[2] = 5.7;
    nodes.push_back(point);
    point[0] = 1.4;
    point[1] = 2.4;
    point[2] = 5.4;
    nodes.push_back(point);
    point[0] = 1.7;
    point[1] = 2.7;
    point[2] = 5.7;
    nodes.push_back(point);
    point[0] = 3.5;
    point[1] = 4.6;
    point[2] = 7.4;
    nodes.push_back(point);
    point[0] = 6.7;
    point[1] = 2.8;
    point[2] = 7.7;
    nodes.push_back(point);

    FE::BasisFunctionSet* basisFunctions = new FE::BasisFunctionSet(3);

    FE::assembleMonomialBasisFunctions3D(*basisFunctions, 3);

    std::vector<std::shared_ptr<const FE::BasisFunctionSet>> vectorOfFunctions(
        1, std::shared_ptr<FE::BasisFunctionSet>(basisFunctions));

    Base::Zone zone ({"dummy", 0});

    Base::Element element(pointIndexes, &vectorOfFunctions, nodes, 3, 14, 18, &zone);
    for (std::size_t i = 0; i < 3; ++i) {
        element.setDefaultBasisFunctionSet(0, i);
    }

    // Face 4 of the cube, local node indices 2,3,6,7 -> points 10,11,14,15
    // (3.5, 4.6, 5.4), (6.7, 2.8, 5.7), (3.5, 4.6, 7.4), (6.7, 2,8, 7.7)
    Base::Face test(&element, 4, Geometry::FaceType::WALL_BC, 3);

    INFO("quadrature rule");
    CHECK((test.getGaussQuadratureRule() != nullptr));

    test.setGaussQuadratureRule(&QuadratureRules::Cn2_3_4::Instance());

    const QuadratureRules::GaussQuadratureRule& rule =
        *test.getGaussQuadratureRule();
    ///\todo figure out why the '::Instance' alters the typeid (it also alters
    /// the hash_code of the typeid)
    INFO("setQuadratureRule");
    CHECK((typeid(rule) == typeid(QuadratureRules::Cn2_3_4::Instance())));

    // check set*BasisFunctionSet without breaking preconditions...

    INFO("getElementPtr");
    CHECK((test.getPtrElementLeft() == &element));

    Geometry::PointReference<2> refPoint;
    Geometry::PointReference<3> point3D;
    for (std::size_t i = 0; i < basisFunctions->size(); ++i) {
        for (refPoint[0] = -1.5; refPoint[0] < 1.51; refPoint[0] += 0.2) {
            for (refPoint[1] = -1.5; refPoint[1] < 1.51; refPoint[1] += 0.2) {
                point3D = test.mapRefFaceToRefElemL((refPoint));
                INFO("basisFunctions");
                CHECK((test.basisFunction(i, (refPoint)) ==
                       (*basisFunctions)[i]->eval((point3D))));
                INFO("basisFunctions");
                CHECK((test.basisFunctionDeriv(i, 0, (refPoint)) ==
                       (*basisFunctions)[i]->evalDeriv0((point3D))));
                INFO("basisFunctions");
                CHECK((test.basisFunctionDeriv(i, 1, (refPoint)) ==
                       (*basisFunctions)[i]->evalDeriv1((point3D))));
                INFO("basisFunctions");
                CHECK((test.basisFunctionDeriv(i, 2, (refPoint)) ==
                       (*basisFunctions)[i]->evalDeriv2((point3D))));
            }
        }
    }

    // Expected diam: sqrt(18.77) ~ 4.3324
    double diam = (nodes[10] - nodes[15]).getCoordinates().l2Norm();
    INFO("getDiameter");
    CHECK(std::abs(test.getDiameter() - diam) < 1e-12);
}
