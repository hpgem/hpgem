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
#include "Geometry/FaceGeometry.h"
#include "Logger.h"

#include "Geometry/ReferenceLine.h"
#include "Geometry/ReferencePoint.h"
#include "Geometry/ReferenceTetrahedron.h"
#include "Geometry/PointPhysical.h"
#include "Geometry/PointReference.h"
#include "Geometry/ElementGeometry.h"
#include "Geometry/Jacobian.h"
#include "Geometry/PhysicalGeometry.h"

#include <cmath>

#include "../catch.hpp"

using namespace hpgem;
TEST_CASE("230FaceGeometry_UnitTest", "[230FaceGeometry_UnitTest]") {

    // dim0

    std::vector<std::size_t> pointIndexes, leftIndices, rightIndices;
    std::vector<Geometry::PointPhysical<1>> nodes1D;

    Geometry::PointReference<1> point1D, compare1D;
    Geometry::PointPhysical<1> point1Dphys, compare1Dphys;
    Geometry::PointReference<0> orig1D;

    pointIndexes.push_back(4);
    pointIndexes.push_back(7);

    for (double i = 0.; i < 1; i += 0.1) {
        point1Dphys[0] = 1. + i / 10.;
        nodes1D.push_back(point1Dphys);
    }

    Geometry::ElementGeometry* left(
        new Geometry::ElementGeometry(pointIndexes, nodes1D));
    pointIndexes[1] = 2;
    Geometry::ElementGeometry* right(
        new Geometry::ElementGeometry(pointIndexes, nodes1D));

    Geometry::FaceGeometry* test(new Geometry::FaceGeometry(left, 0, right, 0));

    leftIndices.resize(1);
    rightIndices.resize(1);

    ///\todo test the normal vector
    ///\todo also test individual functions

    leftIndices =
        test->getElementGLeft()->getPhysicalGeometry()->getLocalFaceNodeIndices(
            test->localFaceNumberLeft());
    rightIndices = test->getPtrElementGRight()
                       ->getPhysicalGeometry()
                       ->getLocalFaceNodeIndices(test->localFaceNumberRight());

    test->initialiseFaceToFaceMapIndex(leftIndices, rightIndices);

    for (std::size_t i = 0;
         i < test->getReferenceGeometry()->getNumberOfNodes(); ++i) {
        orig1D = test->getReferenceGeometry()->getReferenceNodeCoordinate(i);
        compare1D = test->getElementGLeft()
                        ->getReferenceGeometry()
                        ->getReferenceNodeCoordinate(leftIndices[i]);
        compare1Dphys =
            test->getElementGLeft()->referenceToPhysical((compare1D));
        point1D = test->mapRefFaceToRefElemL((orig1D));
        point1Dphys = test->referenceToPhysical((orig1D));
        INFO("getElementGLeft or localFaceNumberLeft or mapRefFaceToRefElemL");
        CHECK((std::abs(compare1D[0] - point1D[0]) < 1e-12));
        INFO("referenceToPhysical");
        CHECK((std::abs(compare1Dphys[0] - point1Dphys[0]) < 1e-12));
        point1D = test->refFaceToRefElemMapL()->transform((orig1D));
        INFO("refFaceToRefElemMapL");
        CHECK((std::abs(compare1D[0] - point1D[0]) < 1e-12));

        compare1D = test->getPtrElementGRight()
                        ->getReferenceGeometry()
                        ->getReferenceNodeCoordinate(rightIndices[i]);
        compare1Dphys =
            test->getPtrElementGRight()->referenceToPhysical((compare1D));
        point1D = test->mapRefFaceToRefElemR((orig1D));
        INFO(
            "getPtrElementGRight or localFaceNumberRight or "
            "mapRefFaceToRefElemR or mapRefFaceToRefFace");
        CHECK(std::abs(compare1D[0] - point1D[0]) < 1e-12);

        INFO("referenceToPhysical");
        CHECK((std::abs(compare1Dphys[0] - point1Dphys[0]) <
               1e-12));  // probably indirectly verified already,
                         // but this is the most important feature
                         // of a face
        point1D = test->refFaceToRefElemMapR()->transform((orig1D));
        INFO("refFaceToRefElemMapR");
        CHECK((std::abs(compare1D[0] - point1D[0]) < 1e-12));
    }

    const Geometry::ReferenceGeometry& refGeo1 = *test->getReferenceGeometry();
    INFO("Face type internal");
    CHECK((test->getFaceType() == Geometry::FaceType::INTERNAL));
    INFO("Reference geometry point");
    CHECK((typeid(refGeo1) == typeid(Geometry::ReferencePoint)));

    // dim 1

    std::vector<Geometry::PointPhysical<2>> nodes2D;

    Geometry::PointReference<2> point2D, compare2D;
    Geometry::PointPhysical<2> point2Dphys, compare2Dphys;
    Geometry::PointReference<1> orig2D;

    for (double i = 0.; i < 1; i += 0.1) {
        point2Dphys[0] = 1. + i / 10.;
        point2Dphys[1] = 2. + i / 10.;
        nodes2D.push_back(point2Dphys);
    }

    point2Dphys[0] = 3.5;
    point2Dphys[1] = 4.6;
    nodes2D.push_back(point2Dphys);
    point2Dphys[0] = 6.7;
    point2Dphys[1] = 2.8;
    nodes2D.push_back(point2Dphys);
    point2Dphys[0] = 1.41;
    point2Dphys[1] = 6.82;
    nodes2D.push_back(point2Dphys);

    pointIndexes.push_back(12);

    delete left;
    delete test;
    delete right;
    left = new Geometry::ElementGeometry(pointIndexes, nodes2D);
    pointIndexes[2] = 10;
    pointIndexes.push_back(11);
    right = new Geometry::ElementGeometry(pointIndexes, nodes2D);

    test = new Geometry::FaceGeometry(left, 0, right, 0);

    leftIndices.resize(2);
    rightIndices.resize(2);

    ///\todo test the normal vector
    ///\todo also test individual functions

    leftIndices =
        test->getElementGLeft()->getPhysicalGeometry()->getLocalFaceNodeIndices(
            test->localFaceNumberLeft());
    rightIndices = test->getPtrElementGRight()
                       ->getPhysicalGeometry()
                       ->getLocalFaceNodeIndices(test->localFaceNumberRight());

    test->initialiseFaceToFaceMapIndex(leftIndices, rightIndices);

    for (std::size_t i = 0;
         i < test->getReferenceGeometry()->getNumberOfNodes(); ++i) {
        orig2D = test->getReferenceGeometry()->getReferenceNodeCoordinate(i);
        compare2D = test->getElementGLeft()
                        ->getReferenceGeometry()
                        ->getReferenceNodeCoordinate(leftIndices[i]);
        compare2Dphys =
            test->getElementGLeft()->referenceToPhysical((compare2D));
        point2D = test->mapRefFaceToRefElemL((orig2D));
        point2Dphys = test->referenceToPhysical((orig2D));
        INFO("getElementGLeft or localFaceNumberLeft or mapRefFaceToRefElemL");
        CHECK((std::abs(compare2D[0] - point2D[0]) < 1e-12));
        INFO("getElementGLeft or localFaceNumberLeft or mapRefFaceToRefElemL");
        CHECK((std::abs(compare2D[1] - point2D[1]) < 1e-12));
        INFO("referenceToPhysical");
        CHECK((std::abs(compare2Dphys[0] - point2Dphys[0]) < 1e-12));
        INFO("referenceToPhysical");
        CHECK((std::abs(compare2Dphys[1] - point2Dphys[1]) < 1e-12));
        point2D = test->refFaceToRefElemMapL()->transform((orig2D));
        INFO("refFaceToRefElemMapL");
        CHECK((std::abs(compare2D[0] - point2D[0]) < 1e-12));
        INFO("refFaceToRefElemMapL");
        CHECK((std::abs(compare2D[1] - point2D[1]) < 1e-12));

        compare2D = test->getPtrElementGRight()
                        ->getReferenceGeometry()
                        ->getReferenceNodeCoordinate(rightIndices[i]);
        compare2Dphys =
            test->getPtrElementGRight()->referenceToPhysical((compare2D));
        point2D = test->mapRefFaceToRefElemR((orig2D));
        INFO(
            "getPtrElementGRight or localFaceNumberRight or "
            "mapRefFaceToRefElemR or mapRefFaceToRefFace");
        CHECK(std::abs(compare2D[0] - point2D[0]) < 1e-12);

        INFO(
            "getPtrElementGRight or localFaceNumberRight or "
            "mapRefFaceToRefElemR or mapRefFaceToRefFace");
        CHECK(std::abs(compare2D[1] - point2D[1]) < 1e-12);

        INFO("referenceToPhysical");
        CHECK((std::abs(compare2Dphys[0] - point2Dphys[0]) <
               1e-12));  // probably indirectly verified already,
                         // but this is the most important feature
                         // of a face
        INFO("referenceToPhysical");
        CHECK((std::abs(compare2Dphys[1] - point2Dphys[1]) < 1e-12));
        point2D = test->refFaceToRefElemMapR()->transform((orig2D));
        INFO("refFaceToRefElemMapR");
        CHECK((std::abs(compare2D[0] - point2D[0]) < 1e-12));
        INFO("refFaceToRefElemMapR");
        CHECK((std::abs(compare2D[1] - point2D[1]) < 1e-12));
    }

    const Geometry::ReferenceGeometry& refGeo2 = *test->getReferenceGeometry();
    INFO("Internal face type");
    CHECK((test->getFaceType() == Geometry::FaceType::INTERNAL));
    INFO("Reference geometry line");
    CHECK((typeid(refGeo2) == typeid(Geometry::ReferenceLine)));

    delete left;
    delete test;
    delete right;
    left = new Geometry::ElementGeometry(pointIndexes, nodes2D);

    test = new Geometry::FaceGeometry(left, 1, Geometry::FaceType::WALL_BC);

    leftIndices.resize(2);
    rightIndices.resize(2);

    ///\todo test the normal vector
    ///\todo also test individual functions

    leftIndices =
        test->getElementGLeft()->getPhysicalGeometry()->getLocalFaceNodeIndices(
            test->localFaceNumberLeft());

    for (std::size_t i = 0;
         i < test->getReferenceGeometry()->getNumberOfNodes(); ++i) {
        orig2D = test->getReferenceGeometry()->getReferenceNodeCoordinate(i);
        compare2D = test->getElementGLeft()
                        ->getReferenceGeometry()
                        ->getReferenceNodeCoordinate(leftIndices[i]);
        compare2Dphys =
            test->getElementGLeft()->referenceToPhysical((compare2D));
        point2D = test->mapRefFaceToRefElemL((orig2D));
        point2Dphys = test->referenceToPhysical((orig2D));
        INFO("getElementGLeft or localFaceNumberLeft or mapRefFaceToRefElemL");
        CHECK((std::abs(compare2D[0] - point2D[0]) < 1e-12));
        INFO("getElementGLeft or localFaceNumberLeft or mapRefFaceToRefElemL");
        CHECK((std::abs(compare2D[1] - point2D[1]) < 1e-12));
        INFO("referenceToPhysical");
        CHECK((std::abs(compare2Dphys[0] - point2Dphys[0]) < 1e-12));
        INFO("referenceToPhysical");
        CHECK((std::abs(compare2Dphys[1] - point2Dphys[1]) < 1e-12));
        point2D = test->refFaceToRefElemMapL()->transform((orig2D));
        INFO("refFaceToRefElemMapL");
        CHECK((std::abs(compare2D[0] - point2D[0]) < 1e-12));
        INFO("refFaceToRefElemMapL");
        CHECK((std::abs(compare2D[1] - point2D[1]) < 1e-12));
    }

    const Geometry::ReferenceGeometry& refGeo3 = *test->getReferenceGeometry();
    INFO("FaceType wall");
    CHECK((test->getFaceType() == Geometry::FaceType::WALL_BC));
    INFO("Reference geometry line");
    CHECK((typeid(refGeo3) == typeid(Geometry::ReferenceLine)));
    INFO("No right element");
    CHECK((test->getPtrElementGRight() == nullptr));
}

TEST_CASE("Mapping left-right", "[230FaceGeometry_UnitTest]") {
    // Test case to check the mapping between the left and right face/element
    // coordinates. This checks that the two ways to map a face coordinate to a
    // physical coordinate commute, that is, that for any reference coordinate
    // the following two mappings give the same physical coordinate:
    //  - ref. face -> ref. left element -> physical left element
    //  - ref. face -> right ref. face -> ref. right element
    //         -> physical right element
    // (Note: This is only the case for non-periodic faces)
    //
    // Note that the face coordinate with respect to the right element may not
    // be the same as that for the left element (which is used as coordinate
    // system for the face.). This is due to the reference face having fixed
    // face ordering, which may not agree with that of the mesh.
    //
    // Testing occurs by iterating over the 6 node configurations for the face
    // of the right element. For each configuration we check for several face
    // reference coordinates that the mappings commute.

    using namespace Geometry;
    const auto& refTet = ReferenceTetrahedron::Instance();
    std::vector<PointPhysical<3>> coords({
        // Three face nodes
        {-.5, -0.25, 0},
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, -1.0},  // Non face node for left tetrahedron
        {1, 1, 1}      // Non face node for the right tetrahedron
    });

    // The indices in coords for the nodes that form the left/right element
    // Note: The order matters, as the i-th entry is for the i-the reference
    // node.
    // Note: Technically the indices are for the nodes, but those are identical
    // to the coords here, as there are no (hypothetical) nodes with multiple
    // coordinates.
    std::vector<std::size_t> leftNodes, rightNodes;
    leftNodes = {0, 1, 2, 3};
    rightNodes = {0, 1, 2, 4};

    // The first three entries in left/rightNodes are shared. These correspond
    // to face 2.
    const std::size_t FACE_NUMBER = 2;

    ElementGeometry left(leftNodes, coords);

    // Storage for the coordinate indices on the face of the left/right element.
    // That is, the i-th entry is the coordinate index of the i-th node on the
    // FACE_NUMBER-face of the left/right element. For the left face this is
    // stable, for the right face this changes as we permute the order of the
    // reference nodes.
    std::vector<std::size_t> leftFaceNodes(3), rightFaceNodes(3);
    // Ordering of the nodes on FACE_NUMBER
    const auto rightFaceOrdering =
        refTet.getCodim1EntityLocalIndices(FACE_NUMBER);
    // Permute over the configurations for the right nodes
    do {
        // Lookup the node indices
        for (std::size_t i = 0; i < rightFaceOrdering.size(); ++i) {
            rightFaceNodes[i] = rightNodes[rightFaceOrdering[i]];
            leftFaceNodes[i] = leftNodes[rightFaceOrdering[i]];
        }
        // Construct right, face geometry from scratch every step to prevent
        // contamination from the previous test.
        ElementGeometry right(rightNodes, coords);
        FaceGeometry face(&left, FACE_NUMBER, &right, FACE_NUMBER);
        face.initialiseFaceToFaceMapIndex(leftFaceNodes, rightFaceNodes);
        const auto& faceGeom = *refTet.getCodim1ReferenceGeometry(FACE_NUMBER);
        // Check the mapping for each of the corners of the face.
        for (std::size_t i = 0; i < faceGeom.getNumberOfNodes(); ++i) {
            const auto& faceCoord =
                faceGeom.getReferenceNodeCoordinate(i).castDimension<2>();
            auto leftCoord =
                left.referenceToPhysical(face.mapRefFaceToRefElemL(faceCoord));
            auto rightCoord =
                right.referenceToPhysical(face.mapRefFaceToRefElemR(faceCoord));
            CHECK((leftCoord - rightCoord).l2NormSquared() < 1e-8);
        }
        // Next permutation of the face nodes. Note: --, to exclude the last
        // rightNode from the permutation as this is the node that is not on the
        // face.
    } while (std::next_permutation(rightNodes.begin(), --rightNodes.end()));
}
