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
#include "Geometry/PointPhysical.h"
#include "Geometry/PointReference.h"
#include "Geometry/ElementGeometry.h"
#include "Geometry/Jacobian.h"
#include "Geometry/PhysicalGeometry.h"

#include <cmath>
using namespace hpgem;
int main() {

    // dim0

    std::vector<std::size_t> pointIndexes, leftIndices, rightIndices;
    std::vector<Geometry::PointPhysical<1> > nodes1D;

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
        logger.assert_always(
            (std::abs(compare1D[0] - point1D[0]) < 1e-12),
            "getElementGLeft or localFaceNumberLeft or mapRefFaceToRefElemL");
        logger.assert_always(
            (std::abs(compare1Dphys[0] - point1Dphys[0]) < 1e-12),
            "referenceToPhysical");
        point1D = test->refFaceToRefElemMapL()->transform((orig1D));
        logger.assert_always((std::abs(compare1D[0] - point1D[0]) < 1e-12),
                             "refFaceToRefElemMapL");

        compare1D = test->getPtrElementGRight()
                        ->getReferenceGeometry()
                        ->getReferenceNodeCoordinate(rightIndices[i]);
        compare1Dphys =
            test->getPtrElementGRight()->referenceToPhysical((compare1D));
        point1D = test->mapRefFaceToRefElemR((orig1D));
        logger.assert_always((std::abs(compare1D[0] - point1D[0]) < 1e-12),
                             "getPtrElementGRight or localFaceNumberRight or "
                             "mapRefFaceToRefElemR or mapRefFaceToRefFace");
        logger.assert_always(
            (std::abs(compare1Dphys[0] - point1Dphys[0]) < 1e-12),
            "referenceToPhysical");  // probably indirectly verified already,
                                     // but this is the most important feature
                                     // of a face
        point1D = test->refFaceToRefElemMapR()->transform((orig1D));
        logger.assert_always((std::abs(compare1D[0] - point1D[0]) < 1e-12),
                             "refFaceToRefElemMapR");
    }

    const Geometry::ReferenceGeometry& refGeo1 = *test->getReferenceGeometry();
    logger.assert_always((test->getFaceType() == Geometry::FaceType::INTERNAL),
                         "Face type internal");
    logger.assert_always((typeid(refGeo1) == typeid(Geometry::ReferencePoint)),
                         "Reference geometry point");

    // dim 1

    std::vector<Geometry::PointPhysical<2> > nodes2D;

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
        logger.assert_always(
            (std::abs(compare2D[0] - point2D[0]) < 1e-12),
            "getElementGLeft or localFaceNumberLeft or mapRefFaceToRefElemL");
        logger.assert_always(
            (std::abs(compare2D[1] - point2D[1]) < 1e-12),
            "getElementGLeft or localFaceNumberLeft or mapRefFaceToRefElemL");
        logger.assert_always(
            (std::abs(compare2Dphys[0] - point2Dphys[0]) < 1e-12),
            "referenceToPhysical");
        logger.assert_always(
            (std::abs(compare2Dphys[1] - point2Dphys[1]) < 1e-12),
            "referenceToPhysical");
        point2D = test->refFaceToRefElemMapL()->transform((orig2D));
        logger.assert_always((std::abs(compare2D[0] - point2D[0]) < 1e-12),
                             "refFaceToRefElemMapL");
        logger.assert_always((std::abs(compare2D[1] - point2D[1]) < 1e-12),
                             "refFaceToRefElemMapL");

        compare2D = test->getPtrElementGRight()
                        ->getReferenceGeometry()
                        ->getReferenceNodeCoordinate(rightIndices[i]);
        compare2Dphys =
            test->getPtrElementGRight()->referenceToPhysical((compare2D));
        point2D = test->mapRefFaceToRefElemR((orig2D));
        logger.assert_always((std::abs(compare2D[0] - point2D[0]) < 1e-12),
                             "getPtrElementGRight or localFaceNumberRight or "
                             "mapRefFaceToRefElemR or mapRefFaceToRefFace");
        logger.assert_always((std::abs(compare2D[1] - point2D[1]) < 1e-12),
                             "getPtrElementGRight or localFaceNumberRight or "
                             "mapRefFaceToRefElemR or mapRefFaceToRefFace");
        logger.assert_always(
            (std::abs(compare2Dphys[0] - point2Dphys[0]) < 1e-12),
            "referenceToPhysical");  // probably indirectly verified already,
                                     // but this is the most important feature
                                     // of a face
        logger.assert_always(
            (std::abs(compare2Dphys[1] - point2Dphys[1]) < 1e-12),
            "referenceToPhysical");
        point2D = test->refFaceToRefElemMapR()->transform((orig2D));
        logger.assert_always((std::abs(compare2D[0] - point2D[0]) < 1e-12),
                             "refFaceToRefElemMapR");
        logger.assert_always((std::abs(compare2D[1] - point2D[1]) < 1e-12),
                             "refFaceToRefElemMapR");
    }

    const Geometry::ReferenceGeometry& refGeo2 = *test->getReferenceGeometry();
    logger.assert_always((test->getFaceType() == Geometry::FaceType::INTERNAL),
                         "Internal face type");
    logger.assert_always((typeid(refGeo2) == typeid(Geometry::ReferenceLine)),
                         "Reference geometry line");

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
        logger.assert_always(
            (std::abs(compare2D[0] - point2D[0]) < 1e-12),
            "getElementGLeft or localFaceNumberLeft or mapRefFaceToRefElemL");
        logger.assert_always(
            (std::abs(compare2D[1] - point2D[1]) < 1e-12),
            "getElementGLeft or localFaceNumberLeft or mapRefFaceToRefElemL");
        logger.assert_always(
            (std::abs(compare2Dphys[0] - point2Dphys[0]) < 1e-12),
            "referenceToPhysical");
        logger.assert_always(
            (std::abs(compare2Dphys[1] - point2Dphys[1]) < 1e-12),
            "referenceToPhysical");
        point2D = test->refFaceToRefElemMapL()->transform((orig2D));
        logger.assert_always((std::abs(compare2D[0] - point2D[0]) < 1e-12),
                             "refFaceToRefElemMapL");
        logger.assert_always((std::abs(compare2D[1] - point2D[1]) < 1e-12),
                             "refFaceToRefElemMapL");
    }

    const Geometry::ReferenceGeometry& refGeo3 = *test->getReferenceGeometry();
    logger.assert_always((test->getFaceType() == Geometry::FaceType::WALL_BC),
                         "FaceType wall");
    logger.assert_always((typeid(refGeo3) == typeid(Geometry::ReferenceLine)),
                         "Reference geometry line");
    logger.assert_always((test->getPtrElementGRight() == nullptr),
                         "No right element");
}
