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
// unit test indicate the culprit class and other 'unit' tests may assume correct
// execution of all prior unit tests
#include "Geometry/Mappings/ConcatenatedMapping.h"
#include "Geometry/Mappings/MappingToRefPointToPoint.h"
#include "Geometry/Mappings/MappingToRefPointToLine.h"
#include "Geometry/Mappings/MappingToRefLineToLine.h"
#include "Geometry/Mappings/MappingToRefLineToSquare.h"
#include "Geometry/Mappings/MappingToRefSquareToSquare.h"
#include "Geometry/Mappings/MappingToRefSquareToCube.h"
#include "Geometry/Mappings/MappingToRefCubeToCube.h"
#include "Geometry/Mappings/MappingToRefCubeToHypercube.h"
#include "Logger.h"

#include "Geometry/ReferencePoint.h"
#include "Geometry/ReferenceLine.h"
#include "Geometry/ReferenceSquare.h"
#include "Geometry/ReferenceCube.h"
#include "Geometry/ReferenceHypercube.h"
#include "Geometry/PointReference.h"
#include "Geometry/Jacobian.h"
#include "LinearAlgebra/MiddleSizeVector.h"
#include <cmath>

int main() {

    Geometry::PointReference<0> point0D, compare0D, orig0D;
    Geometry::PointReference<1> point1D, compare1D, orig1D;
    Geometry::PointReference<2> point2D, compare2D, orig2D;
    Geometry::PointReference<3> point3D, compare3D, orig3D;
    Geometry::PointReference<4> point4D, compare4D, orig4D;
    Geometry::ReferenceGeometry* source;
    Geometry::ReferenceGeometry* target;

    Geometry::Jacobian<0, 1> jac;

    std::vector<std::size_t> nodesAfterTransformation(1);

    Geometry::ConcatenatedMapping* test;

    // 0->0->1

    source = &Geometry::ReferencePoint::Instance();
    target = &Geometry::ReferenceLine::Instance();
    test = new Geometry::ConcatenatedMapping(
        Geometry::MappingToRefPointToPoint::Instance(),
        Geometry::MappingToRefPointToLine0::Instance());

    point1D = test->transform((orig0D));
    logger.assert_always((source->isInternalPoint((orig0D)) ==
                          target->isInternalPoint((point1D))),
                         "transform");

    jac = test->calcJacobian((orig0D));

    for (std::size_t i = 0; i < source->getNumberOfNodes(); ++i) {
        orig0D = source->getReferenceNodeCoordinate(i);
        compare1D =
            target->getReferenceNodeCoordinate(nodesAfterTransformation[i]);
        point1D = test->transform((orig0D));
    }

    logger.assert_always((test->getTargetDimension() == 1),
                         "getTargetDimension");

    // 1->1->2

    source = &Geometry::ReferenceLine::Instance();
    target = &Geometry::ReferenceSquare::Instance();
    Geometry::Jacobian<1, 2> jac2;
    delete test;
    test = new Geometry::ConcatenatedMapping(
        Geometry::MappingToRefLineToLine1::Instance(),
        Geometry::MappingToRefLineToSquare3::Instance());
    nodesAfterTransformation.resize(2);
    nodesAfterTransformation[0] = 3;
    nodesAfterTransformation[1] = 2;

    for (orig1D[0] = -2.8189; orig1D[0] < 3.141; orig1D[0] += 0.1) {
        point2D = test->transform((orig1D));
        logger.assert_always((source->isInternalPoint((orig1D)) ==
                              target->isInternalPoint((point2D))),
                             "transform");

        orig1D[0] += -1.e-8;
        compare2D = test->transform((orig1D));
        orig1D[0] += 2.e-8;
        point2D = test->transform((orig1D));

        orig1D[0] += -1e-8;
        jac2 = test->calcJacobian((orig1D));
        logger.assert_always(
            (std::abs(jac2[0] - 5.e7 * (point2D[0] - compare2D[0])) < 1e-5),
            "jacobian");  // estimate is a bit rough, but should work for most
                          // mappings
        logger.assert_always(
            (std::abs(jac2[1] - 5.e7 * (point2D[1] - compare2D[1])) < 1e-5),
            "jacobian");  // implementations are very strongly recommended to be
                          // more accurate
    }

    for (std::size_t i = 0; i < source->getNumberOfNodes(); ++i) {
        orig1D = source->getReferenceNodeCoordinate(i);
        compare2D =
            target->getReferenceNodeCoordinate(nodesAfterTransformation[i]);
        point2D = test->transform((orig1D));
        logger.assert_always((std::abs(point2D[0] - compare2D[0]) < 1e-12),
                             "transform");
        logger.assert_always((std::abs(point2D[1] - compare2D[1]) < 1e-12),
                             "transform");
    }

    logger.assert_always((test->getTargetDimension() == 2),
                         "getTargetDimension");

    // 2->2->3

    source = &Geometry::ReferenceSquare::Instance();
    target = &Geometry::ReferenceCube::Instance();
    Geometry::Jacobian<2, 3> jac3;
    delete test;
    test = new Geometry::ConcatenatedMapping(
        Geometry::MappingToRefSquareToSquare3::Instance(),
        Geometry::MappingToRefSquareToCube2::Instance());
    nodesAfterTransformation.resize(4);
    nodesAfterTransformation[0] = 4;
    nodesAfterTransformation[1] = 0;
    nodesAfterTransformation[2] = 6;
    nodesAfterTransformation[3] = 2;

    for (orig2D[0] = -1.51; orig2D[0] < 1.51; orig2D[0] += 0.2) {
        for (orig2D[1] = -1.51; orig2D[1] < 1.51; orig2D[1] += 0.2) {
            point3D = test->transform((orig2D));
            logger.assert_always((source->isInternalPoint((orig2D)) ==
                                  target->isInternalPoint((point3D))),
                                 "transform");

            orig2D[0] += -1.e-8;
            compare3D = test->transform((orig2D));
            orig2D[0] += 2.e-8;
            point3D = test->transform((orig2D));

            orig2D[0] += -1e-8;
            jac3 = test->calcJacobian((orig2D));
            logger.assert_always(
                (std::abs(jac3[0] - 5.e7 * (point3D[0] - compare3D[0])) < 1e-5),
                "jacobian");  // estimate is a bit rough, but should work for
                              // most mappings
            logger.assert_always(
                (std::abs(jac3[1] - 5.e7 * (point3D[1] - compare3D[1])) < 1e-5),
                "jacobian");  // implementations are very strongly recommended
                              // to be more accurate
            logger.assert_always(
                (std::abs(jac3[2] - 5.e7 * (point3D[2] - compare3D[2])) < 1e-5),
                "jacobian");

            orig2D[1] += -1.e-8;
            compare3D = test->transform((orig2D));
            orig2D[1] += 2.e-8;
            point3D = test->transform((orig2D));

            orig2D[1] += -1e-8;
            jac3 = test->calcJacobian((orig2D));
            logger.assert_always(
                (std::abs(jac3[3] - 5.e7 * (point3D[0] - compare3D[0])) < 1e-5),
                "jacobian");
            logger.assert_always(
                (std::abs(jac3[4] - 5.e7 * (point3D[1] - compare3D[1])) < 1e-5),
                "jacobian");
            logger.assert_always(
                (std::abs(jac3[5] - 5.e7 * (point3D[2] - compare3D[2])) < 1e-5),
                "jacobian");
        }
    }

    for (std::size_t i = 0; i < source->getNumberOfNodes(); ++i) {
        orig2D = source->getReferenceNodeCoordinate(i);
        compare3D =
            target->getReferenceNodeCoordinate(nodesAfterTransformation[i]);
        point3D = test->transform((orig2D));
        logger.assert_always((std::abs(point3D[0] - compare3D[0]) < 1e-12),
                             "transform");
        logger.assert_always((std::abs(point3D[1] - compare3D[1]) < 1e-12),
                             "transform");
        logger.assert_always((std::abs(point3D[2] - compare3D[2]) < 1e-12),
                             "transform");
    }

    logger.assert_always((test->getTargetDimension() == 3),
                         "getTargetDimension");

    // 3->3->4

    source = &Geometry::ReferenceCube::Instance();
    target = &Geometry::ReferenceHypercube::Instance();
    Geometry::Jacobian<3, 4> jac4;
    delete test;
    test = new Geometry::ConcatenatedMapping(
        Geometry::MappingToRefCubeToCube3::Instance(),
        Geometry::MappingToRefCubeToHypercube7::Instance());
    nodesAfterTransformation.resize(8);
    nodesAfterTransformation[0] = 12;
    nodesAfterTransformation[1] = 13;
    nodesAfterTransformation[2] = 8;
    nodesAfterTransformation[3] = 9;
    nodesAfterTransformation[4] = 14;
    nodesAfterTransformation[5] = 15;
    nodesAfterTransformation[6] = 10;
    nodesAfterTransformation[7] = 11;

    for (orig3D[0] = -1.51; orig3D[0] < 1.51; orig3D[0] += 0.6) {
        for (orig3D[1] = -1.51; orig3D[1] < 1.51; orig3D[1] += 0.7) {
            for (orig3D[2] = -1.51; orig3D[2] < 1.51; orig3D[2] += 0.8) {
                point4D = test->transform((orig3D));
                logger.assert_always((source->isInternalPoint((orig3D)) ==
                                      target->isInternalPoint((point4D))),
                                     "transform");

                orig3D[0] += -1.e-8;
                compare4D = test->transform((orig3D));
                orig3D[0] += 2.e-8;
                point4D = test->transform((orig3D));

                orig3D[0] += -1e-8;
                jac4 = test->calcJacobian((orig3D));
                logger.assert_always(
                    (std::abs(jac4[0] - 5.e7 * (point4D[0] - compare4D[0])) <
                     1e-5),
                    "jacobian");  // estimate is a bit rough, but should work
                                  // for most mappings
                logger.assert_always(
                    (std::abs(jac4[1] - 5.e7 * (point4D[1] - compare4D[1])) <
                     1e-5),
                    "jacobian");  // implementations are very strongly
                                  // recommended to be more accurate
                logger.assert_always(
                    (std::abs(jac4[2] - 5.e7 * (point4D[2] - compare4D[2])) <
                     1e-5),
                    "jacobian");
                logger.assert_always(
                    (std::abs(jac4[3] - 5.e7 * (point4D[3] - compare4D[3])) <
                     1e-5),
                    "jacobian");

                orig3D[1] += -1.e-8;
                compare4D = test->transform((orig3D));
                orig3D[1] += 2.e-8;
                point4D = test->transform((orig3D));

                orig3D[1] += -1e-8;
                jac4 = test->calcJacobian((orig3D));
                logger.assert_always(
                    (std::abs(jac4[4] - 5.e7 * (point4D[0] - compare4D[0])) <
                     1e-5),
                    "jacobian");
                logger.assert_always(
                    (std::abs(jac4[5] - 5.e7 * (point4D[1] - compare4D[1])) <
                     1e-5),
                    "jacobian");
                logger.assert_always(
                    (std::abs(jac4[6] - 5.e7 * (point4D[2] - compare4D[2])) <
                     1e-5),
                    "jacobian");
                logger.assert_always(
                    (std::abs(jac4[7] - 5.e7 * (point4D[3] - compare4D[3])) <
                     1e-5),
                    "jacobian");

                orig3D[2] += -1.e-8;
                compare4D = test->transform((orig3D));
                orig3D[2] += 2.e-8;
                point4D = test->transform((orig3D));

                orig3D[2] += -1e-8;
                jac4 = test->calcJacobian((orig3D));
                logger.assert_always(
                    (std::abs(jac4[8] - 5.e7 * (point4D[0] - compare4D[0])) <
                     1e-5),
                    "jacobian");
                logger.assert_always(
                    (std::abs(jac4[9] - 5.e7 * (point4D[1] - compare4D[1])) <
                     1e-5),
                    "jacobian");
                logger.assert_always(
                    (std::abs(jac4[10] - 5.e7 * (point4D[2] - compare4D[2])) <
                     1e-5),
                    "jacobian");
                logger.assert_always(
                    (std::abs(jac4[11] - 5.e7 * (point4D[3] - compare4D[3])) <
                     1e-5),
                    "jacobian");
            }
        }
    }

    for (std::size_t i = 0; i < source->getNumberOfNodes(); ++i) {
        orig3D = source->getReferenceNodeCoordinate(i);
        compare4D =
            target->getReferenceNodeCoordinate(nodesAfterTransformation[i]);
        point4D = test->transform((orig3D));
        logger.assert_always((std::abs(point4D[0] - compare4D[0]) < 1e-12),
                             "transform");
        logger.assert_always((std::abs(point4D[1] - compare4D[1]) < 1e-12),
                             "transform");
        logger.assert_always((std::abs(point4D[2] - compare4D[2]) < 1e-12),
                             "transform");
        logger.assert_always((std::abs(point4D[3] - compare4D[3]) < 1e-12),
                             "transform");
    }

    logger.assert_always((test->getTargetDimension() == 4),
                         "getTargetDimension");
    return 0;
}
