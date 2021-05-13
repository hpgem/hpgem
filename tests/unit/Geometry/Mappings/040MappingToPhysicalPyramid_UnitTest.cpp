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
#include "Geometry/Mappings/MappingToPhysPyramid.h"
#include "Logger.h"

#include "Geometry/ReferencePyramid.h"
#include "Geometry/PhysicalGeometry.h"
#include "Geometry/PointPhysical.h"
#include "Geometry/Jacobian.h"
#include "Geometry/PointReference.h"
#include <cmath>
#include <Base/L2Norm.h>

#include "../catch.hpp"

// transformations should map internal points to internal points, external
// points to external points and nodes to nodes so construct the physical
// geometries such that this can be checked :(
using namespace hpgem;
bool isInternal3D(const Geometry::PointPhysical<3>& p) {
    return (p[1] - p[0]) > (1.) && (p[1] + p[0]) < (7. - p[2] * 0.65) &&
           (2.4 * p[0] - 0.8 * p[1]) > (1.44 + 0.12 * p[2]) &&
           (-1.5 * p[0] + 1.1 * p[1]) > (.42) && p[2] > 0;
}

TEST_CASE("040MappingToPhysicalPyramid_UnitTest",
          "[040MappingToPhysicalPyramid_UnitTest]") {

    std::vector<std::size_t> pointIndexes;

    pointIndexes.push_back(27);
    pointIndexes.push_back(4);
    pointIndexes.push_back(7);
    pointIndexes.push_back(12);

    std::vector<Geometry::PointPhysical<3> > nodes3D;

    Geometry::PointPhysical<3> point3D, compare3D;
    Geometry::PointReference<3> refPoint3D;

    pointIndexes.push_back(13);

    for (double i = 0.; i < 10; ++i) {
        point3D[0] = 1. + i / 10.;
        point3D[1] = 2. + i / 10.;
        point3D[2] = 0.;
        nodes3D.push_back(point3D);
    }
    for (double i = 0.; i < 10; ++i) {
        point3D[0] = 2. + i / 10.;
        point3D[1] = 5. - i / 10.;
        point3D[2] = 0.;
        nodes3D.push_back(point3D);
    }
    for (double i = 0.; i < 10; ++i) {
        point3D[0] = 1. + i / 10.;
        point3D[1] = 2. + i / 10.;
        point3D[2] = 4.;
        nodes3D.push_back(point3D);
    }
    for (double i = 0.; i < 10; ++i) {
        point3D[0] = 2. + i / 10.;
        point3D[1] = 5. - i / 10.;
        point3D[2] = 4.;
        nodes3D.push_back(point3D);
    }

    Geometry::ReferencePyramid& rGeom3D =
        Geometry::ReferencePyramid::Instance();

    Geometry::PhysicalGeometry<3> oops3D {pointIndexes, nodes3D, &rGeom3D};
    pointIndexes[4] = 18;
    Geometry::PhysicalGeometry<3> pGeom3D {pointIndexes, nodes3D, &rGeom3D};

    Geometry::MappingToPhysPyramid mapping3D(&pGeom3D), reinit3D(&oops3D);
    reinit3D.reinit();

    Geometry::Jacobian<3, 3> jac3D;

    for (refPoint3D[0] = -1.5189; refPoint3D[0] < 1.541; refPoint3D[0] += 0.4) {
        for (refPoint3D[1] = -1.5188; refPoint3D[1] < 1.541;
             refPoint3D[1] += 0.4) {
            for (refPoint3D[2] = -1.5188; refPoint3D[2] < 1.541;
                 refPoint3D[2] += 0.4) {
                point3D = mapping3D.transform((refPoint3D));
                if (rGeom3D.isInternalPoint(
                        (refPoint3D))) {  // not perfect, but the degenerate
                                          // cube face makes the mapping less
                                          // linear than desired
                    INFO("transform");
                    CHECK(isInternal3D(point3D));
                }
                point3D = reinit3D.transform((refPoint3D));
                if (rGeom3D.isInternalPoint((refPoint3D))) {
                    INFO("reinit");
                    CHECK((isInternal3D(point3D)));
                }

                refPoint3D[0] += -1.e-8;
                compare3D = mapping3D.transform((refPoint3D));
                refPoint3D[0] += 2.e-8;
                point3D = mapping3D.transform((refPoint3D));

                refPoint3D[0] += -1e-8;
                jac3D = mapping3D.calcJacobian((refPoint3D));
                INFO("jacobian");
                CHECK((std::abs(jac3D[0] - 5.e7 * (point3D[0] - compare3D[0])) <
                       1e-5));  // for most mappings
                INFO("jacobian");
                CHECK((std::abs(jac3D[1] - 5.e7 * (point3D[1] - compare3D[1])) <
                       1e-5));  // be more accurate
                INFO("jacobian");
                CHECK((std::abs(jac3D[2] - 5.e7 * (point3D[2] - compare3D[2])) <
                       1e-5));
                refPoint3D[1] += -1.e-8;
                compare3D = mapping3D.transform((refPoint3D));
                refPoint3D[1] += 2.e-8;
                point3D = mapping3D.transform((refPoint3D));

                refPoint3D[1] += -1e-8;
                jac3D = mapping3D.calcJacobian((refPoint3D));
                INFO("jacobian");
                CHECK((std::abs(jac3D[3] - 5.e7 * (point3D[0] - compare3D[0])) <
                       1e-5));
                INFO("jacobian");
                CHECK((std::abs(jac3D[4] - 5.e7 * (point3D[1] - compare3D[1])) <
                       1e-5));
                INFO("jacobian");
                CHECK((std::abs(jac3D[5] - 5.e7 * (point3D[2] - compare3D[2])) <
                       1e-5));
                refPoint3D[2] += -1.e-8;
                compare3D = mapping3D.transform((refPoint3D));
                refPoint3D[2] += 2.e-8;
                point3D = mapping3D.transform((refPoint3D));

                refPoint3D[2] += -1e-8;
                jac3D = mapping3D.calcJacobian((refPoint3D));
                INFO("jacobian");
                CHECK((std::abs(jac3D[6] - 5.e7 * (point3D[0] - compare3D[0])) <
                       1e-5));
                INFO("jacobian");
                CHECK((std::abs(jac3D[7] - 5.e7 * (point3D[1] - compare3D[1])) <
                       1e-5));
                INFO("jacobian");
                CHECK((std::abs(jac3D[8] - 5.e7 * (point3D[2] - compare3D[2])) <
                       1e-5));
                refPoint3D[2] += 1e-8;
                // the reference point is outside the square
                //(but its image after transformation and inverse transformation
                // might not be; due to nonlinearities) or they are inside and
                // on
                // the same location
                double dist = Base::L2Norm(refPoint3D -
                                           mapping3D.inverseTransform(point3D));

                INFO("inverse transformation, (distance is "
                     << dist << ", point is " << refPoint3D << "/"
                     << mapping3D.inverseTransform(point3D) << ")");

                bool check =
                    (!rGeom3D.isInternalPoint(refPoint3D)) || dist < 1e-12;
                CHECK(check);
            }
        }
    }

    for (std::size_t i = 0; i < rGeom3D.getNumberOfNodes(); ++i) {
        refPoint3D = rGeom3D.getReferenceNodeCoordinate(i);
        compare3D = pGeom3D.getLocalNodeCoordinates(i);
        point3D = mapping3D.transform((refPoint3D));
        INFO("transform");
        CHECK(std::abs(point3D[0] - compare3D[0]) < 1e-12);
        CHECK(std::abs(point3D[1] - compare3D[1]) < 1e-12);
        CHECK(std::abs(point3D[2] - compare3D[2]) < 1e-12);
    }

    INFO("getTargetDimension");
    CHECK((mapping3D.getTargetDimension() == 3));
}
