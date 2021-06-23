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
#include "Geometry/ElementGeometry.h"
#include <Logger.h>

#include "Geometry/ReferenceCube.h"
#include "Geometry/ReferenceHypercube.h"
#include "Geometry/ReferenceLine.h"
#include "Geometry/ReferencePyramid.h"
#include "Geometry/ReferenceSquare.h"
#include "Geometry/ReferenceTetrahedron.h"
#include "Geometry/ReferenceTriangle.h"
#include "Geometry/ReferenceTriangularPrism.h"
#include "Geometry/Mappings/MappingToPhysHypercubeLinear.h"
#include "Geometry/Mappings/MappingToPhysPyramid.h"
#include "Geometry/Mappings/MappingToPhysSimplexLinear.h"
#include "Geometry/Mappings/MappingToPhysTriangularPrism.h"

#include <cmath>
#include <typeinfo>

#include "../catch.hpp"

using namespace hpgem;
TEST_CASE("220ElementGeometry_UnitTest", "[220ElementGeometry_UnitTest]") {

    // dim 1

    std::vector<std::size_t> pointIndexes;
    std::vector<Geometry::PointPhysical<1> > nodes1D;

    Geometry::PointPhysical<1> point1D, compare1D;
    Geometry::PointReference<1> orig1D;

    Geometry::Jacobian<1, 1> jac, jaccompare;

    pointIndexes.push_back(4);
    pointIndexes.push_back(7);

    for (double i = 0.; i < 1; i += 0.1) {
        point1D[0] = 1. + i / 10.;
        nodes1D.push_back(point1D);
    }

    Geometry::ElementGeometry* test =
        new Geometry::ElementGeometry(pointIndexes, nodes1D);
    const Geometry::MappingReferenceToPhysicalBase& refMap1 =
        *test->getReferenceToPhysicalMap();
    const Geometry::ReferenceGeometry& refGeo1 = *test->getReferenceGeometry();
    INFO("getReferenceToPhysicalMap");
    CHECK(
        (typeid(Geometry::MappingToPhysHypercubeLinear<1>) == typeid(refMap1)));
    INFO("getReferenceGeometry");
    CHECK((typeid(Geometry::ReferenceLine) == typeid(refGeo1)));
    INFO("getNumberOfNodes");
    CHECK((test->getNumberOfNodes() == 2));

    for (orig1D[0] = -1.51; orig1D[0] < 1.51; orig1D[0] += 0.1) {
        compare1D = test->getReferenceToPhysicalMap()->transform((orig1D));
        point1D = test->referenceToPhysical((orig1D));
        INFO("referenceToPhysical");
        CHECK((std::abs(point1D[0] - compare1D[0]) < 1e-12));

        jaccompare = test->getReferenceToPhysicalMap()->calcJacobian((orig1D));
        jac = test->calcJacobian((orig1D));
        INFO("calcJacobian");
        CHECK((std::abs(jac[0] - jaccompare[0]) < 1e-12));
    }
    logger(INFO, "%", *test);

    // dim2
    std::vector<Geometry::PointPhysical<2> > nodes2D;

    Geometry::PointPhysical<2> point2D, compare2D;
    Geometry::PointReference<2> orig2D;
    Geometry::Jacobian<2, 2> jac2, jaccompare2;

    pointIndexes.push_back(10);

    for (double i = 0.; i < 1; i += 0.1) {
        point2D[0] = 1. + i;
        point2D[1] = 2. + i;
        nodes2D.push_back(point2D);
    }

    point2D[0] = 3.5;
    point2D[1] = 4.6;
    nodes2D.push_back(point2D);
    point2D[0] = 6.7;
    point2D[1] = 2.8;
    nodes2D.push_back(point2D);

    delete test;
    test = new Geometry::ElementGeometry(pointIndexes, nodes2D);

    const Geometry::MappingReferenceToPhysicalBase& refMap2 =
        *test->getReferenceToPhysicalMap();
    const Geometry::ReferenceGeometry& refGeo2 = *test->getReferenceGeometry();
    INFO("getReferenceToPhysicalMap");
    CHECK((typeid(Geometry::MappingToPhysSimplexLinear<2>) == typeid(refMap2)));
    INFO("getReferenceGeometry");
    CHECK((typeid(Geometry::ReferenceTriangle) == typeid(refGeo2)));
    INFO("getNumberOfNodes");
    CHECK((test->getNumberOfNodes() == 3));

    for (orig2D[0] = -1.51; orig2D[0] < 1.51; orig2D[0] += 0.2) {
        for (orig2D[1] = -1.511; orig2D[1] < 1.51; orig2D[1] += 0.2) {
            compare2D = test->getReferenceToPhysicalMap()->transform((orig2D));
            point2D = test->referenceToPhysical((orig2D));
            INFO("referenceToPhysical");
            CHECK((std::abs(point2D[0] - compare2D[0]) < 1e-12));
            INFO("referenceToPhysical");
            CHECK((std::abs(point2D[1] - compare2D[1]) < 1e-12));

            jaccompare2 =
                test->getReferenceToPhysicalMap()->calcJacobian((orig2D));
            jac2 = test->calcJacobian((orig2D));
            INFO("calcJacobian");
            CHECK((std::abs(jac2[0] - jaccompare2[0]) < 1e-12));
            INFO("calcJacobian");
            CHECK((std::abs(jac2[1] - jaccompare2[1]) < 1e-12));
            INFO("calcJacobian");
            CHECK((std::abs(jac2[2] - jaccompare2[2]) < 1e-12));
            INFO("calcJacobian");
            CHECK((std::abs(jac2[3] - jaccompare2[3]) < 1e-12));
        }
    }
    logger(INFO, "%", *test);

    pointIndexes.push_back(11);

    delete test;
    test = new Geometry::ElementGeometry(pointIndexes, nodes2D);

    const Geometry::MappingReferenceToPhysicalBase& refMap3 =
        *test->getReferenceToPhysicalMap();
    const Geometry::ReferenceGeometry& refGeo3 = *test->getReferenceGeometry();
    INFO("getReferenceToPhysicalMap");
    CHECK(
        (typeid(Geometry::MappingToPhysHypercubeLinear<2>) == typeid(refMap3)));
    INFO("getReferenceGeometry");
    CHECK((typeid(Geometry::ReferenceSquare) == typeid(refGeo3)));
    INFO("getNumberOfNodes");
    CHECK((test->getNumberOfNodes() == 4));

    for (orig2D[0] = -1.51; orig2D[0] < 1.51; orig2D[0] += 0.2) {
        for (orig2D[1] = -1.511; orig2D[1] < 1.51; orig2D[1] += 0.2) {
            compare2D = test->getReferenceToPhysicalMap()->transform((orig2D));
            point2D = test->referenceToPhysical((orig2D));
            INFO("referenceToPhysical");
            CHECK((std::abs(point2D[0] - compare2D[0]) < 1e-12));
            INFO("referenceToPhysical");
            CHECK((std::abs(point2D[1] - compare2D[1]) < 1e-12));

            jaccompare2 =
                test->getReferenceToPhysicalMap()->calcJacobian((orig2D));
            jac2 = test->calcJacobian((orig2D));
            INFO("calcJacobian");
            CHECK((std::abs(jac2[0] - jaccompare2[0]) < 1e-12));
            INFO("calcJacobian");
            CHECK((std::abs(jac2[1] - jaccompare2[1]) < 1e-12));
            INFO("calcJacobian");
            CHECK((std::abs(jac2[2] - jaccompare2[2]) < 1e-12));
            INFO("calcJacobian");
            CHECK((std::abs(jac2[3] - jaccompare2[3]) < 1e-12));
        }
    }
    logger(INFO, "%", *test);

    // dim 3

    std::vector<Geometry::PointPhysical<3> > nodes3D;

    Geometry::PointPhysical<3> point3D, compare3D;
    Geometry::PointReference<3> orig3D;
    Geometry::Jacobian<3, 3> jac3, jaccompare3;

    for (double i = 0.; i < 1; i += 0.1) {
        point3D[0] = 1. + i / 10.;
        point3D[1] = 2. + i / 10.;
        point3D[2] = 3. + i / 10.;
        nodes3D.push_back(point3D);
    }

    point3D[0] = 3.5;
    point3D[1] = 4.6;
    point3D[2] = 5.4;
    nodes3D.push_back(point3D);
    point3D[0] = 6.7;
    point3D[1] = 2.8;
    point3D[2] = 5.7;
    nodes3D.push_back(point3D);
    point3D[0] = 1.4;
    point3D[1] = 2.4;
    point3D[2] = 5.4;
    nodes3D.push_back(point3D);
    point3D[0] = 1.7;
    point3D[1] = 2.7;
    point3D[2] = 5.7;
    nodes3D.push_back(point3D);
    point3D[0] = 3.5;
    point3D[1] = 4.6;
    point3D[2] = 7.4;
    nodes3D.push_back(point3D);
    point3D[0] = 6.7;
    point3D[1] = 2.8;
    point3D[2] = 7.7;
    nodes3D.push_back(point3D);

    delete test;
    test = new Geometry::ElementGeometry(pointIndexes, nodes3D);

    const Geometry::MappingReferenceToPhysicalBase& refMap4 =
        *test->getReferenceToPhysicalMap();
    const Geometry::ReferenceGeometry& refGeo4 = *test->getReferenceGeometry();
    INFO("getReferenceToPhysicalMap");
    CHECK((typeid(Geometry::MappingToPhysSimplexLinear<3>) == typeid(refMap4)));
    INFO("getReferenceGeometry");
    CHECK((typeid(Geometry::ReferenceTetrahedron) == typeid(refGeo4)));
    INFO("getNumberOfNodes");
    CHECK((test->getNumberOfNodes() == 4));

    for (orig3D[0] = -1.51; orig3D[0] < 1.51; orig3D[0] += 0.3) {
        for (orig3D[1] = -1.511; orig3D[1] < 1.51; orig3D[1] += 0.35) {
            for (orig3D[2] = -1.512; orig3D[2] < 1.51; orig3D[2] += 0.3) {
                compare3D =
                    test->getReferenceToPhysicalMap()->transform((orig3D));
                point3D = test->referenceToPhysical((orig3D));
                INFO("referenceToPhysical");
                CHECK((std::abs(point3D[0] - compare3D[0]) < 1e-12));
                INFO("referenceToPhysical");
                CHECK((std::abs(point3D[1] - compare3D[1]) < 1e-12));
                INFO("referenceToPhysical");
                CHECK((std::abs(point3D[2] - compare3D[2]) < 1e-12));

                jaccompare3 =
                    test->getReferenceToPhysicalMap()->calcJacobian((orig3D));
                jac3 = test->calcJacobian((orig3D));
                INFO("calcJacobian");
                CHECK((std::abs(jac3[0] - jaccompare3[0]) < 1e-12));
                INFO("calcJacobian");
                CHECK((std::abs(jac3[1] - jaccompare3[1]) < 1e-12));
                INFO("calcJacobian");
                CHECK((std::abs(jac3[2] - jaccompare3[2]) < 1e-12));
                INFO("calcJacobian");
                CHECK((std::abs(jac3[3] - jaccompare3[3]) < 1e-12));
                INFO("calcJacobian");
                CHECK((std::abs(jac3[4] - jaccompare3[4]) < 1e-12));
                INFO("calcJacobian");
                CHECK((std::abs(jac3[5] - jaccompare3[5]) < 1e-12));
                INFO("calcJacobian");
                CHECK((std::abs(jac3[6] - jaccompare3[6]) < 1e-12));
                INFO("calcJacobian");
                CHECK((std::abs(jac3[7] - jaccompare3[7]) < 1e-12));
                INFO("calcJacobian");
                CHECK((std::abs(jac3[8] - jaccompare3[8]) < 1e-12));
            }
        }
    }
    logger(INFO, "%", *test);

    pointIndexes.push_back(12);

    delete test;
    test = new Geometry::ElementGeometry(pointIndexes, nodes3D);

    const Geometry::MappingReferenceToPhysicalBase& refMap5 =
        *test->getReferenceToPhysicalMap();
    const Geometry::ReferenceGeometry& refGeo5 = *test->getReferenceGeometry();
    INFO("getReferenceToPhysicalMap");
    CHECK((typeid(Geometry::MappingToPhysPyramid) == typeid(refMap5)));
    INFO("getReferenceGeometry");
    CHECK((typeid(Geometry::ReferencePyramid) == typeid(refGeo5)));
    INFO("getNumberOfNodes");
    CHECK((test->getNumberOfNodes() == 5));

    for (orig3D[0] = -1.51; orig3D[0] < 1.51; orig3D[0] += 0.3) {
        for (orig3D[1] = -1.511; orig3D[1] < 1.51; orig3D[1] += 0.35) {
            for (orig3D[2] = -1.512; orig3D[2] < 1.51; orig3D[2] += 0.3) {
                compare3D =
                    test->getReferenceToPhysicalMap()->transform((orig3D));
                point3D = test->referenceToPhysical((orig3D));
                INFO("referenceToPhysical");
                CHECK((std::abs(point3D[0] - compare3D[0]) < 1e-12));
                INFO("referenceToPhysical");
                CHECK((std::abs(point3D[1] - compare3D[1]) < 1e-12));
                INFO("referenceToPhysical");
                CHECK((std::abs(point3D[2] - compare3D[2]) < 1e-12));

                jaccompare3 =
                    test->getReferenceToPhysicalMap()->calcJacobian((orig3D));
                jac3 = test->calcJacobian((orig3D));
                INFO("calcJacobian");
                CHECK((std::abs(jac3[0] - jaccompare3[0]) < 1e-12));
                INFO("calcJacobian");
                CHECK((std::abs(jac3[1] - jaccompare3[1]) < 1e-12));
                INFO("calcJacobian");
                CHECK((std::abs(jac3[2] - jaccompare3[2]) < 1e-12));
                INFO("calcJacobian");
                CHECK((std::abs(jac3[3] - jaccompare3[3]) < 1e-12));
                INFO("calcJacobian");
                CHECK((std::abs(jac3[4] - jaccompare3[4]) < 1e-12));
                INFO("calcJacobian");
                CHECK((std::abs(jac3[5] - jaccompare3[5]) < 1e-12));
                INFO("calcJacobian");
                CHECK((std::abs(jac3[6] - jaccompare3[6]) < 1e-12));
                INFO("calcJacobian");
                CHECK((std::abs(jac3[7] - jaccompare3[7]) < 1e-12));
                INFO("calcJacobian");
                CHECK((std::abs(jac3[8] - jaccompare3[8]) < 1e-12));
            }
        }
    }
    logger(INFO, "%", *test);

    pointIndexes.push_back(13);

    delete test;
    test = new Geometry::ElementGeometry(pointIndexes, nodes3D);

    const Geometry::MappingReferenceToPhysicalBase& refMap6 =
        *test->getReferenceToPhysicalMap();
    const Geometry::ReferenceGeometry& refGeo6 = *test->getReferenceGeometry();
    INFO("getReferenceToPhysicalMap");
    CHECK((typeid(Geometry::MappingToPhysTriangularPrism) == typeid(refMap6)));
    INFO("getReferenceGeometry");
    CHECK((typeid(Geometry::ReferenceTriangularPrism) == typeid(refGeo6)));
    INFO("getNumberOfNodes");
    CHECK((test->getNumberOfNodes() == 6));

    for (orig3D[0] = -1.51; orig3D[0] < 1.51; orig3D[0] += 0.3) {
        for (orig3D[1] = -1.511; orig3D[1] < 1.51; orig3D[1] += 0.35) {
            for (orig3D[2] = -1.512; orig3D[2] < 1.51; orig3D[2] += 0.3) {
                compare3D =
                    test->getReferenceToPhysicalMap()->transform((orig3D));
                point3D = test->referenceToPhysical((orig3D));
                INFO("referenceToPhysical");
                CHECK((std::abs(point3D[0] - compare3D[0]) < 1e-12));
                INFO("referenceToPhysical");
                CHECK((std::abs(point3D[1] - compare3D[1]) < 1e-12));
                INFO("referenceToPhysical");
                CHECK((std::abs(point3D[2] - compare3D[2]) < 1e-12));

                jaccompare3 =
                    test->getReferenceToPhysicalMap()->calcJacobian((orig3D));
                jac3 = test->calcJacobian((orig3D));
                INFO("calcJacobian");
                CHECK((std::abs(jac3[0] - jaccompare3[0]) < 1e-12));
                INFO("calcJacobian");
                CHECK((std::abs(jac3[1] - jaccompare3[1]) < 1e-12));
                INFO("calcJacobian");
                CHECK((std::abs(jac3[2] - jaccompare3[2]) < 1e-12));
                INFO("calcJacobian");
                CHECK((std::abs(jac3[3] - jaccompare3[3]) < 1e-12));
                INFO("calcJacobian");
                CHECK((std::abs(jac3[4] - jaccompare3[4]) < 1e-12));
                INFO("calcJacobian");
                CHECK((std::abs(jac3[5] - jaccompare3[5]) < 1e-12));
                INFO("calcJacobian");
                CHECK((std::abs(jac3[6] - jaccompare3[6]) < 1e-12));
                INFO("calcJacobian");
                CHECK((std::abs(jac3[7] - jaccompare3[7]) < 1e-12));
                INFO("calcJacobian");
                CHECK((std::abs(jac3[8] - jaccompare3[8]) < 1e-12));
            }
        }
    }
    logger(INFO, "%", *test);

    pointIndexes.push_back(14);
    pointIndexes.push_back(15);

    delete test;
    test = new Geometry::ElementGeometry(pointIndexes, nodes3D);

    const Geometry::MappingReferenceToPhysicalBase& refMap7 =
        *test->getReferenceToPhysicalMap();
    const Geometry::ReferenceGeometry& refGeo7 = *test->getReferenceGeometry();
    INFO("getReferenceToPhysicalMap");
    CHECK(
        (typeid(Geometry::MappingToPhysHypercubeLinear<3>) == typeid(refMap7)));
    INFO("getReferenceGeometry");
    CHECK((typeid(Geometry::ReferenceCube) == typeid(refGeo7)));
    INFO("getNumberOfNodes");
    CHECK((test->getNumberOfNodes() == 8));

    for (orig3D[0] = -1.51; orig3D[0] < 1.51; orig3D[0] += 0.3) {
        for (orig3D[1] = -1.511; orig3D[1] < 1.51; orig3D[1] += 0.35) {
            for (orig3D[2] = -1.512; orig3D[2] < 1.51; orig3D[2] += 0.3) {
                compare3D =
                    test->getReferenceToPhysicalMap()->transform((orig3D));
                point3D = test->referenceToPhysical((orig3D));
                INFO("referenceToPhysical");
                CHECK((std::abs(point3D[0] - compare3D[0]) < 1e-12));
                INFO("referenceToPhysical");
                CHECK((std::abs(point3D[1] - compare3D[1]) < 1e-12));
                INFO("referenceToPhysical");
                CHECK((std::abs(point3D[2] - compare3D[2]) < 1e-12));

                jaccompare3 =
                    test->getReferenceToPhysicalMap()->calcJacobian((orig3D));
                jac3 = test->calcJacobian((orig3D));
                INFO("calcJacobian");
                CHECK((std::abs(jac3[0] - jaccompare3[0]) < 1e-12));
                INFO("calcJacobian");
                CHECK((std::abs(jac3[1] - jaccompare3[1]) < 1e-12));
                INFO("calcJacobian");
                CHECK((std::abs(jac3[2] - jaccompare3[2]) < 1e-12));
                INFO("calcJacobian");
                CHECK((std::abs(jac3[3] - jaccompare3[3]) < 1e-12));
                INFO("calcJacobian");
                CHECK((std::abs(jac3[4] - jaccompare3[4]) < 1e-12));
                INFO("calcJacobian");
                CHECK((std::abs(jac3[5] - jaccompare3[5]) < 1e-12));
                INFO("calcJacobian");
                CHECK((std::abs(jac3[6] - jaccompare3[6]) < 1e-12));
                INFO("calcJacobian");
                CHECK((std::abs(jac3[7] - jaccompare3[7]) < 1e-12));
                INFO("calcJacobian");
                CHECK((std::abs(jac3[8] - jaccompare3[8]) < 1e-12));
            }
        }
    }
    logger(INFO, "%", *test);

    std::vector<Geometry::PointPhysical<4> > nodes4D;

    Geometry::PointPhysical<4> point4D, compare4D;
    Geometry::PointReference<4> orig4D;
    Geometry::Jacobian<4, 4> jac4, jaccompare4;

    pointIndexes.push_back(16);
    pointIndexes.push_back(17);
    pointIndexes.push_back(18);
    pointIndexes.push_back(19);
    pointIndexes.push_back(20);
    pointIndexes.push_back(21);
    pointIndexes.push_back(22);
    pointIndexes.push_back(23);

    for (double i = 0.; i < 1; i += 0.1) {
        point4D[0] = 1. + i;
        point4D[1] = 2. + i;
        point4D[2] = 3. + i;
        point4D[3] = 1.;
        nodes4D.push_back(point4D);
    }

    point4D[0] = 3.5;
    point4D[1] = 4.6;
    point4D[2] = 5.4;
    point4D[3] = 1.;
    nodes4D.push_back(point4D);
    point4D[0] = 6.7;
    point4D[1] = 2.8;
    point4D[2] = 5.7;
    point4D[3] = 1.;
    nodes4D.push_back(point4D);
    point4D[0] = 1.4;
    point4D[1] = 2.4;
    point4D[2] = 5.4;
    point4D[3] = 2.;
    nodes4D.push_back(point4D);
    point4D[0] = 1.7;
    point4D[1] = 2.7;
    point4D[2] = 5.7;
    point4D[3] = 2.;
    nodes4D.push_back(point4D);
    point4D[0] = 3.5;
    point4D[1] = 4.6;
    point4D[2] = 7.4;
    point4D[3] = 2.;
    nodes4D.push_back(point4D);
    point4D[0] = 6.7;
    point4D[1] = 2.8;
    point4D[2] = 7.7;
    point4D[3] = 2.;
    nodes4D.push_back(point4D);
    point4D[0] = 1.4;
    point4D[1] = 2.4;
    point4D[2] = 5.4;
    point4D[3] = 3.;
    nodes4D.push_back(point4D);
    point4D[0] = 1.7;
    point4D[1] = 2.7;
    point4D[2] = 5.7;
    point4D[3] = 3.;
    nodes4D.push_back(point4D);
    point4D[0] = 3.5;
    point4D[1] = 4.6;
    point4D[2] = 5.4;
    point4D[3] = 3.;
    nodes4D.push_back(point4D);
    point4D[0] = 6.7;
    point4D[1] = 2.8;
    point4D[2] = 5.7;
    point4D[3] = 3.;
    nodes4D.push_back(point4D);
    point4D[0] = 1.4;
    point4D[1] = 2.4;
    point4D[2] = 5.4;
    point4D[3] = 4.;
    nodes4D.push_back(point4D);
    point4D[0] = 1.7;
    point4D[1] = 2.7;
    point4D[2] = 5.7;
    point4D[3] = 4.;
    nodes4D.push_back(point4D);
    point4D[0] = 3.5;
    point4D[1] = 4.6;
    point4D[2] = 7.4;
    point4D[3] = 4.;
    nodes4D.push_back(point4D);
    point4D[0] = 6.7;
    point4D[1] = 2.8;
    point4D[2] = 7.7;
    point4D[3] = 4.;
    nodes4D.push_back(point4D);

    delete test;
    test = new Geometry::ElementGeometry(pointIndexes, nodes4D);

    const Geometry::MappingReferenceToPhysicalBase& refMap8 =
        *test->getReferenceToPhysicalMap();
    const Geometry::ReferenceGeometry& refGeo8 = *test->getReferenceGeometry();
    INFO("getReferenceToPhysicalMap");
    CHECK(
        (typeid(Geometry::MappingToPhysHypercubeLinear<4>) == typeid(refMap8)));
    INFO("getReferenceGeometry");
    CHECK((typeid(Geometry::ReferenceHypercube) == typeid(refGeo8)));
    INFO("getNumberOfNodes");
    CHECK((test->getNumberOfNodes() == 16));

    for (orig4D[0] = -1.5189; orig4D[0] < 1.541; orig4D[0] += 0.4) {
        for (orig4D[1] = -1.5189; orig4D[1] < 1.541; orig4D[1] += 0.4) {
            for (orig4D[2] = -1.5189; orig4D[2] < 1.541; orig4D[2] += 0.5) {
                for (orig4D[3] = -1.5189; orig4D[3] < 1.541; orig4D[3] += 0.5) {
                    compare4D =
                        test->getReferenceToPhysicalMap()->transform((orig4D));
                    point4D = test->referenceToPhysical((orig4D));
                    INFO("referenceToPhysical");
                    CHECK((std::abs(point4D[0] - compare4D[0]) < 1e-12));
                    INFO("referenceToPhysical");
                    CHECK((std::abs(point4D[1] - compare4D[1]) < 1e-12));
                    INFO("referenceToPhysical");
                    CHECK((std::abs(point4D[2] - compare4D[2]) < 1e-12));
                    INFO("referenceToPhysical");
                    CHECK((std::abs(point4D[3] - compare4D[3]) < 1e-12));

                    jaccompare4 =
                        test->getReferenceToPhysicalMap()->calcJacobian(
                            (orig4D));
                    jac4 = test->calcJacobian((orig4D));
                    INFO("calcJacobian");
                    CHECK((std::abs(jac4[0] - jaccompare4[0]) < 1e-12));
                    INFO("calcJacobian");
                    CHECK((std::abs(jac4[1] - jaccompare4[1]) < 1e-12));
                    INFO("calcJacobian");
                    CHECK((std::abs(jac4[2] - jaccompare4[2]) < 1e-12));
                    INFO("calcJacobian");
                    CHECK((std::abs(jac4[3] - jaccompare4[3]) < 1e-12));
                    INFO("calcJacobian");
                    CHECK((std::abs(jac4[4] - jaccompare4[4]) < 1e-12));
                    INFO("calcJacobian");
                    CHECK((std::abs(jac4[5] - jaccompare4[5]) < 1e-12));
                    INFO("calcJacobian");
                    CHECK((std::abs(jac4[6] - jaccompare4[6]) < 1e-12));
                    INFO("calcJacobian");
                    CHECK((std::abs(jac4[7] - jaccompare4[7]) < 1e-12));
                    INFO("calcJacobian");
                    CHECK((std::abs(jac4[8] - jaccompare4[8]) < 1e-12));
                    INFO("calcJacobian");
                    CHECK((std::abs(jac4[9] - jaccompare4[9]) < 1e-12));
                    INFO("calcJacobian");
                    CHECK((std::abs(jac4[10] - jaccompare4[10]) < 1e-12));
                    INFO("calcJacobian");
                    CHECK((std::abs(jac4[11] - jaccompare4[11]) < 1e-12));
                    INFO("calcJacobian");
                    CHECK((std::abs(jac4[12] - jaccompare4[12]) < 1e-12));
                    INFO("calcJacobian");
                    CHECK((std::abs(jac4[13] - jaccompare4[13]) < 1e-12));
                    INFO("calcJacobian");
                    CHECK((std::abs(jac4[14] - jaccompare4[14]) < 1e-12));
                    INFO("calcJacobian");
                    CHECK((std::abs(jac4[15] - jaccompare4[15]) < 1e-12));
                }
            }
        }
    }
    logger(INFO, "%", *test);
    delete test;
}
