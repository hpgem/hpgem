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
#include "Geometry/Mappings/MappingToRefPointToLine.h"
#include "Logger.h"

#include "Geometry/ReferenceLine.h"
#include "Geometry/ReferencePoint.h"
#include "Geometry/PointReference.h"
#include "Geometry/Jacobian.h"
#include "LinearAlgebra/MiddleSizeVector.h"
#include <cmath>
using namespace hpgem;
int main() {

    Geometry::PointReference<0> refPoint;
    Geometry::PointReference<1> point, compare;

    Geometry::ReferenceLine& eGeom = Geometry::ReferenceLine::Instance();
    Geometry::ReferencePoint& fGeom = Geometry::ReferencePoint::Instance();

    Geometry::Jacobian<0, 1> jac;

    std::vector<std::size_t> nodesAfterTransformation(1);

    const Geometry::MappingReferenceToReference<1>* test =
        &Geometry::MappingToRefPointToLine0::Instance();
    nodesAfterTransformation[0] = 0;

    point = test->transform((refPoint));
    logger.assert_always(
        (fGeom.isInternalPoint((refPoint)) == eGeom.isInternalPoint((point))),
        "transform");

    jac = test->calcJacobian((refPoint));

    for (std::size_t i = 0; i < fGeom.getNumberOfNodes(); ++i) {
        refPoint = fGeom.getReferenceNodeCoordinate(i);
        compare = eGeom.getReferenceNodeCoordinate(nodesAfterTransformation[i]);
        point = test->transform((refPoint));
        logger.assert_always((std::abs(point[0] - compare[0]) < 1e-12),
                             "transform");
    }

    logger.assert_always((test->getTargetDimension() == 1),
                         "getTargetDimension");

    test = &Geometry::MappingToRefPointToLine1::Instance();
    nodesAfterTransformation[0] = 1;

    point = test->transform((refPoint));
    logger.assert_always(
        (fGeom.isInternalPoint((refPoint)) == eGeom.isInternalPoint((point))),
        "transform");

    jac = test->calcJacobian((refPoint));

    for (std::size_t i = 0; i < fGeom.getNumberOfNodes(); ++i) {
        refPoint = fGeom.getReferenceNodeCoordinate(i);
        compare = eGeom.getReferenceNodeCoordinate(nodesAfterTransformation[i]);
        point = test->transform((refPoint));
    }

    logger.assert_always((test->getTargetDimension() == 1),
                         "getTargetDimension");

    return 0;
}
