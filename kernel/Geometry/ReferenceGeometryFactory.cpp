/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2021, University of Twente
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
#include "ReferenceGeometryFactory.h"

#include "ReferenceCube.h"
#include "ReferenceHypercube.h"
#include "ReferenceCurvilinearLine.h"
#include "ReferenceCurvilinearTetrahedron.h"
#include "ReferenceCurvilinearTriangle.h"
#include "ReferenceLine.h"
#include "ReferencePoint.h"
#include "ReferencePyramid.h"
#include "ReferenceSquare.h"
#include "ReferenceTetrahedron.h"
#include "ReferenceTriangle.h"
#include "ReferenceTriangularPrism.h"

namespace hpgem {
namespace Geometry {

ReferenceGeometry& ReferenceGeometryFactory::getGeometry(
    std::size_t dimension, std::size_t numberOfPoints) {
    // Split out decision logic
    switch (dimension) {
        case 0:
            return getGeometry0(numberOfPoints);

        case 1:
            return getGeometry1(numberOfPoints);
        case 2:
            return getGeometry2(numberOfPoints);
        case 3:
            return getGeometry3(numberOfPoints);
        case 4:
            return getGeometry4(numberOfPoints);
        default:
            logger.fail("No shapes defined for dimension % (% points)",
                        dimension, numberOfPoints);
    }
}

inline ReferenceGeometry& ReferenceGeometryFactory::getGeometry0(
    std::size_t numberOfPoints) {
    if (numberOfPoints == 1) {
        return ReferencePoint::Instance();
    }
    logger.fail("A zero-th dimensional shape must have 1 point");
}

ReferenceGeometry& ReferenceGeometryFactory::getGeometry1(
    std::size_t numberOfPoints) {
    // Dimension 1 is completely decidable
    if (numberOfPoints == 2) {
        return ReferenceLine::Instance();
    } else {
        logger.assert_always(numberOfPoints > 2, "Too few points for a Line");
        return ReferenceCurvilinearLine::getReferenceLagrangeLine(
            numberOfPoints - 1);
    }
}

ReferenceGeometry& ReferenceGeometryFactory::getGeometry2(
    std::size_t numberOfPoints) {
    switch (numberOfPoints) {
        case 3:
            return ReferenceTriangle::Instance();
        case 4:
            return ReferenceSquare::Instance();
        default:
            // Handled below
            break;
    }
    ReferenceGeometry*& geometry = cached2DGeometries_[numberOfPoints];
    if (geometry == nullptr) {
        int order =
            ReferenceCurvilinearTriangle::getOrderFromPoints(numberOfPoints);
        logger.assert_always(order > 0, "No 2D shape with % points",
                             numberOfPoints);
        geometry =
            &ReferenceCurvilinearTriangle::getReferenceLagrangeTriangle(order);
    }
    return *geometry;
}
ReferenceGeometry& ReferenceGeometryFactory::getGeometry3(
    std::size_t numberOfPoints) {

    switch (numberOfPoints) {
        case 4:
            return ReferenceTetrahedron::Instance();
        case 5:
            return ReferencePyramid::Instance();
        case 6:
            return ReferenceTriangularPrism::Instance();
        case 8:
            return ReferenceCube::Instance();
        default:
            // Handled below
            break;
    }
    ReferenceGeometry*& geometry = cached3DGeometries_[numberOfPoints];
    if (geometry == nullptr) {
        int order =
            ReferenceCurvilinearTetrahedron::getOrderFromPoints(numberOfPoints);
        logger.assert_always(order > 0, "No 3D shape with % points",
                             numberOfPoints);
        geometry = &ReferenceCurvilinearTetrahedron::
                       getReferenceCurvilinearTetrahedron(order);
    }
    return *geometry;
}
ReferenceGeometry& ReferenceGeometryFactory::getGeometry4(
    std::size_t numberOfPoints) {
    if (numberOfPoints == 16) {
        return ReferenceHypercube::Instance();
    } else {
        logger.fail("No 4D shape with % points known", numberOfPoints);
    }
}

}  // namespace Geometry
}  // namespace hpgem
