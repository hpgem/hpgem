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
#ifndef HPGEM_REFERENCECURVILINEARELEMENTCHECKS_H
#define HPGEM_REFERENCECURVILINEARELEMENTCHECKS_H

#include "Geometry/ReferenceCurvilinearElement.h"
#include "Geometry/Mappings/MappingReferenceToReference.h"
#include "../catch.hpp"
#include "typeinfo"

namespace hpgem {

/// Check
template <std::size_t d>
void testReferenceGeometryOfCodimReferencePoints(
    const Geometry::ReferenceCurvilinearElement<d>& geom) {
    static_assert(d > 0 && d <= 2, "Not implemented for this dimension");
    std::size_t nent;
    if (d == 1) {
        nent = geom.getNumberOfCodim1Entities();
    } else if (d == 2) {
        nent = geom.getNumberOfCodim2Entities();
    }
    for (std::size_t i = 0; i < nent; ++i) {
        const Geometry::ReferenceGeometry* pointlike;
        if (d == 1) {
            pointlike = geom.getCodim1ReferenceGeometry(i);
        } else if (d == 2) {
            pointlike = geom.getCodim2ReferenceGeometry(i);
        }
        INFO("Check codim " << d << " entity " << i << "to be a point.");
        // Stringent test, it should not even be a subclass.
        CHECK(typeid(Geometry::ReferencePoint) == typeid(*pointlike));
    }
}

template <std::size_t d>
void testAllReferencePointsAreInside(
    const Geometry::ReferenceCurvilinearElement<d>& geom) {
    std::size_t numPoints = geom.getNumberOfNodes();

    auto* baseGeom = geom.getBaseGeometry();
    REQUIRE(baseGeom != nullptr);

    INFO("All reference points should be internal");
    for (std::size_t i = 0; i < numPoints; ++i) {
        const Geometry::PointReference<d> p =
            geom.getReferenceNodeCoordinate(i);
        // Internal by both base and curvilinear geometry
        CHECK(baseGeom->isInternalPoint(p));
        CHECK(geom.isInternalPoint(p));
    }
}

template <std::size_t d>
void testCodim1CorrespondenceWithBaseGeometry(
    const Geometry::ReferenceCurvilinearElement<d>& geom) {
    static_assert(d >= 1, "Codim 1 entities only relevant for dim >= 1");
    std::size_t numCodim1 = geom.getNumberOfCodim1Entities();
    Geometry::ReferenceGeometry* baseGeom = geom.getBaseGeometry();

    REQUIRE(numCodim1 == baseGeom->getNumberOfCodim1Entities());

    for (std::size_t c1 = 0; c1 < numCodim1; ++c1) {
        const auto* mapping = baseGeom->getCodim1MappingPtr(c1);
        const auto* c1Geom = geom.getCodim1ReferenceGeometry(c1);
        INFO("Codim-1 entity: " + std::to_string(c1));
        REQUIRE(c1Geom->getGeometryType() ==
                baseGeom->getCodim1ReferenceGeometry(c1)->getGeometryType());
        auto pointIndices = geom.getCodim1EntityLocalIndices(c1);
        for (std::size_t p = 0; p < c1Geom->getNumberOfNodes(); ++p) {
            Geometry::PointReference<d - 1> pface =
                c1Geom->getReferenceNodeCoordinate(p);
            Geometry::PointReference<d> mappedpface = mapping->transform(pface);
            Geometry::PointReference<d> localp =
                geom.getReferenceNodeCoordinate(pointIndices[p]);
            CHECK(localp == mappedpface);
        }
    }
}

template <std::size_t d>
void testCodim2CorrespondenceWithBaseGeometry(
    const Geometry::ReferenceCurvilinearElement<d>& geom) {
    static_assert(d > 2, "Codim 2 entities only relevant for d > 2");
    std::size_t numCodim2 = geom.getNumberOfCodim2Entities();
    Geometry::ReferenceGeometry* baseGeom = geom.getBaseGeometry();

    using Geometry::PointReference;

    REQUIRE(numCodim2 == baseGeom->getNumberOfCodim2Entities());
    for (std::size_t c2 = 0; c2 < numCodim2; ++c2) {
        const auto* mapping = baseGeom->getCodim2MappingPtr(c2);
        const auto* c2Geom = geom.getCodim2ReferenceGeometry(c2);
        INFO("Codim-2 entity: " + std::to_string(c2));
        REQUIRE(c2Geom->getGeometryType() ==
                baseGeom->getCodim2ReferenceGeometry(c2)->getGeometryType());
        auto pointIndices = baseGeom->getCodim2EntityLocalIndices(c2);
        for (std::size_t p = 0; p < c2Geom->getNumberOfNodes(); ++p) {
            PointReference<d - 2> pridge =
                c2Geom->getReferenceNodeCoordinate(p);
            PointReference<d> mappedPoint = mapping->transform(pridge);
            PointReference<d> localPoint =
                geom.getReferenceNodeCoordinate(pointIndices[p]);
            CHECK(localPoint == mappedPoint);
        }
    }
}

}  // namespace hpgem

#endif  // HPGEM_REFERENCECURVILINEARELEMENTCHECKS_H
