/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2022, University of Twente
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

#include "../../catch.hpp"
#include "MappingTestHelpers.h"

#include <Geometry/Mappings/MappingToPhysTetrahedronQuadratic.h>
#include <LinearAlgebra/SmallMatrix.h>
#include <Geometry/ReferenceCurvilinearTetrahedron.h>

using namespace hpgem;
using namespace hpgem::LinearAlgebra;
using namespace hpgem::Geometry;

using PRef = PointReference<3>;
using PPhys = PointPhysical<3>;

TEST_CASE("Affine mapped coordinates", "[MappingToPhysTetrahedronQuadratic]") {
    // x' = 0.1 x + 0.6 y - 0.3 z + 1.5
    // y' = 0.9 x - 0.1 y + 0.7 z + 0.3
    // z' = -0.3 x + 0  y + 0.5 z - 3.5
    SmallMatrix<3, 3> transform = {
        {0.1, 0.6, -0.3}, {0.9, -0.1, 0.7}, {-0.3, 0.0, 0.5}};
    SmallVector<3> offset({1.5, 0.3, -3.5});

    const ReferenceGeometry& geom =
        ReferenceCurvilinearTetrahedron::getReferenceCurvilinearTetrahedron(2);

    std::vector<std::size_t> nodeIndices(geom.getNumberOfNodes());
    std::vector<PPhys> pphyss(geom.getNumberOfNodes());

    for (std::size_t i = 0; i < geom.getNumberOfNodes(); ++i) {
        // Map the reference points
        const PRef& pref = geom.getReferenceNodeCoordinate(i);
        SmallVector<3> coords = pref.getCoordinates();
        pphyss[i] = transform * coords + offset;
        nodeIndices[i] = i;
    }

    PhysicalGeometry<3> physGeom(nodeIndices, pphyss, &geom);
    MappingToPhysTetrahedronQuadratic mapping(&physGeom);

    INFO("Check reference point mapping");
    for (std::size_t i = 0; i < geom.getNumberOfNodes(); ++i) {
        auto pphys = mapping.transform(geom.getReferenceNodeCoordinate(i));
        REQUIRE((pphys - pphyss[i]).l2NormSquared() < 1e-16);

        auto pref = mapping.inverseTransform(pphys);
        REQUIRE((pref - geom.getReferenceNodeCoordinate(i)).l2NormSquared() <
                1e-16);
    }
    INFO("Check Jacobian")
    // Some varying points inside
    testJacobian({0.1, 0.1, 0.1}, mapping);
    testJacobian({0.1, 0.1, 0.8}, mapping);
    testJacobian({0.33333, 0.33333, 0.33333}, mapping);
}

PPhys transformNonAffine(const PRef& p) {
    // Arbitrarily chosen quadratic transformation.
    // Values chosen so that the transformation is predominantly linear.
    // x = 2y + 0.1 x^2 - 0.3 xz + 0.1
    // y = 2z - 0.1 y^2 - 0.2 yz - 1.5
    // z = 2x + 0.2 y^2 - 0.1 xy + 3
    double x = p[0];
    double y = p[1];
    double z = p[2];
    PPhys result{2 * y + 0.1 * x * x - 0.3 * x * z + 0.1,
                 2 * z - 0.1 * y * y - 0.2 * y * z - 1.5,
                 2 * x + 0.2 * y * y - 0.1 * x * y + 3.0};
    return result;
}

TEST_CASE("Non affine mapping", "MappingToPhysTetrahedronQuadratic") {
    // Use a non affine transformation of physical space from the function
    // transformNonAffine

    const ReferenceGeometry& geom =
        ReferenceCurvilinearTetrahedron::getReferenceCurvilinearTetrahedron(2);

    std::vector<std::size_t> nodeIndices(geom.getNumberOfNodes());
    std::vector<PPhys> pphyss(geom.getNumberOfNodes());

    for (std::size_t i = 0; i < geom.getNumberOfNodes(); ++i) {
        const PRef& pref = geom.getReferenceNodeCoordinate(i);
        pphyss[i] = transformNonAffine(pref);
        nodeIndices[i] = i;
    }

    PhysicalGeometry<3> physGeom(nodeIndices, pphyss, &geom);
    MappingToPhysTetrahedronQuadratic mapping(&physGeom);

    INFO("Check reference point mapping")
    for (std::size_t i = 0; i < geom.getNumberOfNodes(); ++i) {
        auto pphys = mapping.transform(geom.getReferenceNodeCoordinate(i));
        REQUIRE((pphyss[i] - pphys).l2NormSquared() < 1e-16);

        auto pref = mapping.inverseTransform(pphyss[i]);
        REQUIRE((pref - geom.getReferenceNodeCoordinate(i)).l2NormSquared() <
                1e-16);
    }

    INFO("Check single point mapping")
    for (const PRef& pref : {PRef{0.1, 0.5, 0.1}, PRef{0.7, 0.1, 0.1},
                             PRef{0.33333, 0.3333, 0.333}}) {
        auto mapped = mapping.transform(pref);
        auto directMapped = transformNonAffine(pref);
        REQUIRE((mapped - directMapped).l2NormSquared() < 1e-16);
        testJacobian(pref, mapping);
    }
}
