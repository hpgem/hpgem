/*
This file forms part of hpGEM. This package has been developed over a number of
years by various people at the University of Twente and a full list of
contributors can be found at http://hpgem.org/about-the-code/team

This code is distributed using BSD 3-Clause License. A copy of which can found
below.


Copyright (c) 2014, University of Twente
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "MappingToPhysTriangularPrism.h"

#include "Geometry/PointPhysical.h"
#include "Geometry/PointReference.h"
#include "Geometry/PhysicalGeometry.h"
#include "Geometry/Jacobian.h"

#include <limits>

namespace hpgem {

namespace Geometry {
MappingToPhysTriangularPrism::MappingToPhysTriangularPrism(
    const PhysicalGeometry<3>* const physicalGeometry)
    : MappingReferenceToPhysical<3>(physicalGeometry) {
    logger.assert_debug(physicalGeometry != nullptr,
                        "Invalid physical geometry passed");
    reinit();
}

PointPhysical<3> MappingToPhysTriangularPrism::transform(
    const PointReference<3>& pR) const {
#ifdef SAVECOEFFS
    return a0 + a1 * xi[0] + a2 * xi[1] + a3 * xi[2] + a4 * xi[0] * xi[2] +
           a5 * xi[1] * xi[2];
#else
    PointPhysical<3> pP;
    const double t1 = pR[0] * pR[2];
    const double t2 = pR[1] * pR[2];
    LinearAlgebra::SmallVector<6> f2;

    f2[0] = 0.5 * (1.0 - pR[0] - pR[1] - pR[2] + t1 + t2);
    f2[1] = 0.5 * (pR[0] - t1);
    f2[2] = 0.5 * (pR[1] - t2);
    f2[3] = 0.5 * (-t2 + 1.0 - pR[0] - t1 - pR[1] + pR[2]);
    f2[4] = 0.5 * (t1 + pR[0]);
    f2[5] = 0.5 * (t2 + pR[1]);

    PointPhysical<3> p;

    for (std::size_t i = 0; i < 6; ++i) {
        p = geometry_->getLocalNodeCoordinates(i);
        pP += f2[i] * p;
    }
    return pP;
#endif
}

PointReference<3> MappingToPhysTriangularPrism::inverseTransform(
    const PointPhysical<3>& pointPhysical) const {
    Geometry::PointReference<3> result;
    Geometry::PointPhysical<3> comparison = transform(result);
    LinearAlgebra::SmallVector<3> correction;
    double error = (pointPhysical - comparison).l2NormSquared();
    std::size_t loop_count{0};
    while (error > 1e-28 && loop_count++ < 100) {
        correction = (pointPhysical - comparison).getCoordinates();
        calcJacobian(result).solve(correction);
        result = PointReference<3>(result + correction);
        comparison = transform(result);
        error = (pointPhysical - comparison).l2NormSquared();
    }
    if (loop_count == 100) {
        return {std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN()};
    }
    return result;
}

Jacobian<3, 3> MappingToPhysTriangularPrism::calcJacobian(
    const PointReference<3>& pR) const {
    Jacobian<3, 3> jacobian;
#ifdef SAVECOEFFS
    Geometry::PointPhysical<3> d_dxi0(a1 + xi[2] * a4);
    Geometry::PointPhysical<3> d_dxi1(a2 + xi[2] * a5);
    Geometry::PointPhysical<3> d_dxi2(a3 + xi[0] * a4 + xi[1] * a5);

    for (int i = 0; i < 3; ++i) {
        jacobian(i, 0) = d_dxi0[i];
        jacobian(i, 1) = d_dxi1[i];
        jacobian(i, 2) = d_dxi2[i];
    }
#else
    LinearAlgebra::SmallVector<6> df_dxi0, df_dxi1, df_dxi2;

    df_dxi0[0] = +0.5 * (-1. + pR[2]);
    df_dxi0[1] = +0.5 * (+1. - pR[2]);
    df_dxi0[2] = +0.0;
    df_dxi0[3] = +0.5 * (-1. - pR[2]);
    df_dxi0[4] = +0.5 * (+1. + pR[2]);
    df_dxi0[5] = +0.0;

    df_dxi1[0] = +0.5 * (-1. + pR[2]);
    df_dxi1[1] = +0.0;
    df_dxi1[2] = +0.5 * (+1. - pR[2]);
    df_dxi1[3] = +0.5 * (-1. - pR[2]);
    df_dxi1[4] = +0.0;
    df_dxi1[5] = +0.5 * (+1. + pR[2]);

    df_dxi2[0] = +0.5 * (-1. + pR[0] + pR[1]);
    df_dxi2[1] = -0.5 * pR[0];
    df_dxi2[2] = -0.5 * pR[1];
    df_dxi2[3] = +0.5 * (+1. - pR[1] - pR[0]);
    df_dxi2[4] = +0.5 * pR[0];
    df_dxi2[5] = +0.5 * pR[1];

    Geometry::PointPhysical<3> d_dxi0;
    Geometry::PointPhysical<3> d_dxi1;
    Geometry::PointPhysical<3> d_dxi2;

    for (std::size_t i = 0; i < 3; ++i) {
        d_dxi0[i] = 0.;
        d_dxi1[i] = 0.;
        d_dxi2[i] = 0.;
    }

    Geometry::PointPhysical<3> p;

    for (std::size_t i = 0; i < 6; ++i) {
        p = geometry_->getLocalNodeCoordinates(i);

        d_dxi0 += df_dxi0[i] * p;
        d_dxi1 += df_dxi1[i] * p;
        d_dxi2 += df_dxi2[i] * p;
    }

    for (std::size_t i = 0; i < 3; ++i) {
        jacobian(i, 0) = d_dxi0[i];
        jacobian(i, 1) = d_dxi1[i];
        jacobian(i, 2) = d_dxi2[i];
    }
#endif
    return jacobian;
}

void MappingToPhysTriangularPrism::reinit() {
#ifdef SAVECOEFFS
    FixedVector<Geometry::PointPhysical<3>, 6> p;

    for (std::size_t i = 0; i < 6; ++i) {
        physicalGeometry->getVertexPoint(i, p[i]);
    }

    // see Maple/triangularprismMapping1.mws
    a0 = 0.5 * (p[0] + p[3]);
    a1 = 0.5 * (-p[0] + p[4] - p[3] + p[1]);
    a2 = 0.5 * (-p[0] + p[5] - p[3] + p[2]);
    a3 = 0.5 * (-p[0] + p[3]);
    a4 = 0.5 * (p[0] + p[4] - p[3] - p[1]);
    a5 = 0.5 * (p[0] + p[5] - p[3] - p[2]);
#endif
}

bool MappingToPhysTriangularPrism::isValidPoint(
    const PointReference<3>& pointReference) const {
    return (0. <= pointReference[0]) && (pointReference[0] <= 1.) &&
           (0. <= pointReference[1]) && (pointReference[1] <= 1.) &&
           (pointReference[0] + pointReference[1] - 1. <= 1.e-16) &&
           (pointReference[2] >= -1.) && (pointReference[2] <= 1.);
}
}  // namespace Geometry

}  // namespace hpgem
