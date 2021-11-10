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
#include "MappingToPhysTriangleQuadratic.h"

#include "Geometry/PhysicalGeometry.h"

namespace hpgem {
namespace Geometry {

/*
 * Implementation notes:
 * Ordering of the points in the LagrangeTriangle: (multiplied by 2)
 * (0,0) -> (1-x-y)(1-2x-2y)
 * (1,0) -> 4x(1-x-y)
 * (2,0) -> x(2x-1)
 * (0,1) -> 4y(1-x-y)
 * (1,1) -> 4xy
 * (0,2) -> y(2y-1)
 *
 *
 */

MappingToPhysTriangleQuadratic::MappingToPhysTriangleQuadratic(
    const PhysicalGeometry<2> *const physicalGeometry)
    : MappingReferenceToPhysical(physicalGeometry) {
    logger.assert_always(physicalGeometry != nullptr,
                         "Invalid physical geometry passed");
    reinit();
}

PointPhysical<2> MappingToPhysTriangleQuadratic::transform(
    const PointReference<2> &pr) const {
    PointPhysical<2> result;
    const double x = pr[0];
    const double y = pr[1];
    const double r = 1 - x - y;

    std::vector<double> weights(6);
    weights[0] = r * (2 * r - 1);
    weights[1] = 4 * x * r;
    weights[2] = x * (2 * x - 1);
    weights[3] = 4 * y * r;
    weights[4] = 4 * x * y;
    weights[5] = y * (2 * y - 1);

    PointPhysical<2> p;
    for (std::size_t i = 0; i < 6; ++i) {
        p = geometry_->getLocalNodeCoordinates(i);
        result += weights[i] * p;
    }
    return result;
}

PointReference<2> MappingToPhysTriangleQuadratic::inverseTransform(
    const PointPhysical<2> &p) const {
    PointReference<2> result = {0.3333, 0.333};  // Centroid
    // Dirty way: Newton iteration
    PointReference<2> correction;
    PointPhysical<2> physError;
    do {
        physError = transform(result);
        physError -= p;
        LinearAlgebra::SmallVector<2> temp = physError.getCoordinates();
        Jacobian<2, 2> jacobian = calcJacobian(result);
        jacobian.solve(temp);
        correction = Geometry::PointReference<2>(temp);
        result -= correction;
    } while (correction.getCoordinates().l2NormSquared() > 1e-10);

    return result;
}
Jacobian<2, 2> MappingToPhysTriangleQuadratic::calcJacobian(
    const PointReference<2> &p) const {
    const double x = p[0];
    const double y = p[1];
    const double r = 1 - x - y;

    std::vector<std::vector<double>> derivs(2);
    derivs[0].resize(6);
    derivs[1].resize(6);

    derivs[0][0] = -4 * r + 1;
    derivs[1][0] = derivs[0][0];

    derivs[0][1] = 4 * (r - x);
    derivs[1][1] = -4 * x;

    derivs[0][2] = 4 * x - 1;
    derivs[1][2] = 0;

    derivs[0][3] = -4 * y;
    derivs[1][3] = 4 * (r - y);

    derivs[0][4] = 4 * y;
    derivs[1][4] = 4 * x;

    derivs[0][5] = 0;
    derivs[1][5] = 4 * y - 1;

    Jacobian<2, 2> result;
    PointPhysical<2> pphys;
    for (std::size_t i = 0; i < 6; ++i) {
        pphys = geometry_->getLocalNodeCoordinates(i);
        for (std::size_t j = 0; j < 2; ++j) {
            for (std::size_t k = 0; k < 2; ++k) {
                result(j, k) += pphys[j] * derivs[k][i];
            }
        }
    }
    return result;
}
void MappingToPhysTriangleQuadratic::reinit() {}
MappingReferenceToPhysicalBase *MappingToPhysTriangleQuadratic::copy() const {
    return new MappingToPhysTriangleQuadratic(*this);
}

}  // namespace Geometry
}  // namespace hpgem