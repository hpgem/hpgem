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
#include "VTKLagrangeTriangle.h"

#include <valarray>

namespace hpgem {
namespace Output {

VTKLagrangeTriangle::VTKLagrangeTriangle(std::size_t order)
    : points_((order + 1) * (order + 2) / 2) {
    std::vector<std::valarray<std::size_t>> bary =
        computeBaryIntegerPoints(order);
    double h = 1.0 / order;
    for (std::size_t i = 0; i < bary.size(); ++i) {
        // Conversion from barycentric to actual coordinates.
        points_[i][0] = h * bary[i][1];
        points_[i][1] = h * bary[i][2];
    }
}

std::vector<std::valarray<std::size_t>>
    VTKLagrangeTriangle::computeBaryIntegerPoints(std::size_t order) {

    // Integer variation on barycenteric coordinates. When multiplying by h =
    // 1/order these are actual barycentric coordinates. While the actual
    // coordiantes of the points don't matter, they will be mapped such that
    // [order, 0, 0] = (0,0), [0, order, 0] = (1,0) and [0, 0, order] = (0,1).
    using IVec = std::valarray<std::size_t>;

    std::size_t numPoints = (order + 1) * (order + 2) / 2;
    std::vector<IVec> points(numPoints);

    std::size_t index = 0;  // Index in the output
    for (std::size_t layer = 0; 3 * layer <= order; ++layer) {
        // Split the points on this layer into offset + localPoint
        // Offset is a vector that ensures that we are away from each corner
        // point by 'layer'-points. The localPoint is a localPoint on the layer,
        // which is ensured by having at least one zero coordinate.
        IVec offset = IVec({layer, layer, layer});
        // sum(offset[i]) = 3*layer, so ew need sum(localPoint[i]) = rOrder
        std::size_t rOrder = order - 3 * layer;

        // Either 1 or 3 corner points
        points[index++] = offset + IVec({rOrder, 0, 0});
        if (rOrder == 0) {
            break;
        }
        points[index++] = offset + IVec({0, rOrder, 0});
        points[index++] = offset + IVec({0, 0, rOrder});

        // For rOrder >= 2 there are rOrder - 1 intermediate points on each edge
        // These are enumerated in counter clockwise order
        for (std::size_t i = 0; i < rOrder - 1; ++i) {
            std::size_t i1 = i + 1;
            std::size_t i2 = rOrder - i - 1;
            points[index + i + 0 * (rOrder - 1)] = offset + IVec({i2, i1, 0});
            points[index + i + 1 * (rOrder - 1)] = offset + IVec({0, i2, i1});
            points[index + i + 2 * (rOrder - 1)] = offset + IVec({i1, 0, i2});
        }
        index += 3 * (rOrder - 1);
    }
    // Safety check
    logger.assert_always(
        index == numPoints,
        "Incorrect number of points for the VTKLangrangeTriangle of order %",
        order);
    return points;
}

}  // namespace Output
}  // namespace hpgem
