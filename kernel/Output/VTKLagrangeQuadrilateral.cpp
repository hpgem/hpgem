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
#include "VTKLagrangeQuadrilateral.h"

#include "Geometry/PointReference.h"

namespace hpgem {

namespace Output {

VTKLagrangeQuadrilateral::VTKLagrangeQuadrilateral(std::size_t order)
    : points_((order + 1) * (order + 1)) {

    double h = 2.0 / order;

    // Corner points
    points_[0] = Geometry::PointReference<2>({-1.0, -1.0});
    points_[1] = Geometry::PointReference<2>({+1.0, -1.0});
    points_[2] = Geometry::PointReference<2>({+1.0, +1.0});
    points_[3] = Geometry::PointReference<2>({-1.0, +1.0});
    // Four sides
    // Each side has order + 1 points, of which the outer two are the corner
    // points
    for (std::size_t i = 0; i < order - 1; ++i) {
        points_[4 + 0 * order + i] =
            Geometry::PointReference<2>({-1.0 + (i + 1) * h, -1.0});
        points_[3 + 1 * order + i] =
            Geometry::PointReference<2>({+1.0, -1.0 + (i + 1) * h});
        points_[2 + 2 * order + i] =
            Geometry::PointReference<2>({-1.0 + (i + 1) * h, +1.0});
        points_[1 + 3 * order + i] =
            Geometry::PointReference<2>({-1.0, -1.0 + (i + 1) * h});
    }
    // Indices 0 - 4*p-1 are now filled with the boundary points, the rest is
    // for the (p-1)^2 inner points
    std::size_t index = 4 * order;
    for (std::size_t i = 0; i < order - 1; ++i) {
        double y = -1.0 + (i + 1) * h;
        for (std::size_t j = 0; j < order - 1; ++j) {
            double x = -1.0 + (j + 1) * h;
            points_[index++] = Geometry::PointReference<2>({x, y});
        }
    }
    logger.assert_debug(index == points_.size(),
                        "Incorrect number of points on VTK quadrilateral");
}

}  // namespace Output

}  // namespace hpgem
