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
#include "VTKLagrangeHexahedron.h"

namespace hpgem {
namespace Output {

VTKLagrangeHexahedron::VTKLagrangeHexahedron(std::size_t order)
    : points_((order + 1) * (order + 1) * (order + 1)) {
    // Corner points
    points_[0] = Geometry::PointReference<3>({-1.0, -1.0, -1.0});
    points_[1] = Geometry::PointReference<3>({+1.0, -1.0, -1.0});
    points_[2] = Geometry::PointReference<3>({+1.0, +1.0, -1.0});
    points_[3] = Geometry::PointReference<3>({-1.0, +1.0, -1.0});
    points_[4] = Geometry::PointReference<3>({-1.0, -1.0, +1.0});
    points_[5] = Geometry::PointReference<3>({+1.0, -1.0, +1.0});
    points_[6] = Geometry::PointReference<3>({+1.0, +1.0, +1.0});
    points_[7] = Geometry::PointReference<3>({-1.0, +1.0, +1.0});

    double h = 2.0 / order;

    // On each edge we now get (order - 1) edge points
    std::size_t offset = 8;
    for (std::size_t i = 0; i < order - 1; ++i) {
        // z=-1.0 plane
        double dx = (i + 1) * h;
        points_[offset + 0 * (order - 1) + i] =
            Geometry::PointReference<3>({-1.0 + dx, -1.0, -1.0});
        points_[offset + 1 * (order - 1) + i] =
            Geometry::PointReference<3>({+1.0, -1.0 + dx, -1.0});
        points_[offset + 2 * (order - 1) + i] =
            Geometry::PointReference<3>({-1.0 + dx, +1.0, -1.0});
        points_[offset + 3 * (order - 1) + i] =
            Geometry::PointReference<3>({-1.0, -1.0 + dx, -1.0});
        // z=+1.0 plane
        points_[offset + 4 * (order - 1) + i] =
            Geometry::PointReference<3>({-1.0 + dx, -1.0, +1.0});
        points_[offset + 5 * (order - 1) + i] =
            Geometry::PointReference<3>({+1.0, -1.0 + dx, +1.0});
        points_[offset + 6 * (order - 1) + i] =
            Geometry::PointReference<3>({-1.0 + dx, +1.0, +1.0});
        points_[offset + 7 * (order - 1) + i] =
            Geometry::PointReference<3>({-1.0, -1.0 + dx, +1.0});
        // Side edges
        points_[offset + 8 * (order - 1) + i] =
            Geometry::PointReference<3>({-1.0, -1.0, -1.0 + dx});
        points_[offset + 9 * (order - 1) + i] =
            Geometry::PointReference<3>({+1.0, -1.0, -1.0 + dx});
        // Note: this is different from the usual order
        points_[offset + 10 * (order - 1) + i] =
            Geometry::PointReference<3>({-1.0, +1.0, -1.0 + dx});
        points_[offset + 11 * (order - 1) + i] =
            Geometry::PointReference<3>({+1.0, +1.0, -1.0 + dx});
    }
    // The 12*(order-1) points have been used for the edges
    offset += 12 * (order - 1);
    // The first 8 + 12*(order-1) points have now been filled
    // The next entries will be the 6 faces, each with (order-1)^2 points on
    // them.
    // They ordering of the faces is x=-1.0, x=1.0, y=-1.0, y=1.0, z=-1.0 and
    // lastly z=1.0. Within the face the ordering of the points is so that the
    // x-coordinate changes fastest, then the y and slowest is z (exclude the
    // one which is kept fixed).
    std::size_t facePoints = (order - 1) * (order - 1);
    for (std::size_t i = 0; i < order - 1; ++i) {
        double slow = -1.0 + (i + 1) * h;
        for (std::size_t j = 0; j < order - 1; ++j) {
            double fast = -1.0 + (j + 1) * h;
            points_[offset + 0 * facePoints + i * (order - 1) + j] =
                Geometry::PointReference<3>({-1.0, fast, slow});
            points_[offset + 1 * facePoints + i * (order - 1) + j] =
                Geometry::PointReference<3>({+1.0, fast, slow});
            points_[offset + 2 * facePoints + i * (order - 1) + j] =
                Geometry::PointReference<3>({fast, -1.0, slow});
            points_[offset + 3 * facePoints + i * (order - 1) + j] =
                Geometry::PointReference<3>({fast, +1.0, slow});
            points_[offset + 4 * facePoints + i * (order - 1) + j] =
                Geometry::PointReference<3>({fast, slow, -1.0});
            points_[offset + 5 * facePoints + i * (order - 1) + j] =
                Geometry::PointReference<3>({fast, slow, +1.0});
        }
    }
    // 6 faces have been added
    offset += 6 * facePoints;
    // Inside points, again so that the x coordinate is fastest and z is slowest
    for (std::size_t zi = 0; zi < order - 1; ++zi) {
        double z = -1.0 + (zi + 1) * h;
        for (std::size_t yi = 0; yi < order - 1; ++yi) {
            double y = -1.0 + (yi + 1) * h;
            for (std::size_t xi = 0; xi < order - 1; ++xi) {
                double x = -1.0 + (xi + 1) * h;
                points_[offset++] = Geometry::PointReference<3>({x, y, z});
            }
        }
    }
    logger.assert_debug(offset == points_.size(),
                        "Incorrect number of VTKLangrangeHexehedron points");
}

}  // namespace Output
}  // namespace hpgem
