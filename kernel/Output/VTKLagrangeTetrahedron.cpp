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
#include "VTKLagrangeTetrahedron.h"

#include <array>
#include <valarray>
#include "LinearAlgebra/SmallVector.h"

#include "VTKLagrangeTriangle.h"

namespace hpgem {
namespace Output {
VTKLagrangeTetrahedron::VTKLagrangeTetrahedron(std::size_t order)
    : points_((order + 1) * (order + 2) * (order + 3) / 6) {
    double h = 1.0 / order;

    // The VTK ordering of the points is similar to that of the
    // LagrangeTriangle. The outer points of the tetrahedron are listed in the
    // ordering:
    //  - Corner vertices
    //  - Edge points
    //  - Face points (using triangle ordering)
    // Either this has listed all the points (then stop) or an inner tetrahedron
    // is revealed, for which the process repeats.

    // The approach to generating these is to directly implement this process.
    // But instead of doing it directly for the reference coordinates, do it for
    // scaled barycentric coordinates. Scaling by a factor of 1/h so that all
    // points are at non negative integer coordinates, with all coordinates
    // summing up to 'order'.
    // In addition to be easier to work with, it has the advantage that all
    // points p in the l-th layer satisfy min_i(p[i]) = l, i.e. the smallest
    // coordinate has value l.

    std::size_t numPoints = points_.size();
    // Scaled barycentric coordinates
    using IVec = std::valarray<std::size_t>;

    std::vector<IVec> baryPoints(numPoints);
    std::size_t index = 0;
    for (std::size_t layer = 0; 4 * layer <= order; ++layer) {
        // We split the coordinates into two parts: offset and localPoint.
        // The offset ensures that we are on the right onion-peeling layer of
        // the tetrahedron by being layer-points away from each corner. The
        // localPoint then results in the local point on this layer, which
        // corresponds to having at least one zero coordinate in localPoint.
        // The sum of coordinates in localPoint is fixed by the barycentric
        // constraint at rOrder.
        IVec offset = IVec({layer, layer, layer, layer});
        std::size_t rOrder = order - 4 * layer;

        // Corner points, one or four, each having three zero local coordinates
        // (or four if there is only one point).
        baryPoints[index++] = offset + IVec({rOrder, 0, 0, 0});
        if (rOrder == 0) {
            break;
        }
        baryPoints[index++] = offset + IVec({0, rOrder, 0, 0});
        baryPoints[index++] = offset + IVec({0, 0, rOrder, 0});
        baryPoints[index++] = offset + IVec({0, 0, 0, rOrder});

        // 6 edges, each contributing rOrder - 1 points
        // Each with two zero local coordinates
        for (std::size_t i = 0; i < rOrder - 1; ++i) {
            std::size_t ito = i + 1;
            std::size_t ifrom = rOrder - ito;
            baryPoints[index + 0 * (rOrder - 1) + i] =
                offset + IVec({ifrom, ito, 0, 0});
            baryPoints[index + 1 * (rOrder - 1) + i] =
                offset + IVec({0, ifrom, ito, 0});
            baryPoints[index + 2 * (rOrder - 1) + i] =
                offset + IVec({ito, 0, ifrom, 0});
            baryPoints[index + 3 * (rOrder - 1) + i] =
                offset + IVec({ifrom, 0, 0, ito});
            baryPoints[index + 4 * (rOrder - 1) + i] =
                offset + IVec({0, ifrom, 0, ito});
            baryPoints[index + 5 * (rOrder - 1) + i] =
                offset + IVec({0, 0, ifrom, ito});
        }
        index += 6 * (rOrder - 1);

        // The points on the triangular faces follow the triangular ordering,
        // they have exactly one zero in the local coordinate.

        // Reuse the ordering as generated by the VTKLagrangeTriangle. Note that
        // these are barycentric coordinates for a triangle and thus have
        // length 3. The missing coordinate will be zero.
        std::vector<IVec> facePoints =
            VTKLagrangeTriangle::computeBaryIntegerPoints(rOrder);
        // The face points generated by VTKLagrangeTriangle will include the
        // points on the vertices and edges of the triangle, which are also the
        // points on the vertices and edges of the tetrahedron. Hence, they need
        // to be skipped.
        std::size_t skip = 3 * rOrder;
        std::size_t toAdd = facePoints.size() - skip;

        // The three component barycentric coordinates on the triangle need to
        // be mapped to the four component face coordinates of the tetrahedron.
        // As noted before, this is done by setting one of the local tetrahedron
        // coordinates to zero and using the other three coordinates as is. This
        // requires four mappings, one for each face. They are specified by
        // three indices i, j and k, denoting that the i-th, j-th and k-th
        // coordinate of the local point on the tetrahedron corresponds to the
        // first, second and third barycentric coordinate on the triangle, the
        // fourth (unlisted) coordinate will be left at zero. Thus ensuring the
        // point is on the layer. For example, for this first face we have
        // (i,j,k) -> (i,j,0,k).
        //
        // The ordering of the faces is different from the illustration in [1]
        // and obtained by looking at the resulting geometry in paraview.
        // Incorrect orderings will give distorted geometry, when looking at the
        // result, as the connectivity between the face points to those on the
        // corners/edges of the tetrahedron is done implicitly through the
        // ordering.
        //
        // Note: Each triple would defines a reference face, and is
        // in counter clock wise order.
        //
        // [1]
        // https://blog.kitware.com/modeling-arbitrary-order-lagrange-finite-elements-in-the-visualization-toolkit/
        std::array<std::array<std::size_t, 3>, 4> faceMapping = {
            std::array<std::size_t, 3>({0, 1, 3}),
            std::array<std::size_t, 3>({2, 3, 1}),
            std::array<std::size_t, 3>({0, 3, 2}),
            std::array<std::size_t, 3>({0, 2, 1}),

        };
        for (std::size_t face = 0; face < 4; ++face) {
            for (std::size_t i = 0; i < toAdd; ++i) {
                // Map the local point
                IVec localPoint({0, 0, 0, 0});
                for (std::size_t j = 0; j < 3; ++j) {
                    localPoint[faceMapping[face][j]] = facePoints[i + skip][j];
                }
                // Add it to the output
                baryPoints[index++] = offset + localPoint;
            }
        }
    }

    // Safety check, did we generate exactly the right amount of points
    logger.assert_always(
        index == numPoints,
        "Incorrect number of points for reference tetrahedron");
    // unscale and convert the scaled barycentric to the reference coordinates.
    for (std::size_t i = 0; i < points_.size(); ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            points_[i][j] = h * baryPoints[i][j + 1];
        }
    }
}

}  // namespace Output
}  // namespace hpgem