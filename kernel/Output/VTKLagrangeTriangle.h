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
#ifndef HPGEM_VTKLAGRANGETRIANGLE_H
#define HPGEM_VTKLAGRANGETRIANGLE_H

#include "VTKElement.h"

#include <valarray>

namespace hpgem {
namespace Output {
class VTKLagrangeTriangle final : public VTKElement<2> {
   public:
    explicit VTKLagrangeTriangle(std::size_t order);

    std::uint8_t vtkId() const final { return 69; }

    const std::vector<Geometry::PointReference<2>>& getPoints() const final {
        return points_;
    }

    /// Compute the scaled barycentric coordinates of the Lagrange points.
    ///
    /// This computes the Lagrange points in VTK ordering for a given
    /// (polynomial) order. The coordinates for the points are given in scaled
    /// barycentric coordinates, that is bary centric coordinates scaled by a
    /// factor order, so that the coordinates are all non negative integers. As
    /// these are barycentric coordinates the sum of the coordinates sum up to
    /// order.
    ///
    /// The VTK ordering for a triangle is recursive. First one lists the outer
    /// corner points. Then in CCW direction the points along the edge. This
    /// either exhausts the points or reveals an inner triangle (possibly a
    /// single point) for which the process can be repeated. The first three
    /// nodes in this process are [order, 0, 0], [0, order, 0] and [0, 0,
    /// order].
    ///
    /// \param order The order for which to generate the Lagrange points.
    /// \return A vector of size (order+1)*(order+2)/2 with the scaled
    /// barycentric coordinates.
    static std::vector<std::valarray<std::size_t>> computeBaryIntegerPoints(
        std::size_t order);

   private:
    std::vector<Geometry::PointReference<2>> points_;
};
}  // namespace Output
}  // namespace hpgem

#endif  // HPGEM_VTKLAGRANGETRIANGLE_H
