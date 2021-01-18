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

#include "Output/VTKLagrangeTetrahedron.h"

#include "../catch.hpp"

TEST_CASE("Tetrahedron points", "[VTKLagrangeTetrahedron]") {
    // The ordering of the Lagrange points in VTKLagrangeTetrahedron is
    // complicated. However, the set of points is easy to describe as the points
    // h*(i,j,k) where h = 1/order and i,j,k are all non negative integers with
    // i + j + k <= order.
    //
    // The idea of this test is to see if all these points are generated in some
    // order (not necessarily the right one for VTK).

    using namespace hpgem;
    using PRef = Geometry::PointReference<3>;

    auto order = GENERATE(1, 2, 3, 4, 5, 6, 7, 8, 9, 10);
    double h = 1.0 / order;

    // GENERATE REFERENCE POINTS //
    ///////////////////////////////
    INFO("Testing order " << order);
    // Tetrahedral number order+1
    std::size_t numPoints = (order + 1) * (order + 2) * (order + 3) / 6;

    std::vector<PRef> orderedPoints(numPoints);
    // Note all <= order
    std::size_t index = 0;
    for (std::size_t i = 0; i <= order; ++i) {
        for (std::size_t j = 0; i + j <= order; ++j) {
            for (std::size_t k = 0; i + j + k <= order; ++k) {
                orderedPoints[index++] = PRef({i * h, j * h, k * h});
            }
        }
    }
    // Check that we did not make a mistake in the generation
    UNSCOPED_INFO("Invalid reference generation");
    REQUIRE(index == numPoints);

    // ACTUAL POINTS //
    ///////////////////

    Output::VTKLagrangeTetrahedron vtkTetrahedron(order);
    std::vector<PRef> vtkPoints(vtkTetrahedron.getPoints());

    UNSCOPED_INFO("Incorrect number of VTK points");
    REQUIRE(vtkPoints.size() == numPoints);

    // COMPARE //
    /////////////

    // Sort, and compare
    std::sort(orderedPoints.begin(), orderedPoints.end());
    std::sort(vtkPoints.begin(), vtkPoints.end());
    for (std::size_t i = 0; i < numPoints; ++i) {
        INFO("Checking point " << i)
        UNSCOPED_INFO("Reference " << orderedPoints[i] << " actual "
                                   << vtkPoints[i]);
        CHECK((vtkPoints[i] - orderedPoints[i]).getCoordinates().l2Norm() <
              1e-8);
    }
}