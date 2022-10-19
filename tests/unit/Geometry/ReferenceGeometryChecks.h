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
#ifndef HPGEM_REFERENCEGEOMETRYCHECKS_H
#define HPGEM_REFERENCEGEOMETRYCHECKS_H

#include "Geometry/ReferenceGeometry.h"
#include "../catch.hpp"
#include <memory>
#include <utility>

namespace hpgem {

template <std::size_t dim>
void testMeasure(const Geometry::ReferenceGeometry& geom) {
    // Test the measure of the geometry
    // Approach:
    // Compute the bounding box of the reference element. Create an
    // equidistant grid of point for each of the dimensions, and
    // compute the fraction of points inside the reference element.
    // The measure should equal this fraction of the measure of the
    // bounding box.

    using namespace hpgem::Geometry;
    using namespace hpgem::LinearAlgebra;
    using PRef = Geometry::PointReference<dim>;
    using Vec = LinearAlgebra::SmallVector<dim>;

    // Compute a bounding box
    PRef minPoint, maxPoint;
    for (std::size_t i = 0; i < dim; ++i) {
        minPoint[i] = std::numeric_limits<double>::infinity();
        maxPoint[i] = -std::numeric_limits<double>::infinity();
    }
    for (std::size_t i = 0; i < geom.getNumberOfNodes(); ++i) {
        PRef refPoint = geom.getReferenceNodeCoordinate(i);
        for (std::size_t j = 0; j < dim; ++j) {
            minPoint[j] = std::min(minPoint[j], refPoint[j]);
            maxPoint[j] = std::max(maxPoint[j], refPoint[j]);
        }
    }
    // Compute bounding box measure and side lengths
    double bbMeasure = 1.0;
    Vec sideLength;
    for (std::size_t i = 0; i < dim; ++i) {
        sideLength[i] = maxPoint[i] - minPoint[i];
        bbMeasure *= sideLength[i];
    }
    // Test how many points of the bounding box are inside
    std::size_t testPointsPerDimension = 100;
    std::size_t totalTestingPoints = 1;
    Vec offset, stride;
    for (std::size_t i = 0; i < dim; ++i) {
        stride[i] = sideLength[i] / testPointsPerDimension;
        offset[i] = minPoint[i] + stride[i] / 2.0;
    }
    // Numerical index of the point to test
    // pref[i] == offset[i] + index[i] * stride[i]
    std::array<std::size_t, dim> index;
    for (std::size_t i = 0; i < dim; ++i) {
        index[i] = 0;
        totalTestingPoints *= testPointsPerDimension;
    }
    PRef pref = PRef(offset);  // Actual point to test

    std::size_t pointsInside = 0;
    while (true) {
        if (geom.isInternalPoint(pref)) {
            pointsInside++;
        }
        // Lexicographical increment
        std::size_t i = 0;
        for (; i < dim; ++i) {
            if (index[i] < testPointsPerDimension - 1) {
                index[i]++;
                pref[i] += stride[i];
                break;
            } else {
                pref[i] = offset[i];
                index[i] = 0;
            }
        }
        if (i == dim) {
            break;
        }
    }
    // Compute measure
    double expectedMeasure = (bbMeasure * pointsInside) / totalTestingPoints;
    // Epsilon reasoning:
    // Draw lines along the coordinate axes through all the test points
    // - There are dim directions
    // - Each line intersects the surface of the reference volume at zero or two
    //   points
    // - The error in both points is approximately
    //   bbMeasure/testPointsPerDimension
    // Additional note: Measure is probably hard coded according to a logical
    // formula. Thus, large errors are far more likely that small rounding
    // errors.
    double eps = 2.0 * dim * bbMeasure / testPointsPerDimension;
    REQUIRE(geom.measure() == Approx(expectedMeasure).epsilon(eps));
}

/**
 * Test the correctness of the codimMappingIndex + codim0MappingPtr pair
 * by trying out each possible node permutation and testing whether the corner
 * nodes of the reference shape are correctly mapped.
 *
 * @tparam dim The dimension of the shape
 * @param geom The reference shape to test
 * @param nextPermutation Function that computes the next valid permutation of
 * the nodes. For a geometry with N nodes this takes a valid permutation of the
 * indices [0..N-1] and should update the vector in place with the next valid
 * permutation. It should return true if this was a new permutation, and false
 * if the input was the last permutation.
 */
template <std::size_t dim>
void testCodim0Mapping(
    const Geometry::ReferenceGeometry& geom,
    const std::function<bool(std::vector<std::size_t>&)>& nextPermutation) {

    // Reference node order
    std::vector<std::size_t> reference;
    std::vector<std::size_t> permutation;
    std::vector<Geometry::PointReference<dim>> points;
    for (std::size_t i = 0; i < geom.getNumberOfNodes(); ++i) {
        reference.push_back(i);
        permutation.push_back(i);
        points.push_back(geom.getReferenceNodeCoordinate(i));
    }
    do {
        std::size_t index = geom.getCodim0MappingIndex(reference, permutation);
        const auto* mapping = geom.getCodim0MappingPtr(index);
        INFO("Checking mapping " << index);
        CHECK(mapping != nullptr);
        for (std::size_t i = 0; i < geom.getNumberOfNodes(); ++i) {
            // Note the ordering here
            Geometry::PointReference<dim> mpoint =
                mapping->transform(points[permutation[i]]);
            auto diff = mpoint - points[i];
            REQUIRE(diff.l2NormSquared() < 1e-24);
        }
    } while (nextPermutation(permutation));
}

}  // namespace hpgem

#endif  // HPGEM_REFERENCEGEOMETRYCHECKS_H
