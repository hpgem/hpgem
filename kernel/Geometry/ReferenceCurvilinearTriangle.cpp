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
#include "ReferenceCurvilinearTriangle.h"

#include <map>

#include "ReferenceCurvilinearElement_Impl.h"

#include "ReferenceCurvilinearLine.h"
#include "ReferenceLine.h"
#include "ReferencePoint.h"
#include "ReferenceTriangle.h"

namespace hpgem {
namespace Geometry {

int ReferenceCurvilinearTriangle::getOrderFromPoints(std::size_t numberOfPoints) {
    // We need to invert N = (order + 2)*(order + 1)/2;
    // Using standard mathematics we get
    // order = (-3 + sqrt(9 + 8*(N-1))/2
    double orderD = -1.5 + std::sqrt(9 + 8*(numberOfPoints - 1))/2;
    // Round and check the result using integer arithmetic
    int order = static_cast<int> (std::lround(orderD));
    std::size_t actualPoints = (order + 2) * (order + 1)/2;
    if (actualPoints == numberOfPoints) {
        return order;
    } else {
        return -1;
    }
}

ReferenceCurvilinearTriangle&
    ReferenceCurvilinearTriangle::getReferenceLagrangeTriangle(std::size_t order) {
    static std::map<std::size_t, ReferenceCurvilinearTriangle*> triangles;
    ReferenceCurvilinearTriangle*& triangle = triangles[order];
    if (triangle == nullptr) {
        triangle = new ReferenceCurvilinearTriangle(order);
    }
    return *triangle;
}

std::string getReferenceLagrangeTriangeName(std::size_t order) {
    std::stringstream name;
    name << "Lagrange-ReferenceTriangle-" << order;
    return name.str();
}

std::vector<Geometry::PointReference<2>>
    ReferenceCurvilinearTriangle::createPoints(std::size_t order) {
    double h = 1.0 / order;
    std::size_t numPoints = (order + 2) * (order + 1) / 2;
    std::vector<PointReference<2>> result;
    result.reserve(numPoints);
    for (std::size_t i = 0; i < order + 1; ++i) {
        for (std::size_t j = 0; i + j < order + 1; ++j) {
            result.emplace_back(LinearAlgebra::SmallVector<2>({j * h, i * h}));
        }
    }
    logger.assert_debug(result.size() == numPoints,
                        "Incorrect number of points constructed.");
    return result;
}

ReferenceCurvilinearTriangle::ReferenceCurvilinearTriangle(std::size_t order)
    : ReferenceCurvilinearElement<2>(
          &ReferenceTriangle::Instance(),
          // Codim 1 3 Lagrange lines
          std::vector<ReferenceGeometry*>(
              3, &ReferenceCurvilinearLine::getReferenceLagrangeLine(order)),
          // Codim 2: points -> thus not included
          std::vector<ReferenceGeometry*>(),
          // The actual Lagrange points used for the triangle
          createPoints(order), getReferenceLagrangeTriangeName(order), order) {}

}  // namespace Geometry
}  // namespace hpgem
