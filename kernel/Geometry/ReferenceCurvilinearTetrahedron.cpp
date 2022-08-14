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
#include "ReferenceCurvilinearTetrahedron.h"

#include <map>

#include "ReferenceCurvilinearElement_Impl.h"

#include "ReferenceCurvilinearLine.h"
#include "ReferenceCurvilinearTriangle.h"
#include "ReferenceLine.h"
#include "ReferencePoint.h"
#include "ReferenceTetrahedron.h"
#include "ReferenceTriangle.h"

namespace hpgem {
namespace Geometry {

int ReferenceCurvilinearTetrahedron::getOrderFromPoints(
    std::size_t numberOfPoints) {
    std::size_t tetraNumber = 1;
    int order = 0;
    while (tetraNumber <= numberOfPoints) {
        if (tetraNumber == numberOfPoints) {
            return order;
        }
        // Add triangular number
        order += 1;
        tetraNumber += ((order + 1) * (order + 2)) / 2;
    }
    return -1;
}

ReferenceCurvilinearTetrahedron&
    ReferenceCurvilinearTetrahedron::getReferenceCurvilinearTetrahedron(
        std::size_t order) {
    static std::map<std::size_t, ReferenceCurvilinearTetrahedron*> tetrahedra;
    auto*& tetrahedron = tetrahedra[order];
    if (tetrahedron == nullptr) {
        tetrahedron = new ReferenceCurvilinearTetrahedron(order);
    }
    return *tetrahedron;
}
namespace {
std::string getReferenceTetrahedronName(std::size_t order) {
    std::stringstream name;
    name << "Curvilinear-ReferenceTetrahedron-" << order;
    return name.str();
}
}  // namespace

std::vector<Geometry::PointReference<3>>
    ReferenceCurvilinearTetrahedron::createPoints(std::size_t order) {
    double h = 1.0 / order;
    std::size_t numberOfPoints = (order + 1) * (order + 2) * (order + 3) / 6;
    std::vector<PointReference<3>> result;
    result.reserve(numberOfPoints);
    for (std::size_t iz = 0; iz < order + 1; ++iz) {
        for (std::size_t iy = 0; iy + iz < order + 1; ++iy) {
            for (std::size_t ix = 0; ix + iy + iz < order + 1; ++ix) {
                result.emplace_back(
                    LinearAlgebra::SmallVector<3>({ix * h, iy * h, iz * h}));
            }
        }
    }
    logger.assert_debug(result.size() == numberOfPoints,
                        "Incorrect number of points");
    return result;
}

ReferenceCurvilinearTetrahedron::ReferenceCurvilinearTetrahedron(
    std::size_t order)
    : ReferenceCurvilinearElement<3>(
          &ReferenceTetrahedron::Instance(),
          // Codim 1 => Triangle sides
          std::vector<ReferenceGeometry*>(
              4, &ReferenceCurvilinearTriangle::getReferenceLagrangeTriangle(
                     order)),
          // Codim 2 => Lines
          std::vector<ReferenceGeometry*>(
              6, &ReferenceCurvilinearLine::getReferenceLagrangeLine(order)),
          // points, name, order
          createPoints(order), getReferenceTetrahedronName(order), order) {}

}  // namespace Geometry
}  // namespace hpgem