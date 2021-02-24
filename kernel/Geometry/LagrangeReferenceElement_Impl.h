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
#include "LagrangeReferenceElement.h"

#include <limits>

#include "Mappings/MappingReferenceToReference.h"

namespace hpgem {
namespace Geometry {

template <std::size_t dim>
LagrangeReferenceElement<dim>::LagrangeReferenceElement(
    ReferenceGeometry* baseGeometry,
    std::vector<ReferenceGeometry*> codim1Geometries,
    std::vector<ReferenceGeometry*> codim2Geometries,
    std::vector<Geometry::PointReference<dim>> points, const std::string& name)
    : baseGeometry_(baseGeometry),
      codim1Geometries_(std::move(codim1Geometries)),
      codim2Geometries_(std::move(codim2Geometries)),
      points_(std::move(points)),
      // Fill these by invalid values
      baseGeometryIndicices_(baseGeometry->getNumberOfNodes(),
                             std::numeric_limits<std::size_t>::max()),
      codim1Indices_(baseGeometry->getNumberOfCodim1Entities()),
      codim2Indices_(baseGeometry->getNumberOfCodim2Entities()),
      ReferenceGeometry(baseGeometry->getGeometryType(), name) {

    // Verify the correct size of codimDGeometries
    if (dim > 1) {
        logger.assert_debug(codim1Geometries_.size() ==
                                baseGeometry->getNumberOfCodim1Entities(),
                            "Incorrect number of codim 1 geometries");
    } else {
        logger.assert_debug(codim1Geometries_.empty(),
                            "Dimension too low for special codim 1 geometries");
    }
    if (dim > 2) {
        logger.assert_debug(codim2Geometries_.size() ==
                                baseGeometry->getNumberOfCodim2Entities(),
                            "Incorrect number of codim 2 geometries");
    } else {
        logger.assert_debug(codim2Geometries_.empty(),
                            "Dimension too low for special codim 2 geometries");
    }

    // Compute the base geometry indices
    for (std::size_t i = 0; i < baseGeometry->getNumberOfNodes(); ++i) {
        const PointReference<dim>& basePoint =
            baseGeometry->getReferenceNodeCoordinate(i);
        baseGeometryIndicices_[i] = findPoint(basePoint);
    }

    computeCodim1Indices();
    computeCodim2Indices();
}

template <std::size_t dim>
std::size_t LagrangeReferenceElement<dim>::findPoint(
    const PointReference<dim>& point) const {
    // Epsilon for comparing point coordinates. Reference points are of order 1,
    // and should match up to small rounding.
    const double EPS = 1e-8;

    for (std::size_t j = 0; j < points_.size(); ++j) {
        if ((point - points_[j]).getCoordinates().l2NormSquared() < EPS * EPS) {
            return j;
        }
    }
    logger.assert_always(false, "Point not found.");
    return std::numeric_limits<std::size_t>::max();
}

// Computing CODIM-D Indices //
///////////////////////////////
//
// The indices for codim D are only needed if dim > D. They are computed via the
// codim-D mapping in the base geometry. Implementing this requires template
// specialization as the (dim - D) < 0, underflows. This would therefore require
// a PointReference of a very high dimension (that would not be used).
//
// The specialization for dim < D is explicitly written without template magic,
// to keep it readable. This is feasible as we only have dim <= 3 (and therefore
// D < 3).

template <>
void LagrangeReferenceElement<0>::computeCodim1Indices() {}
template <>
void LagrangeReferenceElement<1>::computeCodim1Indices() {}

template <std::size_t dim>
void LagrangeReferenceElement<dim>::computeCodim1Indices() {
    for (std::size_t i = 0; i < getNumberOfCodim1Entities(); ++i) {
        const auto* codim1Mapping = getCodim1MappingPtr(i);
        const auto* codim1Geom = getCodim1ReferenceGeometry(i);
        codim1Indices_[i].resize(codim1Geom->getNumberOfNodes());
        for (std::size_t nodeId = 0; nodeId < codim1Geom->getNumberOfNodes();
             ++nodeId) {
            // Map the node from the codim1 geometry to this geometry
            const PointReference<dim - 1>& codimPoint =
                codim1Geom->getReferenceNodeCoordinate(nodeId);
            const PointReference<dim> point =
                codim1Mapping->transform(codimPoint);
            codim1Indices_[i][nodeId] = findPoint(point);
        }
    }
}

// Nothing to compute for dim=[0,1,2]
template <>
void LagrangeReferenceElement<0>::computeCodim2Indices() {}
template <>
void LagrangeReferenceElement<1>::computeCodim2Indices() {}
template <>
void LagrangeReferenceElement<2>::computeCodim2Indices() {}

template <std::size_t dim>
void LagrangeReferenceElement<dim>::computeCodim2Indices() {
    // Derive the codim2Indices
    for (std::size_t i = 0; i < getNumberOfCodim2Entities(); ++i) {
        const auto* codim2Mapping = getCodim2MappingPtr(i);
        const auto* codim2Geom = getCodim2ReferenceGeometry(i);
        codim2Indices_[i].resize(codim2Geom->getNumberOfNodes());
        for (std::size_t nodeId = 0; nodeId < codim2Geom->getNumberOfNodes();
             ++nodeId) {
            // Map the node from the codim1 geometry to this geometry
            const PointReference<dim - 2>& codimPoint =
                codim2Geom->getReferenceNodeCoordinate(nodeId);
            const PointReference<dim> point =
                codim2Mapping->transform(codimPoint);
            codim2Indices_[i][nodeId] = findPoint(point);
        }
    }
}

}  // namespace Geometry
}  // namespace hpgem
