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
#ifndef HPGEM_LAGRANGEREFERENCEELEMENT_H
#define HPGEM_LAGRANGEREFERENCEELEMENT_H

#include "ReferenceGeometry.h"
#include "ReferencePoint.h"

namespace hpgem {
namespace Geometry {

template <std::size_t dim>
class LagrangeReferenceElement : public ReferenceGeometry {
   public:
    LagrangeReferenceElement(ReferenceGeometry* baseGeometry,
                             std::vector<ReferenceGeometry*> codim1Geometries,
                             std::vector<ReferenceGeometry*> codim2Geometries,
                             std::vector<Geometry::PointReference<dim>> points,
                             const std::string& name);

    bool isInternalPoint(const PointReference<0>& point) const final {
        return baseGeometry_->isInternalPoint(point);
    }
    bool isInternalPoint(const PointReference<1>& point) const final {
        return baseGeometry_->isInternalPoint(point);
    }
    bool isInternalPoint(const PointReference<2>& point) const final {
        return baseGeometry_->isInternalPoint(point);
    }
    bool isInternalPoint(const PointReference<3>& point) const final {
        return baseGeometry_->isInternalPoint(point);
    }
    bool isInternalPoint(const PointReference<4>& point) const final {
        return baseGeometry_->isInternalPoint(point);
    }

    const PointReferenceBase& getCenter() const final {
        return baseGeometry_->getCenter();
    }

    std::size_t getNumberOfNodes() const final { return points_.size(); }

    const PointReferenceBase& getReferenceNodeCoordinate(
        const std::size_t& localIndex) const final {
        logger.assert_debug(localIndex < points_.size(),
                            "Index % larger than the number of nodes %",
                            localIndex, points_.size());
        return points_[localIndex];
    }

    // Mapping codims

    std::size_t getCodim0MappingIndex(
        const std::vector<std::size_t>& v1,
        const std::vector<std::size_t>& v2) const final {
        // Compute the mapping by using the base geometry
        // As v1 and v2 also include all the extra (Lagrange) nodes, we need to
        // select the nodes that correspond to the base geometry.
        std::size_t baseNumberOfNodes = baseGeometry_->getNumberOfNodes();
        std::vector<std::size_t> v1selection(baseNumberOfNodes);
        std::vector<std::size_t> v2selection(baseNumberOfNodes);
        for (std::size_t i = 0; i < baseNumberOfNodes; ++i) {
            v1selection[i] = v1[baseGeometryIndicices_[i]];
            v2selection[i] = v2[baseGeometryIndicices_[i]];
        }
        return baseGeometry_->getCodim0MappingIndex(v1selection, v2selection);
    }

    const MappingReferenceToReference<0>* getCodim0MappingPtr(
        const std::size_t index) const final {
        return baseGeometry_->getCodim0MappingPtr(index);
    }

    std::size_t getNumberOfCodim1Entities() const final {
        if (dim == 1) {
            return points_.size();
        } else {
            return baseGeometry_->getNumberOfCodim1Entities();
        }
    }

    std::vector<std::size_t> getCodim1EntityLocalIndices(
        const std::size_t index) const final {
        logger.assert_debug(index < getNumberOfCodim1Entities(),
                            "Too large face index %", index);
        if (dim > 1) {
            return codim1Indices_[index];
        } else if (dim == 1) {
            return {index};  // Implicit conversion
        } else {
            logger.assert_always(false,
                                 "Dimension too low for codim 1 entities");
            return {};
        }
    }

    const MappingReferenceToReference<1>* getCodim1MappingPtr(
        const std::size_t faceIndex) const final {
        return baseGeometry_->getCodim1MappingPtr(faceIndex);
    }

    const ReferenceGeometry* getCodim1ReferenceGeometry(
        const std::size_t index) const final {
        logger.assert_debug(index < getNumberOfCodim1Entities(),
                            "Too large face number %", index);
        if (dim == 1) {
            return &ReferencePoint::Instance();
        } else {
            return codim1Geometries_[index];
        }
    }

    std::size_t getNumberOfCodim2Entities() const final {
        if (dim == 2) {
            return points_.size();
        } else {
            return baseGeometry_->getNumberOfCodim3Entities();
        }
    }

    std::vector<std::size_t> getCodim2EntityLocalIndices(
        const std::size_t index) const final {
        logger.assert_debug(index < getNumberOfCodim2Entities(),
                            "Too large codim 2 index %", index);
        if (dim > 2) {
            return codim2Indices_[index];
        } else if (dim == 2) {
            return {index};
        } else {
            logger.assert_always(false,
                                 "Dimension too low for codim 2 entities");
            return {};
        }
    }

    const MappingReferenceToReference<2>* getCodim2MappingPtr(
        const std::size_t index) const final {
        return baseGeometry_->getCodim2MappingPtr(index);
    }

    const ReferenceGeometry* getCodim2ReferenceGeometry(
        const std::size_t index) const final {
        logger.assert_debug(index < getNumberOfCodim2Entities(),
                            "Too large codim 2 index %", index);
        if (dim == 2) {
            return &ReferencePoint::Instance();
        } else {
            return codim2Geometries_[index];
        }
    }

    std::size_t getNumberOfCodim3Entities() const final {
        if (dim == 3) {
            return points_.size();
        } else {
            return baseGeometry_->getNumberOfCodim3Entities();
        }
    }

    std::vector<std::size_t> getCodim3EntityLocalIndices(
        const std::size_t index) const final {
        logger.assert_debug(index < getNumberOfCodim3Entities(),
                            "Too large codim 3 index %", index);
        if (dim == 3) {
            return {index};
        } else if (dim > 3) {
            // Don't support 4D
            logger.assert_always(
                false, "Codim 3 not implemented for object of dim %", dim);
            return {};
        } else {
            logger.assert_always(false,
                                 "Dimension too low for codim 3 entities");
            return {};
        }
    }

   private:
    /// The base (non Lagrange) geometry which this extends upon
    ReferenceGeometry* baseGeometry_;
    /// Reference Lagrange geometries for the codim1 entities
    std::vector<ReferenceGeometry*> codim1Geometries_;
    /// Reference Lagrange geometries for the codim2 entities
    std::vector<ReferenceGeometry*> codim2Geometries_;

    /// Lagrange points, in arbitrary order
    std::vector<Geometry::PointReference<dim>> points_;

    /// Indices in points_ of the nodes of the baseGeometry
    /// This should satisfy:
    /// baseGeometry->getReferenceNodeCoordinate(i) ==
    ///     getReferenceNodeCoordinate(baseGeometryIndices[i])
    std::vector<std::size_t> baseGeometryIndicices_;
    /// For each codim 1 entity, the mapping between the reference points on
    /// that entity to to the corresponding reference points on this element.
    /// Only available if dim > 1
    ///
    /// i.e. The following identity should hold:
    /// getReferenceNodeCoordinate(codim1Indices[i][n])
    ///  == getCodim1MappingPtr(i).transform(
    ///     getCodim1ReferenceGeometry(i)->getReferenceNodeCoordinate(n))
    std::vector<std::vector<std::size_t>> codim1Indices_;
    /// Same as the codim1Indices, but for the codim2 entities.
    /// Only available if dim > 2
    std::vector<std::vector<std::size_t>> codim2Indices_;

    /// Find the index of a point in the point list
    /// \param point The point to find
    /// \return The index of that point in points_
    std::size_t findPoint(const PointReference<dim>& point) const;

    void computeCodim1Indices();
    void computeCodim2Indices();
};
}  // namespace Geometry
}  // namespace hpgem

#endif  // HPGEM_LAGRANGEREFERENCEELEMENT_H
