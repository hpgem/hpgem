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
#ifndef HPGEM_REFERENCECURVILINEARELEMENT_H
#define HPGEM_REFERENCECURVILINEARELEMENT_H

#include "ReferenceGeometry.h"
#include "ReferencePoint.h"

namespace hpgem {
namespace Geometry {

/// \brief Dimensionless base class for ReferenceCurvilinearElement
class ReferenceCurvilinearElementBase : public ReferenceGeometry {
   protected:
    ReferenceCurvilinearElementBase(const ReferenceGeometryType& geo,
                                    std::string name, std::size_t order)
        : ReferenceGeometry(geo, std::move(name)), order_(order){};

   public:
    /// The linear version of this curvilinear element
    /// \return
    virtual ReferenceGeometry* getBaseGeometry() const = 0;
    /// The polynomial order of this curvilinear element
    std::size_t getOrder() const { return order_; }

    double measure() const final { return getBaseGeometry()->measure(); }

   private:
    std::size_t order_;
};

/// Parent class for Curvilinear elements
///
/// Each curvilinear element is defined by
///  - A linear base geometry (e.g. ReferenceTriangle)
///  - A reference geometry for the boundary
///  - A set of reference points
///
/// Note: For curvilinear elements there are usually multiple nodes on each
/// side. Hence, for dim-dimensional curvilinear element we have that the number
/// of codim-dim entities is smaller than the number of nodes. As the codim-dim
/// entities are the nodes of the base geometry.
///
/// Implementation:
///  - The mappings of the linear element are reused
///  - The mappings are used to compute the which reference points correspond to
///    what reference points on the faces and other codim-X entities.
/// \tparam dim
template <std::size_t dim>
class ReferenceCurvilinearElement : public ReferenceCurvilinearElementBase {
   public:
    /// Create a curvilinear element based on a linear element
    /// \param baseGeometry The geometry of the base element
    /// \param codim1Geometries
    ///   For dim > 1, the description of the codim 1 entities on the boundary.
    ///   These should follow the same order as the ones in the base element,
    ///   but correspond to the curvilinear nature of the boundary.
    /// \param codim2Geometries
    ///  Same as codim1Geometries for codim 2, for dim > 2.
    /// \param points The reference points, these should include:
    ///   - The reference points of the base geometry
    ///   - The reference points on each codim1/codim2 geometry mapped by the
    ///     corresponding mapping.
    /// \param name A name for this element (as required by ReferenceGeometry)
    /// \param order The order of the element
    ReferenceCurvilinearElement(
        ReferenceGeometry* baseGeometry,
        std::vector<ReferenceGeometry*> codim1Geometries,
        std::vector<ReferenceGeometry*> codim2Geometries,
        std::vector<Geometry::PointReference<dim>> points,
        const std::string& name, std::size_t order);

    ReferenceGeometry* getBaseGeometry() const final { return baseGeometry_; }

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

    const PointReference<dim>& getCenter() const final {
        return baseGeometry_->getCenter();
    }

    std::size_t getNumberOfNodes() const final { return points_.size(); }

    const PointReference<dim>& getReferenceNodeCoordinate(
        const std::size_t& localIndex) const final {
        logger.assert_debug(localIndex < points_.size(),
                            "Index % larger than the number of nodes %",
                            localIndex, points_.size());
        return points_[localIndex];
    }

    // Mapping codims

    // CODIM 0 //
    // ------- //

    std::size_t getCodim0MappingIndex(
        const std::vector<std::size_t>& v1,
        const std::vector<std::size_t>& v2) const final;

    const MappingReferenceToReference<0>* getCodim0MappingPtr(
        const std::size_t index) const final {
        return baseGeometry_->getCodim0MappingPtr(index);
    }

    // CODIM 1 //
    // ------- //

    std::size_t getNumberOfCodim1Entities() const final {
        return baseGeometry_->getNumberOfCodim1Entities();
    }

    std::vector<std::size_t> getCodim1EntityLocalIndices(
        const std::size_t index) const final;

    const BoundaryFaceMapping* getCodim1MappingPtr(
        const std::size_t faceIndex) const final {
        return baseGeometry_->getCodim1MappingPtr(faceIndex);
    }

    const ReferenceGeometry* getCodim1ReferenceGeometry(
        const std::size_t index) const final;

    // CODIM 2 //
    // ------- //

    std::size_t getNumberOfCodim2Entities() const final {
        return baseGeometry_->getNumberOfCodim2Entities();
    }

    std::vector<std::size_t> getCodim2EntityLocalIndices(
        const std::size_t index) const final;

    const MappingReferenceToReference<2>* getCodim2MappingPtr(
        const std::size_t index) const final {
        return baseGeometry_->getCodim2MappingPtr(index);
    }

    const ReferenceGeometry* getCodim2ReferenceGeometry(
        const std::size_t index) const final;

    // CODIM 3 //
    // ------- //

    std::size_t getNumberOfCodim3Entities() const final {
        return baseGeometry_->getNumberOfCodim3Entities();
    }

    std::vector<std::size_t> getCodim3EntityLocalIndices(
        const std::size_t index) const final;

   private:
    /// The base linear geometry which this extends upon
    ReferenceGeometry* baseGeometry_;
    /// Reference Curvilinear geometries for the codim1 entities
    std::vector<ReferenceGeometry*> codim1Geometries_;
    /// Reference Curvilinear geometries for the codim2 entities
    std::vector<ReferenceGeometry*> codim2Geometries_;

    /// Reference points that are used to define the mapping between reference
    /// and physical elements. The following assumptions are made:
    ///  - Each reference point is unique
    ///  - The reference points of the baseGeometry are a subset
    ///  - The reference points of the curvilinear boundary geometries map onto
    ///    reference points of this element.
    std::vector<Geometry::PointReference<dim>> points_;

    /// Indices in points_ of the nodes of the baseGeometry
    /// This should satisfy:
    /// baseGeometry->getReferenceNodeCoordinate(i) ==
    ///     getReferenceNodeCoordinate(baseGeometryIndices[i])
    std::vector<std::size_t> baseGeometryIndicices_;

    /// The mapping of the reference points of each codim-1 element, to the
    /// corresponding reference point on this element. Only available if
    /// dim > 1.
    ///
    /// The mapping is completely determined by the identity:
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

    /// Compute codim1Indices_
    void computeCodim1Indices();
    /// Compute codim2Indices_
    void computeCodim2Indices();
};
}  // namespace Geometry
}  // namespace hpgem

#endif  // HPGEM_REFERENCECURVILINEARELEMENT_H
