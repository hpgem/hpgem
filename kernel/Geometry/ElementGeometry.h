/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2014, University of Twente
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

#ifndef HPGEM_KERNEL_ELEMENTGEOMETRY_H
#define HPGEM_KERNEL_ELEMENTGEOMETRY_H

#include <vector>
#include <iostream>
#include "Point.h"
#include "PointPhysical.h"
#include "Jacobian.h"
#include "Mappings/MappingReferenceToPhysical.h"

#include "ReferenceGeometryFactory.h"
#include "Mappings/RefinementMapping.h"
#include "Mappings/MappingReferenceToPhysical.h"
#include "Mappings/MappingToPhysHypercubeLinear.h"
#include "Mappings/MappingToPhysSimplexLinear.h"
#include "Mappings/MappingToPhysPyramid.h"
#include "Mappings/MappingToPhysTetrahedronQuadratic.h"
#include "Mappings/MappingToPhysTriangularPrism.h"
#include "Mappings/MappingToPhysTriangleQuadratic.h"

#include "PointReference.h"

namespace hpgem {

namespace Geometry {
template <std::size_t DIM>
class PointReference;
template <std::size_t DIM>
class PointPhysical;
class MappingReferenceToPhysicalBase;
class PhysicalGeometryBase;
class ReferenceGeometry;
class RefinementGeometry;
template <std::size_t dimFrom, std::size_t dimTo>
class Jacobian;

class ElementGeometry {
   public:
    /// New style constructor with one less pass
    template <std::size_t DIM>
    ElementGeometry(const std::vector<std::size_t>& globalNodeIndexes,
                    std::vector<PointPhysical<DIM> >& nodes);

    /// Copy constructor
    ElementGeometry(const ElementGeometry& other);

    virtual ~ElementGeometry();

    /// Returns a pointer to the referenceToPhysicalMapping
    const MappingReferenceToPhysicalBase* getReferenceToPhysicalMap() const;
    MappingReferenceToPhysicalBase* getReferenceToPhysicalMap();

    /// Returns a pointer to the physicalGeometry object.
    const PhysicalGeometryBase* getPhysicalGeometry() const;

    /// Returns a pointer to the physicalGeometry object.
    PhysicalGeometryBase* getPhysicalGeometry();

    /// Gets the number of nodes associated with this element.
    std::size_t getNumberOfNodes() const;

    ///\deprecated Does not follow naming conventions, use getNumberOfNodes()
    /// instead.
    std::size_t getNrOfNodes() const;

    /// Returns a pointer to the referenceGeometry object.
    const ReferenceGeometry* getReferenceGeometry() const;

    ReferenceGeometry* getReferenceGeometry();

    /// Returns a pointer to the refinementMapping object.
    const RefinementMapping* getRefinementMap() const;

    /// This method gets a PointReference, which specifies a coordinate in the
    /// ReferenceGeometry, and returns a PointPhysical which is the
    /// corresponding point in the PhysicalGeometry, given the mapping.
    template <std::size_t DIM>
    PointPhysical<DIM> referenceToPhysical(
        const PointReference<DIM>& pointReference) const;

    /// This routine is the inverse of referenceToPhysical. For elements where
    /// the mapping from the reference Element is nonlinear (square, cube,
    /// triangularPrism and pyramid) there might be multiple valid reference
    /// points some of which might be interior to the reference element. In this
    /// case a reference point is selected depending on implementation details.
    template <std::size_t DIM>
    PointReference<DIM> physicalToReference(
        const PointPhysical<DIM>& pointPhysical) const;

    /// This method gets a PointReference and returns the corresponding jacobian
    /// of the referenceToPhysicalMapping.
    template <std::size_t DIM>
    Jacobian<DIM, DIM> calcJacobian(
        const PointReference<DIM>& pointReference) const;

    void enableRefinement();

    /// Output operator.
    friend std::ostream& operator<<(std::ostream& os,
                                    const ElementGeometry& elementGeometry);

   private:
    template <std::size_t DIM>
    static ReferenceGeometry* createReferenceGeometry(std::size_t size);

    template <std::size_t DIM>
    static PhysicalGeometry<DIM>* createPhysicalGeometry(
        const std::vector<std::size_t>& globalNodeIndexes,
        std::vector<PointPhysical<DIM> >& nodes,
        const ReferenceGeometry* const geo);

    template <std::size_t DIM>
    static MappingReferenceToPhysicalBase* createMappings(
        std::size_t size, const PhysicalGeometry<DIM>* const pGeo);

   protected:
    /// The corresponding referenceGeometry object, for integration.
    ReferenceGeometry* const referenceGeometry_;

    /// The physicalGeometry object contains pointers to the actual physical
    /// points, and a container of global node indexes.
    PhysicalGeometryBase* physicalGeometry_;

    /// The referenceToPhysicalMapping relates the coordinates of the reference
    /// object to the physical object; basically a matrix transformation.
    MappingReferenceToPhysicalBase* referenceToPhysicalMapping_;

    /// The corresponding refinementGeometry object
    RefinementMapping* refinementMap_;
};

/// This method gets a PointReference, which specifies a coordinate in the
/// ReferenceGeometry, and returns a PointPhysical which is the corresponding
/// point in the PhysicalGeometry, given the mapping.
template <std::size_t DIM>
PointPhysical<DIM> ElementGeometry::referenceToPhysical(
    const PointReference<DIM>& pointReference) const {
    return referenceToPhysicalMapping_->castDimension<DIM>().transform(
        pointReference);
}

template <std::size_t DIM>
PointReference<DIM> ElementGeometry::physicalToReference(
    const PointPhysical<DIM>& pointPhysical) const {
    return referenceToPhysicalMapping_->castDimension<DIM>().inverseTransform(
        pointPhysical);
}

/// This method gets a PointReference and returns the corresponding Jacobian of
/// the referenceToPhysicalMapping.
template <std::size_t DIM>
Jacobian<DIM, DIM> ElementGeometry::calcJacobian(
    const PointReference<DIM>& pointReference) const {
    return referenceToPhysicalMapping_->castDimension<DIM>().calcJacobian(
        pointReference);
}

/// Create the reference element for the given number of nodes. Since this
/// method is templated on the dimension, there is always a unique reference
/// geometry for the given number of nodes. This method then returns a pointer
/// to the relevant ReferenceGeometry, for example the reference triangle if the
/// given number of nodes equals 3.
template <std::size_t DIM>
ReferenceGeometry* ElementGeometry::createReferenceGeometry(std::size_t size) {
    return &ReferenceGeometryFactory::Instance().getGeometry(DIM, size);
}

template <std::size_t DIM>
PhysicalGeometry<DIM>* ElementGeometry::createPhysicalGeometry(
    const std::vector<std::size_t>& globalNodeIndexes,
    std::vector<PointPhysical<DIM> >& nodes,
    const ReferenceGeometry* const geo) {
    logger.assert_debug(geo != nullptr, "Invalid reference geometry passed");
    return new PhysicalGeometry<DIM>(globalNodeIndexes, nodes, geo);
}

template <std::size_t DIM>
MappingReferenceToPhysicalBase* ElementGeometry::createMappings(
    std::size_t size, const PhysicalGeometry<DIM>* const pGeo) {
    logger(ERROR, "DIM may range from 1 to 4, but it seems to be %", DIM);
    return nullptr;
}

template <>
inline MappingReferenceToPhysicalBase* ElementGeometry::createMappings(
    std::size_t size, const PhysicalGeometry<1>* const pGeo) {
    logger.assert_debug(pGeo != nullptr, "Invalid physical geometry passed");
    logger.assert_debug(size == 2, "1D can only map to a line");
    logger(VERBOSE, "ElementGeometry created a mapping for a line.");
    return new Geometry::MappingToPhysHypercubeLinear<1>(pGeo);
}

template <>
inline MappingReferenceToPhysicalBase* ElementGeometry::createMappings(
    std::size_t size, const PhysicalGeometry<2>* const pGeo) {
    logger.assert_debug(pGeo != nullptr, "Invalid physical geometry passed");
    switch (size) {
        case 3:
            logger(VERBOSE,
                   "ElementGeometry created a mapping for a triangle.");
            return new Geometry::MappingToPhysSimplexLinear<2>(pGeo);
        case 4:
            logger(VERBOSE, "ElementGeometry created a mapping for a square.");
            return new Geometry::MappingToPhysHypercubeLinear<2>(pGeo);
        case 6:
            logger(
                VERBOSE,
                "ElementGeometry created a mapping for a quadratic triangle.");
            return new Geometry::MappingToPhysTriangleQuadratic(pGeo);
    }
    logger(FATAL, "No know entities contain this many nodes. \n");
    return nullptr;
}

template <>
inline MappingReferenceToPhysicalBase* ElementGeometry::createMappings(
    std::size_t size, const PhysicalGeometry<3>* const pGeo) {
    logger.assert_debug(pGeo != nullptr, "Invalid physical geometry passed");
    switch (size) {
        case 4:
            logger(VERBOSE,
                   "ElementGeometry created a mapping for a tetrahedron.");
            return new Geometry::MappingToPhysSimplexLinear<3>(pGeo);
        case 5:
            logger(VERBOSE, "ElementGeometry created a mapping for a pyramid.");
            return new Geometry::MappingToPhysPyramid(pGeo);
        case 6:
            logger(VERBOSE,
                   "ElementGeometry created a mapping for a triangular prism.");
            return new Geometry::MappingToPhysTriangularPrism(pGeo);
        case 8:
            logger(VERBOSE, "ElementGeometry created a mapping for a cube.");
            return new Geometry::MappingToPhysHypercubeLinear<3>(pGeo);
        case 10:
            logger(VERBOSE,
                   "ElementGeometry created a mapping for a quadratic "
                   "tetrahedron.");
            return new Geometry::MappingToPhysTetrahedronQuadratic(pGeo);
    }
    logger(FATAL, "No know entities contain this many nodes. \n");
    return nullptr;
}

template <>
inline MappingReferenceToPhysicalBase* ElementGeometry::createMappings(
    std::size_t size, const PhysicalGeometry<4>* const pGeo) {
    logger.assert_debug(pGeo != nullptr, "Invalid physical geometry passed");
    logger.assert_debug(size == 16, "4D can only map to a hypercube");
    logger(VERBOSE, "ElementGeometry created a mapping for a hypercube.");
    return new Geometry::MappingToPhysHypercubeLinear<4>(pGeo);
}

template <std::size_t DIM>
ElementGeometry::ElementGeometry(
    const std::vector<std::size_t>& globalNodeIndexes,
    std::vector<PointPhysical<DIM> >& nodes)
    : referenceGeometry_(ElementGeometry::createReferenceGeometry<DIM>(
          globalNodeIndexes.size())),
      physicalGeometry_(ElementGeometry::createPhysicalGeometry(
          globalNodeIndexes, nodes, referenceGeometry_)),
      referenceToPhysicalMapping_(ElementGeometry::createMappings<DIM>(
          globalNodeIndexes.size(),
          static_cast<PhysicalGeometry<DIM>*>(physicalGeometry_))),
      refinementMap_(nullptr)  // refinement is turned off by default, to enable
                               // it one needs to call enableRefinement
{}
}  // namespace Geometry
}  // namespace hpgem

#endif  // HPGEM_KERNEL_ELEMENTGEOMETRY_H
