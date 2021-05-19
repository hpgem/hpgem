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

#include <ElementGeometry.h>

#include "PointReference.h"

namespace hpgem {

namespace Geometry {
class ElementGeometry;

std::ostream& operator<<(std::ostream& os,
                         const ElementGeometry& elementGeometry) {
    os << "PhysicalGeometry={";
    for (std::size_t i = 0;
         i < elementGeometry.physicalGeometry_->getNumberOfNodes(); i++) {
        os << (elementGeometry.physicalGeometry_)->getNodeIndex(i) << " ";
    }
    os << '}' << std::endl;
    os << "Corners: " << std::endl;
    for (std::size_t index :
         elementGeometry.physicalGeometry_->getNodeIndexes()) {
        os << elementGeometry.physicalGeometry_->getGlobalNodeCoordinates(index)
           << std::endl;
    }
    return os;
}

/// Copy constructor that actually makes a deep copy
ElementGeometry::ElementGeometry(const ElementGeometry& other)
    : referenceGeometry_(other.referenceGeometry_),
      refinementMap_(other.refinementMap_)  // refinement is turned off by
                                            // default, to  enable it one needs
                                            // to call enableRefinement
{
    // we have to un-hide the template arguments to make a copy :(
    switch (referenceGeometry_->getGeometryType()) {
        case ReferenceGeometryType::LINE:
            physicalGeometry_ = createPhysicalGeometry<1>(
                other.physicalGeometry_->getNodeIndexes(),
                static_cast<PhysicalGeometry<1>*>(other.physicalGeometry_)
                    ->getNodeCoordinates(),
                referenceGeometry_);
            referenceToPhysicalMapping_ = createMappings<1>(
                other.physicalGeometry_->getNodeIndexes().size(),
                static_cast<PhysicalGeometry<1>*>(physicalGeometry_));
            break;
        case ReferenceGeometryType::SQUARE:
        case ReferenceGeometryType::TRIANGLE:
            physicalGeometry_ = createPhysicalGeometry<2>(
                other.physicalGeometry_->getNodeIndexes(),
                static_cast<PhysicalGeometry<2>*>(other.physicalGeometry_)
                    ->getNodeCoordinates(),
                referenceGeometry_);
            referenceToPhysicalMapping_ = createMappings<2>(
                other.physicalGeometry_->getNodeIndexes().size(),
                static_cast<PhysicalGeometry<2>*>(physicalGeometry_));
            break;
        case ReferenceGeometryType::TETRAHEDRON:
        case ReferenceGeometryType::CUBE:
        case ReferenceGeometryType::TRIANGULARPRISM:
        case ReferenceGeometryType::PYRAMID:
            physicalGeometry_ = createPhysicalGeometry<3>(
                other.physicalGeometry_->getNodeIndexes(),
                static_cast<PhysicalGeometry<3>*>(other.physicalGeometry_)
                    ->getNodeCoordinates(),
                referenceGeometry_);
            referenceToPhysicalMapping_ = createMappings<3>(
                other.physicalGeometry_->getNodeIndexes().size(),
                static_cast<PhysicalGeometry<3>*>(physicalGeometry_));
            break;
        case ReferenceGeometryType::HYPERCUBE:
            physicalGeometry_ = createPhysicalGeometry<4>(
                other.physicalGeometry_->getNodeIndexes(),
                static_cast<PhysicalGeometry<4>*>(other.physicalGeometry_)
                    ->getNodeCoordinates(),
                referenceGeometry_);
            referenceToPhysicalMapping_ = createMappings<4>(
                other.physicalGeometry_->getNodeIndexes().size(),
                static_cast<PhysicalGeometry<4>*>(physicalGeometry_));
            break;
        case ReferenceGeometryType::POINT:
            logger(ERROR, "-1 dimensional faces are not allowed!");
    }
}

ElementGeometry::~ElementGeometry() {
    delete physicalGeometry_;
    delete referenceToPhysicalMapping_;
}

/// Returns a pointer to the referenceToPhysicalMapping

const MappingReferenceToPhysical* ElementGeometry::getReferenceToPhysicalMap()
    const {
    return referenceToPhysicalMapping_;
}

MappingReferenceToPhysical* ElementGeometry::getReferenceToPhysicalMap() {
    return referenceToPhysicalMapping_;
}

/// Returns a pointer to the physicalGeometry object.

const PhysicalGeometryBase* ElementGeometry::getPhysicalGeometry() const {
    return physicalGeometry_;
}

PhysicalGeometryBase* ElementGeometry::getPhysicalGeometry() {
    return physicalGeometry_;
}

/// Returns a pointer to the referenceGeometry object.

const ReferenceGeometry* ElementGeometry::getReferenceGeometry() const {
    return referenceGeometry_;
}

ReferenceGeometry* ElementGeometry::getReferenceGeometry() {
    return referenceGeometry_;
}

/// Returns a pointer to the refinementGeometry object.

const RefinementMapping* ElementGeometry::getRefinementMap() const {
    return refinementMap_;
}

std::size_t ElementGeometry::getNumberOfNodes() const {
    return physicalGeometry_->getNumberOfNodes();
}

std::size_t ElementGeometry::getNrOfNodes() const { return getNumberOfNodes(); }

}  // namespace Geometry

}  // namespace hpgem
