/*
 This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
 http://hpgem.org/about-the-code/team
 
 This code is distributed using BSD 3-Clause License. A copy of which can found below.
 
 
 Copyright (c) 2014, University of Twente
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <climits>
#include <limits>

#include "FaceGeometry.h"

#include "Geometry/PhysicalGeometry.h"
#include "Geometry/ElementGeometry.h"
#include "Geometry/GlobalNamespaceGeometry.h"
#include "Geometry/ReferenceGeometry.h"
#include "Geometry/PointPhysical.h"
#include "Geometry/Mappings/MappingReferenceToReference.h"
#include "Geometry/PointReference.h"
#include "Geometry/Jacobian.h"
#include "Geometry/Mappings/OutwardNormalVectorSign.h"
#include "Geometry/Mappings/ConcatenatedMapping.h"

#include "Logger.h"

namespace Geometry
{
    
    class FaceGeometry;
    
    ///Constructor: define the left and right element, the relative position of 
    ///this face on the elements. Note that this constructor is only for interior
    ///faces. This constructor will not select the correct face to face mapping
    ///because the information needed to do so can be readily generated in the
    ///constructor of Face. (currently Face is the only subclass and FaceGeometry is not
    ///directly used.) If a user want to directly use this class himself, he should be
    ///aware that this constructor assumes he will call initialiseFaceToFaceMapIndex
    ///before actually using the FaceGeometry.
    FaceGeometry::FaceGeometry(ElementGeometry* ptrElemL, const LocalFaceNumberType& localFaceNumberL, ElementGeometry* ptrElemR, const LocalFaceNumberType& localFaceNumberR)
            : leftElementGeom_(ptrElemL), rightElementGeom_(ptrElemR), localFaceNumberLeft_(localFaceNumberL), localFaceNumberRight_(localFaceNumberR), faceToFaceMapIndex_(Geometry::MAXSIZET), faceType_(FaceType::INTERNAL)
    {
        logger.assert(ptrElemL!=nullptr, "Invalid main element passed");
        logger.assert(ptrElemR!=nullptr, "This constructor is intended for internal faces");
    }
    
    //! Constructor for boundary faces.
    FaceGeometry::FaceGeometry(ElementGeometry* ptrElemL, const LocalFaceNumberType& localFaceNumberL, const FaceType& boundaryLabel)
            : leftElementGeom_(ptrElemL), rightElementGeom_(nullptr), localFaceNumberLeft_(localFaceNumberL), localFaceNumberRight_(Geometry::MAXSIZET), faceToFaceMapIndex_(0), faceType_(boundaryLabel)
    {
        logger.assert(ptrElemL!=nullptr, "Invalid main element passed");
    }
    
    FaceGeometry::FaceGeometry(const FaceGeometry& other, 
                               ElementGeometry* ptrElemL, const LocalFaceNumberType& localFaceNumberL, 
                               ElementGeometry* ptrElemRight, const LocalFaceNumberType& localFaceNumberR) 
    {
        logger.assert(ptrElemL!=nullptr, "Invalid main element passed");
        leftElementGeom_ = ptrElemL;
        rightElementGeom_ = ptrElemRight;
        localFaceNumberLeft_ = localFaceNumberL;
        localFaceNumberRight_ = localFaceNumberR;
        faceToFaceMapIndex_ = other.faceToFaceMapIndex_;
        faceType_ = other.faceType_;        
    }

    
    /// The referenceGeometry is returned. It is taken from left element, to always ensure its existence
    // should use the same interface as in the Element!
    const ReferenceGeometry* FaceGeometry::getReferenceGeometry() const
    {
        return leftElementGeom_->getReferenceGeometry()->getCodim1ReferenceGeometry(localFaceNumberLeft_);
    }
    
    //! Return a Mapping
    FaceGeometry::RefFaceToRefElementMappingPtr FaceGeometry::refFaceToRefElemMapL() const
    {
        //we abuse shared_ptr here to provide common behaviour with refFaceToRefElemMapR, but we don't actually want to delete the mapping when the shared_ptr is destroyed
        return RefFaceToRefElementMappingPtr(leftElementGeom_->getReferenceGeometry()->getCodim1MappingPtr(localFaceNumberLeft_), [](const MappingReferenceToReference<1>*){});
    }
    
    //! Return a mapping to the right reference element.
    typename FaceGeometry::RefFaceToRefElementMappingPtr FaceGeometry::refFaceToRefElemMapR() const
    {
        const MappingReferenceToReference<0> * const m1Ptr = this->getReferenceGeometry()->getCodim0MappingPtr(faceToFaceMapIndex_);
        const MappingReferenceToReference<1> * const m2Ptr = rightElementGeom_->getReferenceGeometry()->getCodim1MappingPtr(localFaceNumberRight_);
        
        return RefFaceToRefElementMappingPtr(new ConcatenatedMapping(*m1Ptr, *m2Ptr));
    }
    
    //finding node numbers here is way to hard, leave that to someplace else
    void FaceGeometry::initialiseFaceToFaceMapIndex(const std::vector<std::size_t>& leftNodes, const std::vector<std::size_t>& rightNodes)
    {
        logger.assert(leftNodes.size()==rightNodes.size(), "Inconsistent amount of nodes for left and right face");
        faceToFaceMapIndex_ = getReferenceGeometry()->getCodim0MappingIndex(leftNodes, rightNodes);
    }
    
    bool FaceGeometry::isInternal() const
    {
        if (faceType_ == FaceType::INTERNAL || faceType_ == FaceType::PERIODIC_BC || faceType_ == FaceType::SUBDOMAIN_BOUNDARY || faceType_ == FaceType::PERIODIC_SUBDOMAIN_BC)
        {
            ///\todo move this assertion to the unit test
            logger.assert(rightElementGeom_ != nullptr, "There is no right element, so no internal face.");
            return true;
        }
        else
        {
            logger.assert(rightElementGeom_ == nullptr, "There is a right element, so no boundary face.");
            return false;
        }
    }

}
