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
    FaceGeometry::FaceGeometry(ElementGeometry* ptrElemL, const LocalFaceNrType& localFaceNumL, ElementGeometry* ptrElemR, const LocalFaceNrType& localFaceNumR)
            : leftElementGeom_(ptrElemL), rightElementGeom_(ptrElemR), localFaceNumberLeft_(localFaceNumL), localFaceNumberRight_(localFaceNumR), faceToFaceMapIndex_(Geometry::MAXSIZET), faceType_(FaceType::INTERNAL)
    {
        logger.assert(ptrElemL!=nullptr, "Invalid main element passed");
        logger.assert(ptrElemR!=nullptr, "This constructor is intended for internal faces");
    }
    
    //! Constructor for boundary faces.
    FaceGeometry::FaceGeometry(ElementGeometry* ptrElemL, const LocalFaceNrType& localFaceNumL, const FaceType& boundaryLabel)
            : leftElementGeom_(ptrElemL), rightElementGeom_(nullptr), localFaceNumberLeft_(localFaceNumL), localFaceNumberRight_(Geometry::MAXSIZET), faceToFaceMapIndex_(0), faceType_(boundaryLabel)
    {
        logger.assert(ptrElemL!=nullptr, "Invalid main element passed");
    }
    
    FaceGeometry::FaceGeometry(const FaceGeometry& other, 
                               ElementGeometry* ptrElemL, const LocalFaceNrType& localFaceNumL, 
                               ElementGeometry* ptrElemRight, const LocalFaceNrType& localFaceNumR) 
    {
        logger.assert(ptrElemL!=nullptr, "Invalid main element passed");
        leftElementGeom_ = ptrElemL;
        rightElementGeom_ = ptrElemRight;
        localFaceNumberLeft_ = localFaceNumL;
        localFaceNumberRight_ = localFaceNumR;
        faceToFaceMapIndex_ = other.faceToFaceMapIndex_;
        faceType_ = other.faceType_;        
    }

    
    /// The referenceGeometry is returned. It is taken from left element, to always ensure its existence
    // should use the same interface as in the Element!
    const ReferenceGeometry* FaceGeometry::getReferenceGeometry() const
    {
        return leftElementGeom_->getReferenceGeometry()->getCodim1ReferenceGeometry(localFaceNumberLeft_);
    }
    
    /*! Map a point in coordinates of the reference geometry of the face to
     *  the reference geometry of the left (L) element. */
    const PointReference& FaceGeometry::mapRefFaceToRefElemL(const ReferencePointT& pRefFace) const
    {
        return leftElementGeom_->getReferenceGeometry()->getCodim1MappingPtr(localFaceNumberLeft_)->transform(pRefFace);
    }
    
    /*! Map a point in coordinates of the reference geometry of the face to
     *  the reference geometry of the right (R) element. */
    const PointReference& FaceGeometry::mapRefFaceToRefElemR(const ReferencePointT& pRefFace) const
    {
        // In the L function we have assumed the point pRefFace to be
        // given in coordinates of the system used by the reference face
        // on the L side. That means that now to make sure that the
        // point transformed from the right side of the face onto the
        // right element is the same as the one on the left side of the
        // face, we have to use the refFace2RefFace mapping.
        
        return rightElementGeom_->getReferenceGeometry()->getCodim1MappingPtr(localFaceNumberRight_)->transform(mapRefFaceToRefFace(pRefFace));
        
    }
    
    /*! Map from reference face coordinates on the left side to those on the
     *  right side. */
    const PointReference& FaceGeometry::mapRefFaceToRefFace(const ReferencePointT& pIn) const
    {
        return getReferenceGeometry()->getCodim0MappingPtr(faceToFaceMapIndex_)->transform(pIn);
    }
    
    //! Return a Mapping
    FaceGeometry::RefFaceToRefElementMappingPtr FaceGeometry::refFaceToRefElemMapL() const
    {
        return RefFaceToRefElementMappingPtr(leftElementGeom_->getReferenceGeometry()->getCodim1MappingPtr(localFaceNumberLeft_));
    }
    
    //! Return a mapping to the right reference element.
    typename FaceGeometry::RefFaceToRefElementMappingPtr FaceGeometry::refFaceToRefElemMapR() const
    {
        const MappingReferenceToReference * const m1Ptr = this->getReferenceGeometry()->getCodim0MappingPtr(faceToFaceMapIndex_);
        const MappingReferenceToReference * const m2Ptr = leftElementGeom_->getReferenceGeometry()->getCodim1MappingPtr(localFaceNumberRight_);
        
        return RefFaceToRefElementMappingPtr(new ConcatenatedMapping(*m1Ptr, *m2Ptr));
    }
    
    /*!Compute the normal vector at a given point (in face coords).
     
     Faces turn up in DG discretizations as domains for some of the
     integrals. For these integrals, the face must be able to provide the
     geometric information. This is completely included in the normal
     vector, which gives
     <UL>
     <LI> the direction of the normal vector oriented from left (L) to
     right (R) element; that way, for a boundary face, it is an external
     normal vector;
     <LI> the transformation of the integration element between reference
     face and physical space: this scalar is given by the 2-norm of the
     returned normal vector. I.e. the length of the vector equals the ratio of the
     physical face and reference face.
     </UL> */
    LinearAlgebra::NumericalVector FaceGeometry::getNormalVector(const ReferencePointT& pRefFace) const
    {
        std::size_t DIM = pRefFace.size() + 1;
        LinearAlgebra::NumericalVector result;
        if (DIM > 1)
        {
            // first Jacobian (mapping reference face -> reference element)
            
            Jacobian j1 = leftElementGeom_->getReferenceGeometry()->getCodim1MappingPtr(localFaceNumberLeft_) // this is the refFace2RefElemMap
            ->calcJacobian(pRefFace);
            
            // second Jacobian (mapping reference element -> phys. element),
            // for this we first need the points coordinates on the
            // reference element
            
            Jacobian j2 = leftElementGeom_->calcJacobian(mapRefFaceToRefElemL(pRefFace));
            
            Jacobian j3 = j2.multiplyJacobiansInto(j1);
            
            result = j3.computeWedgeStuffVector();
            
            double det = j2.determinant();
            
            double sign = OutwardNormalVectorSign(leftElementGeom_->getReferenceGeometry()->getCodim1MappingPtr(localFaceNumberLeft_));
            result *= ((det > 0) ? 1 : -1) * sign;
        }
        else
        { //if DIM==1
          //for one dimension the fancy wedge stuff wont work
          //but we know the left point has outward pointing vector -1
          //and the right point has outward pointing vector 1
          //so this is the same as the physical point
            const PointReference& pRefElement = mapRefFaceToRefElemL(pRefFace);
            
            //but if someone has mirrored the physical line
            //we also have to also mirror the normal vector
            //the face cant be made larger or smaller so
            //the vector should have length one
            
            Jacobian j = leftElementGeom_->calcJacobian(mapRefFaceToRefElemL(pRefFace));
            int sgn = (j[0] > 0) ? 1 : -1;
            result.resize(DIM);
            result[0] = pRefElement[0] * sgn;
        }
        return result;
        
    }
    
    PointPhysical FaceGeometry::referenceToPhysical(const PointReference& p) const
    {
        return getElementGLeft()->referenceToPhysical(mapRefFaceToRefElemL(p));
    }
    
    //finding node numbers here is way to hard, leave that to someplace else
    void FaceGeometry::initialiseFaceToFaceMapIndex(const std::vector<std::size_t>& leftNodes, const std::vector<std::size_t>& rightNodes)
    {
        logger.assert(leftNodes.size()==rightNodes.size(), "Inconsistent amount of nodes for left and right face");
        faceToFaceMapIndex_ = getReferenceGeometry()->getCodim0MappingIndex(leftNodes, rightNodes);
    }
    
    bool FaceGeometry::isInternal() const
    {
        if (faceType_ == FaceType::INTERNAL || faceType_ == FaceType::PERIODIC_BC || faceType_ == FaceType::SUBDOMAIN_BOUNDARY)
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
