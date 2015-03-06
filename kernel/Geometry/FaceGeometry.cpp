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

#include "FaceGeometry.hpp"

#include "Geometry/PhysicalGeometry.hpp"
#include "Geometry/ElementGeometry.hpp"
#include "Geometry/GlobalNamespaceGeometry.hpp"
#include "Geometry/ReferenceGeometry.hpp"
#include "Geometry/PointPhysical.hpp"
#include "Geometry/Mappings/MappingReferenceToReference.hpp"
#include "Geometry/PointReference.hpp"
#include "Geometry/Jacobian.hpp"
#include "Geometry/Mappings/OutwardNormalVectorSign.hpp"
#include "Geometry/Mappings/ConcatenatedMapping.hpp"

#include "Logger.h"

namespace Geometry
{

    class FaceGeometry;
    
    ///Constructor: define the left and right element, the relative position of 
    ///this face on the elements. Note that this constructor is only for interior
    ///faces. 
    ///\todo Find out what faceToFaceMapIndex_ does and why it takes UINT_MAX
    FaceGeometry::FaceGeometry(ElementGeometryT* ptrElemL,
                               const LocalFaceNrType& localFaceNumL,
                               ElementGeometryT* ptrElemR,
                               const LocalFaceNrType& localFaceNumR) :
    leftElementGeom_(ptrElemL),
    localFaceNumberLeft_(localFaceNumL),
    rightElementGeom_(ptrElemR),
    localFaceNumberRight_(localFaceNumR),
    faceType_(FaceType::INTERNAL),
    faceToFaceMapIndex_(Geometry::MAXSIZET) { }


    //! Constructor for boundary faces.
    FaceGeometry::FaceGeometry(ElementGeometryT* ptrElemL, 
                               const LocalFaceNrType& localFaceNumL, 
                               const FaceType& boundaryLabel) :
    leftElementGeom_(ptrElemL),
    rightElementGeom_(nullptr),
    localFaceNumberLeft_(localFaceNumL),
    localFaceNumberRight_(Geometry::MAXSIZET),
    faceToFaceMapIndex_(0),
    faceType_(boundaryLabel) { }
    
    /// The referenceGeometry is returned. It is taken from left element, to always ensure its existence
    // should use the same interface as in the Element!
    const ReferenceGeometry* FaceGeometry::getReferenceGeometry()const
    {
        return leftElementGeom_->getReferenceGeometry()->getCodim1ReferenceGeometry(localFaceNumberLeft_);
    }

    /*! Map a point in coordinates of the reference geometry of the face to
     *  the reference geometry of the left (L) element. */
    PointReference FaceGeometry::mapRefFaceToRefElemL(const ReferencePointT& pRefFace) const
    {
        return leftElementGeom_->getReferenceGeometry()->getCodim1MappingPtr(localFaceNumberLeft_)->transform(pRefFace);
    }
    
    /*! Map a point in coordinates of the reference geometry of the face to
     *  the reference geometry of the right (R) element. */
    PointReference FaceGeometry::mapRefFaceToRefElemR(const ReferencePointT& pRefFace) const
    {
        // In the L function we have assumed the point pRefFace to be
        // given in coordinates of the system used by the reference face
        // on the L side. That means that now to make sure that the
        // point transformed from the right side of the face onto the
        // right element is the same as the one on the left side of the
        // face, we have to use the refFace2RefFace mapping.

        PointReference pOtherSide = mapRefFaceToRefFace(pRefFace);
        return rightElementGeom_->getReferenceGeometry()->getCodim1MappingPtr(localFaceNumberRight_)->transform(pOtherSide);

    }
    
    /*! Map from reference face coordinates on the left side to those on the
     *  right side. */
    PointReference FaceGeometry::mapRefFaceToRefFace(const ReferencePointT& pIn) const
    {
        return getReferenceGeometry()->getCodim0MappingPtr(faceToFaceMapIndex_)->transform(pIn);
    }

    //! Return a Mapping
    FaceGeometry::RefFaceToRefElementMapping FaceGeometry::refFaceToRefElemMapL() const
    {
        return RefFaceToRefElementMapping(leftElementGeom_->getReferenceGeometry()->getCodim1MappingPtr(localFaceNumberLeft_));
    }

    //! Return a mapping to the right reference element.
    typename FaceGeometry::RefFaceToRefElementMapping
    FaceGeometry::refFaceToRefElemMapR() const
    {
        const MappingReferenceToReference * const m1Ptr = this->getReferenceGeometry()->getCodim0MappingPtr(faceToFaceMapIndex_);
        const MappingReferenceToReference * const m2Ptr = leftElementGeom_->getReferenceGeometry()->getCodim1MappingPtr(localFaceNumberRight_);

        return RefFaceToRefElementMapping(new ConcatenatedMapping(*m1Ptr, *m2Ptr));
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
        std::size_t DIM = pRefFace.size()+1;
        LinearAlgebra::NumericalVector result;
        if (DIM > 1)
        {
            // first Jacobian (mapping reference face -> reference element)

            Jacobian j1 = leftElementGeom_->getReferenceGeometry()->getCodim1MappingPtr(localFaceNumberLeft_) // this is the refFace2RefElemMap
                ->calcJacobian(pRefFace);

            // second Jacobian (mapping reference element -> phys. element),
            // for this we first need the points coordinates on the
            // reference element
            PointReference pRefElement = mapRefFaceToRefElemL(pRefFace);

            Jacobian j2 = leftElementGeom_->calcJacobian(pRefElement);

            Jacobian j3 = j2.multiplyJacobiansInto(j1);

            result = j3.computeWedgeStuffVector();

            double det = j2.determinant();

            double sign = OutwardNormalVectorSign(leftElementGeom_->getReferenceGeometry()->getCodim1MappingPtr(localFaceNumberLeft_));
            result *= ((det > 0)?1:-1) * sign;
        }
        else
        {   //if DIM==1
            //for one dimension the fancy wedge stuff wont work
            //but we know the left point has outward pointing vector -1
            //and the right point has outward pointing vector 1
            //so this is the same as the physical point
            PointReference pRefElement = mapRefFaceToRefElemL(pRefFace);

            //but if someone has mirrored the physical line
            //we also have to also mirror the normal vector
            //the face cant be made larger or smaller so
            //the vector should have length one

            Jacobian j = leftElementGeom_->calcJacobian(pRefElement);
            std::size_t sgn = (j[0] > 0)? 1 : -1;
            result.resize(DIM);
            result[0] = pRefElement[0] * sgn;
        }
        return result;

    }
    
    PointPhysical FaceGeometry::referenceToPhysical(const PointReference& p) const
    {
        PointReference pElement = mapRefFaceToRefElemL(p);
        return getElementGLeft()->referenceToPhysical(pElement);
    }

    //finding node numbers here is way to hard (see also the 288 lines of commented out constructor), leave that to someplace else
    void FaceGeometry::initialiseFaceToFaceMapIndex(const std::vector<std::size_t>& leftVertices, const std::vector<std::size_t>& rightVertices) {
        faceToFaceMapIndex_=getReferenceGeometry()->getCodim0MappingIndex(leftVertices,rightVertices);
    }
    
    bool FaceGeometry::isInternal() const
    {
        if(faceType_ == FaceType::INTERNAL ||
           faceType_ == FaceType::PERIODIC_BC ||
           faceType_ == FaceType::SUBDOMAIN_BOUNDARY)
        {
            ///\todo move this to the unit test
            logger.assert(rightElementGeom_!=nullptr, "There is no right element, so no internal face.");
            return true;
        }
        else
        {
            logger.assert(rightElementGeom_==nullptr, "There is a right element, so no boundary face.");
            return false;
        }
    }
    
    ///(@tito) why ask for a face when you only copy fields of the faceGeometry???
    /*void FaceGeometry::copyFromParent(const FaceGeometryT& fa)
    {
        faceToFaceMapIndex_ = fa.faceToFaceMapIndex_;
        faceToFaceMapMatrix_ = fa.faceToFaceMapMatrix_;
    }*/
    
    /*void FaceGeometry::invertFaceToFaceMapMatrix()
    {
        // invert the matrix
        double a = faceToFaceMapMatrix_(0, 0);
        double b = faceToFaceMapMatrix_(0, 1);
        double c = faceToFaceMapMatrix_(0, 2);
        double d = faceToFaceMapMatrix_(1, 0);
        double e = faceToFaceMapMatrix_(1, 1);
        double f = faceToFaceMapMatrix_(1, 2);
        double g = faceToFaceMapMatrix_(2, 0);
        double h = faceToFaceMapMatrix_(2, 1);
        double k = faceToFaceMapMatrix_(2, 2);

        double det = a * (e * k - f * h) - b * (k * d - f * g) + c * (d * h - e * g);

        faceToFaceMapMatrix_(0, 0) = (e * k - f * h) / det;
        faceToFaceMapMatrix_(1, 0) = (f * g - d * k) / det;
        faceToFaceMapMatrix_(2, 0) = (d * h - e * g) / det;
        faceToFaceMapMatrix_(0, 1) = (c * h - b * k) / det;
        faceToFaceMapMatrix_(1, 1) = (a * k - c * g) / det;
        faceToFaceMapMatrix_(2, 1) = (g * b - a * h) / det;
        faceToFaceMapMatrix_(0, 2) = (b * f - c * e) / det;
        faceToFaceMapMatrix_(1, 2) = (c * d - a * f) / det;
        faceToFaceMapMatrix_(2, 2) = (a * e - b * d) / det;
    }*/
    
    //!Print the matrix that maps the physical face to the reference face.
    /*void FaceGeometry::printRefMatrix() const
    {
        std::cout << "(" << faceToFaceMapIndex_ << ") refMatrix: " << faceToFaceMapMatrix_ << "\n";
    }*/
    
};
