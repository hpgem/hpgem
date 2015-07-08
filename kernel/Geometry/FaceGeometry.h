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

#ifndef ____FaceGeometry__
#define ____FaceGeometry__

#include <set>
#include <vector>
#include <memory>

#include <Logger.h>
#include "LinearAlgebra/MiddleSizeMatrix.h"
#include "ReferencePoint.h"
#include "PointReference.h"
#include "PointPhysical.h"
#include "ElementGeometry.h"
#include "Mappings/MappingReferenceToReference.h"
#include "Jacobian.h"
#include "Mappings/OutwardNormalVectorSign.h"

//--------------------------------------------------------------------------------------------------
//
//
// Example Internal Face Illustration:
//
//     f------e---d
//     | E1  /    |    Physical Elements E1 & E2
//     |    /  E2 |    with nodes a,b,..,f
//     a---b------c    How the elements are connected can be found here.
//
//       ^     ^        Mappings G1 & G2
//    G1 |     | G2     G1:(1,2,3,4)->(a,b,e,f)
//                      G2:(1,2,3,4)->(b,c,d,e)
//  4----3     4----3
//  |    |     |    |  Reference Element K
//  |    |     |    |  with nodes 1,2,3,4
//  1----2     1----2  (Basis)Functions are defined on K.
//
//       ^     ^        Mappings H1 & H2
//    H1 |     | H2     H1:(n1,n2)->(2,3), H2:(n1,n2)->(1,4)
//
//      n1-----n2      Reference Face S with nodes n1,n2
//                     Face Integration rules are defined on S.
//
// To evaluate an integral over the Face b-e, this integral is rewritten as an
// integral over the reference Face S and approximated by some Gauss
// rule. Evaluations are needed in both G1^(-1)(E1) and G2^(-1)(E2) since in DG
// functions are not defined uniquely on element Faces. These evaluations
// however need to be matched (since for example 2-3 may be attached to 1-4 or
// to 4-1), which can be done by calculating H1^(-1)(G1^(-1)(G2(H2)))! In fact
// this is a mapping of the reference Face S onto itself, and can be seen as a
// rotation or mirroring T. This Rotation (matrix) only changes when the
// connectivity changes. To avoid having to calculate this over and over, it is
// useful to store the rotation somewhere. Furthermore, since this is reference
// element related, all possible rotations only need to be calculated once for
// every possible Face.
// (One other point is that when dealing with degenerate elements, the mappings
// G1^(-1) and G2^(-1) are not defined well. However this exception should be
// handled in the code.)
//--------------------------------------------------------------------------------------------------

namespace Base
{
    
    class Face;
}

namespace LinearAlgebra
{
    class MiddleSizeMatrix;
    class MiddleSizeVector;
}

namespace Geometry
{
    template<std::size_t DIM>
    class PointPhysical;
    class PointPhysicalBase;
    template<std::size_t DIM>
    class PointReference;
    class PointReferenceBase;
    template<int codim>
    class MappingReferenceToReference;
    class ElementGeometry;
    class ReferenceGeometry;
    
    ///FaceType provides a classification of faces. The following are distinguished
    ///WALL_BC: Domain boundary (default set by internal mesh generators)
    ///OPEN_BC: Domain boundary (can be set manually to signal a different boundary type)
    ///INTERNAL: an internal face with nothing special going on
    ///SUBDOMAIN_BOUNDARY: an internal face where the left element is not on the same processor as the right element
    ///PERIODIC_BC: an internal face where the left element and the right element do not agree on the physical coordinates (usually caused by connecting the mesh across a periodic boundary) (used e.g. when splitting the face to create two extra coordinates instead of one)
    ///PERIODIC_SUBDOMAIN_BC: a combination of SUBDOMAIN_BOUNDARY and PERIODIC_BC
    enum class FaceType
    {
        OPEN_BC, WALL_BC, PERIODIC_BC, INTERNAL, SUBDOMAIN_BOUNDARY, PERIODIC_SUBDOMAIN_BC
    };
    
    //For sake of consistency, placed here.
    inline std::ostream & operator<<(std::ostream& out, FaceType ft)
    {
        switch (ft)
        {
            case FaceType::OPEN_BC:
                out << "Open boundary condition";
                break;
            case FaceType::WALL_BC:
                out << "Wall boundary condition";
                break;
            case FaceType::INTERNAL:
                out << "Internal";
                break;
            case FaceType::PERIODIC_BC:
                out << "Periodic boundary condition";
                break;
            case FaceType::PERIODIC_SUBDOMAIN_BC:
                out << "Periodic subdomain boundary";
                break;
            case FaceType::SUBDOMAIN_BOUNDARY:
                out << "Subdomain boundary";
                break;
        }
        return out;
    }
    
    /*!
     \brief Class to represent the geometry and topology of a face.
     
     Faces are the objects of codimension 1 that bound elements. They are
     part of the mesh in that they provide some of the topological
     information: For DG methods, the faces hold the information about
     neighboring element(s) (2 neighbors which are arbitrarily called left
     (L) and right (R) on an internal face, one neighboring element on a
     boundary face, always called left (L)). Additionally they contain
     information about their own geometry, so that it is possible to
     integrate on them.
     
     A boundary flag should be assigned to each face (even internal ones: in
     two-fluid problems also internal faces may need special treatment).
     
     Note the material on our discussions whether internal and boundary faces
     should be different polymorphic classes. In the current implementation
     they are NOT, mainly to avoid the complicated (multiple) derivations
     that had to be made by the user to implement their own face
     classes. Instead, to quickly decide whether a face is a boundary one,
     test whether elementPtrR() == 0.
     
     Face does not need a reference to its own reference geometry, since it
     can easily access it via the reference geometry of its left neighbour
     and its local face number in that one. */

    class FaceGeometry
    {
    public:
        using MatrixT = LinearAlgebra::MiddleSizeMatrix;
        using SetOfGlobalNodes = std::set<std::size_t>;
        using VectorOfLocalNodes = std::vector<std::size_t>;
        using LocalFaceNrType = std::size_t;
        using RefFaceToRefElementMappingPtr = std::shared_ptr<const MappingReferenceToReference<1> >;
        
        using ReferenceFaceGeometryT = ReferenceGeometry;

    public:
        ///Constructor for interior faces.
        //constructor will not initialize faceToFaceMapIndex, because it doesnt know how the elements are connected
        FaceGeometry(ElementGeometry* ptrElemL, const LocalFaceNrType& localFaceNumL, ElementGeometry* ptrElemRight, const LocalFaceNrType& localFaceNumR);

        /// Constructor for boundary faces.
        FaceGeometry(ElementGeometry* ptrElemL, const LocalFaceNrType&localFaceNumL, const FaceType& boundaryLabel);
        
        /// Copy constructor with new elements, for both internal and boundary faces.
        FaceGeometry(const FaceGeometry& other, ElementGeometry* ptrElemL, const LocalFaceNrType& localFaceNumL, ElementGeometry* ptrElemRight, const LocalFaceNrType& localFaceNumR);
        
        /// Don't use this copy constructor, but use the one with new elements instead
        FaceGeometry(const FaceGeometry &other) = delete;

        virtual ~FaceGeometry() = default;
        
        /// Return the pointer to the left element geometry.
        const ElementGeometry* getElementGLeft() const
        {
            return leftElementGeom_;
        }
        
        /// Return the pointer to the right element geometry, nullptr if inexistent for boundaries.
        const ElementGeometry* getPtrElementGRight() const
        {
            return rightElementGeom_;
        }
        
        /// Return local face number of the face in the left element.
        std::size_t localFaceNumberLeft() const
        {
            return localFaceNumberLeft_;
        }
        /// Return local face number of the face in the right element.
        std::size_t localFaceNumberRight() const
        {
            return localFaceNumberRight_;
        }
        
        FaceType getFaceType() const
        {
            return faceType_;
        }
        
        void setFaceType(const FaceType& newFace)
        {
            if (isInternal())
            {
                //isInternal checks that faceType_ is appropriate for an internal face, so we should update faceType_ first to make the assertion work
                faceType_ = newFace;
                logger.assert(isInternal(), "This face should be internal (%)", newFace);
            }
            else
            {
                faceType_ = newFace;
                logger.assert(!isInternal(), "This face should not be internal (%)", newFace);
            }
        }
        
        std::size_t getFaceToFaceMapIndex() const
        {
            return faceToFaceMapIndex_;
        }
        
        const ReferenceFaceGeometryT* getReferenceGeometry() const;

        /** \brief Map a point in coordinates of the reference geometry of the face to
         *  the reference geometry of the left (L) element. */
        template<std::size_t DIM>
        const PointReference<DIM + 1>& mapRefFaceToRefElemL(const PointReference<DIM>& pRefFace) const;

        /** \brief Map a point in coordinates of the reference geometry of the face to
         *  the reference geometry of the right (R) element. */
        template<std::size_t DIM>
        const PointReference<DIM + 1>& mapRefFaceToRefElemR(const PointReference<DIM>& pRefFace) const;

        /** \brief Map from reference face coordinates on the left side to those on the
         *  right side. */
        template<std::size_t DIM>
        const PointReference<DIM>& mapRefFaceToRefFace(const PointReference<DIM>& pIn) const;
        
        /// Get the normal vector corresponding to a given RefPoint
        template<std::size_t DIM>
        LinearAlgebra::SmallVector<DIM + 1> getNormalVector(const PointReference<DIM>& pRefFace) const;

        /// \brief Return the mapping of the reference face to the left reference element
        RefFaceToRefElementMappingPtr refFaceToRefElemMapL() const;

        /// \brief Return the mapping of the reference face to the right reference element.
        RefFaceToRefElementMappingPtr refFaceToRefElemMapR() const;

        /// \brief Returns the physical point corresponding to the given face reference point.
        template<std::size_t DIM>
        PointPhysical<DIM + 1> referenceToPhysical(const Geometry::PointReference<DIM>& pointReference) const;

        ///\brief set up the faceToFaceMapIndex based on node connectivity information instead of node location
        void initialiseFaceToFaceMapIndex(const std::vector<std::size_t>& leftNodes, const std::vector<std::size_t>& rightNodes);

        void copyFromParent(const FaceGeometry& fa);

        void invertFaceToFaceMapMatrix();

        void recalculateRefinementMatrix(MatrixT& Lmat, MatrixT& Rmat);

        void printRefMatrix() const;

        /// \brief Returns true if the face is internal and false otherwise.
        bool isInternal() const;

    protected:
        
        /// pointer to the left element geometry
        const ElementGeometry* leftElementGeom_;
        
        /// pointer to the right element geometry
        const ElementGeometry* rightElementGeom_;

        /// face index of the left element corresponding to this face
        LocalFaceNrType localFaceNumberLeft_;
        
        /// face index of the right element corresponding to this face
        LocalFaceNrType localFaceNumberRight_;

        /// Index corresponding to a face-to-face mapping. 
        std::size_t faceToFaceMapIndex_;
        
        /// Type of face (internal, boundary, ...)
        FaceType faceType_;
        
    };

    /*! Map a point in coordinates of the reference geometry of the face to
     *  the reference geometry of the left (L) element. */
    template<std::size_t DIM>
    const PointReference<DIM + 1>& FaceGeometry::mapRefFaceToRefElemL(const PointReference<DIM>& pRefFace) const
    {
        return leftElementGeom_->getReferenceGeometry()->getCodim1MappingPtr(localFaceNumberLeft_)->transform(pRefFace);
    }

    /*! Map a point in coordinates of the reference geometry of the face to
     *  the reference geometry of the right (R) element. */
    template<std::size_t DIM>
    const PointReference<DIM + 1>& FaceGeometry::mapRefFaceToRefElemR(const PointReference<DIM>& pRefFace) const
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
    template<std::size_t DIM>
    const PointReference<DIM>& FaceGeometry::mapRefFaceToRefFace(const PointReference<DIM>& pIn) const
    {
        return getReferenceGeometry()->getCodim0MappingPtr(faceToFaceMapIndex_)->transform(pIn);
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
    template<std::size_t DIM>
    LinearAlgebra::SmallVector<DIM + 1> FaceGeometry::getNormalVector(const PointReference<DIM>& pRefFace) const
    {
        LinearAlgebra::SmallVector<DIM + 1> result;
        if (DIM > 0)
        {
            // first Jacobian (mapping reference face -> reference element)

            Jacobian<DIM, DIM + 1> j1 = leftElementGeom_->getReferenceGeometry()->getCodim1MappingPtr(localFaceNumberLeft_) // this is the refFace2RefElemMap
            ->calcJacobian(pRefFace);

            // second Jacobian (mapping reference element -> phys. element),
            // for this we first need the points coordinates on the
            // reference element

            Jacobian<DIM + 1, DIM + 1> j2 = leftElementGeom_->calcJacobian(mapRefFaceToRefElemL(pRefFace));

            Jacobian<DIM, DIM + 1> j3 = j2.multiplyJacobiansInto(j1);

            result = j3.computeWedgeStuffVector();

            double det = j2.determinant();

            double sign = OutwardNormalVectorSign(leftElementGeom_->getReferenceGeometry()->getCodim1MappingPtr(localFaceNumberLeft_));
            result *= ((det > 0) ? 1 : -1) * sign;
        }
        else
        { //if DIM==0
          //for one dimension the fancy wedge stuff wont work
          //but we know the left point has outward pointing vector -1
          //and the right point has outward pointing vector 1
          //so this is the same as the physical point
            const PointReference<1>& pRefElement = mapRefFaceToRefElemL(pRefFace);

            //but if someone has mirrored the physical line
            //we also have to also mirror the normal vector
            //the face cant be made larger or smaller so
            //the vector should have length one

            Jacobian<DIM + 1, DIM + 1> j = leftElementGeom_->calcJacobian(mapRefFaceToRefElemL(pRefFace));
            int sgn = (j[0] > 0) ? 1 : -1;
            result[0] = pRefElement[0] * sgn;
        }
        return result;

    }

    template<std::size_t DIM>
    PointPhysical<DIM + 1> FaceGeometry::referenceToPhysical(const PointReference<DIM>& p) const
    {
        return getElementGLeft()->referenceToPhysical(mapRefFaceToRefElemL(p));
    }
}
#endif /* defined(____FaceGeometry__) */
