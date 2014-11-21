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

///\BUG resolves field has incomplete type
#include "LinearAlgebra/Matrix.hpp"

//--------------------------------------------------------------------------------------------------
//
// SHOULD THIS GO INTO THE DOXYGEN DOCUMENTATION?
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


namespace Base{

	class Face;
}

namespace LinearAlgebra {
	class Matrix;
	class NumericalVector;
}

namespace Geometry
{
    
	class PointPhysical;
	class PointReference;
	class MappingReferenceToReference;
	class ElementGeometry;
	class ReferenceGeometry;

    enum  FaceType { OPEN_BC, WALL_BC, PERIODIC_BC, INTERNAL};

    
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
        typedef PointPhysical                                              PointPhysicalT;
        typedef LinearAlgebra::Matrix                                           MatrixT;
        //typedef std::vector<PointPhysicalT* >                                     ListOfPointPtrsT;
        //typedef Geometry::MappingReferenceToPhysical                       MappingReferenceToPhysicalT;
        typedef std::set<unsigned int>                                          SetOfGlobalNodes;
        typedef std::vector<unsigned int>                                       VectorOfLocalNodes;
        typedef PointPhysical                                              PhysicalPointT;
        typedef PointPhysical                                            PhysicalPointOnTheFaceT;
        typedef PointReference                                             ReferencePointT;
        typedef PointReference                                           ReferencePointOnTheFaceT;
        typedef unsigned int                                                    LocalFaceNrType;
        typedef std::auto_ptr<const MappingReferenceToReference >   RefFaceToRefElementMapping;// doing new later and passing, should handle its own deletion.
        																				//the ref->ref mappings are all singletons and should never be deleted (bug?) -FB
        typedef ReferenceGeometry                                        ReferenceFaceGeometryT;

        typedef ElementGeometry                                            ElementGeometryT;
        //typedef PhysicalGeometry                                           PhysicalGeometryT;

        ///\BUG some of the previous developpers think this should be Base::Face...
        typedef FaceGeometry                                                      FaceT;


    public:
            //FaceGeometry();
        //constructor will not initialize faceToFaceMapIndex, because it doesnt know how the elements are connected
        FaceGeometry(ElementGeometryT* ptrElemL, const LocalFaceNrType&  localFaceNumL, ElementGeometryT* ptrElemRight, const LocalFaceNrType&  localFaceNumR);

            //! Ctor for boundary faces.
        FaceGeometry(ElementGeometryT* ptrElemL, const LocalFaceNrType&localFaceNumL, const FaceType& boundaryLabel);


        virtual ~FaceGeometry(){}



            // Sets.
        /*virtual void                        setPtrElementLeft(const ElementGeometryT* value)   {leftElementGeom_ = value;}
        virtual void                        setPtrElementRight(const ElementGeometryT* value)  {rightElementGeom_ = value;}
        virtual void                        setFaceType(const FaceType& value)           {faceType_ = value;}
        virtual void                        setFaceToFaceTransformation(int value)       {faceToFaceMapIndex_ = value;}*/

            // Gets.
            /// Return the pointer to the left element.
        virtual const ElementGeometryT*           getElementGLeft()       const     {return leftElementGeom_;}
            /// Return the pointer to the right element, NULL if inexistent for boundaries.
        virtual const ElementGeometryT*           getPtrElementGRight()   const     {return rightElementGeom_;}

            /// Return local face number of the face in the left element.
        virtual unsigned int                localFaceNumberLeft()  const      {return localFaceNumberLeft_;}
            /// Return local face number of the face in the right element.
        virtual unsigned int                localFaceNumberRight()   const    {return localFaceNumberRight_;}

        virtual FaceType                    getFaceType() const          {return faceType_;}
        
        virtual void                        setFaceType(const FaceType& newFace) {
            if(faceType_==FaceType::INTERNAL){
                throw "It is currently impossible to alter the type of internal faces";
            }else if(newFace==FaceType::INTERNAL){
                throw "Boundary faces cannot also be internal faces";
            }else{
                faceType_=newFace;
            }
        }

        virtual int                         getFaceToFaceMapIndex()const {return faceToFaceMapIndex_;}

        virtual const ReferenceFaceGeometryT*       getReferenceGeometry() const;


        /*! Map a point in coordinates of the reference geometry of the face to
         *  the reference geometry of the left (L) element. */
        virtual void    mapRefFaceToRefElemL(const ReferencePointOnTheFaceT& pRefFace, ReferencePointT& pRefEl) const;

        /*! Map a point in coordinates of the reference geometry of the face to
         *  the reference geometry of the right (R) element. */
        virtual void    mapRefFaceToRefElemR(const ReferencePointOnTheFaceT& pRefFace, ReferencePointT& pRefEl) const;
        
        /*! Map from reference face coordinates on the left side to those on the
         *  right side. */
        virtual void    mapRefFaceToRefFace(const ReferencePointOnTheFaceT& pIn, ReferencePointOnTheFaceT& pOut) const;
            /// Get a normal at a given RefPoint
        virtual void    getNormalVector(const ReferencePointOnTheFaceT& pRefFace, LinearAlgebra::NumericalVector& v) const;

            //! Return a Mapping (not pointer or reference! Ok, wrapped by auto_ptr) /bug why?
        virtual         RefFaceToRefElementMapping refFaceToRefElemMapL() const;
        
            //! Return a mapping to the right reference element.
        virtual         RefFaceToRefElementMapping refFaceToRefElemMapR() const;
        

        virtual void            referenceToPhysical(const Geometry::PointReference& pointReference, PointPhysicalT& pointPhysical)const;

        
        ///\brief set up the faceToFaceMapIndex based on vertex connectivity information instead of node location
        void initialiseFaceToFaceMapIndex(const std::vector<unsigned int>& leftVertices, const std::vector<unsigned int>& rightVertices);

//-MTJ-start--------------

        void            copyFromParent(const FaceT& fa);
        
        void            invertFaceToFaceMapMatrix();

        // Q_new = R * Q_old * L
        void            recalculateRefinementMatrix(MatrixT& Lmat, MatrixT& Rmat);
        
        void            printRefMatrix() const;
        
    protected:
        MatrixT                     faceToFaceMapMatrix_;
//-MTJ-end--------------

    protected:

        ///\brief default constructor - for use with wrapper classes
        FaceGeometry():rightElementGeom_(NULL),leftElementGeom_(NULL){}
        
    protected:
        const ElementGeometryT*           rightElementGeom_;
        const ElementGeometryT*           leftElementGeom_;
        
        LocalFaceNrType             localFaceNumberLeft_;
        LocalFaceNrType             localFaceNumberRight_;
        
        LocalFaceNrType             faceToFaceMapIndex_;
        FaceType                    faceType_;
        
    };
};
#endif /* defined(____FaceGeometry__) */
