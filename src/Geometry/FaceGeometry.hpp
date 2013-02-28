//
//  FaceGeometry.h
//  
//
//  Created by Shavarsh Nurijanyan on 2/10/13.
//
//

#ifndef ____FaceGeometry__
#define ____FaceGeometry__

#include "../LinearAlgebra/Matrix.hpp"
#include "ElementGeometry.hpp"
#include "Mappings/MappingReferenceToReference.hpp"
#include "Mappings/ConcatenatedMapping.hpp"
#include "Mappings/OutwardNormalVectorSign.hpp"
#include "Geometry/GlobalNamespaceGeometry.hpp"
#include "Geometry/PhysicalGeometry.hpp"

#include <set>
#include <vector>
#include <memory>
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

namespace Geometry
{
    
    enum  FaceType { OPEN_BC, WALL_BC, PERIODIC_BC, INTERNAL};

    template <unsigned int DIM>
    class Face;
    
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

  template <unsigned int DIM>
    class FaceGeometry
    {
    public:
        typedef PointPhysical<DIM>                                              PointPhysicalT;
        typedef LinearAlgebra::Matrix                                           MatrixT;
        typedef std::list<PointPhysicalT* >                                     ListOfPointPtrsT;
        typedef Geometry::MappingReferenceToPhysical<DIM>                       MappingReferenceToPhysicalT;
        typedef std::set<unsigned int>                                          SetOfGlobalNodes;
        typedef std::vector<unsigned int>                                       VectorOfLocalNodes;
        typedef PointPhysical<DIM>                                              PhysicalPointT;
        typedef PointPhysical<DIM-1>                                            PhysicalPointOnTheFaceT;
        typedef PointReference<DIM>                                             ReferencePointT;
        typedef PointReference<DIM-1>                                           ReferencePointOnTheFaceT;
        typedef unsigned int                                                    LocalFaceNrType;
        typedef std::auto_ptr<const MappingReferenceToReference<DIM-1, DIM> >   RefFaceToRefElementMapping;// doing new later and passing, should handle its own deletion.
        
        typedef       ReferenceGeometry<DIM-1>                                  ReferenceFaceGeometryT;
        
        typedef       ElementGeometry<DIM>                                      ElementGeometryT;
        typedef       PhysicalGeometry<DIM>                                     PhysicalGeometryT;
        

    public:
            //FaceGeometry();
        FaceGeometry(ElementGeometryT* ptrElemL, const LocalFaceNrType&  localFaceNumL, ElementGeometryT* ptrElemRight, const LocalFaceNrType&  localFaceNumR);
        
            //! Ctor for boundary faces.
        FaceGeometry(ElementGeometryT* ptrElemL, const LocalFaceNrType&localFaceNumL, const FaceType& boundaryLabel);
        
        
        virtual ~FaceGeometry(){}
        
        
        
            // Sets.
        virtual void setPtrElementLeft(ElementGeometryT* value)   {leftElementGeom_ = value;}
        virtual void setPtrElementRight(ElementGeometryT* value)  {rightElementGeom_ = value;}
        virtual void setFaceType(const FaceType& value)         {faceType_ = value;}
        virtual void setFaceToFaceTransformation(int value)     {faceToFaceMapIndex_ = value;}
        
            // Gets.
            /// Return the pointer to the left element.
        virtual ElementGeometryT*          getElementGLeft()            {return leftElementGeom_;}
            /// Return the pointer to the right element, NULL if inexistent for boundaries.
        virtual ElementGeometryT*          getPtrElementGRight()        {return rightElementGeom_;}
        
            /// Return local face number of the face in the left element.
        virtual unsigned int              localFaceNumberLeft()        {return localFaceNumberLeft_;}
            /// Return local face number of the face in the right element.
        virtual unsigned int              localFaceNumberRight()       {return localFaceNumberRight_;}
        
        virtual FaceType                  getFaceType() const          {return faceType_;}
        
        virtual int                       getFaceToFaceMapIndex()const {return faceToFaceMapIndex_;}
        
        const ReferenceFaceGeometryT*       getReferenceGeometry() const;  // ElementGeometry use this interface!

        
        /*! Map a point in coordinates of the reference geometry of the face to
         *  the reference geometry of the left (L) element. */
        virtual void mapRefFaceToRefElemL(const ReferencePointOnTheFaceT& pRefFace, ReferencePointT& pRefEl) const;
        
        /*! Map a point in coordinates of the reference geometry of the face to
         *  the reference geometry of the right (R) element. */
        virtual void mapRefFaceToRefElemR(const ReferencePointOnTheFaceT& pRefFace, ReferencePointT& pRefEl) const;
        
        /*! Map from reference face coordinates on the left side to those on the
         *  right side. */
        virtual void mapRefFaceToRefFace(const ReferencePointOnTheFaceT& pIn, ReferencePointOnTheFaceT& pOut) const;
            /// Get a normal at a given RefPoint
        virtual void getNormalVector(const ReferencePointOnTheFaceT& pRefFace, PointPhysicalT& v) const;

            //! Return a Mapping (not pointer or reference! Ok, wrapped by auto_ptr) /bug why?
        virtual RefFaceToRefElementMapping refFaceToRefElemMapL() const;
        
            //! Return a mapping to the right reference element.
        virtual RefFaceToRefElementMapping refFaceToRefElemMapR() const;
        
//-MTJ-start--------------
// #define MTJ
#ifdef MTJ
        void copyFromParent(const Face<DIM>& fa);
        
        void invertFaceToFaceMapMatrix();

        // Q_new = R * Q_old * L
        void recalculateRefinementMatrix(MatrixT& Lmat, MatrixT& Rmat);
        
        void printRefMatrix() const;
        
    protected:
        MatrixT                 faceToFaceMapMatrix_;
#endif
//-MTJ-end--------------

    protected:
        ElementGeometryT*        rightElementGeom_;
        ElementGeometryT*        leftElementGeom_;
        
        LocalFaceNrType         localFaceNumberLeft_;
        LocalFaceNrType         localFaceNumberRight_;
        
        LocalFaceNrType         faceToFaceMapIndex_;
        FaceType                faceType_;
        
    };
};
#include"FaceGeometry_Impl.hpp"
#endif /* defined(____FaceGeometry__) */