//
//  ReferenceTriangularPrism.cpp
//  
//
//  Created by Shavarsh Nurijanyan on 2/5/13.
//
//

#include "ReferenceTriangularPrism.hpp"

namespace Geometry
{
    int ReferenceTriangularPrism::localNodeIndexes_[5][4] =
    {
        { 0, 2, 1 },
        { 3, 4, 5 },
        { 2, 0, 5, 3 },
        { 0, 1, 3, 4 },
        { 1, 2, 4, 5 }
    };

    int ReferenceTriangularPrism::localNodesOnEdge_[9][2] =
    {
        { 0, 1 },
        { 0, 2 },
        { 1, 2 },
        { 3, 4 },
        { 3, 5 },
        { 4, 5 },
        { 0, 3 },
        { 1, 4 },
        { 2, 5 }
    };

    ReferenceTriangularPrism::ReferenceTriangularPrism():
        ReferenceGeometry<ThreeD>(ThreeD+3, TRIANGULARPRISM),
        referenceGeometryCodim1TrianglePtr_(&ReferenceTriangle::Instance()),
        referenceGeometryCodim1SquarePtr_(&ReferenceSquare::Instance()),
        referenceGeometryCodim2Ptr_(&ReferenceLine::Instance())
    {
        PointReferenceT p1, p2, p3, p4, p5, p6;
        
        p1[0] = +0.0; p1[1] = +0.0; p1[2] = -1.0;
        p2[0] = +1.0; p2[1] = +0.0; p2[2] = -1.0;
        p3[0] = +0.0; p3[1] = +1.0; p3[2] = -1.0;
        p4[0] = +0.0; p4[1] = +0.0; p4[2] = +1.0;
        p5[0] = +1.0; p5[1] = +0.0; p5[2] = +1.0;
        p6[0] = +0.0; p6[1] = +1.0; p6[2] = +1.0;
        
        points_[0] = p1;
        points_[1] = p2;
        points_[2] = p3;
        points_[3] = p4;
        points_[4] = p5;
        points_[5] = p6;

        /// Mappings between triangular prisms are not implemented
        mappingsTriangularPrismToTriangularPrism_[0] = 0;

        mappingsFaceToTriangularPrism_[0] = &MappingToRefFaceToTriangularPrism0::Instance();
        mappingsFaceToTriangularPrism_[1] = &MappingToRefFaceToTriangularPrism1::Instance();
        mappingsFaceToTriangularPrism_[2] = &MappingToRefFaceToTriangularPrism2::Instance();
        mappingsFaceToTriangularPrism_[3] = &MappingToRefFaceToTriangularPrism3::Instance();
        mappingsFaceToTriangularPrism_[4] = &MappingToRefFaceToTriangularPrism4::Instance();
    }
    
    ReferenceTriangularPrism::ReferenceTriangularPrism(const ReferenceTriangularPrism& copy):
        ReferenceGeometry<ThreeD>(copy),
        referenceGeometryCodim1TrianglePtr_(copy.referenceGeometryCodim1TrianglePtr_),
        referenceGeometryCodim1SquarePtr_(copy.referenceGeometryCodim1SquarePtr_),
        referenceGeometryCodim2Ptr_(copy.referenceGeometryCodim2Ptr_)
    { }
    
    bool ReferenceTriangularPrism::isInternalPoint(const PointReferenceT& p) const
    {
        return ((-1. <= p[2]) && (1. >= p[2]) && (p[0] >= 0.) && (p[0] <= 1.) && (p[1] >= 0.) && (p[1] <= 1. - p[0]));
    }
    
    void ReferenceTriangularPrism::getCenter(PointReferenceT& p) const
    {
        p[0] = 1. / 3.;
        p[1] = 1. / 3.;
        p[2] = 0.;
    }
    
    void ReferenceTriangularPrism::getNode(const IndexT& i, PointReferenceT& point) const
    {
        point = points_[i];
    }
    
    std::ostream& operator<<(std::ostream& os, const ReferenceTriangularPrism& prism)
    {
        os <<prism.getName()<<" = ( ";
        ReferenceTriangularPrism::const_iterator it = prism.points_.begin();
        ReferenceTriangularPrism::const_iterator end = prism.points_.end();

        for ( ; it != end; ++it)
        {
            os << (*it) << '\t';
        }
        os <<')'<<std::endl;

        return os;
    }

    // ================================== Codimension 0 ============================================

    int ReferenceTriangularPrism::
    getCodim0MappingIndex(const ListOfIndexesT& list1, const ListOfIndexesT& list2) const
    {
        /// TODO: Implement tetrahedron to tetrahedron mappings.
        throw "ReferenceTriangularPrism::getCodim0MappingIndex: T.p to t.p mappings do not exist";
    }

    const MappingReferenceToReference<3, 3>*
    ReferenceTriangularPrism::getCodim0MappingPtr(const IndexT i) const
    {
        /// TODO: Implement tetrahedron to tetrahedron mappings.
        throw "ReferenceTetrahedron::getCodim0MappingPtr: T.p to T.p mappings do not exist";
    }

    // ================================== Codimension 1 ============================================

    void ReferenceTriangularPrism::
    getCodim1EntityLocalIndices(const IndexT faceIndex, ListOfIndexesT& faceNodesLocal) const
    {
        if (faceIndex < 2)
        {
            faceNodesLocal.resize(3); // triangles
            faceNodesLocal[0] = (IndexT) localNodeIndexes_[faceIndex][0];
            faceNodesLocal[1] = (IndexT) localNodeIndexes_[faceIndex][1];
            faceNodesLocal[2] = (IndexT) localNodeIndexes_[faceIndex][2];
        }
        else if (faceIndex < 5)
        {
            faceNodesLocal.resize(4); // squares
            faceNodesLocal[0] = (IndexT) localNodeIndexes_[faceIndex][0];
            faceNodesLocal[1] = (IndexT) localNodeIndexes_[faceIndex][1];
            faceNodesLocal[2] = (IndexT) localNodeIndexes_[faceIndex][2];
            faceNodesLocal[3] = (IndexT) localNodeIndexes_[faceIndex][3];
        }
        else
        {
            throw "ReferenceTriangularPrism::getCodim1EntityLocalIndices: Index out of range. T.p has 5 faces.";
        }
    }

    const ReferenceGeometry<2>*
    ReferenceTriangularPrism::getCodim1ReferenceGeometry(const IndexT faceIndex) const
    {
        if (faceIndex < 2)
        {
            return referenceGeometryCodim1TrianglePtr_;
        }
        else if (faceIndex < 5)
        {
            return referenceGeometryCodim1SquarePtr_;
        }
        else
        {
            throw "ReferenceTriangularPrism::getCodim1ReferenceGeometry: Index out of range. T.p has 5 faces.";
        }
    }

    const MappingReferenceToReference<2, 3>*
    ReferenceTriangularPrism::getCodim1MappingPtr(const IndexT faceIndex) const
    {
        if (faceIndex < 5)
        {
            return mappingsFaceToTriangularPrism_[faceIndex];
        }
        else
        {
            throw "ERROR: Asked for a square point index larger than 3. There are only 4 nodes in a square!";
        }
    }

    // ================================== Codimension 2 ============================================

    void ReferenceTriangularPrism::
    getCodim2EntityLocalIndices(const IndexT edgeIndex, ListOfIndexesT& edgeNodesLocal) const
    {
        if (edgeIndex < 9)
        {
            edgeNodesLocal.resize(2); // 2 nodes per edge
            edgeNodesLocal[0] = (IndexT) localNodesOnEdge_[edgeIndex][0];
            edgeNodesLocal[1] = (IndexT) localNodesOnEdge_[edgeIndex][1];
        }
        else
        {
            throw "ReferenceTriangularPrism::getCodim2EntityLocalIndices Index out of range. T.p has only 9 edges.";
        }
    }

    const ReferenceGeometry<1>*
    ReferenceTriangularPrism::getCodim2ReferenceGeometry(const IndexT edgeIndex) const
    {
        if (edgeIndex < 9)
        {
            return referenceGeometryCodim2Ptr_;
        }
        else
        {
            throw "ReferenceTriangularPrism::getCodim2ReferenceGeometry Index out of range. T.p has only 9 edges.";
        }
    }

    const MappingReferenceToReference<1, 3>*
    ReferenceTriangularPrism::getCodim2MappingPtr(const IndexT faceIndex) const
    {
        /// TODO: Implement line to t.p. mappings.
        throw "ReferenceTriangularPrism::getCodim2MappingPtr: Line to TP mappings do not exist";
    }

    // ================================== Codimension 3 ============================================

    void ReferenceTriangularPrism::
    getCodim3EntityLocalIndices(const IndexT nodeIndex, ListOfIndexesT& nodeNodesLocal) const
    {
        if (nodeIndex < 6)
        {
            nodeNodesLocal.resize(1); // 2 nodes per edge
            nodeNodesLocal[0] = nodeIndex;
        }
        else
        {
            throw "ReferenceTriangularPrism::Index out of range. TP has only 6 nodes.";
        }
    }

    // ================================== Quadrature rules =====================================

    /// Add a quadrature rule into the list of valid quadrature rules for this geometry.
    void ReferenceTriangularPrism::addGaussQuadratureRule(QuadratureRules::GaussQuadratureRule<3>* const qr) 
    {
        std::list<QuadratureRules::GaussQuadratureRule<3>*>::iterator it = lstGaussQuadratureRules_.begin();
        while (it != lstGaussQuadratureRules_.end())
        {
          if ((*it)->order() < qr->order()) ++it;
          else break;
        }
        lstGaussQuadratureRules_.insert(it,qr);
    }

    /// Get a valid quadrature for this geometry.
    QuadratureRules::GaussQuadratureRule<3>* const ReferenceTriangularPrism::getGaussQuadratureRule(int order) const 
    {
        for (std::list<QuadratureRules::GaussQuadratureRule<3>*>::const_iterator it = lstGaussQuadratureRules_.begin();
              it != lstGaussQuadratureRules_.end(); ++it)
          if ((*it)->order() >= order) return *it;

        return NULL;
    }
            
    // =============================== Refinement mappings =====================================
    
    void ReferenceTriangularPrism::refinementTransform(int refineType, int subElementIdx, 
                  const PointReferenceT& p, PointReferenceT& pMap) const 
    {
        switch (refineType)
        {
          case 0:
          case 1:
          case 2:
          case 3:
            break;

          case 4:
            switch (subElementIdx)
            {
              case 0:
                pMap[0]= (p[0])/2.;
                pMap[1]= (p[1])/2.;
                pMap[2]=p[2];
                break;
              case 1:
                pMap[0]= (p[0]+1)/2.;
                pMap[1]= (p[1])/2.;
                pMap[2]=p[2];
                break;
              case 2:
                pMap[0]= (p[0])/2.;
                pMap[1]= (p[1]+1)/2.;
                pMap[2]=p[2];
                break;
              case 3:
                pMap[0]= (1.-p[1])/2.;
                pMap[1]= (1.-p[0])/2.;
                pMap[2]=p[2];
                break;
            }
            break;
            
          case 5:       //????????????????????????
          case 6:       //????????????????????????
            switch (subElementIdx)
            {
              case 0:
                pMap[0]= (1. + p[0])/4.;
                pMap[1]= (3. - p[0] +3.*p[1] - p[0]*p[1])/8.;
                pMap[2]=p[2];
                break;
              case 1:
                pMap[0]= (p[0]+1)/2.;
                pMap[1]= (p[1])/2.;
                pMap[2]=p[2];
                break;
            }
            break;
            
          case 7:       //????????????????????????
            switch (subElementIdx)
            {
              case 0:
                pMap[0]= (3. - p[1] +3.*p[0] - p[0]*p[1])/8.;
                pMap[1]= (1. + p[1])/4.;
                pMap[2]=p[2];
                break;
              case 1:
                pMap[0]= (p[0])/2.;
                pMap[1]= (p[1]+1)/2.;
                pMap[2]=p[2];
                break;
            }
            break;
            
          case 20:
            // We use this for data transfer between a coarse element, which is
            // a hexahedron, and its subelements, which are three triangular-prisms.
            // The hexahedron itself is a subelement of a coarser triangular-prism
            // which element is refined in x_0 direction.
            switch (subElementIdx)
            {
              case 0:
                pMap[0]= -1. + p[0]*( 1. - (-1.)) + p[1]*(-1. - (-1.));
                pMap[1]= -1. + p[0]*(-1. - (-1.)) + p[1]*( 0. - (-1.));
                pMap[2]=p[2];
                break;
              case 1:
                pMap[0]= -1. + p[0]*( 1. - (-1.)) + p[1]*(-1. - (-1.));
                pMap[1]=  0. + p[0]*( 1. - ( 0.)) + p[1]*( 1. - ( 0.));
                pMap[2]=p[2];
                break;
              case 2:
                pMap[0]= 1. + p[0]*(-1. - ( 1.)) + p[1]*( 1. - ( 1.));
                pMap[1]= 1. + p[0]*( 0. - ( 1.)) + p[1]*(-1. - ( 1.));
                pMap[2]=p[2];
                break;
            }
            break;

          case 21:
            // We use this for data transfer between a coarse element, which is
            // a hexahedron, and its subelements, which are three triangular-prisms.
            // The hexahedron itself is a subelement of a coarser triangular-prism
            // element which is refined in x_1 direction.
            switch (subElementIdx)
            {
              case 0:
                pMap[0]= -1. + p[0]*( 0. - (-1.)) + p[1]*(-1. - (-1.));
                pMap[1]= -1. + p[0]*(-1. - (-1.)) + p[1]*( 1. - (-1.));
                pMap[2]=p[2];
                break;
              case 1:
                pMap[0]=  0. + p[0]*( 1. - ( 0.)) + p[1]*( 1. - ( 0.));
                pMap[1]= -1. + p[0]*(-1. - (-1.)) + p[1]*( 1. - (-1.));
                pMap[2]=p[2];
                break;
              case 2:
                pMap[0]=  1. + p[0]*(-1. - ( 1.)) + p[1]*( 0. - ( 1.));
                pMap[1]=  1. + p[0]*( 1. - ( 1.)) + p[1]*(-1. - ( 1.));
                pMap[2]=p[2];
                break;
            }
            break;

          case 22:
            // We use this for data transfer between a coarse element, which is
            // a hexahedron, and its triangular subelement. The coarse element 
            // is split into two hexahedrons and one triangular-prism.
            // The parent hexahedron itself is a subsubelement of a coarser 
            // triangular-prism element which is refined in x_0 direction.
            switch (subElementIdx)
            {
              case 2:
                pMap[0]=  1.    + p[0]*(-1. - ( 1.   )) + p[1]*( 1.    - ( 1.   ));
                pMap[1]=  1./3. + p[0]*( 0. - ( 1./3.)) + p[1]*(-1./3. - ( 1./3.));
                pMap[2]=p[2];
                break;
                
              case 0:
              case 1:
                std::cout << "gotcha!  this must be an error!  subelement-0 and -1 are hexahedrons\n";
                break;
            }
            break;

          case 23:
            // We use this for data transfer between a coarse element, which is
            // a hexahedron, and its triangular subelement. The coarse element 
            // is split into two hexahedrons and one triangular-prism.
            // The parent hexahedron itself is a subsubelement of a coarser 
            // triangular-prism element which is refined in x_1 direction.
            switch (subElementIdx)
            {
              case 2:
                pMap[0]=  1./3. + p[0]*( -1./3. - ( 1./3.)) + p[1]*( 0. - ( 1./3.));
                pMap[1]=  1.    + p[0]*(  1.    - ( 1.   )) + p[1]*(-1. - ( 1.   ));
                pMap[2]=p[2];
                break;
                
              case 0:
              case 1:
                std::cout << "gotcha!  this must be an error!  subelement-0 and -1 are hexahedrons\n";
                break;
            }
            break;

          case 24:
            // We use this for data transfer between a coarse element, which is
            // a hexahedron, and its triangular subelement. The coarse element 
            // is split into one hexahedron and two triangular-prisms.
            // The parent hexahedron itself is a subsubelement of a coarser 
            // triangular-prism element which is refined in x_0 direction.
            switch (subElementIdx)
            {
              case 0:
                pMap[0]= -1. + p[0]*( 1. - (-1.)) + p[1]*(-1.    - (-1.));
                pMap[1]= -1. + p[0]*(-1. - (-1.)) + p[1]*(-1./3. - (-1.));
                pMap[2]= p[2];
                break;

              case 1:
                pMap[0]= -1.    + p[0]*(1. - (-1.))   + p[1]*(-1. - (-1.));
                pMap[1]=  1./3. + p[0]*(1. - (1./3.)) + p[1]*( 1. - (1./3.));
                pMap[2]= p[2];
                break;

              case 2:
                std::cout << "gotcha!  this must be an error!  subelement-2 is a hexahedron\n";
                break;
            }
            break;

          case 25:
            // We use this for data transfer between a coarse element, which is
            // a hexahedron, and its triangular subelement. The coarse element 
            // is split into one hexahedron and two triangular-prisms.
            // The parent hexahedron itself is a subsubelement of a coarser 
            // triangular-prism element which is refined in x_1 direction.
            switch (subElementIdx)
            {
              case 0:
                pMap[0]= -1. + p[0]*(-1./3. - (-1.)) + p[1]*(-1. - (-1.));
                pMap[1]= -1. + p[0]*(-1.    - (-1.)) + p[1]*( 1. - (-1.));
                pMap[2]= p[2];
                break;

              case 1:
                pMap[0]=  1./3. + p[0]*( 1. - (1./3.)) + p[1]*(1. - (1./3.));
                pMap[1]= -1.    + p[0]*(-1. - (-1.))   + p[1]*(1. - (-1.));
                pMap[2]= p[2];
                break;

              case 2:
                std::cout << "gotcha!  this must be an error!  subelement-2 is a hexahedron\n";
                break;
            }
            break;

          default:
            pMap = p;
            break;
        }
    }  // end of refinementTransform

    void ReferenceTriangularPrism::getRefinementMappingMatrixL(int refineType, int subElementIdx, 
                LinearAlgebra::Matrix& Q) const 
    {
        Q.resize(4,4);
        Q = 0.;
        Q(3,3) = 1.;
        switch (refineType)
        {
          case 0:
            switch (subElementIdx)
            {
              case 0:
                Q(0,0) = .5;  Q(0,3) = -.5;
                Q(1,1) = 1.;
                Q(2,2) = 1.;
                break;
              case 1:
                Q(0,0) = .5;  Q(0,3) =  .5;
                Q(1,1) = 1.;
                Q(2,2) = 1.;
                break;
            }
            break;

          case 1:
            switch (subElementIdx)
            {
              case 0:
                Q(0,0) = 1.;
                Q(1,1) = .5;  Q(1,3) = -.5;
                Q(2,2) = 1.;
                break;
              case 1:
                Q(0,0) = 1.;
                Q(1,1) = .5;  Q(1,3) =  .5;
                Q(2,2) = 1.;
                break;
            }
            break;

          case 2:
            switch (subElementIdx)
            {
              case 0:
                Q(0,0) = 1.;
                Q(1,1) = 1.;
                Q(2,2) = .5;  Q(2,3) = -.5;
                break;
              case 1:
                Q(0,0) = 1.;
                Q(1,1) = 1.;
                Q(2,2) = .5;  Q(2,3) =  .5;
                break;
            }
            break;

          case 3:
            switch (subElementIdx)
            {
              case 0:
                Q(0,0) = .5;  Q(0,3) = -.5;
                Q(1,1) = .5;  Q(1,3) = -.5;
                Q(2,2) = 1.;
                break;
              case 1:
                Q(0,0) = .5;  Q(0,3) =  .5;
                Q(1,1) = .5;  Q(1,3) = -.5;
                Q(2,2) = 1.;
                break;
              case 2:
                Q(0,0) = .5;  Q(0,3) = -.5;
                Q(1,1) = .5;  Q(1,3) =  .5;
                Q(2,2) = 1.;
                break;
              case 3:
                Q(0,0) = .5;  Q(0,3) =  .5;
                Q(1,1) = .5;  Q(1,3) =  .5;
                Q(2,2) = 1.;
                break;
            }
            break;

          case 4:
            switch (subElementIdx)
            {
              case 0:
                Q(0,0) = .5;  Q(0,3) = -.5;
                Q(1,1) = 1.;
                Q(2,2) = .5;  Q(2,3) = -.5;
                break;
              case 1:
                Q(0,0) = .5;  Q(0,3) =  .5;
                Q(1,1) = 1.;
                Q(2,2) = .5;  Q(2,3) = -.5;
                break;
              case 2:
                Q(0,0) = .5;  Q(0,3) = -.5;
                Q(1,1) = 1.;
                Q(2,2) = .5;  Q(2,3) =  .5;
                break;
              case 3:
                Q(0,0) = .5;  Q(0,3) =  .5;
                Q(1,1) = 1.;
                Q(2,2) = .5;  Q(2,3) =  .5;
                break;
            }
            break;

          case 5:
            switch (subElementIdx)
            {
              case 0:
                Q(0,0) = 1.;
                Q(1,1) = .5;  Q(1,3) = -.5;
                Q(2,2) = .5;  Q(2,3) = -.5;
                break;
              case 1:
                Q(0,0) = 1.;
                Q(1,1) = .5;  Q(1,3) =  .5;
                Q(2,2) = .5;  Q(2,3) = -.5;
                break;
              case 2:
                Q(0,0) = 1.;
                Q(1,1) = .5;  Q(1,3) = -.5;
                Q(2,2) = .5;  Q(2,3) =  .5;
                break;
              case 3:
                Q(0,0) = 1.;
                Q(1,1) = .5;  Q(1,3) =  .5;
                Q(2,2) = .5;  Q(2,3) =  .5;
                break;
            }
            break;

          case 6:
            switch (subElementIdx)
            {
              case 0:
                Q(0,0) = .5;  Q(0,3) = -.5;
                Q(1,1) = .5;  Q(1,3) = -.5;
                Q(2,2) = .5;  Q(2,3) = -.5;
                break;
              case 1:
                Q(0,0) = .5;  Q(0,3) =  .5;
                Q(1,1) = .5;  Q(1,3) = -.5;
                Q(2,2) = .5;  Q(2,3) = -.5;
                break;
              case 2:
                Q(0,0) = .5;  Q(0,3) = -.5;
                Q(1,1) = .5;  Q(1,3) =  .5;
                Q(2,2) = .5;  Q(2,3) = -.5;
                break;
              case 3:
                Q(0,0) = .5;  Q(0,3) =  .5;
                Q(1,1) = .5;  Q(1,3) =  .5;
                Q(2,2) = .5;  Q(2,3) = -.5;
                break;
              case 4:
                Q(0,0) = .5;  Q(0,3) = -.5;
                Q(1,1) = .5;  Q(1,3) = -.5;
                Q(2,2) = .5;  Q(2,3) =  .5;
                break;
              case 5:
                Q(0,0) = .5;  Q(0,3) =  .5;
                Q(1,1) = .5;  Q(1,3) = -.5;
                Q(2,2) = .5;  Q(2,3) =  .5;
                break;
              case 6:
                Q(0,0) = .5;  Q(0,3) = -.5;
                Q(1,1) = .5;  Q(1,3) =  .5;
                Q(2,2) = .5;  Q(2,3) =  .5;
                break;
              case 7:
                Q(0,0) = .5;  Q(0,3) =  .5;
                Q(1,1) = .5;  Q(1,3) =  .5;
                Q(2,2) = .5;  Q(2,3) =  .5;
                break;
            }
            break;

          default:
            Q(0,0) = 1.;
            Q(1,1) = 1.;
            Q(2,2) = 1.;
            break;
        }
    }  // end of getRefinementMappingMatrixL

    void ReferenceTriangularPrism::getRefinementMappingMatrixR(int refineType, int subElementIdx, 
                LinearAlgebra::Matrix& Q) const 
    {
        Q.resize(4,4);
        Q = 0.;
        Q(3,3) = 1.;
        switch (refineType)
        {
          case 0:
            switch (subElementIdx)
            {
              case 0:
                Q(0,0) = 2.;  Q(0,3) =  1.;
                Q(1,1) = 1.;
                Q(2,2) = 1.;
                break;
              case 1:
                Q(0,0) = 2.;  Q(0,3) = -1.;
                Q(1,1) = 1.;
                Q(2,2) = 1.;
                break;
            }
            break;

          case 1:
            switch (subElementIdx)
            {
              case 0:
                Q(0,0) = 1.;
                Q(1,1) = 2.;  Q(1,3) =  1.;
                Q(2,2) = 1.;
                break;
              case 1:
                Q(0,0) = 1.;
                Q(1,1) = .5;  Q(1,3) = -1.;
                Q(2,2) = 1.;
                break;
            }
            break;

          case 2:
            switch (subElementIdx)
            {
              case 0:
                Q(0,0) = 1.;
                Q(1,1) = 1.;
                Q(2,2) = 2.;  Q(2,3) =  1.;
                break;
              case 1:
                Q(0,0) = 1.;
                Q(1,1) = 1.;
                Q(2,2) = 2.;  Q(2,3) = -1.;
                break;
            }
            break;

          case 3:
            switch (subElementIdx)
            {
              case 0:
                Q(0,0) = 2.;  Q(0,3) =  1.;
                Q(1,1) = 2.;  Q(1,3) =  1.;
                Q(2,2) = 1.;
                break;
              case 1:
                Q(0,0) = 2.;  Q(0,3) = -1.;
                Q(1,1) = 2.;  Q(1,3) =  1.;
                Q(2,2) = 1.;
                break;
              case 2:
                Q(0,0) = 2.;  Q(0,3) =  1.;
                Q(1,1) = 2.;  Q(1,3) = -1.;
                Q(2,2) = 1.;
                break;
              case 3:
                Q(0,0) = 2.;  Q(0,3) = -1.;
                Q(1,1) = 2.;  Q(1,3) = -1.;
                Q(2,2) = 1.;
                break;
            }
            break;

          case 4:
            switch (subElementIdx)
            {
              case 0:
                Q(0,0) = 2.;  Q(0,3) =  1.;
                Q(1,1) = 1.;
                Q(2,2) = 2.;  Q(2,3) =  1.;
                break;
              case 1:
                Q(0,0) = 2.;  Q(0,3) = -1.;
                Q(1,1) = 1.;
                Q(2,2) = 2.;  Q(2,3) =  1.;
                break;
              case 2:
                Q(0,0) = 2.;  Q(0,3) =  1.;
                Q(1,1) = 1.;
                Q(2,2) = 2.;  Q(2,3) = -1.;
                break;
              case 3:
                Q(0,0) = 2.;  Q(0,3) = -1.;
                Q(1,1) = 1.;
                Q(2,2) = 2.;  Q(2,3) = -1.;
                break;
            }
            break;

          case 5:
            switch (subElementIdx)
            {
              case 0:
                Q(0,0) = 1.;
                Q(1,1) = 2.;  Q(1,3) =  1.;
                Q(2,2) = 2.;  Q(2,3) =  1.;
                break;
              case 1:
                Q(0,0) = 1.;
                Q(1,1) = 2.;  Q(1,3) = -1.;
                Q(2,2) = 2.;  Q(2,3) =  1.;
                break;
              case 2:
                Q(0,0) = 1.;
                Q(1,1) = 2.;  Q(1,3) =  1.;
                Q(2,2) = 2.;  Q(2,3) = -1.;
                break;
              case 3:
                Q(0,0) = 1.;
                Q(1,1) = 2.;  Q(1,3) = -1.;
                Q(2,2) = 2.;  Q(2,3) = -1.;
                break;
            }
            break;

          case 6:
            switch (subElementIdx)
            {
              case 0:
                Q(0,0) = 2.;  Q(0,3) =  1.;
                Q(1,1) = 2.;  Q(1,3) =  1.;
                Q(2,2) = 2.;  Q(2,3) =  1.;
                break;
              case 1:
                Q(0,0) = 2.;  Q(0,3) = -1.;
                Q(1,1) = 2.;  Q(1,3) =  1.;
                Q(2,2) = 2.;  Q(2,3) =  1.;
                break;
              case 2:
                Q(0,0) = 2.;  Q(0,3) =  1.;
                Q(1,1) = 2.;  Q(1,3) = -1.;
                Q(2,2) = 2.;  Q(2,3) =  1.;
                break;
              case 3:
                Q(0,0) = 2.;  Q(0,3) = -1.;
                Q(1,1) = 2.;  Q(1,3) = -1.;
                Q(2,2) = 2.;  Q(2,3) =  1.;
                break;
              case 4:
                Q(0,0) = 2.;  Q(0,3) =  1.;
                Q(1,1) = 2.;  Q(1,3) =  1.;
                Q(2,2) = 2.;  Q(2,3) = -1.;
                break;
              case 5:
                Q(0,0) = 2.;  Q(0,3) = -1.;
                Q(1,1) = 2.;  Q(1,3) =  1.;
                Q(2,2) = 2.;  Q(2,3) = -1.;
                break;
              case 6:
                Q(0,0) = 2.;  Q(0,3) =  1.;
                Q(1,1) = 2.;  Q(1,3) = -1.;
                Q(2,2) = 2.;  Q(2,3) = -1.;
                break;
              case 7:
                Q(0,0) = 2.;  Q(0,3) = -1.;
                Q(1,1) = 2.;  Q(1,3) = -1.;
                Q(2,2) = 2.;  Q(2,3) = -1.;
                break;
            }
            break;

          default:
            Q(0,0) = 1.;
            Q(1,1) = 1.;
            Q(2,2) = 1.;
            break;
        }
    }  // end of getRefinementMappingMatrixR

    void ReferenceTriangularPrism::getCodim1RefinementMappingMatrixL(int refineType, DimT subElementIdx, 
                              DimT faLocalIndex, LinearAlgebra::Matrix& Q) const 
    {
        int faRefinementType(-1);
        DimT subFaceIndex(0);

        Q.resize(3,3);
        Q = 0.;
        Q(2,2) = 1.;

        switch (refineType)
        {
          case 0:
              switch (faLocalIndex)
              {
                case 2:
                  faRefinementType = 0;
                  subFaceIndex = ((subElementIdx+1) % 2);
                  break;

                default:
                  faRefinementType = -1;
                  subFaceIndex = 100;
                  break;
              }
              getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixL(faRefinementType,subFaceIndex,Q);
              break;

          case 1:
              switch (faLocalIndex)
              {
                case 3:
                  faRefinementType = 0;
                  subFaceIndex = subElementIdx % 2;
                  break;

                default:
                  faRefinementType = -1;
                  subFaceIndex = 100;
                  break;
              }
              getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixL(faRefinementType,subFaceIndex,Q);
              break;

          case 2:
              switch (faLocalIndex)
              {
                case 4:
                  faRefinementType = 0;
                  subFaceIndex = subElementIdx % 2;
                  break;
                  
                default:
                  faRefinementType = -1;
                  subFaceIndex = 100;
                  break;
              }
              getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixL(faRefinementType,subFaceIndex,Q);
              break;

          case 3:
              switch (faLocalIndex)
              {
                case 2:
                case 3:
                case 4:
                  faRefinementType = 1;
                  subFaceIndex = subElementIdx;
                  break;
                  
                default:
                  faRefinementType = -1;
                  subFaceIndex = 100;
                  break;
              }
              getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixL(faRefinementType,subFaceIndex,Q);
              break;

          case 4:
          case 5:
              switch (faLocalIndex)
              {
                case 2:
                  faRefinementType = 0;
                  subFaceIndex = ((subElementIdx/2)+1) % 2;
                  break;
                  
                case 3:
                  faRefinementType = 0;
                  subFaceIndex = subElementIdx % 2;
                  break;

                case 4:
                  faRefinementType = 0;
                  subFaceIndex = (subElementIdx+1) % 2;
                  break;
                  
                default:
                  faRefinementType = -1;
                  subFaceIndex = 100;
                  break;
              }
              getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixL(faRefinementType,subFaceIndex,Q);
              break;

          case 6:
              switch (faLocalIndex)
              {
                case 2:
                  Q(0,0) = -1.;
                  Q(1,1) = 1.;
                  break;
                  
                case 3:
                  faRefinementType = 0;
                  subFaceIndex = subElementIdx;
                  getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixL(faRefinementType,subFaceIndex,Q);
                  break;
                  
                case 4:
                  faRefinementType = 0;
                  subFaceIndex = (subElementIdx+1) % 2;
                  getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixL(faRefinementType,subFaceIndex,Q);
                  if (subElementIdx==0) 
                  {
                      Q(0,0) = -Q(0,0);
                  }
                  break;
                  
                default:
                  Q(0,0) = 1.;
                  Q(1,1) = 1.;
                  break;
              }
              break;

          case 7:
              switch (faLocalIndex)
              {
                case 2:
                  faRefinementType = 0;
                  subFaceIndex = (subElementIdx+1) % 2;
                  getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixL(faRefinementType,subFaceIndex,Q);
                  if (subElementIdx==0) 
                  {
                      Q(0,0) = -Q(0,0);
                  }
                  break;
                  
                case 4:
                  faRefinementType = 0;
                  subFaceIndex = subElementIdx;
                  getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixL(faRefinementType,subFaceIndex,Q);
                  break;
                  
                default:
                  Q(0,0) = 1.;
                  Q(1,1) = 1.;
                  break;
              }
              break;

          default:
              Q(0,0) = 1.;
              Q(1,1) = 1.;
              break;
        }
    }  // end of getCodim1RefinementMappingMatrixL

    void ReferenceTriangularPrism::getCodim1RefinementMappingMatrixR(int refineType, DimT subElementIdx, 
                              DimT faLocalIndex, LinearAlgebra::Matrix& Q) const 
    {
        int faRefinementType(-1);
        DimT subFaceIndex(0);

        Q.resize(3,3);
        Q = 0.;
        Q(2,2) = 1.;

        switch (refineType)
        {
          case 0:
              switch (faLocalIndex)
              {
                case 2:
                  faRefinementType = 0;
                  subFaceIndex = ((subElementIdx+1) % 2);
                  break;

                default:
                  faRefinementType = -1;
                  subFaceIndex = 100;
                  break;
              }
              getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixR(faRefinementType,subFaceIndex,Q);
              break;

          case 1:
              switch (faLocalIndex)
              {
                case 3:
                  faRefinementType = 0;
                  subFaceIndex = subElementIdx % 2;
                  break;

                default:
                  faRefinementType = -1;
                  subFaceIndex = 100;
                  break;
              }
              getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixR(faRefinementType,subFaceIndex,Q);
              break;

          case 2:
              switch (faLocalIndex)
              {
                case 4:
                  faRefinementType = 0;
                  subFaceIndex = subElementIdx % 2;
                  break;
                  
                default:
                  faRefinementType = -1;
                  subFaceIndex = 100;
                  break;
              }
              getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixR(faRefinementType,subFaceIndex,Q);
              break;

          case 3:
              switch (faLocalIndex)
              {
                case 2:
                case 3:
                case 4:
                  faRefinementType = 1;
                  subFaceIndex = subElementIdx;
                  break;
                  
                default:
                  faRefinementType = -1;
                  subFaceIndex = 100;
                  break;
              }
              getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixR(faRefinementType,subFaceIndex,Q);
              break;

          case 4:
          case 5:
              switch (faLocalIndex)
              {
                case 2:
                  faRefinementType = 0;
                  subFaceIndex = ((subElementIdx/2)+1) % 2;
                  break;
                  
                case 3:
                  faRefinementType = 0;
                  subFaceIndex = subElementIdx % 2;
                  break;

                case 4:
                  faRefinementType = 0;
                  subFaceIndex = (subElementIdx+1) % 2;
                  break;
                  
                default:
                  faRefinementType = -1;
                  subFaceIndex = 100;
                  break;
              }
              getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixR(faRefinementType,subFaceIndex,Q);
              break;

          case 6:
              switch (faLocalIndex)
              {
                case 2:
                  Q(0,0) = -1.;
                  Q(1,1) = 1.;
                  break;

                case 3:
                  faRefinementType = 0;
                  subFaceIndex = subElementIdx;
                  getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixR(faRefinementType,subFaceIndex,Q);
                  break;
                  
                case 4:
                  faRefinementType = 0;
                  subFaceIndex = (subElementIdx+1) % 2;
                  getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixR(faRefinementType,subFaceIndex,Q);
                  if (subElementIdx==0) 
                  {
                      Q(0,0) = -Q(0,0);
                      Q(0,1) = -Q(0,1);
                      Q(0,2) = -Q(0,2);
                  }
                  break;
                  
                default:
                  Q(0,0) = 1.;
                  Q(1,1) = 1.;
                  break;
              }
              break;

          case 7:
              switch (faLocalIndex)
              {
                case 2:
                  faRefinementType = 0;
                  subFaceIndex = (subElementIdx+1) % 2;
                  getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixR(faRefinementType,subFaceIndex,Q);
                  if (subElementIdx==0) 
                  {
                      Q(0,0) = -Q(0,0);
                      Q(0,1) = -Q(0,1);
                      Q(0,2) = -Q(0,2);
                  }
                  break;
                  
                case 4:
                  faRefinementType = 0;
                  subFaceIndex = subElementIdx;
                  getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixR(faRefinementType,subFaceIndex,Q);
                  break;
                  
                default:
                  Q(0,0) = 1.;
                  Q(1,1) = 1.;
                  break;
              }
              break;

          default:
              Q(0,0) = 1.;
              Q(1,1) = 1.;
              break;
        }
    }  // end of getCodim1RefinementMappingMatrixR

};
