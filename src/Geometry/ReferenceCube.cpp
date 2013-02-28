#include "ReferenceLine.hpp"
#include "ReferenceSquare.hpp"
#include "ReferenceCube.hpp"
#include "Mappings/MappingReferenceToReference.hpp"
#include "Mappings/MappingToRefSquareToCube.hpp"
#include "Mappings/MappingToRefCubeToCube.hpp"

namespace Geometry
{

    int ReferenceCube::localNodeIndexes_[6][4] =
    {
        { 0, 1, 2, 3 },
        { 0, 1, 4, 5 },
        { 0, 2, 4, 6 },
        { 1, 3, 5, 7 },
        { 2, 3, 6, 7 },
        { 4, 5, 6, 7 }
    };

    int ReferenceCube::localNodesOnEdge_[12][2] =
    {
        { 0 , 1 },
        { 2 , 3 },
        { 4 , 5 },
        { 6 , 7 },
        { 0 , 2 },
        { 1 , 3 },
        { 4 , 6 },
        { 5 , 7 },
        { 0 , 4 },
        { 1 , 5 },
        { 2 , 6 },
        { 3 , 7 },
    };

    ReferenceCube::ReferenceCube():
        ReferenceGeometry<ThreeD>(ThreeD+5, CUBE),
        referenceGeometryCodim1Ptr_(&ReferenceSquare::Instance()),
        referenceGeometryCodim2Ptr_(&ReferenceLine::Instance())
    {
        PointReferenceT p1,p2,p3,p4,p5,p6,p7,p8;
        
        p1[0] = -1.0; p1[1] = -1.0; p1[2] = -1.0;
        p2[0] = +1.0; p2[1] = -1.0; p2[2] = -1.0;
        p3[0] = -1.0; p3[1] = +1.0; p3[2] = -1.0;
        p4[0] = +1.0; p4[1] = +1.0; p4[2] = -1.0;
        p5[0] = -1.0; p5[1] = -1.0; p5[2] = +1.0;
        p6[0] = +1.0; p6[1] = -1.0; p6[2] = +1.0;
        p7[0] = -1.0; p7[1] = +1.0; p7[2] = +1.0;
        p8[0] = +1.0; p8[1] = +1.0; p8[2] = +1.0;
        
        points_[0] = p1;
        points_[1] = p2;
        points_[2] = p3;
        points_[3] = p4;
        points_[4] = p5;
        points_[5] = p6;
        points_[6] = p7;
        points_[7] = p8;

        mappingsSquareToCube_[0] = &MappingToRefSquareToCube0::Instance();
        mappingsSquareToCube_[1] = &MappingToRefSquareToCube1::Instance();
        mappingsSquareToCube_[2] = &MappingToRefSquareToCube2::Instance();
        mappingsSquareToCube_[3] = &MappingToRefSquareToCube3::Instance();
        mappingsSquareToCube_[4] = &MappingToRefSquareToCube4::Instance();
        mappingsSquareToCube_[5] = &MappingToRefSquareToCube5::Instance();

        mappingsCubeToCube_[0] = &MappingToRefCubeToCube0::Instance();
        mappingsCubeToCube_[1] = &MappingToRefCubeToCube1::Instance();
        mappingsCubeToCube_[2] = &MappingToRefCubeToCube2::Instance();
        mappingsCubeToCube_[3] = &MappingToRefCubeToCube3::Instance();
        mappingsCubeToCube_[4] = &MappingToRefCubeToCube4::Instance();
        mappingsCubeToCube_[5] = &MappingToRefCubeToCube5::Instance();
        mappingsCubeToCube_[6] = &MappingToRefCubeToCube6::Instance();
        mappingsCubeToCube_[7] = &MappingToRefCubeToCube7::Instance();
    }

    void ReferenceCube::getCenter(PointReferenceT& p) const
    {
        p[2] = p[1] = p[0] = 0.;
    }

    void ReferenceCube::getNode(const IndexT& i, PointReferenceT& point) const
    {
        if (i < 8)
        {
            point = points_[i];
        }
        else
        {
            throw "ERROR: Index is greater than the number of points";
        }
    }

    bool ReferenceCube::isInternalPoint(const PointReferenceT& p) const
    {
        return ((p[0] >= -1.) && (p[0] <= 1.) && (p[1] >= -1.) && (p[1] <= 1.) && (p[2] >= -1.) && (p[2] <= 1.));
    }
    
    std::ostream& operator<<(std::ostream& os, const ReferenceCube& cube)
    {
        os << cube.getName()<<" =( ";
        ReferenceCube::const_iterator it = cube.points_.begin();
        ReferenceCube::const_iterator end = cube.points_.end();

        for ( ; it != end; ++it)
        {
            os << (*it) << '\t';
        }
        os <<')'<<std::endl;

        return os;
    }

    // ================================== Codimension 0 ============================================

    int ReferenceCube::getCodim0MappingIndex(const ListOfIndexesT& list1, const ListOfIndexesT& list2) const
    {
        if (list1.size() == 8 && list2.size() == 8)
        {
            if ( (list1[0] == list2[0]) && (list1[4]==list2[4]) )
            {
                if ( (list1[1] == list2[1]) )
                    return 0;
                else
                    return 7;
            }
            else if ( (list1[0] == list2[1]) && (list1[4]==list2[5]) )
            {
                if (list1[1] == list2[0])
                    return 5;
                else
                    return 1;
            }
            else if ( (list1[0] == list2[2]) && (list1[4]==list2[6]) )
            {
                if (list1[2] == list2[0])
                    return 4;
                else
                    return 3;
            }
            else if ( (list1[0] == list2[3]) && (list1[4]==list2[7]) )
            {
                if ( (list1[3] == list2[0]) ) // (list1(3)==list2(0)) // Holds for both 2 and 6!
                    return 6;
                else
                    return 2;
            }
        }
        else
        {
            throw "ERROR: Number of node indexes was different than 8 -> not a cube.";
        }
        return -1;
    }

    const MappingReferenceToReference<3, 3>*
    ReferenceCube::getCodim0MappingPtr(const IndexT i) const
    {
        if (i < 8)
        {
            return mappingsCubeToCube_[i];
        }
        else
        {
            throw "ERROR: Cube50";
        }
    }

    // ================================== Codimension 1 ============================================

    const MappingReferenceToReference<2, 3>*
    ReferenceCube::getCodim1MappingPtr(const IndexT faceIndex) const
    {
        if (faceIndex < 6)
        {
            return mappingsSquareToCube_[faceIndex];
        }
        else
        {
            throw "ERROR: Cube100";
        }
    }

    const ReferenceGeometry<2>*
    ReferenceCube::getCodim1ReferenceGeometry(const IndexT e) const
    {
        if (e < 8)
        {
            return referenceGeometryCodim1Ptr_;
        }
        else
        {
            throw "ERROR: Cube150";
        }
    }

    void ReferenceCube::getCodim1EntityLocalIndices(const IndexT i, ListOfIndexesT& faceNodesLocal) const
    {
        if (i < 6)
        {
            faceNodesLocal.resize(4); // 2 nodes per edge
            faceNodesLocal[0] = localNodeIndexes_[i][0];
            faceNodesLocal[1] = localNodeIndexes_[i][1];
            faceNodesLocal[2] = localNodeIndexes_[i][2];
            faceNodesLocal[3] = localNodeIndexes_[i][3];
        }
        else
        {
            throw "ERROR: Cube75";
        }

    }

    // ================================== Codimension 2 ============================================

    const MappingReferenceToReference<1, 3>*
    ReferenceCube::getCodim2MappingPtr(const IndexT lineIndex) const
    {
        return 0;
    }

    const ReferenceGeometry<1>*
    ReferenceCube::getCodim2ReferenceGeometry(const IndexT e) const
    {
        if (e < 12)
        {
            return referenceGeometryCodim2Ptr_;
        }
        else
        {
            throw "ERROR: Cube150";
        }
    }

    void ReferenceCube::getCodim2EntityLocalIndices(const IndexT i, ListOfIndexesT& edgeNodesLocal) const
    {
        if (i < 12)
        {
            edgeNodesLocal.resize(2); // 2 nodes per edge
            edgeNodesLocal[0] = localNodesOnEdge_[i][0];
            edgeNodesLocal[1] = localNodesOnEdge_[i][1];
        }
        else
        {
            throw "ERROR: Cube200";
        }

    }
    
    // ================================== Quadrature rules =====================================

    /// Add a quadrature rule into the list of valid quadrature rules for this geometry.
    void ReferenceCube::addGaussQuadratureRule(QuadratureRules::GaussQuadratureRule<3>* const qr)
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
    QuadratureRules::GaussQuadratureRule<3>* const ReferenceCube::getGaussQuadratureRule(int order) const
    {
        for (std::list<QuadratureRules::GaussQuadratureRule<3>*>::const_iterator it = lstGaussQuadratureRules_.begin();
              it != lstGaussQuadratureRules_.end(); ++it)
          if ((*it)->order() >= order) return *it;

        return NULL;
    }

    // =============================== Refinement mappings =====================================
    
    void ReferenceCube::refinementTransform(int refineType, int subElementIdx, 
                  const PointReferenceT& p, PointReferenceT& pMap) const 
    {
        switch (refineType)
        {
          case 0:
            switch (subElementIdx)
            {
              case 0:
                  pMap[0]= (p[0]-1)/2.;
                  pMap[1]=p[1];
                  pMap[2]=p[2];
                  break;
              case 1:
                  pMap[0]= (p[0]+1)/2.;
                  pMap[1]=p[1];
                  pMap[2]=p[2];
                  break;
            }
            break;
            
          case 1:
            switch (subElementIdx)
            {
              case 0:
                  pMap[0]=p[0];
                  pMap[1]= (p[1]-1)/2.;
                  pMap[2]=p[2];
                  break;
              case 1:
                  pMap[0]=p[0];
                  pMap[1]= (p[1]+1)/2.;
                  pMap[2]=p[2];
                  break;
            }
            break;
            
          case 2:
            switch (subElementIdx)
            {
              case 0:
                  pMap[0]=p[0];
                  pMap[1]=p[1];
                  pMap[2]= (p[2]-1)/2.;
                  break;
              case 1:
                  pMap[0]=p[0];
                  pMap[1]=p[1];
                  pMap[2]= (p[2]+1)/2.;
                  break;
            }
            break;
            
          case 3:
            switch (subElementIdx)
            {
              case 0:
                  pMap[0]= (p[0]-1)/2.;
                  pMap[1]= (p[1]-1)/2.;
                  pMap[2]=p[2];
                  break;
              case 1:
                  pMap[0]= (p[0]+1)/2.;
                  pMap[1]= (p[1]-1)/2.;
                  pMap[2]=p[2];
                  break;
              case 2:
                  pMap[0]= (p[0]-1)/2.;
                  pMap[1]= (p[1]+1)/2.;
                  pMap[2]=p[2];
                  break;
              case 3:
                  pMap[0]= (p[0]+1)/2.;
                  pMap[1]= (p[1]+1)/2.;
                  pMap[2]=p[2];
                  break;
            }
            break;

          case 4:
            switch (subElementIdx)
            {
              case 0:
                  pMap[0]= (p[0]-1)/2.;
                  pMap[1]=p[1];
                  pMap[2]= (p[2]-1)/2.;
                  break;
              case 1:
                  pMap[0]= (p[0]+1)/2.;
                  pMap[1]=p[1];
                  pMap[2]= (p[2]-1)/2.;
                  break;
              case 2:
                  pMap[0]= (p[0]-1)/2.;
                  pMap[1]=p[1];
                  pMap[2]= (p[2]+1)/2.;
                  break;
              case 3:
                  pMap[0]= (p[0]+1)/2.;
                  pMap[1]=p[1];
                  pMap[2]= (p[2]+1)/2.;
                  break;
            }
            break;

          case 5:
            switch (subElementIdx)
            {
              case 0:
                  pMap[0]=p[0];
                  pMap[1]= (p[1]-1)/2.;
                  pMap[2]= (p[2]-1)/2.;
                  break;
              case 1:
                  pMap[0]=p[0];
                  pMap[1]= (p[1]+1)/2.;
                  pMap[2]= (p[2]-1)/2.;
                  break;
              case 2:
                  pMap[0]=p[0];
                  pMap[1]= (p[1]-1)/2.;
                  pMap[2]= (p[2]+1)/2.;
                  break;
              case 3:
                  pMap[0]=p[0];
                  pMap[1]= (p[1]+1)/2.;
                  pMap[2]= (p[2]+1)/2.;
                  break;
            }
            break;

          case 6:
            switch (subElementIdx)
            {
              case 0:
                  pMap[0]= (p[0]-1)/2.;
                  pMap[1]= (p[1]-1)/2.;
                  pMap[2]= (p[2]-1)/2.;
                  break;
              case 1:
                  pMap[0]= (p[0]+1)/2.;
                  pMap[1]= (p[1]-1)/2.;
                  pMap[2]= (p[2]-1)/2.;
                  break;
              case 2:
                  pMap[0]= (p[0]-1)/2.;
                  pMap[1]= (p[1]+1)/2.;
                  pMap[2]= (p[2]-1)/2.;
                  break;
              case 3:
                  pMap[0]= (p[0]+1)/2.;
                  pMap[1]= (p[1]+1)/2.;
                  pMap[2]= (p[2]-1)/2.;
                  break;
              case 4:
                  pMap[0]= (p[0]-1)/2.;
                  pMap[1]= (p[1]-1)/2.;
                  pMap[2]= (p[2]+1)/2.;
                  break;
              case 5:
                  pMap[0]= (p[0]+1)/2.;
                  pMap[1]= (p[1]-1)/2.;
                  pMap[2]= (p[2]+1)/2.;
                  break;
              case 6:
                  pMap[0]= (p[0]-1)/2.;
                  pMap[1]= (p[1]+1)/2.;
                  pMap[2]= (p[2]+1)/2.;
                  break;
              case 7:
                  pMap[0]= (p[0]+1)/2.;
                  pMap[1]= (p[1]+1)/2.;
                  pMap[2]= (p[2]+1)/2.;
                  break;
            }
            break;

          case 22:
            // We use this for data transfer between a coarse element, which is
            // a hexahedron, and its hexahedron subelement. The coarse element 
            // is split into two hexahedrons and one triangular-prism.
            // The parent hexahedron itself is a subsubelement of a coarser 
            // triangular-prism element which is refined in x_0 direction.
            switch (subElementIdx)
            {
              case 0:
                  pMap[0]= p[0];
                  pMap[1]= (-7. - p[0] + p[1]*(5. - p[0]))/12.;
                  pMap[2]= p[2];
                  break;

              case 1:
                  pMap[0]= p[0];
                  pMap[1]= (7. + p[0] + p[1]*(5. - p[0]))/12.;
                  pMap[2]= p[2];
                  break;

              case 2:
                  std::cout << "gotcha!  this must be an error!  subelement-2 is a triangular-prism\n";
                  break;
            }
            break;

          case 23:
            // We use this for data transfer between a coarse element, which is
            // a hexahedron, and its hexahedron subelement. The coarse element 
            // is split into two hexahedrons and one triangular-prism.
            // The parent hexahedron itself is a subsubelement of a coarser 
            // triangular-prism element which is refined in x_1 direction.
            switch (subElementIdx)
            {
              case 0:
                  pMap[0]= (-7. + 5.*p[0] + p[1]*(-1. - p[0]))/12.;
                  pMap[1]= p[1];
                  pMap[2]= p[2];
                  break;

              case 1:
                  pMap[0]= (7. + 5.*p[0] + p[1]*(1. - p[0]))/12.;
                  pMap[1]= p[1];
                  pMap[2]= p[2];
                  break;

              case 2:
                  std::cout << "gotcha!  this must be an error!  subelement-2 is a triangular-prism\n";
                  break;
            }
            break;

          case 24:
            // We use this for data transfer between a coarse element, which is
            // a hexahedron, and its hexahedron subelement. The coarse element 
            // is split into one hexahedron and two triangular-prisms.
            // The parent hexahedron itself is a subsubelement of a coarser 
            // triangular-prism element which is refined in x_0 direction.
            switch (subElementIdx)
            {
              case 0:
                  pMap[0]= p[0];
                  pMap[1]= p[1]*(2. + p[0])/3.;
                  pMap[2]= p[2];
                  break;
            }
            break;

          case 25:
            // We use this for data transfer between a coarse element, which is
            // a hexahedron, and its hexahedron subelement. The coarse element 
            // is split into one hexahedron and two triangular-prisms.
            // The parent hexahedron itself is a subsubelement of a coarser 
            // triangular-prism element which is refined in x_1 direction.
            switch (subElementIdx)
            {
              case 0:
                  pMap[0]= p[0]*(2. + p[1])/3.;
                  pMap[1]= p[1];
                  pMap[2]= p[2];
                  break;
            }
            break;

          default:
            pMap = p;
            break;
        }
    }  // end of refinementTransform

    void ReferenceCube::getRefinementMappingMatrixL(int refineType, int subElementIdx, 
                LinearAlgebra::Matrix& Q) const 
    {
        Q.resize(4,4);
        Q = 0.;
        Q(3,3) = 1.;
        switch (refineType)
        {
          case -2:
          case -1:
            Q(0,0) = 1.;
            Q(1,1) = 1.;
            Q(2,2) = 1.;
            break;
            
          case 0:
          case 1:
          case 2:
          case 3:
            Q(0,0) = 1.;
            Q(1,1) = 1.;
            Q(2,2) = 1.;
            break;
            
          case 4:
            switch (subElementIdx)
            {
              case 0:
                Q(0,0) = .5;
                Q(1,1) = .5;
                Q(2,2) = 1.;
                break;
              case 1:
                Q(0,0) = .5;   Q(0,2) = .5;
                Q(1,1) = .5;
                Q(2,2) = 1.;
                break;
              case 2:
                Q(0,0) = .5;
                Q(1,1) = .5;   Q(1,2) = .5;
                Q(2,2) = 1.;
                break;
              case 3:
                Q(0,1) = -.5;  Q(0,2) = .5;
                Q(1,0) = -.5;  Q(1,2) = .5;
                Q(2,2) = 1.;
                break;
            }
            break;

          case 5:       //????????????????????????
          case 6:       //????????????????????????
            switch (subElementIdx)
            {
              case 0:
                Q(0,0) = .5;   Q(0,2) = .5;
                Q(1,1) = .5;
                Q(2,2) = 1.;
                break;
              case 1:
                Q(0,0) = .5;   Q(0,2) = .5;
                Q(1,1) = .5;
                Q(2,2) = 1.;
                break;
            }
            break;

          case 7:       //????????????????????????
            break;
        }
    }  // end of getRefinementMappingMatrixL

    void ReferenceCube::getRefinementMappingMatrixR(int refineType, int subElementIdx, 
                LinearAlgebra::Matrix& Q) const 
    {
        Q.resize(4,4);
        Q = 0.;
        Q(3,3) = 1.;
        switch (refineType)
        {
          case -2:
          case -1:
            Q(0,0) = 1.;
            Q(1,1) = 1.;
            Q(2,2) = 1.;
            break;
            
          case 0:
          case 1:
          case 2:
          case 3:
          case 4:
          case 5:       //????????????????????????
          case 6:       //????????????????????????
          case 7:       //????????????????????????
            switch (subElementIdx)
            {
              case 0:
                Q(0,0) = 1.;
                Q(1,1) = 1.;
                Q(2,2) = 1.;
                break;
              case 1:
                Q(0,0) = 1.;
                Q(1,1) = 1.;
                Q(2,2) = 1.;
                break;
              case 2:
                Q(0,0) = 1.;
                Q(1,1) = 1.;
                Q(2,2) = 1.;
                break;
              case 3:
                Q(0,0) = 1.;
                Q(1,1) = 1.;
                Q(2,2) = 1.;
                break;
            }
            break;
        }
    }  // end of getRefinementMappingMatrixR

    void ReferenceCube::getCodim1RefinementMappingMatrixL(int refineType, DimT subElementIdx, 
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
                case 3:
                  faRefinementType = -1;
                  subFaceIndex = 100;
                  break;
                  
                default:
                  faRefinementType = 0;
                  subFaceIndex = subElementIdx;
                  break;
              }
              getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixL(faRefinementType,subFaceIndex,Q);
              break;

          case 1:
              switch (faLocalIndex)
              {
                case 1:
                case 4:
                  faRefinementType = -1;
                  subFaceIndex = 100;
                  break;
                  
                case 0:
                case 5:
                  faRefinementType = 1;
                  subFaceIndex = subElementIdx;
                  break;
                  
                case 2:
                case 3:
                  faRefinementType = 0;
                  subFaceIndex = subElementIdx;
                  break;
              }
              getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixL(faRefinementType,subFaceIndex,Q);
              break;

          case 2:
              switch (faLocalIndex)
              {
                case 0:
                case 5:
                  faRefinementType = -1;
                  subFaceIndex = 100;
                  break;
                  
                default:
                  faRefinementType = 1;
                  subFaceIndex = subElementIdx;
                  break;
              }
              getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixL(faRefinementType,subFaceIndex,Q);
              break;

          case 3:
              switch (faLocalIndex)
              {
                case 0:
                case 5:
                  faRefinementType = 2;
                  subFaceIndex = subElementIdx;
                  break;
                  
                case 1:
                case 4:
                  faRefinementType = 0;
                  subFaceIndex = subElementIdx % 2;
                  break;
                  
                case 2:
                case 3:
                  faRefinementType = 0;
                  subFaceIndex = subElementIdx / 2;
                  break;
              }
              getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixL(faRefinementType,subFaceIndex,Q);
              break;

          case 4:
              switch (faLocalIndex)
              {
                case 1:
                case 4:
                  faRefinementType = 2;
                  subFaceIndex = subElementIdx;
                  break;
                  
                case 0:
                case 5:
                  faRefinementType = 0;
                  subFaceIndex = subElementIdx % 2;
                  break;
                  
                case 2:
                case 3:
                  faRefinementType = 1;
                  subFaceIndex = subElementIdx / 2;
                  break;
              }
              getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixL(faRefinementType,subFaceIndex,Q);
              break;

          case 5:
              switch (faLocalIndex)
              {
                case 2:
                case 3:
                  faRefinementType = 2;
                  subFaceIndex = subElementIdx;
                  break;
                  
                case 1:
                case 4:
                  faRefinementType = 1;
                  subFaceIndex = subElementIdx % 2;
                  break;
                  
                case 0:
                case 5:
                  faRefinementType = 1;
                  subFaceIndex = subElementIdx / 2;
                  break;
              }
              getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixL(faRefinementType,subFaceIndex,Q);
              break;

          case 6:
              faRefinementType = 2;
              switch (faLocalIndex)
              {
                case 0:
                case 5:
                  subFaceIndex = subElementIdx % 4;
                  break;
                  
                case 1:
                case 4:
                  subFaceIndex = 2*(subElementIdx / 4) + (subElementIdx % 2);
                  break;
                  
                case 2:
                case 3:
                  subFaceIndex = subElementIdx / 2;
                  break;
              }
              getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixL(faRefinementType,subFaceIndex,Q);
              break;

          default:
              Q(0,0) = 1.;
              Q(1,1) = 1.;
              break;
        }
    }  // end of getCodim1RefinementMappingMatrixL

    void ReferenceCube::getCodim1RefinementMappingMatrixR(int refineType, DimT subElementIdx, 
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
                case 3:
                  faRefinementType = -1;
                  subFaceIndex = 100;
                  break;
                  
                default:
                  faRefinementType = 0;
                  subFaceIndex = subElementIdx;
                  break;
              }
              getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixR(faRefinementType,subFaceIndex,Q);
              break;

          case 1:
              switch (faLocalIndex)
              {
                case 1:
                case 4:
                  faRefinementType = -1;
                  subFaceIndex = 100;
                  break;
                  
                case 0:
                case 5:
                  faRefinementType = 1;
                  subFaceIndex = subElementIdx;
                  break;
                  
                case 2:
                case 3:
                  faRefinementType = 0;
                  subFaceIndex = subElementIdx;
                  break;
              }
              getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixR(faRefinementType,subFaceIndex,Q);
              break;

          case 2:
              switch (faLocalIndex)
              {
                case 0:
                case 5:
                  faRefinementType = -1;
                  subFaceIndex = 100;
                  break;
                  
                default:
                  faRefinementType = 1;
                  subFaceIndex = subElementIdx;
                  break;
              }
              getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixR(faRefinementType,subFaceIndex,Q);
              break;

          case 3:
              switch (faLocalIndex)
              {
                case 0:
                case 5:
                  faRefinementType = 2;
                  subFaceIndex = subElementIdx;
                  break;
                  
                case 1:
                case 4:
                  faRefinementType = 0;
                  subFaceIndex = subElementIdx % 2;
                  break;
                  
                case 2:
                case 3:
                  faRefinementType = 0;
                  subFaceIndex = subElementIdx / 2;
                  break;
              }
              getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixR(faRefinementType,subFaceIndex,Q);
              break;

          case 4:
              switch (faLocalIndex)
              {
                case 1:
                case 4:
                  faRefinementType = 2;
                  subFaceIndex = subElementIdx;
                  break;
                  
                case 0:
                case 5:
                  faRefinementType = 0;
                  subFaceIndex = subElementIdx % 2;
                  break;
                  
                case 2:
                case 3:
                  faRefinementType = 1;
                  subFaceIndex = subElementIdx / 2;
                  break;
              }
              getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixR(faRefinementType,subFaceIndex,Q);
              break;

          case 5:
              switch (faLocalIndex)
              {
                case 2:
                case 3:
                  faRefinementType = 2;
                  subFaceIndex = subElementIdx;
                  break;
                  
                case 1:
                case 4:
                  faRefinementType = 1;
                  subFaceIndex = subElementIdx % 2;
                  break;
                  
                case 0:
                case 5:
                  faRefinementType = 1;
                  subFaceIndex = subElementIdx / 2;
                  break;
              }
              getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixR(faRefinementType,subFaceIndex,Q);
              break;

          case 6:
              faRefinementType = 2;
              switch (faLocalIndex)
              {
                case 0:
                case 5:
                  subFaceIndex = subElementIdx % 4;
                  break;
                  
                case 1:
                case 4:
                  subFaceIndex = 2*(subElementIdx / 4) + (subElementIdx % 2);
                  break;
                  
                case 2:
                case 3:
                  subFaceIndex = subElementIdx / 2;
                  break;
              }
              getCodim1ReferenceGeometry(faLocalIndex)->getRefinementMappingMatrixR(faRefinementType,subFaceIndex,Q);
              break;
              
          default:
              Q(0,0) = 1.;
              Q(1,1) = 1.;
              break;
        }
    }  // end of getCodim1RefinementMappingMatrixR

    
};
