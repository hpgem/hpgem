#include <iostream>

#include "Base/TestErrorDebug.hpp"
#include "Geometry/RefinementTriangle.hpp"
#include "LinearAlgebra/Matrix.hpp"


namespace Geometry
{
    unsigned int RefinementTriangle::nrOfNewNodes(int refineType) const 
    { 
        switch (refineType)
        {
          case 0:   // r0: face-0 refinement, split into 2 subelements
          case 1:   // r1: face-1 refinement, split into 2 subelements
          case 2:   // r2: face-2 refinement, split into 2 subelements
            return 1;
            
          case 3:   // r3: all-face refinement, split into 4 subelements
            return 3;
        }
        
        return 0;
    }
    
    void RefinementTriangle::getAllNodes(int refineType, VectorOfPointPhysicalsT& nodes) const 
    {
        // get all element's nodes
        nodes.clear();
        PointPhysicalT p;
        for (unsigned int i=0; i<referenceGeometry_->getNumberOfNodes(); ++i)
        {
            physicalGeometry_->getLocalNodeCoordinates(i, p);
            nodes.push_back(p);
        }

        // add new nodes
        switch (refineType)
        {
          case 0:  // r0: face-0 refinement, split into 2 subelements
            nodes.push_back(0.5 * (nodes[0]+nodes[1]));  // between 0-1
            break;

          case 1:  // r1: face-1 refinement, split into 2 subelements
            nodes.push_back(0.5 * (nodes[0]+nodes[2]));  // between 0-2
            break;

          case 2:  // r2: face-2 refinement, split into 2 subelements
            nodes.push_back(0.5 * (nodes[1]+nodes[2]));  // between 1-2
            break;

          case 3:  // r3: all-faces refinement, split into 4 subelements
            nodes.push_back(0.5 * (nodes[0]+nodes[1]));  // between 0-1
            nodes.push_back(0.5 * (nodes[0]+nodes[2]));  // between 0-2
            nodes.push_back(0.5 * (nodes[1]+nodes[2]));  // between 1-2
            break;
        }
    }

    unsigned int RefinementTriangle::nrOfSubElements(int refineType) const 
    { 
        switch (refineType)
        {
          case 0:   // r0: face-0 refinement, split into 2 subelements
          case 1:   // r1: face-1 refinement, split into 2 subelements
          case 2:   // r2: face-2 refinement, split into 2 subelements
            return 2;

          case 3:   // r3: all-face refinement, split into 4 subelements
            return 4;
        }

        return 0;
    }

    void RefinementTriangle::subElementLocalNodeIndices(int refineType, unsigned int iSubElement, VectorOfIndicesT& LocalNodeIdx) const 
    {
        TestErrorDebug((iSubElement<nrOfSubElements(refineType)),
                        "RefinementQuadrilateral: invalid sub-element index while getting its local node indices!");

        LocalNodeIdx.clear();
        switch (refineType)
        {
            case 0:    // r0: face-0 refinement, split into 2 subelements
              switch (iSubElement)
              {
                case 0:
                  {
                    unsigned int nodes[] = { 0, 3, 2 };
                    for (unsigned int i=0; i<3; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 1:
                  {
                    unsigned int nodes[] = { 3, 1, 2 };
                    for (unsigned int i=0; i<3; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
              }
              break;
              
            case 1:    // r1: face-1 refinement, split into 2 subelements
              switch (iSubElement)
              {
                case 0:
                  // left
                  {
                    unsigned int nodes[] = { 0, 1, 3 };
                    for (unsigned int i=0; i<3; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 1:
                  // right
                  {
                    unsigned int nodes[] = { 3, 1, 2 };
                    for (unsigned int i=0; i<3; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
              }
              break;
              
            case 2:    // r2: face-2 refinement, split into 2 subelements
              switch (iSubElement)
              {
                case 0:
                  // bottom left
                  {
                    unsigned int nodes[] = { 0, 1, 3 };
                    for (unsigned int i=0; i<3; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 1:
                  // bottom right
                  {
                    unsigned int nodes[] = { 0, 3, 2 };
                    for (unsigned int i=0; i<3; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
              }
              break;
                  
            case 3:    // r3: all-face refinement, split into 4 subelements
              switch (iSubElement)
              {
                case 0:
                  // bottom left
                  {
                    unsigned int nodes[] = { 0, 3, 4 };
                    for (unsigned int i=0; i<3; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 1:
                  // bottom right
                  {
                    unsigned int nodes[] = { 3, 1, 5 };
                    for (unsigned int i=0; i<3; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 2:
                  // top left
                  {
                    unsigned int nodes[] = { 4, 5, 2 };
                    for (unsigned int i=0; i<3; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 3:
                  // top right
                  {
                    unsigned int nodes[] = { 5, 4, 3 };
                    for (unsigned int i=0; i<3; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
              }
              break;
        }

    }
    
    void RefinementTriangle::adjacentSubElementsPairs(int refineType,
                    VectorOfIndicesT& elemIdx1, VectorOfIndicesT& localFaceIdx1,
                    VectorOfIndicesT& elemIdx2, VectorOfIndicesT& localFaceIdx2) const 
    {
        elemIdx1.clear();
        elemIdx2.clear();
        localFaceIdx1.clear();
        localFaceIdx2.clear();
        
        switch (refineType) 
        {
          case 0:
              // the first elements pair
              elemIdx1.push_back(1);   localFaceIdx1.push_back(1);  // elem-1, local face-1
              elemIdx2.push_back(0);   localFaceIdx2.push_back(2);  // elem-0, local face-2
              break;
              
          case 1:
              // the first elements pair
              elemIdx1.push_back(1);   localFaceIdx1.push_back(0);  // elem-1, local face-0
              elemIdx2.push_back(0);   localFaceIdx2.push_back(2);  // elem-0, local face-2
              break;
              
          case 2:
              // the first elements pair
              elemIdx1.push_back(1);   localFaceIdx1.push_back(0);  // elem-1, local face-0
              elemIdx2.push_back(0);   localFaceIdx2.push_back(1);  // elem-0, local face-1
              break;
              
          case 3:
              // the first elements pair
              elemIdx1.push_back(3);   localFaceIdx1.push_back(2);  // elem-3, local face-2
              elemIdx2.push_back(0);   localFaceIdx2.push_back(2);  // elem-0, local face-2
              
              // the second elements pair
              elemIdx1.push_back(3);   localFaceIdx1.push_back(1);  // elem-3, local face-1
              elemIdx2.push_back(1);   localFaceIdx2.push_back(1);  // elem-1, local face-1
              
              // the third elements pair
              elemIdx1.push_back(3);   localFaceIdx1.push_back(0);  // elem-3, local face-0
              elemIdx2.push_back(2);   localFaceIdx2.push_back(0);  // elem-2, local face-0
              
              break;
        } // end of switch
    }

    unsigned int RefinementTriangle::nrOfSubElementsOnFace(int refineType, unsigned int faLocalIndex) const 
    { 
        switch (refineType)
        {
          case 0:   // r0 refinement
            switch (faLocalIndex)
              {
                case 0:
                  return 2;

                case 1:
                case 2:
                  return 1;
              }
              break;
            
          case 1:   // r1 refinement
            switch (faLocalIndex)
              {
                case 1:
                  return 2;

                case 0:
                case 2:
                  return 1;
              }
              break;
            
          case 2:   // r2 refinement
            switch (faLocalIndex)
              {
                case 2:
                  return 2;

                case 0:
                case 1:
                  return 1;
              }
              break;
            
          case 3:   // r3 refinement
              return 2;
        }
        
        return 0;
    }

    void RefinementTriangle::subElementsOnFace(int refineType, unsigned int faLocalIndex, VectorOfIndicesT& localSubElemIdx) const 
    {
        localSubElemIdx.clear();
        switch (refineType)
        {
          case 0:   // r0 refinement
            switch (faLocalIndex)
              {
                case 0:
                  localSubElemIdx.push_back(0);
                  localSubElemIdx.push_back(1);
                  break;

                case 1:
                  localSubElemIdx.push_back(0);
                  break;

                case 2:
                  localSubElemIdx.push_back(1);
                  break;
              }
              break;
            
          case 1:   // r1 refinement
            switch (faLocalIndex)
              {
                case 0:
                  localSubElemIdx.push_back(0);
                  break;

                case 1:
                  localSubElemIdx.push_back(0);
                  localSubElemIdx.push_back(1);
                  break;

                case 2:
                  localSubElemIdx.push_back(1);
                  break;
              }
              break;
            
          case 2:   // r2 refinement
            switch (faLocalIndex)
              {
                case 0:
                  localSubElemIdx.push_back(0);
                  break;

                case 1:
                  localSubElemIdx.push_back(1);
                  break;

                case 2:
                  localSubElemIdx.push_back(0);
                  localSubElemIdx.push_back(1);
                  break;
              }
              break;
            
          case 3:   // r3 refinement
            switch (faLocalIndex)
              {
                case 0:
                  localSubElemIdx.push_back(0);
                  localSubElemIdx.push_back(1);
                  break;

                case 1:
                  localSubElemIdx.push_back(0);
                  localSubElemIdx.push_back(2);
                  break;

                case 2:
                  localSubElemIdx.push_back(1);
                  localSubElemIdx.push_back(2);
                  break;
              }
              break;
            
        }
    }
    
    unsigned int RefinementTriangle::getLocalSubFaceNr(int refineType, unsigned int localFaceNr, unsigned int subElementIdx) const 
    { 
      return localFaceNr;
    }

}  // end of namespace Geometry
