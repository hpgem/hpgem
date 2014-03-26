#include <iostream>

#include "Base/TestErrorDebug.hpp"
#include "Geometry/RefinementLine.hpp"
#include "LinearAlgebra/Matrix.hpp"


namespace Geometry
{
    unsigned int RefinementLine::nrOfNewNodes(int refineType) const 
    { 
        if (refineType==0)
            return 1;
        else
            return 0;
    }
    
    void RefinementLine::getAllNodes(int refineType, VectorOfPointPhysicalsT& nodes) const 
    {
        // get all element's nodes
        nodes.clear();
        PointPhysicalT p(1);
        for (unsigned int i=0; i<referenceGeometry_->getNumberOfNodes(); ++i)
        {
            physicalGeometry_->getLocalNodeCoordinates(i, p);
            nodes.push_back(p);
        }

        // add new nodes
        if (refineType==0)
        {
            // r0: split into 2 subelements
            nodes.push_back(0.5 * (nodes[0]+nodes[1]));  // between 0-1
        }
    }

    unsigned int RefinementLine::nrOfSubElements(int refineType) const 
    { 
        if (refineType==0)
            return 2;
        else
            return 0;
    }

    void RefinementLine::subElementLocalNodeIndices(int refineType, unsigned int iSubElement, VectorOfIndicesT& LocalNodeIdx) const 
    {
        TestErrorDebug((iSubElement<nrOfSubElements(refineType)),
                        "RefinementQuadrilateral: invalid sub-element index while getting its local node indices!");

        LocalNodeIdx.clear();
        if (refineType==0)
        {
            switch (iSubElement)
            {
              case 0:
                {
                  unsigned int nodes[] = { 0, 2 };
                  for (unsigned int i=0; i<2; ++i)
                    LocalNodeIdx.push_back(nodes[i]);
                }
                break;
                
              case 1:
                {
                  unsigned int nodes[] = { 2, 1 };
                  for (unsigned int i=0; i<2; ++i)
                    LocalNodeIdx.push_back(nodes[i]);
                }
                break;
            }
        }  // end if

    }
    
    void RefinementLine::adjacentSubElementsPairs(int refineType,
                    VectorOfIndicesT& elemIdx1, VectorOfIndicesT& localFaceIdx1,
                    VectorOfIndicesT& elemIdx2, VectorOfIndicesT& localFaceIdx2) const 
    {
        elemIdx1.clear();
        elemIdx2.clear();
        localFaceIdx1.clear();
        localFaceIdx2.clear();
        
        if (refineType==0)
        {
              // the first elements pair
              elemIdx1.push_back(1);   localFaceIdx1.push_back(0);  // elem-1, local face-0
              elemIdx2.push_back(0);   localFaceIdx2.push_back(1);  // elem-0, local face-1
        } // end if
    }

    unsigned int RefinementLine::nrOfSubElementsOnFace(int refineType, unsigned int faLocalIndex) const 
    { 
        if (refineType==0)
            return 1;
        else
            return 0;
    }

    void RefinementLine::subElementsOnFace(int refineType, unsigned int faLocalIndex, VectorOfIndicesT& localSubElemIdx) const 
    {
        localSubElemIdx.clear();
        if (refineType==0)
        {
            switch (faLocalIndex)
              {
                case 0:
                  localSubElemIdx.push_back(0);
                  break;

                case 1:
                  localSubElemIdx.push_back(1);
                  break;
              }
        } // end if
    }
    
    unsigned int RefinementLine::getLocalSubFaceNr(int refineType, unsigned int localFaceNr, unsigned int subElementIdx) const 
    { 
        return localFaceNr;
    }

}  // end of namespace Geometry
