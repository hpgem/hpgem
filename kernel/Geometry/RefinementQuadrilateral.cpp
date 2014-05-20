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

#include <iostream>

#include "Base/TestErrorDebug.hpp"
#include "Geometry/RefinementQuadrilateral.hpp"

#include "PointPhysical.hpp"
#include "PointReference.hpp"
#include "PhysicalGeometry.hpp"
#include "ReferenceGeometry.hpp"

namespace Geometry
{
    unsigned int RefinementQuadrilateral::nrOfNewNodes(int refineType) const 
    { 
        if (refineType == 2) 
          return 5;
        else 
          return 2;
    }
    
    void RefinementQuadrilateral::getAllNodes(int refineType, VectorOfPointPhysicalsT& nodes) const 
    {
        // get all element's nodes
        nodes.clear();
        PointPhysicalT p(2);
        for (unsigned int i=0; i<referenceGeometry_->getNumberOfNodes(); ++i)
        {
            physicalGeometry_->getLocalNodeCoordinates(i, p);
            nodes.push_back(p);
        }

        // add new nodes
        switch (refineType)
        {
          case 0:  // x-refinement
            nodes.push_back(0.5 * (nodes[0]+nodes[1]));  // between 0-1
            nodes.push_back(0.5 * (nodes[2]+nodes[3]));  // between 2-3
            break;
            
          case 1:  // y-refinement
            nodes.push_back(0.5 * (nodes[0]+nodes[2]));  // between 0-2
            nodes.push_back(0.5 * (nodes[1]+nodes[3]));  // between 1-3
            break;
            
            
          case 2:  // xy-refinement
            nodes.push_back( 0.5 * (nodes[0]+nodes[1]));    // between 0-1
            nodes.push_back( 0.5 * (nodes[0]+nodes[2]));    // between 0-2
            nodes.push_back( 0.5 * (nodes[1]+nodes[3]));    // between 1-3
            nodes.push_back( 0.5 * (nodes[2]+nodes[3]));    // between 2-3
            nodes.push_back(0.25 * (nodes[0]+nodes[1]+nodes[2]+nodes[3]));  // center of those four
            break;
        }
    }
    
    unsigned int RefinementQuadrilateral::nrOfSubElements(int refineType) const 
    { 
        if (refineType == 2) 
          return 4;
        else 
          return 2;
    }

    void RefinementQuadrilateral::subElementLocalNodeIndices(int refineType, unsigned int iSubElement, VectorOfIndicesT& LocalNodeIdx) const 
    {
        TestErrorDebug((iSubElement<nrOfSubElements(refineType)),
                        "RefinementQuadrilateral: invalid sub-element index while getting its local node indices!");

        LocalNodeIdx.clear();
        switch (refineType)
        {
            case 0:    // x-refinement
              switch (iSubElement)
              {
                case 0:
                  // left
                  {
                    unsigned int nodes[] = { 0, 4, 2, 5 };
                    for (unsigned int i=0; i<4; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 1:
                  // right
                  {
                    unsigned int nodes[] = { 4, 1, 5, 3 };
                    for (unsigned int i=0; i<4; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
              }
              break;
              
            case 1:    // y-refinement
              switch (iSubElement)
              {
                case 0:
                  // left
                  {
                    unsigned int nodes[] = { 0, 1, 4, 5 };
                    for (unsigned int i=0; i<4; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 1:
                  // right
                  {
                    unsigned int nodes[] = { 4, 5, 2, 3 };
                    for (unsigned int i=0; i<4; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
              }
              break;
              
            case 2:    // xy-refinement
              switch (iSubElement)
              {
                case 0:
                  // bottom left
                  {
                    unsigned int nodes[] = { 0, 4, 5, 8 };
                    for (unsigned int i=0; i<4; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 1:
                  // bottom right
                  {
                    unsigned int nodes[] = { 4, 1, 8, 6 };
                    for (unsigned int i=0; i<4; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 2:
                  // top left
                  {
                    unsigned int nodes[] = { 5, 8, 2, 7 };
                    for (unsigned int i=0; i<4; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 3:
                  // top right
                  {
                    unsigned int nodes[] = { 8, 6, 7, 3 };
                    for (unsigned int i=0; i<4; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
              }
              break;
        }

    }
    
    void RefinementQuadrilateral::adjacentSubElementsPairs(int refineType,
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
              elemIdx2.push_back(0);   localFaceIdx2.push_back(3);  // elem-0, local face-3
              
              break;
              
          case 2:
              // the first elements pair
              elemIdx1.push_back(1);   localFaceIdx1.push_back(1);  // elem-1, local face-1
              elemIdx2.push_back(0);   localFaceIdx2.push_back(2);  // elem-0, local face-2
              
              // the second elements pair
              elemIdx1.push_back(2);   localFaceIdx1.push_back(0);  // elem-2, local face-2
              elemIdx2.push_back(0);   localFaceIdx2.push_back(3);  // elem-0, local face-0
              
              // the third elements pair
              elemIdx1.push_back(3);   localFaceIdx1.push_back(0);  // elem-3, local face-3
              elemIdx2.push_back(1);   localFaceIdx2.push_back(3);  // elem-1, local face-1
              
              // the fourth elements pair
              elemIdx1.push_back(3);   localFaceIdx1.push_back(1);  // elem-3, local face-1
              elemIdx2.push_back(2);   localFaceIdx2.push_back(2);  // elem-2, local face-2
              
              break;
        } // end of switch
    }

    unsigned int RefinementQuadrilateral::nrOfSubElementsOnFace(int refineType, unsigned int faLocalIndex) const 
    { 
        switch (refineType)
        {
            case 0:    // x-refinement
              switch (faLocalIndex)
              {
                case 0:
                case 3:
                  return 2;

                case 1:
                case 2:
                  return 1;
              }
              break;
              
            case 1:    // y-refinement
              switch (faLocalIndex)
              {
                case 1:
                case 2:
                  return 2;

                case 0:
                case 3:
                  return 1;
              }
              break;
              
            case 2:    // xy-refinement
              return 2;
        }
    }

    void RefinementQuadrilateral::subElementsOnFace(int refineType, unsigned int faLocalIndex, VectorOfIndicesT& localSubElemIdx) const 
    {
        localSubElemIdx.clear();
        switch (refineType)
        {
            case 0:    // x-refinement
              switch (faLocalIndex)
              {
                case 0:
                case 3:
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
              
            case 1:    // y-refinement
              switch (faLocalIndex)
              {
                case 1:
                case 2:
                  localSubElemIdx.push_back(0);
                  localSubElemIdx.push_back(1);
                  break;

                case 0:
                  localSubElemIdx.push_back(0);
                  break;

                case 3:
                  localSubElemIdx.push_back(1);
                  break;
              }
              break;
              
            case 2:    // xy-refinement
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
                  localSubElemIdx.push_back(3);
                  break;

                case 3:
                  localSubElemIdx.push_back(2);
                  localSubElemIdx.push_back(3);
                  break;
              }
              break;
        }
    }
    
    unsigned int RefinementQuadrilateral::getLocalSubFaceNr(int refineType, unsigned int localFaceNr, unsigned int subElementIdx) const 
    { 
        return localFaceNr;
    }
}  // end of namespace Geometry
