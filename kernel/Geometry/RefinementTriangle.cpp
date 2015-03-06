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

#include "Logger.h"
#include "Geometry/RefinementTriangle.hpp"
#include "LinearAlgebra/Matrix.hpp"

#include "PointPhysical.hpp"
#include "PointReference.hpp"
#include "ReferenceGeometry.hpp"
#include "PhysicalGeometry.hpp"

namespace Geometry
{
    std::size_t RefinementTriangle::nrOfNewNodes(int refineType) const
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
        PointPhysicalT p(2);
        for (std::size_t i=0; i<referenceGeometry_->getNumberOfNodes(); ++i)
        {
            p = physicalGeometry_->getLocalNodeCoordinates(i);
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

    std::size_t RefinementTriangle::nrOfSubElements(int refineType) const
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

    void RefinementTriangle::subElementLocalNodeIndices(int refineType, std::size_t iSubElement, VectorOfIndicesT& LocalNodeIdx) const
    {
        logger.assert((iSubElement<nrOfSubElements(refineType)),
                        "RefinementQuadrilateral: invalid sub-element index while getting its local node indices!");

        LocalNodeIdx.clear();
        switch (refineType)
        {
            case 0:    // r0: face-0 refinement, split into 2 subelements
              switch (iSubElement)
              {
                case 0:
                  {
                    std::size_t nodes[] = { 0, 3, 2 };
                    for (std::size_t i=0; i<3; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 1:
                  {
                    std::size_t nodes[] = { 3, 1, 2 };
                    for (std::size_t i=0; i<3; ++i)
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
                    std::size_t nodes[] = { 0, 1, 3 };
                    for (std::size_t i=0; i<3; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 1:
                  // right
                  {
                    std::size_t nodes[] = { 3, 1, 2 };
                    for (std::size_t i=0; i<3; ++i)
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
                    std::size_t nodes[] = { 0, 1, 3 };
                    for (std::size_t i=0; i<3; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 1:
                  // bottom right
                  {
                    std::size_t nodes[] = { 0, 3, 2 };
                    for (std::size_t i=0; i<3; ++i)
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
                    std::size_t nodes[] = { 0, 3, 4 };
                    for (std::size_t i=0; i<3; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 1:
                  // bottom right
                  {
                    std::size_t nodes[] = { 3, 1, 5 };
                    for (std::size_t i=0; i<3; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 2:
                  // top left
                  {
                    std::size_t nodes[] = { 4, 5, 2 };
                    for (std::size_t i=0; i<3; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 3:
                  // top right
                  {
                    std::size_t nodes[] = { 5, 4, 3 };
                    for (std::size_t i=0; i<3; ++i)
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

    std::size_t RefinementTriangle::nrOfSubElementsOnFace(int refineType, std::size_t faLocalIndex) const
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

    void RefinementTriangle::subElementsOnFace(int refineType, std::size_t faLocalIndex, VectorOfIndicesT& localSubElemIdx) const
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
    
    std::size_t RefinementTriangle::getLocalSubFaceNr(int refineType, std::size_t localFaceNr, std::size_t subElementIdx) const
    { 
      return localFaceNr;
    }

}  // end of namespace Geometry
