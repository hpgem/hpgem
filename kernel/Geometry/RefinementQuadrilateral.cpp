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
#include "Geometry/RefinementQuadrilateral.h"

#include "PointPhysical.h"
#include "PointReference.h"
#include "PhysicalGeometry.h"
#include "ReferenceGeometry.h"

namespace Geometry
{
    std::size_t RefinementQuadrilateral::getNumberOfNewNodes(std::size_t refineType) const
    {
        if (refineType == 2)
            return 5;
        else
            return 2;
    }
    
    std::vector<const PointPhysicalBase*> RefinementQuadrilateral::getAllNodes(std::size_t refineType) const
    {
        // get all element's nodes
        std::vector<const PointPhysicalBase*> nodes;
        for (std::size_t i = 0; i < referenceGeometry_->getNumberOfNodes(); ++i)
        {
            nodes.push_back(&physicalGeometry_->getLocalNodeCoordinates(i));
        }
        
        // add new nodes
        switch (refineType)
        {
            case 0: // x-refinement
                nodes.push_back(new PointPhysical<2>(0.5 * (*static_cast<const PointPhysical<2>*>(nodes[0]) + *static_cast<const PointPhysical<2>*>(nodes[1])))); // between 0-1
                nodes.push_back(new PointPhysical<2>(0.5 * (*static_cast<const PointPhysical<2>*>(nodes[2]) + *static_cast<const PointPhysical<2>*>(nodes[3])))); // between 2-3
                break;
                
            case 1: // y-refinement
                nodes.push_back(new PointPhysical<2>(0.5 * (*static_cast<const PointPhysical<2>*>(nodes[0]) + *static_cast<const PointPhysical<2>*>(nodes[2])))); // between 0-2
                nodes.push_back(new PointPhysical<2>(0.5 * (*static_cast<const PointPhysical<2>*>(nodes[1]) + *static_cast<const PointPhysical<2>*>(nodes[3])))); // between 1-3
                break;
                
            case 2: // xy-refinement
                nodes.push_back(new PointPhysical<2>(0.5 * (*static_cast<const PointPhysical<2>*>(nodes[0]) + *static_cast<const PointPhysical<2>*>(nodes[1])))); // between 0-1
                nodes.push_back(new PointPhysical<2>(0.5 * (*static_cast<const PointPhysical<2>*>(nodes[2]) + *static_cast<const PointPhysical<2>*>(nodes[3])))); // between 2-3
                nodes.push_back(new PointPhysical<2>(0.5 * (*static_cast<const PointPhysical<2>*>(nodes[0]) + *static_cast<const PointPhysical<2>*>(nodes[2])))); // between 0-2
                nodes.push_back(new PointPhysical<2>(0.5 * (*static_cast<const PointPhysical<2>*>(nodes[1]) + *static_cast<const PointPhysical<2>*>(nodes[3])))); // between 1-3
                nodes.push_back(new PointPhysical<2>(0.5 * (*static_cast<const PointPhysical<2>*>(nodes[4]) + *static_cast<const PointPhysical<2>*>(nodes[5])))); // middle
                break;
        }
        return nodes;
    }
    
    std::size_t RefinementQuadrilateral::getNumberOfSubElements(std::size_t refineType) const
    {
        if (refineType == 2)
            return 4;
        else
            return 2;
    }
    
    void RefinementQuadrilateral::subElementLocalNodeIndices(std::size_t refineType, std::size_t iSubElement, VectorOfIndicesT& LocalNodeIdx) const
    {
        logger.assert((iSubElement < getNumberOfSubElements(refineType)), "RefinementQuadrilateral: invalid sub-element index while getting its local node indices!");
        
        LocalNodeIdx.clear();
        switch (refineType)
        {
            case 0: // x-refinement
                switch (iSubElement)
                {
                    case 0:
                        // left
                    {
                        std::size_t nodes[] = {0, 4, 2, 5};
                        for (std::size_t i = 0; i < 4; ++i)
                            LocalNodeIdx.push_back(nodes[i]);
                    }
                        break;
                        
                    case 1:
                        // right
                    {
                        std::size_t nodes[] = {4, 1, 5, 3};
                        for (std::size_t i = 0; i < 4; ++i)
                            LocalNodeIdx.push_back(nodes[i]);
                    }
                        break;
                }
                break;
                
            case 1: // y-refinement
                switch (iSubElement)
                {
                    case 0:
                        // left
                    {
                        std::size_t nodes[] = {0, 1, 4, 5};
                        for (std::size_t i = 0; i < 4; ++i)
                            LocalNodeIdx.push_back(nodes[i]);
                    }
                        break;
                        
                    case 1:
                        // right
                    {
                        std::size_t nodes[] = {4, 5, 2, 3};
                        for (std::size_t i = 0; i < 4; ++i)
                            LocalNodeIdx.push_back(nodes[i]);
                    }
                        break;
                }
                break;
                
            case 2: // xy-refinement
                switch (iSubElement)
                {
                    case 0:
                        // bottom left
                    {
                        std::size_t nodes[] = {0, 4, 5, 8};
                        for (std::size_t i = 0; i < 4; ++i)
                            LocalNodeIdx.push_back(nodes[i]);
                    }
                        break;
                        
                    case 1:
                        // bottom right
                    {
                        std::size_t nodes[] = {4, 1, 8, 6};
                        for (std::size_t i = 0; i < 4; ++i)
                            LocalNodeIdx.push_back(nodes[i]);
                    }
                        break;
                        
                    case 2:
                        // top left
                    {
                        std::size_t nodes[] = {5, 8, 2, 7};
                        for (std::size_t i = 0; i < 4; ++i)
                            LocalNodeIdx.push_back(nodes[i]);
                    }
                        break;
                        
                    case 3:
                        // top right
                    {
                        std::size_t nodes[] = {8, 6, 7, 3};
                        for (std::size_t i = 0; i < 4; ++i)
                            LocalNodeIdx.push_back(nodes[i]);
                    }
                        break;
                }
                break;
        }
        
    }
    
    void RefinementQuadrilateral::adjacentSubElementsPairs(std::size_t refineType, VectorOfIndicesT& elemIdx1, VectorOfIndicesT& localFaceIdx1, VectorOfIndicesT& elemIdx2, VectorOfIndicesT& localFaceIdx2) const
    {
        elemIdx1.clear();
        elemIdx2.clear();
        localFaceIdx1.clear();
        localFaceIdx2.clear();
        
        switch (refineType)
        {
            case 0:
                // the first elements pair
                elemIdx1.push_back(1);
                localFaceIdx1.push_back(1); // elem-1, local face-1
                elemIdx2.push_back(0);
                localFaceIdx2.push_back(2); // elem-0, local face-2
                        
                break;
                
            case 1:
                // the first elements pair
                elemIdx1.push_back(1);
                localFaceIdx1.push_back(0); // elem-1, local face-0
                elemIdx2.push_back(0);
                localFaceIdx2.push_back(3); // elem-0, local face-3
                        
                break;
                
            case 2:
                // the first elements pair
                elemIdx1.push_back(1);
                localFaceIdx1.push_back(1); // elem-1, local face-1
                elemIdx2.push_back(0);
                localFaceIdx2.push_back(2); // elem-0, local face-2
                        
                // the second elements pair
                elemIdx1.push_back(2);
                localFaceIdx1.push_back(0); // elem-2, local face-2
                elemIdx2.push_back(0);
                localFaceIdx2.push_back(3); // elem-0, local face-0
                        
                // the third elements pair
                elemIdx1.push_back(3);
                localFaceIdx1.push_back(0); // elem-3, local face-3
                elemIdx2.push_back(1);
                localFaceIdx2.push_back(3); // elem-1, local face-1
                        
                // the fourth elements pair
                elemIdx1.push_back(3);
                localFaceIdx1.push_back(1); // elem-3, local face-1
                elemIdx2.push_back(2);
                localFaceIdx2.push_back(2); // elem-2, local face-2
                        
                break;
        } // end of switch
    }
    
    std::size_t RefinementQuadrilateral::getNumberOfSubElementsOnFace(std::size_t refineType, std::size_t faLocalIndex) const
    {
        switch (refineType)
        {
            case 0: // x-refinement
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
                
            case 1: // y-refinement
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
                
            case 2: // xy-refinement
                return 2;
            default:
                logger(WARN, "refineType % not implemented", refineType);
        }
        
        return 1;
    }
    
    void RefinementQuadrilateral::subElementsOnFace(std::size_t refineType, std::size_t faLocalIndex, VectorOfIndicesT& localSubElemIdx) const
    {
        localSubElemIdx.clear();
        switch (refineType)
        {
            case 0: // x-refinement
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
                
            case 1: // y-refinement
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
                
            case 2: // xy-refinement
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
    
    std::size_t RefinementQuadrilateral::getLocalSubFaceNumber(std::size_t refineType, std::size_t localFaceNr, std::size_t subElementIdx) const
    {
        return localFaceNr;
    }
} // end of namespace Geometry
