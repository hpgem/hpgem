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
#include "Geometry/RefinementLine.h"
#include "LinearAlgebra/MiddleSizeMatrix.h"

#include "PointPhysical.h"
#include "ReferenceLine.h"
#include "PhysicalGeometry.h"
#include "PointReference.h"

namespace Geometry
{
    std::size_t RefinementLine::getNumberOfNewNodes(std::size_t refineType) const
    {
        if (refineType == 0)
            return 1;
        else
            return 0;
    }
    
    std::vector<const PointPhysicalBase*> RefinementLine::getAllNodes(std::size_t refineType) const
    {
        // get all element's nodes
        std::vector<const PointPhysicalBase*> nodes;
        for (std::size_t i = 0; i < referenceGeometry_->getNumberOfNodes(); ++i)
        {
            nodes.push_back(&physicalGeometry_->getLocalNodeCoordinates(i));
        }
        
        // add new nodes
        if (refineType == 0)
        {
            // r0: split into 2 subelements
            nodes.push_back(new PointPhysical<1>(0.5 * (*static_cast<const PointPhysical<1>*>(nodes[0]) + *static_cast<const PointPhysical<1>*>(nodes[1])))); // between 0-1
        }
        return nodes;
    }
    
    std::size_t RefinementLine::getNumberOfSubElements(std::size_t refineType) const
    {
        if (refineType == 0)
            return 2;
        else
            return 0;
    }
    
    void RefinementLine::subElementLocalNodeIndices(std::size_t refineType, std::size_t iSubElement, VectorOfIndicesT& LocalNodeIdx) const
    {
        logger.assert((iSubElement < getNumberOfSubElements(refineType)), "RefinementQuadrilateral: invalid sub-element index while getting its local node indices!");
        
        LocalNodeIdx.clear();
        if (refineType == 0)
        {
            switch (iSubElement)
            {
                case 0:
                {
                    std::size_t nodes[] = {0, 2};
                    for (std::size_t i = 0; i < 2; ++i)
                        LocalNodeIdx.push_back(nodes[i]);
                }
                    break;
                    
                case 1:
                {
                    std::size_t nodes[] = {2, 1};
                    for (std::size_t i = 0; i < 2; ++i)
                        LocalNodeIdx.push_back(nodes[i]);
                }
                    break;
            }
        } // end if
        
    }
    
    void RefinementLine::adjacentSubElementsPairs(std::size_t refineType, VectorOfIndicesT& elemIdx1, VectorOfIndicesT& localFaceIdx1, VectorOfIndicesT& elemIdx2, VectorOfIndicesT& localFaceIdx2) const
    {
        elemIdx1.clear();
        elemIdx2.clear();
        localFaceIdx1.clear();
        localFaceIdx2.clear();
        
        if (refineType == 0)
        {
            // the first elements pair
            elemIdx1.push_back(1);
            localFaceIdx1.push_back(0); // elem-1, local face-0
            elemIdx2.push_back(0);
            localFaceIdx2.push_back(1); // elem-0, local face-1
        } // end if
    }
    
    std::size_t RefinementLine::getNumberOfSubElementsOnFace(std::size_t refineType, std::size_t faLocalIndex) const
    {
        if (refineType == 0)
            return 1;
        else
            return 0;
    }
    
    void RefinementLine::subElementsOnFace(std::size_t refineType, std::size_t faLocalIndex, VectorOfIndicesT& localSubElemIdx) const
    {
        localSubElemIdx.clear();
        if (refineType == 0)
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
    
    std::size_t RefinementLine::getLocalSubFaceNumber(std::size_t refineType, std::size_t localFaceNr, std::size_t subElementIdx) const
    {
        return localFaceNr;
    }

} // end of namespace Geometry
