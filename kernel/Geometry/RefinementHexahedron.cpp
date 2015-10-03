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
#include "Geometry/RefinementHexahedron.h"

#include "PointPhysical.h"
#include "PointReference.h"
#include "ReferenceGeometry.h"
#include "PhysicalGeometry.h"

namespace Geometry
{
    std::size_t RefinementHexahedron::nrOfNewNodes(std::size_t refineType) const
    {
        switch (refineType)
        {
            case 0: // x-refinement
            case 1: // y-refinement
            case 2: // z-refinement
                return 4;
                
            case 3: // xy-refinement
            case 4: // xz-refinement
            case 5: // yz-refinement
                return 10;
                
            case 6: // xyz-refinement
                return 19;
        }
        
        return 0;
    }
    
    std::vector<const PointPhysicalBase*> RefinementHexahedron::getAllNodes(std::size_t refineType) const
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
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[0]) + *static_cast<const PointPhysical<3>*>(nodes[1])))); // between 0-1
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[2]) + *static_cast<const PointPhysical<3>*>(nodes[3])))); // between 2-3
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[4]) + *static_cast<const PointPhysical<3>*>(nodes[5])))); // between 4-5
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[6]) + *static_cast<const PointPhysical<3>*>(nodes[7])))); // between 6-7
                break;
                
            case 1: // y-refinement
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[0]) + *static_cast<const PointPhysical<3>*>(nodes[2])))); // between 0-2
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[1]) + *static_cast<const PointPhysical<3>*>(nodes[3])))); // between 1-3
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[4]) + *static_cast<const PointPhysical<3>*>(nodes[6])))); // between 4-6
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[5]) + *static_cast<const PointPhysical<3>*>(nodes[7])))); // between 5-7
                break;
                
            case 2: // z-refinement
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[0]) + *static_cast<const PointPhysical<3>*>(nodes[4])))); // between 0-4
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[1]) + *static_cast<const PointPhysical<3>*>(nodes[5])))); // between 1-5
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[2]) + *static_cast<const PointPhysical<3>*>(nodes[6])))); // between 2-6
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[3]) + *static_cast<const PointPhysical<3>*>(nodes[7])))); // between 3-7
                break;
                
            case 3: // xy-refinement
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[0]) + *static_cast<const PointPhysical<3>*>(nodes[1]))));
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[2]) + *static_cast<const PointPhysical<3>*>(nodes[3]))));
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[0]) + *static_cast<const PointPhysical<3>*>(nodes[2]))));
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[1]) + *static_cast<const PointPhysical<3>*>(nodes[3]))));
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[4]) + *static_cast<const PointPhysical<3>*>(nodes[5]))));
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[6]) + *static_cast<const PointPhysical<3>*>(nodes[7]))));
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[4]) + *static_cast<const PointPhysical<3>*>(nodes[6]))));
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[5]) + *static_cast<const PointPhysical<3>*>(nodes[7]))));
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[8]) + *static_cast<const PointPhysical<3>*>(nodes[9]))));
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[12]) + *static_cast<const PointPhysical<3>*>(nodes[13]))));
                break;
                
            case 4: // xz-refinement
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[0]) + *static_cast<const PointPhysical<3>*>(nodes[1]))));
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[4]) + *static_cast<const PointPhysical<3>*>(nodes[5]))));
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[0]) + *static_cast<const PointPhysical<3>*>(nodes[4]))));
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[1]) + *static_cast<const PointPhysical<3>*>(nodes[5]))));
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[2]) + *static_cast<const PointPhysical<3>*>(nodes[3]))));
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[6]) + *static_cast<const PointPhysical<3>*>(nodes[7]))));
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[2]) + *static_cast<const PointPhysical<3>*>(nodes[6]))));
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[3]) + *static_cast<const PointPhysical<3>*>(nodes[7]))));
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[8]) + *static_cast<const PointPhysical<3>*>(nodes[9]))));
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[12]) + *static_cast<const PointPhysical<3>*>(nodes[13]))));
                break;
                
            case 5: // yz-refinement
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[0]) + *static_cast<const PointPhysical<3>*>(nodes[2]))));
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[4]) + *static_cast<const PointPhysical<3>*>(nodes[6]))));
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[0]) + *static_cast<const PointPhysical<3>*>(nodes[4]))));
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[2]) + *static_cast<const PointPhysical<3>*>(nodes[6]))));
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[1]) + *static_cast<const PointPhysical<3>*>(nodes[3]))));
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[5]) + *static_cast<const PointPhysical<3>*>(nodes[7]))));
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[1]) + *static_cast<const PointPhysical<3>*>(nodes[5]))));
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[3]) + *static_cast<const PointPhysical<3>*>(nodes[7]))));
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[8]) + *static_cast<const PointPhysical<3>*>(nodes[9]))));
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[12]) + *static_cast<const PointPhysical<3>*>(nodes[13]))));
                break;
                
            case 6: // xyz-refinement
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[0]) + *static_cast<const PointPhysical<3>*>(nodes[1]))));
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[2]) + *static_cast<const PointPhysical<3>*>(nodes[3]))));
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[4]) + *static_cast<const PointPhysical<3>*>(nodes[5]))));
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[6]) + *static_cast<const PointPhysical<3>*>(nodes[7]))));
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[0]) + *static_cast<const PointPhysical<3>*>(nodes[2]))));
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[4]) + *static_cast<const PointPhysical<3>*>(nodes[6]))));
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[1]) + *static_cast<const PointPhysical<3>*>(nodes[3]))));
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[5]) + *static_cast<const PointPhysical<3>*>(nodes[7]))));
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[0]) + *static_cast<const PointPhysical<3>*>(nodes[4]))));
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[1]) + *static_cast<const PointPhysical<3>*>(nodes[5]))));
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[2]) + *static_cast<const PointPhysical<3>*>(nodes[6]))));
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[3]) + *static_cast<const PointPhysical<3>*>(nodes[7]))));

                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[8]) + *static_cast<const PointPhysical<3>*>(nodes[9]))));
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[10]) + *static_cast<const PointPhysical<3>*>(nodes[11]))));
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[12]) + *static_cast<const PointPhysical<3>*>(nodes[13]))));
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[14]) + *static_cast<const PointPhysical<3>*>(nodes[15]))));
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[16]) + *static_cast<const PointPhysical<3>*>(nodes[17]))));
                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[18]) + *static_cast<const PointPhysical<3>*>(nodes[19]))));

                nodes.push_back(new PointPhysical<3>(0.5 * (*static_cast<const PointPhysical<3>*>(nodes[20]) + *static_cast<const PointPhysical<3>*>(nodes[21]))));
                break;
                
            default:
                break;
        }
        return nodes;
    }
    
    std::size_t RefinementHexahedron::nrOfSubElements(std::size_t refineType) const
    {
        std::size_t nrSubElements(0);
        
        switch (refineType)
        {
            case 0: // x-refinement
            case 1: // y-refinement
            case 2: // z-refinement
                nrSubElements = 2;
                break;
                
            case 3: // xy-refinement
            case 4: // xz-refinement
            case 5: // yz-refinement
                nrSubElements = 4;
                break;
                
            case 6: // xyz-refinement
                nrSubElements = 8;
                break;
                
            default: // not refined
                nrSubElements = 1;
                break;
                
        }
        
        return nrSubElements;
    }
    
    void RefinementHexahedron::subElementLocalNodeIndices(std::size_t refineType, std::size_t iSubElement, VectorOfIndicesT& LocalNodeIdx) const
    {
        logger.assert((iSubElement < nrOfSubElements(refineType)), "RefinementHexahedron: invalid sub-element index while getting its local node indices!");
        
        LocalNodeIdx.clear();
        switch (refineType)
        {
            case 0: // x-refinement
                switch (iSubElement)
                {
                    case 0:
                        // left
                    {
                        std::size_t nodes[] = {0, 8, 2, 9, 4, 10, 6, 11};
                        for (std::size_t i = 0; i < 8; ++i)
                            LocalNodeIdx.push_back(nodes[i]);
                    }
                        break;
                        
                    case 1:
                        // right
                    {
                        std::size_t nodes[] = {8, 1, 9, 3, 10, 5, 11, 7};
                        for (std::size_t i = 0; i < 8; ++i)
                            LocalNodeIdx.push_back(nodes[i]);
                    }
                        break;
                }
                break;
                
            case 1: // y-refinement
                switch (iSubElement)
                {
                    case 0:
                        // front
                    {
                        std::size_t nodes[] = {0, 1, 8, 9, 4, 5, 10, 11};
                        for (std::size_t i = 0; i < 8; ++i)
                            LocalNodeIdx.push_back(nodes[i]);
                    }
                        break;
                        
                    case 1:
                        // rear
                    {
                        std::size_t nodes[] = {8, 9, 2, 3, 10, 11, 6, 7};
                        for (std::size_t i = 0; i < 8; ++i)
                            LocalNodeIdx.push_back(nodes[i]);
                    }
                        break;
                }
                break;
                
            case 2: // z-refinement
                switch (iSubElement)
                {
                    case 0:
                        // bottom
                    {
                        std::size_t nodes[] = {0, 1, 2, 3, 8, 9, 10, 11};
                        for (std::size_t i = 0; i < 8; ++i)
                            LocalNodeIdx.push_back(nodes[i]);
                    }
                        break;
                        
                    case 1:
                        // top
                    {
                        std::size_t nodes[] = {8, 9, 10, 11, 4, 5, 6, 7};
                        for (std::size_t i = 0; i < 8; ++i)
                            LocalNodeIdx.push_back(nodes[i]);
                    }
                        break;
                }
                break;
                
            case 3: // xy-refinement
                switch (iSubElement)
                {
                    case 0:
                        // left-front
                    {
                        std::size_t nodes[] = {0, 8, 9, 12, 4, 13, 14, 17};
                        for (std::size_t i = 0; i < 8; ++i)
                            LocalNodeIdx.push_back(nodes[i]);
                    }
                        break;
                        
                    case 1:
                        // right-front
                    {
                        std::size_t nodes[] = {8, 1, 12, 10, 13, 5, 17, 15};
                        for (std::size_t i = 0; i < 8; ++i)
                            LocalNodeIdx.push_back(nodes[i]);
                    }
                        break;
                        
                    case 2:
                        // left-rear
                    {
                        std::size_t nodes[] = {9, 12, 2, 11, 14, 17, 6, 16};
                        for (std::size_t i = 0; i < 8; ++i)
                            LocalNodeIdx.push_back(nodes[i]);
                    }
                        break;
                        
                    case 3:
                        // right-rear
                    {
                        std::size_t nodes[] = {12, 10, 11, 3, 17, 15, 16, 7};
                        for (std::size_t i = 0; i < 8; ++i)
                            LocalNodeIdx.push_back(nodes[i]);
                    }
                        break;
                }
                break;
                
            case 4: // xz-refinement
                switch (iSubElement)
                {
                    case 0:
                        // bottom-left
                    {
                        std::size_t nodes[] = {0, 8, 2, 9, 10, 14, 12, 15};
                        for (std::size_t i = 0; i < 8; ++i)
                            LocalNodeIdx.push_back(nodes[i]);
                    }
                        break;
                        
                    case 1:
                        // bottom-right
                    {
                        std::size_t nodes[] = {8, 1, 9, 3, 14, 11, 15, 13};
                        for (std::size_t i = 0; i < 8; ++i)
                            LocalNodeIdx.push_back(nodes[i]);
                    }
                        break;
                        
                    case 2:
                        // top-left
                    {
                        std::size_t nodes[] = {10, 14, 12, 15, 4, 16, 6, 17};
                        for (std::size_t i = 0; i < 8; ++i)
                            LocalNodeIdx.push_back(nodes[i]);
                    }
                        break;
                        
                    case 3:
                        // top-right
                    {
                        std::size_t nodes[] = {14, 11, 15, 13, 16, 5, 17, 7};
                        for (std::size_t i = 0; i < 8; ++i)
                            LocalNodeIdx.push_back(nodes[i]);
                    }
                        break;
                }
                break;
                
            case 5: // yz-refinement
                switch (iSubElement)
                {
                    case 0:
                        // bottom-front
                    {
                        std::size_t nodes[] = {0, 1, 8, 9, 10, 11, 14, 15};
                        for (std::size_t i = 0; i < 8; ++i)
                            LocalNodeIdx.push_back(nodes[i]);
                    }
                        break;
                        
                    case 1:
                        // bottom-rear
                    {
                        std::size_t nodes[] = {8, 9, 2, 3, 14, 15, 12, 13};
                        for (std::size_t i = 0; i < 8; ++i)
                            LocalNodeIdx.push_back(nodes[i]);
                    }
                        break;
                        
                    case 2:
                        // top-front
                    {
                        std::size_t nodes[] = {10, 11, 14, 15, 4, 5, 16, 17};
                        for (std::size_t i = 0; i < 8; ++i)
                            LocalNodeIdx.push_back(nodes[i]);
                    }
                        break;
                        
                    case 3:
                        // top-rear
                    {
                        std::size_t nodes[] = {14, 15, 12, 13, 16, 17, 6, 7};
                        for (std::size_t i = 0; i < 8; ++i)
                            LocalNodeIdx.push_back(nodes[i]);
                    }
                        break;
                }
                break;
                
            case 6: // xyz-refinement
                switch (iSubElement)
                {
                    case 0:
                        // bottom-left-front
                    {
                        std::size_t nodes[] = {0, 8, 9, 12, 13, 17, 18, 21};
                        for (std::size_t i = 0; i < 8; ++i)
                            LocalNodeIdx.push_back(nodes[i]);
                    }
                        break;
                        
                    case 1:
                        // bottom-right-front
                    {
                        std::size_t nodes[] = {8, 1, 12, 10, 17, 14, 21, 19};
                        for (std::size_t i = 0; i < 8; ++i)
                            LocalNodeIdx.push_back(nodes[i]);
                    }
                        break;
                        
                    case 2:
                        // bottom-left-rear
                    {
                        std::size_t nodes[] = {9, 12, 2, 11, 18, 21, 15, 20};
                        for (std::size_t i = 0; i < 8; ++i)
                            LocalNodeIdx.push_back(nodes[i]);
                    }
                        break;
                        
                    case 3:
                        // bottom-right-rear
                    {
                        std::size_t nodes[] = {12, 10, 11, 3, 21, 19, 20, 16};
                        for (std::size_t i = 0; i < 8; ++i)
                            LocalNodeIdx.push_back(nodes[i]);
                    }
                        break;
                        
                    case 4:
                        // top-left-front
                    {
                        std::size_t nodes[] = {13, 17, 18, 21, 4, 22, 23, 26};
                        for (std::size_t i = 0; i < 8; ++i)
                            LocalNodeIdx.push_back(nodes[i]);
                    }
                        break;
                        
                    case 5:
                        // top-right-front
                    {
                        std::size_t nodes[] = {17, 14, 21, 19, 22, 5, 26, 24};
                        for (std::size_t i = 0; i < 8; ++i)
                            LocalNodeIdx.push_back(nodes[i]);
                    }
                        break;
                        
                    case 6:
                        // top-left-rear
                    {
                        std::size_t nodes[] = {18, 21, 15, 20, 23, 26, 6, 25};
                        for (std::size_t i = 0; i < 8; ++i)
                            LocalNodeIdx.push_back(nodes[i]);
                    }
                        break;
                        
                    case 7:
                        // top-right-rear
                    {
                        std::size_t nodes[] = {21, 19, 20, 16, 26, 24, 25, 7};
                        for (std::size_t i = 0; i < 8; ++i)
                            LocalNodeIdx.push_back(nodes[i]);
                    }
                        break;
                }
                break;
                
            default: // not refined
                break;
        }
    } // end of subElementLocalNodeIndices
    
    void RefinementHexahedron::adjacentSubElementsPairs(std::size_t refineType, VectorOfIndicesT& elemIdx1, VectorOfIndicesT& localFaceIdx1, VectorOfIndicesT& elemIdx2, VectorOfIndicesT& localFaceIdx2) const
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
                localFaceIdx1.push_back(2); // elem-1, local face-2
                elemIdx2.push_back(0);
                localFaceIdx2.push_back(3); // elem-0, local face-3
                        
                break;
                
            case 1:
                // the first elements pair
                elemIdx1.push_back(1);
                localFaceIdx1.push_back(1); // elem-1, local face-1
                elemIdx2.push_back(0);
                localFaceIdx2.push_back(4); // elem-0, local face-4
                        
                break;
                
            case 2:
                // the first elements pair
                elemIdx1.push_back(1);
                localFaceIdx1.push_back(0); // elem-1, local face-0
                elemIdx2.push_back(0);
                localFaceIdx2.push_back(5); // elem-0, local face-5
                        
                break;
                
            case 3:
                // the first elements pair
                elemIdx1.push_back(1);
                localFaceIdx1.push_back(2); // elem-1, local face-2
                elemIdx2.push_back(0);
                localFaceIdx2.push_back(3); // elem-0, local face-3
                        
                // the second elements pair
                elemIdx1.push_back(2);
                localFaceIdx1.push_back(1); // elem-2, local face-1
                elemIdx2.push_back(0);
                localFaceIdx2.push_back(4); // elem-0, local face-4
                        
                // the third elements pair
                elemIdx1.push_back(3);
                localFaceIdx1.push_back(1); // elem-3, local face-1
                elemIdx2.push_back(1);
                localFaceIdx2.push_back(4); // elem-1, local face-4
                        
                // the fourth elements pair
                elemIdx1.push_back(3);
                localFaceIdx1.push_back(2); // elem-3, local face-2
                elemIdx2.push_back(2);
                localFaceIdx2.push_back(3); // elem-2, local face-3
                        
                break;
                
            case 4:
                // the first elements pair
                elemIdx1.push_back(1);
                localFaceIdx1.push_back(2); // elem-1, local face-2
                elemIdx2.push_back(0);
                localFaceIdx2.push_back(3); // elem-0, local face-3
                        
                // the second elements pair
                elemIdx1.push_back(2);
                localFaceIdx1.push_back(0); // elem-2, local face-0
                elemIdx2.push_back(0);
                localFaceIdx2.push_back(5); // elem-0, local face-5
                        
                // the third elements pair
                elemIdx1.push_back(3);
                localFaceIdx1.push_back(0); // elem-3, local face-0
                elemIdx2.push_back(1);
                localFaceIdx2.push_back(5); // elem-1, local face-5
                        
                // the fourth elements pair
                elemIdx1.push_back(3);
                localFaceIdx1.push_back(2); // elem-3, local face-2
                elemIdx2.push_back(2);
                localFaceIdx2.push_back(3); // elem-2, local face-3
                        
                break;
                
            case 5:
                // the first elements pair
                elemIdx1.push_back(1);
                localFaceIdx1.push_back(1); // elem-1, local face-1
                elemIdx2.push_back(0);
                localFaceIdx2.push_back(4); // elem-0, local face-4
                        
                // the second elements pair
                elemIdx1.push_back(2);
                localFaceIdx1.push_back(0); // elem-2, local face-0
                elemIdx2.push_back(0);
                localFaceIdx2.push_back(5); // elem-0, local face-5
                        
                // the third elements pair
                elemIdx1.push_back(3);
                localFaceIdx1.push_back(0); // elem-3, local face-0
                elemIdx2.push_back(1);
                localFaceIdx2.push_back(5); // elem-1, local face-5
                        
                // the fourth elements pair
                elemIdx1.push_back(3);
                localFaceIdx1.push_back(1); // elem-3, local face-1
                elemIdx2.push_back(2);
                localFaceIdx2.push_back(4); // elem-2, local face-4
                        
                break;
                
            case 6:
                // the first elements pair
                elemIdx1.push_back(1);
                localFaceIdx1.push_back(2); // elem-1, local face-1
                elemIdx2.push_back(0);
                localFaceIdx2.push_back(3); // elem-0, local face-4
                        
                // the second elements pair
                elemIdx1.push_back(2);
                localFaceIdx1.push_back(1); // elem-2, local face-0
                elemIdx2.push_back(0);
                localFaceIdx2.push_back(4); // elem-0, local face-5
                        
                // the third elements pair
                elemIdx1.push_back(3);
                localFaceIdx1.push_back(1); // elem-3, local face-0
                elemIdx2.push_back(1);
                localFaceIdx2.push_back(4); // elem-1, local face-5
                        
                // the fourth elements pair
                elemIdx1.push_back(3);
                localFaceIdx1.push_back(2); // elem-3, local face-1
                elemIdx2.push_back(2);
                localFaceIdx2.push_back(3); // elem-2, local face-4
                        
                // the fifth elements pair
                elemIdx1.push_back(4);
                localFaceIdx1.push_back(0); // elem-1, local face-1
                elemIdx2.push_back(0);
                localFaceIdx2.push_back(5); // elem-0, local face-4
                        
                // the sixth elements pair
                elemIdx1.push_back(5);
                localFaceIdx1.push_back(0); // elem-2, local face-0
                elemIdx2.push_back(1);
                localFaceIdx2.push_back(5); // elem-0, local face-5
                        
                // the seventh elements pair
                elemIdx1.push_back(6);
                localFaceIdx1.push_back(0); // elem-3, local face-0
                elemIdx2.push_back(2);
                localFaceIdx2.push_back(5); // elem-1, local face-5
                        
                // the eighth elements pair
                elemIdx1.push_back(7);
                localFaceIdx1.push_back(0); // elem-3, local face-1
                elemIdx2.push_back(3);
                localFaceIdx2.push_back(5); // elem-2, local face-4
                        
                // the nineth elements pair
                elemIdx1.push_back(5);
                localFaceIdx1.push_back(2); // elem-1, local face-1
                elemIdx2.push_back(4);
                localFaceIdx2.push_back(3); // elem-0, local face-4
                        
                // the tenth elements pair
                elemIdx1.push_back(6);
                localFaceIdx1.push_back(1); // elem-2, local face-0
                elemIdx2.push_back(4);
                localFaceIdx2.push_back(4); // elem-0, local face-5
                        
                // the eleventh  elements pair
                elemIdx1.push_back(7);
                localFaceIdx1.push_back(1); // elem-3, local face-0
                elemIdx2.push_back(5);
                localFaceIdx2.push_back(4); // elem-1, local face-5
                        
                // the twelveth elements pair
                elemIdx1.push_back(7);
                localFaceIdx1.push_back(2); // elem-3, local face-1
                elemIdx2.push_back(6);
                localFaceIdx2.push_back(3); // elem-2, local face-4
                        
                break;
                
        } // end of switch
    }
    
    std::size_t RefinementHexahedron::nrOfSubElementsOnFace(std::size_t refineType, std::size_t faLocalIndex) const
    {
        switch (refineType)
        {
            case 0:
                switch (faLocalIndex)
                {
                    case 0:
                    case 1:
                    case 4:
                    case 5:
                        return 2;
                        
                    case 2:
                    case 3:
                        return 1;
                }
                break;
                
            case 1:
                switch (faLocalIndex)
                {
                    case 0:
                    case 2:
                    case 3:
                    case 5:
                        return 2;
                        
                    case 1:
                    case 4:
                        return 1;
                }
                break;
                
            case 2:
                switch (faLocalIndex)
                {
                    case 0:
                    case 5:
                        return 1;
                        
                    case 1:
                    case 2:
                    case 3:
                    case 4:
                        return 2;
                }
                break;
                
            case 3:
                switch (faLocalIndex)
                {
                    case 0:
                    case 5:
                        return 4;
                        
                    case 1:
                    case 2:
                    case 3:
                    case 4:
                        return 2;
                }
                break;
                
            case 4:
                switch (faLocalIndex)
                {
                    case 0:
                    case 2:
                    case 3:
                    case 5:
                        return 2;
                        
                    case 1:
                    case 4:
                        return 4;
                }
                break;
                
            case 5:
                switch (faLocalIndex)
                {
                    case 0:
                    case 1:
                    case 4:
                    case 5:
                        return 2;
                        
                    case 2:
                    case 3:
                        return 4;
                }
                break;
                
            case 6:
                return 4;
        } // end switch refinement
        
        return 1;
    }
    
    void RefinementHexahedron::subElementsOnFace(std::size_t refineType, std::size_t faLocalIndex, VectorOfIndicesT& localSubElemIdx) const
    {
        localSubElemIdx.clear();
        switch (refineType)
        {
            case 0:
                switch (faLocalIndex)
                {
                    case 0:
                    case 1:
                    case 4:
                    case 5:
                        localSubElemIdx.push_back(0);
                        localSubElemIdx.push_back(1);
                        break;
                        
                    case 2:
                        localSubElemIdx.push_back(0);
                        break;
                        
                    case 3:
                        localSubElemIdx.push_back(1);
                        break;
                }
                break;
                
            case 1:
                switch (faLocalIndex)
                {
                    case 0:
                    case 2:
                    case 3:
                    case 5:
                        localSubElemIdx.push_back(0);
                        localSubElemIdx.push_back(1);
                        break;
                        
                    case 1:
                        localSubElemIdx.push_back(0);
                        break;
                        
                    case 4:
                        localSubElemIdx.push_back(1);
                        break;
                }
                break;
                
            case 2:
                switch (faLocalIndex)
                {
                    case 0:
                        localSubElemIdx.push_back(0);
                        break;
                        
                    case 1:
                    case 2:
                    case 3:
                    case 4:
                        localSubElemIdx.push_back(0);
                        localSubElemIdx.push_back(1);
                        break;
                        
                    case 5:
                        localSubElemIdx.push_back(1);
                        break;
                }
                break;
                
            case 3:
                switch (faLocalIndex)
                {
                    case 0:
                        localSubElemIdx.push_back(0);
                        localSubElemIdx.push_back(1);
                        localSubElemIdx.push_back(2);
                        localSubElemIdx.push_back(3);
                        break;
                        
                    case 1:
                        localSubElemIdx.push_back(0);
                        localSubElemIdx.push_back(1);
                        break;
                        
                    case 2:
                        localSubElemIdx.push_back(0);
                        localSubElemIdx.push_back(2);
                        break;
                        
                    case 3:
                        localSubElemIdx.push_back(1);
                        localSubElemIdx.push_back(3);
                        break;
                        
                    case 4:
                        localSubElemIdx.push_back(2);
                        localSubElemIdx.push_back(3);
                        break;
                        
                    case 5:
                        localSubElemIdx.push_back(0);
                        localSubElemIdx.push_back(1);
                        localSubElemIdx.push_back(2);
                        localSubElemIdx.push_back(3);
                        break;
                }
                break;
                
            case 4:
                switch (faLocalIndex)
                {
                    case 0:
                        localSubElemIdx.push_back(0);
                        localSubElemIdx.push_back(1);
                        break;
                        
                    case 1:
                        localSubElemIdx.push_back(0);
                        localSubElemIdx.push_back(1);
                        localSubElemIdx.push_back(2);
                        localSubElemIdx.push_back(3);
                        break;
                        
                    case 2:
                        localSubElemIdx.push_back(0);
                        localSubElemIdx.push_back(2);
                        break;
                        
                    case 3:
                        localSubElemIdx.push_back(1);
                        localSubElemIdx.push_back(3);
                        break;
                        
                    case 4:
                        localSubElemIdx.push_back(0);
                        localSubElemIdx.push_back(1);
                        localSubElemIdx.push_back(2);
                        localSubElemIdx.push_back(3);
                        break;
                        
                    case 5:
                        localSubElemIdx.push_back(2);
                        localSubElemIdx.push_back(3);
                        break;
                }
                break;
                
            case 5:
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
                        localSubElemIdx.push_back(0);
                        localSubElemIdx.push_back(1);
                        localSubElemIdx.push_back(2);
                        localSubElemIdx.push_back(3);
                        break;
                        
                    case 3:
                        localSubElemIdx.push_back(0);
                        localSubElemIdx.push_back(1);
                        localSubElemIdx.push_back(2);
                        localSubElemIdx.push_back(3);
                        break;
                        
                    case 4:
                        localSubElemIdx.push_back(1);
                        localSubElemIdx.push_back(3);
                        break;
                        
                    case 5:
                        localSubElemIdx.push_back(2);
                        localSubElemIdx.push_back(3);
                        break;
                }
                break;
                
            case 6:
                switch (faLocalIndex)
                {
                    case 0:
                        localSubElemIdx.push_back(0);
                        localSubElemIdx.push_back(1);
                        localSubElemIdx.push_back(2);
                        localSubElemIdx.push_back(3);
                        break;
                        
                    case 1:
                        localSubElemIdx.push_back(0);
                        localSubElemIdx.push_back(1);
                        localSubElemIdx.push_back(4);
                        localSubElemIdx.push_back(5);
                        break;
                        
                    case 2:
                        localSubElemIdx.push_back(0);
                        localSubElemIdx.push_back(2);
                        localSubElemIdx.push_back(4);
                        localSubElemIdx.push_back(6);
                        break;
                        
                    case 3:
                        localSubElemIdx.push_back(1);
                        localSubElemIdx.push_back(3);
                        localSubElemIdx.push_back(5);
                        localSubElemIdx.push_back(7);
                        break;
                        
                    case 4:
                        localSubElemIdx.push_back(2);
                        localSubElemIdx.push_back(3);
                        localSubElemIdx.push_back(6);
                        localSubElemIdx.push_back(7);
                        break;
                        
                    case 5:
                        localSubElemIdx.push_back(4);
                        localSubElemIdx.push_back(5);
                        localSubElemIdx.push_back(6);
                        localSubElemIdx.push_back(7);
                        break;
                }
                break;
                
            default:
                localSubElemIdx.push_back(0);
                break;
        } // end switch refinement
    }
    
    std::size_t RefinementHexahedron::getLocalSubFaceNr(std::size_t refineType, std::size_t localFaceNr, std::size_t subElementIdx) const
    {
        return localFaceNr;
    }

} // end of namespace Geometry
