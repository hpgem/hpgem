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
#include "Geometry/RefinementHexahedron.hpp"

#include "PointPhysical.hpp"
#include "PointReference.hpp"
#include "ReferenceGeometry.hpp"
#include "PhysicalGeometry.hpp"

namespace Geometry
{
    unsigned int RefinementHexahedron::nrOfNewNodes(int refineType) const 
    { 
      switch (refineType)
      {
        case 0:    // x-refinement
        case 1:    // y-refinement
        case 2:    // z-refinement
          return 4;

        case 3:    // xy-refinement
        case 4:    // xz-refinement
        case 5:    // yz-refinement
          return 10;

        case 6:    // xyz-refinement
          return 19;
      }
      
      return 0;
    }

    void RefinementHexahedron::getAllNodes(int refineType, VectorOfPointPhysicalsT& nodes) const 
    {
        // get all element's nodes
        nodes.clear();
        PointPhysicalT p(3);
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
            nodes.push_back(0.5 * (nodes[4]+nodes[5]));  // between 4-5
            nodes.push_back(0.5 * (nodes[6]+nodes[7]));  // between 6-7
            break;
            
          case 1:  // y-refinement
            nodes.push_back(0.5 * (nodes[0]+nodes[2]));  // between 0-2
            nodes.push_back(0.5 * (nodes[1]+nodes[3]));  // between 1-3
            nodes.push_back(0.5 * (nodes[4]+nodes[6]));  // between 4-6
            nodes.push_back(0.5 * (nodes[5]+nodes[7]));  // between 5-7
            break;
            
          case 2:  // z-refinement
            nodes.push_back(0.5 * (nodes[0]+nodes[4]));  // between 0-4
            nodes.push_back(0.5 * (nodes[1]+nodes[5]));  // between 1-5
            nodes.push_back(0.5 * (nodes[2]+nodes[6]));  // between 2-6
            nodes.push_back(0.5 * (nodes[3]+nodes[7]));  // between 3-7
            break;

          case 3:    // xy-refinement
            nodes.push_back(0.5 * (nodes[0]+nodes[1]));    // between 0-1
            nodes.push_back(0.5 * (nodes[0]+nodes[2]));    // between 0-2
            nodes.push_back(0.5 * (nodes[1]+nodes[3]));    // between 1-3
            nodes.push_back(0.5 * (nodes[2]+nodes[3]));    // between 2-3
            nodes.push_back(0.5 * (nodes[9]+nodes[10]));   // between 9-10
            nodes.push_back(0.5 * (nodes[4]+nodes[5]));    // between 4-5
            nodes.push_back(0.5 * (nodes[4]+nodes[6]));    // between 4-6
            nodes.push_back(0.5 * (nodes[5]+nodes[7]));    // between 5-7
            nodes.push_back(0.5 * (nodes[6]+nodes[7]));    // between 6-7
            nodes.push_back(0.5 * (nodes[14]+nodes[15]));  // between 14-15
            break;
            
          case 4:    // xz-refinement
            
            nodes.push_back(0.5 * (nodes[0]+nodes[1]));    // between 0-1
            nodes.push_back(0.5 * (nodes[2]+nodes[3]));    // between 2-3
            nodes.push_back(0.5 * (nodes[0]+nodes[4]));    // between 0-4
            nodes.push_back(0.5 * (nodes[1]+nodes[5]));    // between 1-5
            nodes.push_back(0.5 * (nodes[2]+nodes[6]));    // between 2-6
            nodes.push_back(0.5 * (nodes[3]+nodes[7]));    // between 3-7
            nodes.push_back(0.5 * (nodes[10]+nodes[11]));  // between 10-11
            nodes.push_back(0.5 * (nodes[12]+nodes[13]));  // between 12-13
            nodes.push_back(0.5 * (nodes[4]+nodes[5]));    // between 4-5
            nodes.push_back(0.5 * (nodes[6]+nodes[7]));    // between 6-7
            break;
            
          case 5:    // yz-refinement
            
            nodes.push_back(0.5 * (nodes[0]+nodes[2]));    // between 0-2
            nodes.push_back(0.5 * (nodes[1]+nodes[3]));    // between 1-3
            nodes.push_back(0.5 * (nodes[0]+nodes[4]));    // between 0-4
            nodes.push_back(0.5 * (nodes[1]+nodes[5]));    // between 1-5
            nodes.push_back(0.5 * (nodes[2]+nodes[6]));    // between 2-6
            nodes.push_back(0.5 * (nodes[3]+nodes[7]));    // between 3-7
            nodes.push_back(0.5 * (nodes[10]+nodes[12]));  // between 10-12
            nodes.push_back(0.5 * (nodes[11]+nodes[13]));  // between 11-13
            nodes.push_back(0.5 * (nodes[4]+nodes[6]));    // between 4-6
            nodes.push_back(0.5 * (nodes[5]+nodes[7]));    // between 5-7
            break;
            
          case 6:    // xyz-refinement
            nodes.push_back(0.5 * (nodes[0]+nodes[1]));    // between 0-1
            nodes.push_back(0.5 * (nodes[0]+nodes[2]));    // between 0-2
            nodes.push_back(0.5 * (nodes[1]+nodes[3]));    // between 1-3
            nodes.push_back(0.5 * (nodes[2]+nodes[3]));    // between 2-3
            nodes.push_back(0.5 * (nodes[9]+nodes[10]));   // between 9-10
            nodes.push_back(0.5 * (nodes[0]+nodes[4]));    // between 0-4
            nodes.push_back(0.5 * (nodes[1]+nodes[5]));    // between 1-5
            nodes.push_back(0.5 * (nodes[2]+nodes[6]));    // between 2-6
            nodes.push_back(0.5 * (nodes[3]+nodes[7]));    // between 3-7
            nodes.push_back(0.5 * (nodes[13]+nodes[14]));  // between 13-14
            nodes.push_back(0.5 * (nodes[13]+nodes[15]));  // between 13-15
            nodes.push_back(0.5 * (nodes[14]+nodes[16]));  // between 14-16
            nodes.push_back(0.5 * (nodes[15]+nodes[16]));  // between 15-16
            nodes.push_back(0.5 * (nodes[18]+nodes[19]));  // between 18-19
            nodes.push_back(0.5 * (nodes[4]+nodes[5]));    // between 4-5
            nodes.push_back(0.5 * (nodes[4]+nodes[6]));    // between 4-6
            nodes.push_back(0.5 * (nodes[5]+nodes[7]));    // between 5-7
            nodes.push_back(0.5 * (nodes[6]+nodes[7]));    // between 6-7
            nodes.push_back(0.5 * (nodes[23]+nodes[24]));  // between 23-24
            break;
            
          default:
            break;
        }
    }
    
    unsigned int RefinementHexahedron::nrOfSubElements(int refineType) const 
    { 
        unsigned int nrSubElements(0);

        switch (refineType)
        {
            case 0:    // x-refinement
            case 1:    // y-refinement
            case 2:    // z-refinement
              nrSubElements = 2;
              break;

            case 3:    // xy-refinement
            case 4:    // xz-refinement
            case 5:    // yz-refinement
              nrSubElements = 4;
              break;

            case 6:    // xyz-refinement
              nrSubElements = 8;
              break;
              
            default:    // not refined
              nrSubElements = 1;
              break;
          
        }

        return nrSubElements;
    }

    void RefinementHexahedron::subElementLocalNodeIndices(int refineType, unsigned int iSubElement, VectorOfIndicesT& LocalNodeIdx) const 
    {
        TestErrorDebug((iSubElement<nrOfSubElements(refineType)),
                        "RefinementHexahedron: invalid sub-element index while getting its local node indices!");

        LocalNodeIdx.clear();
        switch (refineType)
        {
            case 0:    // x-refinement
              switch (iSubElement)
              {
                case 0:
                  // left
                  {
                    unsigned int nodes[] = { 0, 8, 2, 9, 4, 10, 6, 11 };
                    for (unsigned int i=0; i<8; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 1:
                  // right
                  {
                    unsigned int nodes[] = { 8, 1, 9, 3, 10, 5, 11, 7 };
                    for (unsigned int i=0; i<8; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
              }
              break;

            case 1:    // y-refinement
              switch (iSubElement)
              {
                case 0:
                  // front
                  {
                    unsigned int nodes[] = { 0, 1, 8, 9, 4, 5, 10, 11 };
                    for (unsigned int i=0; i<8; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 1:
                  // rear
                  {
                    unsigned int nodes[] = { 8, 9, 2, 3, 10, 11, 6, 7 };
                    for (unsigned int i=0; i<8; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
              }
              break;

            case 2:    // z-refinement
              switch (iSubElement)
              {
                case 0:
                  // bottom
                  {
                    unsigned int nodes[] = { 0, 1, 2, 3, 8, 9, 10, 11 };
                    for (unsigned int i=0; i<8; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 1:
                  // top
                  {
                    unsigned int nodes[] = { 8, 9, 10, 11, 4, 5, 6, 7 };
                    for (unsigned int i=0; i<8; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
              }
              break;

            case 3:    // xy-refinement
              switch (iSubElement)
              {
                case 0:
                  // left-front
                  {
                    unsigned int nodes[] = { 0, 8, 9, 12, 4, 13, 14, 17 };
                    for (unsigned int i=0; i<8; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 1:
                  // right-front
                  {
                    unsigned int nodes[] = { 8, 1, 12, 10, 13, 5, 17, 15 };
                    for (unsigned int i=0; i<8; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 2:
                  // left-rear
                  {
                    unsigned int nodes[] = { 9, 12, 2, 11, 14, 17, 6, 16 };
                    for (unsigned int i=0; i<8; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 3:
                  // right-rear
                  {
                    unsigned int nodes[] = { 12, 10, 11, 3, 17, 15, 16, 7 };
                    for (unsigned int i=0; i<8; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
              }
              break;

            case 4:    // xz-refinement
              switch (iSubElement)
              {
                case 0:
                  // bottom-left
                  {
                    unsigned int nodes[] = { 0, 8, 2, 9, 10, 14, 12, 15 };
                    for (unsigned int i=0; i<8; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 1:
                  // bottom-right
                  {
                    unsigned int nodes[] = { 8, 1, 9, 3, 14, 11, 15, 13 };
                    for (unsigned int i=0; i<8; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 2:
                  // top-left
                  {
                    unsigned int nodes[] = { 10, 14, 12, 15, 4, 16, 6, 17 };
                    for (unsigned int i=0; i<8; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 3:
                  // top-right
                  {
                    unsigned int nodes[] = { 14, 11, 15, 13, 16, 5, 17, 7 };
                    for (unsigned int i=0; i<8; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
              }
              break;

            case 5:    // yz-refinement
              switch (iSubElement)
              {
                case 0:
                  // bottom-front
                  {
                    unsigned int nodes[] = { 0, 1, 8, 9, 10, 11, 14, 15 };
                    for (unsigned int i=0; i<8; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 1:
                  // bottom-rear
                  {
                    unsigned int nodes[] = { 8, 9, 2, 3, 14, 15, 12, 13 };
                    for (unsigned int i=0; i<8; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 2:
                  // top-front
                  {
                    unsigned int nodes[] = { 10, 11, 14, 15, 4, 5, 16, 17 };
                    for (unsigned int i=0; i<8; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 3:
                  // top-rear
                  {
                    unsigned int nodes[] = { 14, 15, 12, 13, 16, 17, 6, 7 };
                    for (unsigned int i=0; i<8; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
              }
              break;

            case 6:    // xyz-refinement
              switch (iSubElement)
              {
                case 0:
                  // bottom-left-front
                  {
                    unsigned int nodes[] = { 0, 8, 9, 12, 13, 17, 18, 21 };
                    for (unsigned int i=0; i<8; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 1:
                  // bottom-right-front
                  {
                    unsigned int nodes[] = { 8, 1, 12, 10, 17, 14, 21, 19 };
                    for (unsigned int i=0; i<8; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 2:
                  // bottom-left-rear
                  {
                    unsigned int nodes[] = { 9, 12, 2, 11, 18, 21, 15, 20 };
                    for (unsigned int i=0; i<8; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 3:
                  // bottom-right-rear
                  {
                    unsigned int nodes[] = { 12, 10, 11, 3, 21, 19, 20, 16 };
                    for (unsigned int i=0; i<8; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 4:
                  // top-left-front
                  {
                    unsigned int nodes[] = { 13, 17, 18, 21, 4, 22, 23, 26 };
                    for (unsigned int i=0; i<8; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 5:
                  // top-right-front
                  {
                    unsigned int nodes[] = { 17, 14, 21, 19, 22, 5, 26, 24 };
                    for (unsigned int i=0; i<8; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 6:
                  // top-left-rear
                  {
                    unsigned int nodes[] = { 18, 21, 15, 20, 23, 26, 6, 25 };
                    for (unsigned int i=0; i<8; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 7:
                  // top-right-rear
                  {
                    unsigned int nodes[] = { 21, 19, 20, 16, 26, 24, 25, 7 };
                    for (unsigned int i=0; i<8; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
              }
              break;

            default:    // not refined
              break;
        }
    } // end of subElementLocalNodeIndices

    void RefinementHexahedron::adjacentSubElementsPairs(int refineType,
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
              elemIdx1.push_back(1);   localFaceIdx1.push_back(2);  // elem-1, local face-2
              elemIdx2.push_back(0);   localFaceIdx2.push_back(3);  // elem-0, local face-3
              
              break;

          case 1:
              // the first elements pair
              elemIdx1.push_back(1);   localFaceIdx1.push_back(1);  // elem-1, local face-1
              elemIdx2.push_back(0);   localFaceIdx2.push_back(4);  // elem-0, local face-4
              
              break;

          case 2:
              // the first elements pair
              elemIdx1.push_back(1);   localFaceIdx1.push_back(0);  // elem-1, local face-0
              elemIdx2.push_back(0);   localFaceIdx2.push_back(5);  // elem-0, local face-5
              
              break;

          case 3:
              // the first elements pair
              elemIdx1.push_back(1);   localFaceIdx1.push_back(2);  // elem-1, local face-2
              elemIdx2.push_back(0);   localFaceIdx2.push_back(3);  // elem-0, local face-3
              
              // the second elements pair
              elemIdx1.push_back(2);   localFaceIdx1.push_back(1);  // elem-2, local face-1
              elemIdx2.push_back(0);   localFaceIdx2.push_back(4);  // elem-0, local face-4
              
              // the third elements pair
              elemIdx1.push_back(3);   localFaceIdx1.push_back(1);  // elem-3, local face-1
              elemIdx2.push_back(1);   localFaceIdx2.push_back(4);  // elem-1, local face-4
              
              // the fourth elements pair
              elemIdx1.push_back(3);   localFaceIdx1.push_back(2);  // elem-3, local face-2
              elemIdx2.push_back(2);   localFaceIdx2.push_back(3);  // elem-2, local face-3
              
              break;
              
          case 4:
              // the first elements pair
              elemIdx1.push_back(1);   localFaceIdx1.push_back(2);  // elem-1, local face-2
              elemIdx2.push_back(0);   localFaceIdx2.push_back(3);  // elem-0, local face-3
              
              // the second elements pair
              elemIdx1.push_back(2);   localFaceIdx1.push_back(0);  // elem-2, local face-0
              elemIdx2.push_back(0);   localFaceIdx2.push_back(5);  // elem-0, local face-5
              
              // the third elements pair
              elemIdx1.push_back(3);   localFaceIdx1.push_back(0);  // elem-3, local face-0
              elemIdx2.push_back(1);   localFaceIdx2.push_back(5);  // elem-1, local face-5
              
              // the fourth elements pair
              elemIdx1.push_back(3);   localFaceIdx1.push_back(2);  // elem-3, local face-2
              elemIdx2.push_back(2);   localFaceIdx2.push_back(3);  // elem-2, local face-3
              
              break;
              
          case 5:
              // the first elements pair
              elemIdx1.push_back(1);   localFaceIdx1.push_back(1);  // elem-1, local face-1
              elemIdx2.push_back(0);   localFaceIdx2.push_back(4);  // elem-0, local face-4
              
              // the second elements pair
              elemIdx1.push_back(2);   localFaceIdx1.push_back(0);  // elem-2, local face-0
              elemIdx2.push_back(0);   localFaceIdx2.push_back(5);  // elem-0, local face-5
              
              // the third elements pair
              elemIdx1.push_back(3);   localFaceIdx1.push_back(0);  // elem-3, local face-0
              elemIdx2.push_back(1);   localFaceIdx2.push_back(5);  // elem-1, local face-5
              
              // the fourth elements pair
              elemIdx1.push_back(3);   localFaceIdx1.push_back(1);  // elem-3, local face-1
              elemIdx2.push_back(2);   localFaceIdx2.push_back(4);  // elem-2, local face-4
              
              break;
              
          case 6:
              // the first elements pair
              elemIdx1.push_back(1);   localFaceIdx1.push_back(2);  // elem-1, local face-1
              elemIdx2.push_back(0);   localFaceIdx2.push_back(3);  // elem-0, local face-4
              
              // the second elements pair
              elemIdx1.push_back(2);   localFaceIdx1.push_back(1);  // elem-2, local face-0
              elemIdx2.push_back(0);   localFaceIdx2.push_back(4);  // elem-0, local face-5
              
              // the third elements pair
              elemIdx1.push_back(3);   localFaceIdx1.push_back(1);  // elem-3, local face-0
              elemIdx2.push_back(1);   localFaceIdx2.push_back(4);  // elem-1, local face-5
              
              // the fourth elements pair
              elemIdx1.push_back(3);   localFaceIdx1.push_back(2);  // elem-3, local face-1
              elemIdx2.push_back(2);   localFaceIdx2.push_back(3);  // elem-2, local face-4
              
              // the fifth elements pair
              elemIdx1.push_back(4);   localFaceIdx1.push_back(0);  // elem-1, local face-1
              elemIdx2.push_back(0);   localFaceIdx2.push_back(5);  // elem-0, local face-4
              
              // the sixth elements pair
              elemIdx1.push_back(5);   localFaceIdx1.push_back(0);  // elem-2, local face-0
              elemIdx2.push_back(1);   localFaceIdx2.push_back(5);  // elem-0, local face-5
              
              // the seventh elements pair
              elemIdx1.push_back(6);   localFaceIdx1.push_back(0);  // elem-3, local face-0
              elemIdx2.push_back(2);   localFaceIdx2.push_back(5);  // elem-1, local face-5
              
              // the eighth elements pair
              elemIdx1.push_back(7);   localFaceIdx1.push_back(0);  // elem-3, local face-1
              elemIdx2.push_back(3);   localFaceIdx2.push_back(5);  // elem-2, local face-4
              
              // the nineth elements pair
              elemIdx1.push_back(5);   localFaceIdx1.push_back(2);  // elem-1, local face-1
              elemIdx2.push_back(4);   localFaceIdx2.push_back(3);  // elem-0, local face-4
              
              // the tenth elements pair
              elemIdx1.push_back(6);   localFaceIdx1.push_back(1);  // elem-2, local face-0
              elemIdx2.push_back(4);   localFaceIdx2.push_back(4);  // elem-0, local face-5
              
              // the eleventh  elements pair
              elemIdx1.push_back(7);   localFaceIdx1.push_back(1);  // elem-3, local face-0
              elemIdx2.push_back(5);   localFaceIdx2.push_back(4);  // elem-1, local face-5
              
              // the twelveth elements pair
              elemIdx1.push_back(7);   localFaceIdx1.push_back(2);  // elem-3, local face-1
              elemIdx2.push_back(6);   localFaceIdx2.push_back(3);  // elem-2, local face-4
              
              break;
              
        } // end of switch
    }

    unsigned int RefinementHexahedron::nrOfSubElementsOnFace(int refineType, unsigned int faLocalIndex) const 
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

    void RefinementHexahedron::subElementsOnFace(int refineType, unsigned int faLocalIndex, VectorOfIndicesT& localSubElemIdx) const 
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
    
    unsigned int RefinementHexahedron::getLocalSubFaceNr(int refineType, unsigned int localFaceNr, unsigned int subElementIdx) const 
    { 
      return localFaceNr; 
    }

}  // end of namespace Geometry
