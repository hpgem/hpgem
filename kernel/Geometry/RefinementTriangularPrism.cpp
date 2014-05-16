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
#include "Geometry/RefinementTriangularPrism.hpp"

#include "PointPhysical.hpp"
#include "PointReference.hpp"
#include "ReferenceGeometry.hpp"
#include "PhysicalGeometry.hpp"

namespace Geometry
{
    unsigned int RefinementTriangularPrism::nrOfNewNodes(int refineType) const 
    { 
        switch (refineType)
        {
          case 0:   // r0 refinement
          case 1:   // r1 refinement
          case 2:   // r2 refinement
              return 2;
              
          case 3:   // r3 refinement
              return 3;
            
          case 4:      // r4 refinement
              return 6;
              
          case 5:      // r5 refinement
              return 8;

          case 6:      // r6 refinement
          case 7:      // r6 refinement
              return 4;
        }
        
        return 0;
    }
    
    void RefinementTriangularPrism::getAllNodes(int refineType, VectorOfPointPhysicalsT& nodes) const 
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
          case 0:  // r0: face-2 refinement, split into 2 subelements
            nodes.push_back(0.5 * (nodes[0]+nodes[2]));  // between 0-2
            nodes.push_back(0.5 * (nodes[3]+nodes[5]));  // between 3-5
            break;
            
          case 1:  // r1: face-3 refinement, split into 2 subelements
            nodes.push_back(0.5 * (nodes[0]+nodes[1]));  // between 0-1
            nodes.push_back(0.5 * (nodes[3]+nodes[4]));  // between 3-4
            break;
            
          case 2:  // r2: face-4 refinement, split into 2 subelements
            nodes.push_back(0.5 * (nodes[1]+nodes[2]));  // between 1-2
            nodes.push_back(0.5 * (nodes[4]+nodes[5]));  // between 4-5
            break;
            
          case 3: // r3: face-234 refinement, split into 2 subelements
            nodes.push_back(0.5 * (nodes[0]+nodes[3]));  // between 0-3
            nodes.push_back(0.5 * (nodes[1]+nodes[4]));  // between 1-4
            nodes.push_back(0.5 * (nodes[2]+nodes[5]));  // between 2-5
            break;
            
          case 4: // r4: all-faces refinement, split into 4 subelements
            nodes.push_back(0.5 * (nodes[0]+nodes[1]));  // between 0-1
            nodes.push_back(0.5 * (nodes[0]+nodes[2]));  // between 0-2
            nodes.push_back(0.5 * (nodes[1]+nodes[2]));  // between 1-2
            
            nodes.push_back(0.5 * (nodes[3]+nodes[4]));  // between 3-4
            nodes.push_back(0.5 * (nodes[3]+nodes[5]));  // between 3-5
            nodes.push_back(0.5 * (nodes[4]+nodes[5]));  // between 4-5
            break;
            
          case 5: // r5: split into 3 hexahedron elements
            nodes.push_back(0.5 * (nodes[0]+nodes[1]));  // between 0-1
            nodes.push_back(0.5 * (nodes[0]+nodes[2]));  // between 0-2
            nodes.push_back(0.5 * (nodes[1]+nodes[2]));  // between 1-2
            nodes.push_back((1./3.) * (nodes[0]+nodes[1]+nodes[2]));  // centroid of triangle 0-1-2

            nodes.push_back(0.5 * (nodes[3]+nodes[4]));  // between 3-4
            nodes.push_back(0.5 * (nodes[3]+nodes[5]));  // between 3-5
            nodes.push_back(0.5 * (nodes[4]+nodes[5]));  // between 4-5
            nodes.push_back((1./3.) * (nodes[3]+nodes[4]+nodes[5]));  // centroid of triangle 3-4-5
            break;
            
          case 6: // r6: split into hexahedron and triangular-prism elements
            nodes.push_back(0.5 * (nodes[0]+nodes[1]));  // between 0-1
            nodes.push_back(0.5 * (nodes[1]+nodes[2]));  // between 1-2
            nodes.push_back(0.5 * (nodes[3]+nodes[4]));  // between 3-4
            nodes.push_back(0.5 * (nodes[4]+nodes[5]));  // between 4-5
            break;
            
          case 7: // r7: split into hexahedron and triangular-prism elements
            nodes.push_back(0.5 * (nodes[0]+nodes[2]));  // between 0-2
            nodes.push_back(0.5 * (nodes[1]+nodes[2]));  // between 1-2
            nodes.push_back(0.5 * (nodes[3]+nodes[5]));  // between 3-5
            nodes.push_back(0.5 * (nodes[4]+nodes[5]));  // between 4-5
            break;
        }
            
    }
    
    unsigned int RefinementTriangularPrism::nrOfSubElements(int refineType) const 
    { 
        switch (refineType)
        {
            case 0: // r0 refinement
            case 1: // r1 refinement
            case 2: // r2 refinement
            case 3: // r3 refinement
            case 6: // r6 refinement
            case 7: // r7 refinement
              return 2;
              
            case 4:    // r4 refinement
              return 4;
            
            case 5:    // r5 refinement
              return 3;
        }

        return 1;
    }

    void RefinementTriangularPrism::subElementLocalNodeIndices(int refineType, unsigned int iSubElement, VectorOfIndicesT& LocalNodeIdx) const 
    {
        TestErrorDebug((iSubElement<nrOfSubElements(refineType)),
                        "RefinementTriangularPrism: invalid sub-element index while getting its local node indices!");

        LocalNodeIdx.clear();
        switch (refineType)
        {
          case 0:  // r0: face-2 refinement, split into 2 subelements
              switch (iSubElement)
              {
                case 0:
                  // angle-1
                  {
                    unsigned int nodes[] = { 0, 1, 6, 3, 4, 7 };
                    for (unsigned int i=0; i<6; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 1:
                  // angle-1
                  {
                    unsigned int nodes[] = { 6, 1, 2, 7, 4, 5 };
                    for (unsigned int i=0; i<6; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
              }
              break;
              
          case 1:  // r1: face-3 refinement, split into 2 subelements
              switch (iSubElement)
              {
                case 0:
                  // angle-2
                  {
                    unsigned int nodes[] = { 0, 6, 2, 3, 7, 5 };
                    for (unsigned int i=0; i<6; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 1:
                  // angle-2
                  {
                    unsigned int nodes[] = { 6, 1, 2, 7, 4, 5 };
                    for (unsigned int i=0; i<6; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
              }
              break;

          case 2:  // r2: face-4 refinement, split into 2 subelements
              switch (iSubElement)
              {
                case 0:
                  // angle-0
                  {
                    unsigned int nodes[] = { 0, 1, 6, 3, 4, 7 };
                    for (unsigned int i=0; i<6; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 1:
                  // angle-0
                  {
                    unsigned int nodes[] = { 0, 6, 2, 3, 7, 5 };
                    for (unsigned int i=0; i<6; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
              }
              break;
              
          case 3: // r3: face-234 refinement, split into 2 subelements
              switch (iSubElement)
              {
                case 0:
                  {
                    unsigned int nodes[] = { 0, 1, 2, 6, 7, 8 };
                    for (unsigned int i=0; i<6; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 1:
                  {
                    unsigned int nodes[] = { 6, 7, 8, 3, 4, 5 };
                    for (unsigned int i=0; i<6; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
              }
              break;
              
          case 4: // r4: all-faces refinement, split into 4 subelements
              switch (iSubElement)
              {
                case 0:
                  // angle-0
                  {
                    unsigned int nodes[] = { 0, 6, 7, 3, 9, 10 };
                    for (unsigned int i=0; i<6; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 1:
                  // angle-1
                  {
                    unsigned int nodes[] = { 6, 1, 8, 9, 4, 11 };
                    for (unsigned int i=0; i<6; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 2:
                  // angle-2
                  {
                    unsigned int nodes[] = { 7, 8, 2, 10, 11, 5 };
                    for (unsigned int i=0; i<6; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 3:
                  // middle
                  {
                    unsigned int nodes[] = { 8, 7, 6, 11, 10, 9 };
                    for (unsigned int i=0; i<6; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
              }
              break;
              
          case 5: // r5: split into 3 hexahedron elements
              switch (iSubElement)
              {
                case 0:
                  // hexahedron-0
                  {
                    unsigned int nodes[] = { 7, 0, 9, 6, 11, 3, 13, 10 };
                    for (unsigned int i=0; i<8; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 1:
                  // hexahedron-1
                  {
                    unsigned int nodes[] = { 6, 1, 9, 8, 10, 4, 13, 12 };
                    for (unsigned int i=0; i<8; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 2:
                  // hexahedron-2
                  {
                    unsigned int nodes[] = { 8, 2, 9, 7, 12, 5, 13, 11 };
                    for (unsigned int i=0; i<8; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
              }
              break;
              
          case 6: // r6: split into hexahedron and triangular-prism elements
              switch (iSubElement)
              {
                case 0:
                  // hexahedron
                  {
                    unsigned int nodes[] = { 0, 6, 2, 7, 3, 8, 5, 9 };
                    for (unsigned int i=0; i<8; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 1:
                  // triangular-prism
                  {
                    unsigned int nodes[] = { 6, 1, 7, 8, 4, 9 };
                    for (unsigned int i=0; i<6; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
              }
              break;
              
          case 7: // r7: split into hexahedron and triangular-prism elements
              switch (iSubElement)
              {
                case 0:
                  // hexahedron
                  {
                    unsigned int nodes[] = { 0, 1, 6, 7, 3, 4, 8, 9 };
                    for (unsigned int i=0; i<8; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
                  
                case 1:
                  // triangular-prism
                  {
                    unsigned int nodes[] = { 6, 7, 2, 8, 9, 5 };
                    for (unsigned int i=0; i<6; ++i)
                      LocalNodeIdx.push_back(nodes[i]);
                  }
                  break;
              }
              break;
        }
        
    }
    
    void RefinementTriangularPrism::adjacentSubElementsPairs(int refineType,
                    VectorOfIndicesT& elemIdx1, VectorOfIndicesT& localFaceIdx1,
                    VectorOfIndicesT& elemIdx2, VectorOfIndicesT& localFaceIdx2) const 
    {
        elemIdx1.clear();
        elemIdx2.clear();
        localFaceIdx1.clear();
        localFaceIdx2.clear();
        
        switch (refineType) 
        {
          case 0:   // r0: face-2 refinement, split into 2 subelements
              elemIdx1.push_back(1);   localFaceIdx1.push_back(3);  // elem-1, local face-3
              elemIdx2.push_back(0);   localFaceIdx2.push_back(4);  // elem-0, local face-4
              break;

          case 1:   // r1: face-3 refinement, split into 2 subelements
              elemIdx1.push_back(1);   localFaceIdx1.push_back(2);  // elem-1, local face-2
              elemIdx2.push_back(0);   localFaceIdx2.push_back(4);  // elem-0, local face-4
              break;

          case 2:   // r2: face-4 refinement, split into 2 subelements
              elemIdx1.push_back(1);   localFaceIdx1.push_back(3);  // elem-1, local face-3
              elemIdx2.push_back(0);   localFaceIdx2.push_back(2);  // elem-0, local face-2
              break;

          case 3:   // r3: face-234 refinement, split into 2 subelements
              elemIdx1.push_back(1);   localFaceIdx1.push_back(0);  // elem-1, local face-0
              elemIdx2.push_back(0);   localFaceIdx2.push_back(1);  // elem-0, local face-1
              break;

          case 4:  // r4: all-faces refinement, split into 4 subelements
              // the first elements pair
              elemIdx1.push_back(3);   localFaceIdx1.push_back(4);  // elem-3, local face-4
              elemIdx2.push_back(0);   localFaceIdx2.push_back(4);  // elem-0, local face-4
              
              // the second elements pair
              elemIdx1.push_back(3);   localFaceIdx1.push_back(2);  // elem-3, local face-2
              elemIdx2.push_back(1);   localFaceIdx2.push_back(2);  // elem-1, local face-2
              
              // the third elements pair
              elemIdx1.push_back(3);   localFaceIdx1.push_back(3);  // elem-3, local face-3
              elemIdx2.push_back(2);   localFaceIdx2.push_back(3);  // elem-2, local face-3
              break;
              
          case 5:  // r5: split into 3 hexahedron elements
              // the first elements pair
              elemIdx1.push_back(1);   localFaceIdx1.push_back(2);  // elem-1, local face-2
              elemIdx2.push_back(0);   localFaceIdx2.push_back(4);  // elem-0, local face-4
              
              // the second elements pair
              elemIdx1.push_back(2);   localFaceIdx1.push_back(4);  // elem-2, local face-4
              elemIdx2.push_back(0);   localFaceIdx2.push_back(2);  // elem-0, local face-2
              
              // the third elements pair
              elemIdx1.push_back(2);   localFaceIdx1.push_back(2);  // elem-2, local face-2
              elemIdx2.push_back(1);   localFaceIdx2.push_back(4);  // elem-1, local face-4
              break;
              
          case 6:  // r6: split into hexahedron and triangular-prism elements
              // the first elements pair
              elemIdx1.push_back(1);   localFaceIdx1.push_back(2);  // elem-1, local face-2
              elemIdx2.push_back(0);   localFaceIdx2.push_back(3);  // elem-0, local face-3
              break;
              
          case 7:  // r7: split into hexahedron and triangular-prism elements
              // the first elements pair
              elemIdx1.push_back(1);   localFaceIdx1.push_back(3);  // elem-1, local face-3
              elemIdx2.push_back(0);   localFaceIdx2.push_back(4);  // elem-0, local face-4
              break;
        } // end of switch
    }

    unsigned int RefinementTriangularPrism::nrOfSubElementsOnFace(int refineType, unsigned int faLocalIndex) const 
    { 
        switch (refineType)
        {
          case 0:  // r0: face-2 refinement, split into 2 subelements
            switch (faLocalIndex)
            {
              case 0:
              case 1:
              case 2:
                return 2;

              case 3:
              case 4:
                return 1;
            }
            break;

          case 1:  // r1: face-3 refinement, split into 2 subelements
            switch (faLocalIndex)
            {
              case 0:
              case 1:
              case 3:
                return 2;

              case 2:
              case 4:
                return 1;
            }
            break;

          case 2:  // r2: face-4 refinement, split into 2 subelements
            switch (faLocalIndex)
            {
              case 0:
              case 1:
              case 4:
                return 2;

              case 2:
              case 3:
                return 1;
            }
            break;

          case 3: // r3: face-234 refinement, split into 2 subelements
            switch (faLocalIndex)
            {
              case 0:
              case 1:
                return 1;

              case 2:
              case 3:
              case 4:
                return 2;
            }
            break;

          case 4: // r4: all-faces refinement, split into 4 subelements
            switch (faLocalIndex)
            {
              case 0:
              case 1:
                return 4;
            
              case 2:
              case 3:
              case 4:
                return 2;
            }
            break;
            
          case 5: // r5: split into 3 hexahedron elements
            switch (faLocalIndex)
            {
              case 0:
              case 1:
                return 3;
            
              case 2:
              case 3:
              case 4:
                return 2;
            }
            break;

          case 6: // r6: split into hexahedron and triangular-prism elements
            switch (faLocalIndex)
            {
              case 0:
              case 1:
              case 3:
              case 4:
                return 2;
            
              case 2:
                return 1;
            }
            break;

          case 7: // r7: split into hexahedron and triangular-prism elements
            switch (faLocalIndex)
            {
              case 0:
              case 1:
              case 2:
              case 4:
                return 2;
            
              case 3:
                return 1;
            }
            break;
        } // end switch refinement
        
        return 1;
    }

    void RefinementTriangularPrism::subElementsOnFace(int refineType, unsigned int faLocalIndex, VectorOfIndicesT& localSubElemIdx) const 
    {
        localSubElemIdx.clear();

        switch (refineType)
        {
          case 0:  // r0: face-2 refinement, split into 2 subelements
            switch (faLocalIndex)
            {
              case 0:
              case 1:
              case 2:
                  localSubElemIdx.push_back(0);
                  localSubElemIdx.push_back(1);
                  break;
            
              case 3:
                  localSubElemIdx.push_back(0);
                  break;
            
              case 4:
                  localSubElemIdx.push_back(1);
                  break;
            }
            break;

          case 1:  // r1: face-3 refinement, split into 2 subelements
            switch (faLocalIndex)
            {
              case 0:
              case 1:
              case 3:
                  localSubElemIdx.push_back(0);
                  localSubElemIdx.push_back(1);
                  break;
            
              case 2:
                  localSubElemIdx.push_back(0);
                  break;
            
              case 4:
                  localSubElemIdx.push_back(1);
                  break;
            }
            break;

          case 2:  // r2: face-4 refinement, split into 2 subelements
            switch (faLocalIndex)
            {
              case 0:
              case 1:
              case 4:
                  localSubElemIdx.push_back(0);
                  localSubElemIdx.push_back(1);
                  break;
            
              case 2:
                  localSubElemIdx.push_back(1);
                  break;
            
              case 3:
                  localSubElemIdx.push_back(0);
                  break;
            }
            break;

          case 3: // r3: face-234 refinement, split into 2 subelements
            switch (faLocalIndex)
            {
              case 0:
                  localSubElemIdx.push_back(0);
                  break;
            
              case 1:
                  localSubElemIdx.push_back(1);
                  break;
            
              case 2:
              case 3:
              case 4:
                  localSubElemIdx.push_back(0);
                  localSubElemIdx.push_back(1);
                  break;
            }
            break;

          case 4: // r4: all-faces refinement, split into 4 subelements
            switch (faLocalIndex)
            {
              case 0:
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
                  localSubElemIdx.push_back(0);
                  localSubElemIdx.push_back(1);
                  break;
            
              case 4:
                  localSubElemIdx.push_back(1);
                  localSubElemIdx.push_back(2);
                  break;
            
            }
            break;

          case 5: // r5: split into 3 hexahedron elements
            switch (faLocalIndex)
            {
              case 0:
              case 1:
                  localSubElemIdx.push_back(0);
                  localSubElemIdx.push_back(1);
                  localSubElemIdx.push_back(2);
                  break;
            
              case 2:
                  localSubElemIdx.push_back(0);
                  localSubElemIdx.push_back(2);
                  break;
            
              case 3:
                  localSubElemIdx.push_back(0);
                  localSubElemIdx.push_back(1);
                  break;
            
              case 4:
                  localSubElemIdx.push_back(1);
                  localSubElemIdx.push_back(2);
                  break;
            
            }
            break;

          case 6: // r6: split into hexahedron and triangular-prism elements
            switch (faLocalIndex)
            {
              case 0:
              case 1:
              case 3:
              case 4:
                  localSubElemIdx.push_back(0);
                  localSubElemIdx.push_back(1);
                  break;
            
              case 2:
                  localSubElemIdx.push_back(0);
                  break;
            }
            break;
            
          case 7: // r7: split into hexahedron and triangular-prism elements
            switch (faLocalIndex)
            {
              case 0:
              case 1:
              case 2:
              case 4:
                  localSubElemIdx.push_back(0);
                  localSubElemIdx.push_back(1);
                  break;
            
              case 3:
                  localSubElemIdx.push_back(0);
                  break;
            }
            break;
        } // end switch refinement
    }
    
    unsigned int RefinementTriangularPrism::getLocalSubFaceNr(int refineType, unsigned int localFaceNr, unsigned int subElementIdx) const 
    { 
        if (refineType == 5)
        {
          // split into hexahedron, different local face numbering
          switch (localFaceNr)
          {
            case 0:
              return 0;
              
            case 1:
              return 5;
              
            case 2:
              switch (subElementIdx)
              {
                case 2:
                  return 3;

                case 0:
                  return 1;
              }
              break;
              
            case 3:
              switch (subElementIdx)
              {
                case 0:
                  return 3;
                  
                case 1:
                  return 1;
              }
              break;
              
            case 4:
              switch (subElementIdx)
              {
                case 1:
                  return 3;
                  
                case 2:
                  return 1;
              }
              break;
          } // end switch
        } // end if
        
        // the same as the parent's faces
        return localFaceNr; 
    }

}  // end of namespace Geometry
