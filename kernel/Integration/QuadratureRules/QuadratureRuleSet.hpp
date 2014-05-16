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

#ifndef QuadratureRuleSet_hpp
#define QuadratureRuleSet_hpp

///\TODO (flagged for deletion) this appears to be a partially finished solution for the quadrature rule problem that is recently resolved. There is no interaction with the rest of the kernel -FB

#include <list>
#include <iostream>

namespace QuadratureRules
{
    class QuadratureRuleSet
    {
      public:
          QuadratureRuleSet() {}
          ~QuadratureRuleSet() {}

          unsigned int size(Geometry::ReferenceLine* refGeo) const;
          unsigned int size(Geometry::ReferenceTriangle* refGeo) const;
          unsigned int size(Geometry::ReferenceSquare* refGeo) const;
          unsigned int size(Geometry::ReferenceTetrahedron* refGeo) const;
          unsigned int size(Geometry::ReferenceTriangularPrism* refGeo) const;
          unsigned int size(Geometry::ReferencePyramid* refGeo) const;
          unsigned int size(Geometry::ReferenceCube* refGeo) const;
          unsigned int size(Geometry::ReferenceHypercube* refGeo) const;
          unsigned int size() const;

          void AddRule(GaussQuadratureRule* qr);

          GaussQuadratureRule* GetRule(Geometry::ReferenceLine* refGeo, int order);
          GaussQuadratureRule* GetRule(Geometry::ReferenceTriangle* refGeo, int order);
          GaussQuadratureRule* GetRule(Geometry::ReferenceSquare* refGeo, int order);
          GaussQuadratureRule* GetRule(Geometry::ReferenceTetrahedron* refGeo, int order);
          GaussQuadratureRule* GetRule(Geometry::ReferencePyramid* refGeo, int order);
          GaussQuadratureRule* GetRule(Geometry::ReferenceTriangularPrism* refGeo, int order);
          GaussQuadratureRule* GetRule(Geometry::ReferenceCube* refGeo, int order);
          GaussQuadratureRule* GetRule(Geometry::ReferenceHypercube* refGeo, int order);
        
      private:
          typedef std::list<GaussQuadratureRule*> listOfQR_1DType;
          typedef std::list<GaussQuadratureRule*> listOfQR_2DType;
          typedef std::list<GaussQuadratureRule*> listOfQR_3DType;
          typedef std::list<GaussQuadratureRule*> listOfQR_4DType;

          listOfQR_1DType listOfQR_Line;
          listOfQR_2DType listOfQR_Triangle;
          listOfQR_2DType listOfQR_Square;
          listOfQR_3DType listOfQR_Tetrahedron;
          listOfQR_3DType listOfQR_TriangularPrism;
          listOfQR_3DType listOfQR_Pyramid;
          listOfQR_3DType listOfQR_Cube;
          listOfQR_4DType listOfQR_Hypercube;
    };
    
};

#endif
