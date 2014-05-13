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
//------------------------------------------------------------------------------

#include <list>
#include <iostream>
using std::cerr;
#include "Integration/QuadratureRules/GaussQuadratureRule.hpp"
#include "Geometry/ReferenceGeometry.hpp"
#include "Geometry/ReferenceLine.hpp"
#include "Geometry/ReferenceTriangle.hpp"
#include "Geometry/ReferenceSquare.hpp"
#include "Geometry/ReferenceTetrahedron.hpp"
#include "Geometry/ReferencePyramid.hpp"
#include "Geometry/ReferenceTriangularPrism.hpp"
#include "Geometry/ReferenceCube.hpp"
#include "Geometry/ReferenceHypercube.hpp"

#include "QuadratureRuleSet.hpp"

namespace QuadratureRules
{
    unsigned int QuadratureRuleSet::size(Geometry::ReferenceLine* refGeo) const
    {
        return listOfQR_Line.size();
    }
    
    unsigned int QuadratureRuleSet::size(Geometry::ReferenceTriangle* refGeo) const
    {
        return listOfQR_Triangle.size();
    }
    
    unsigned int QuadratureRuleSet::size(Geometry::ReferenceSquare* refGeo) const
    {
        return listOfQR_Square.size();
    }

    unsigned int QuadratureRuleSet::size(Geometry::ReferenceTetrahedron* refGeo) const
    {
        return listOfQR_Tetrahedron.size();
    }
    
    unsigned int QuadratureRuleSet::size(Geometry::ReferenceTriangularPrism* refGeo) const
    {
        return listOfQR_TriangularPrism.size();
    }
    
    unsigned int QuadratureRuleSet::size(Geometry::ReferencePyramid* refGeo) const
    {
        return listOfQR_Pyramid.size();
    }
    
    unsigned int QuadratureRuleSet::size(Geometry::ReferenceCube* refGeo) const
    {
        return listOfQR_Cube.size();
    }
    
    unsigned int QuadratureRuleSet::size(Geometry::ReferenceHypercube* refGeo) const
    {
        return listOfQR_Hypercube.size();
    }

    unsigned int QuadratureRuleSet::size() const
    {
        return        listOfQR_Line.size() + listOfQR_Triangle.size() + listOfQR_Square.size() + 
                listOfQR_Tetrahedron.size() + listOfQR_TriangularPrism.size() + 
                listOfQR_Pyramid.size() + listOfQR_Cube.size() + listOfQR_Hypercube.size();
    }
    
    void QuadratureRuleSet::AddRule(GaussQuadratureRule* qr)
    {
    	if(qr->forReferenceGeometry()==&Geometry::ReferenceLine::Instance()){
			listOfQR_1DType::iterator it = listOfQR_Line.begin();
			while (it != listOfQR_Line.end())
			{
			  if (((GaussQuadratureRule*)*it)->order() < qr->order()) ++it;
			  else break;
			}
			listOfQR_Line.insert(it,qr);
    	}else if(qr->forReferenceGeometry()==&Geometry::ReferenceTriangle::Instance()){
			listOfQR_2DType::iterator it = listOfQR_Triangle.begin();
			while (it != listOfQR_Triangle.end())
			{
			if (((GaussQuadratureRule*)*it)->order() < qr->order()) ++it;
			else break;
			}
			listOfQR_Triangle.insert(it,qr);
    	}else if(qr->forReferenceGeometry()== &Geometry::ReferenceSquare::Instance()){
			listOfQR_2DType::iterator it = listOfQR_Square.begin();
			while (it != listOfQR_Square.end())
			{
			if (((GaussQuadratureRule*)*it)->order() < qr->order()) ++it;
			else break;
			}
			listOfQR_Square.insert(it,qr);
    	}else if(qr->forReferenceGeometry()== &Geometry::ReferenceTetrahedron::Instance()){
            listOfQR_3DType::iterator it = listOfQR_Tetrahedron.begin();
			while (it != listOfQR_Tetrahedron.end())
			{
			if (((GaussQuadratureRule*)*it)->order() < qr->order()) ++it;
			else break;
			}
			listOfQR_Tetrahedron.insert(it,qr);
		}else if(qr->forReferenceGeometry()== &Geometry::ReferencePyramid::Instance()){
			listOfQR_3DType::iterator it = listOfQR_Pyramid.begin();
			while (it != listOfQR_Pyramid.end())
			{
			if (((GaussQuadratureRule*)*it)->order() < qr->order()) ++it;
			else break;
			}
			listOfQR_Pyramid.insert(it,qr);
		}else if(qr->forReferenceGeometry()== &Geometry::ReferenceTriangularPrism::Instance()){
			listOfQR_3DType::iterator it = listOfQR_TriangularPrism.begin();
			while (it != listOfQR_TriangularPrism.end())
			{
			if (((GaussQuadratureRule*)*it)->order() < qr->order()) ++it;
			else break;
			}
			listOfQR_TriangularPrism.insert(it,qr);
		}else if(qr->forReferenceGeometry()== &Geometry::ReferenceCube::Instance()){
			listOfQR_3DType::iterator it = listOfQR_Cube.begin();
			while (it != listOfQR_Cube.end())
			{
			if (((GaussQuadratureRule*)*it)->order() < qr->order()) ++it;
			else break;
			}
			listOfQR_Cube.insert(it,qr);
		}else if(qr->forReferenceGeometry()== &Geometry::ReferenceHypercube::Instance()){
			listOfQR_4DType::iterator it = listOfQR_Hypercube.begin();
			while (it != listOfQR_Hypercube.end())
			{
			  if (((GaussQuadratureRule*)*it)->order() < qr->order()) ++it;
			  else break;
			}
			listOfQR_Hypercube.insert(it,qr);
			}else throw "QuadratureRuleSet::AddRule: unknown reference geometry!";
    }

    GaussQuadratureRule* QuadratureRuleSet::GetRule(Geometry::ReferenceLine* refGeo, int order)
    {
        for (listOfQR_1DType::iterator it = listOfQR_Line.begin(); 
              it != listOfQR_Line.end(); ++it)
          if (((GaussQuadratureRule*)*it)->order() >= order) return *it;

        return NULL;
    }
    
    GaussQuadratureRule* QuadratureRuleSet::GetRule(Geometry::ReferenceTriangle* refGeo, int order)
    {
        for (listOfQR_2DType::iterator it = listOfQR_Triangle.begin(); 
              it != listOfQR_Triangle.end(); ++it)
          if (((GaussQuadratureRule*)*it)->order() >= order) return *it;

        return NULL;
    }

    GaussQuadratureRule* QuadratureRuleSet::GetRule(Geometry::ReferenceSquare* refGeo, int order)
    {
        for (listOfQR_2DType::iterator it = listOfQR_Square.begin();
              it != listOfQR_Square.end(); ++it)
          if (((GaussQuadratureRule*)*it)->order() >= order) return *it;

        return NULL;
    }


    GaussQuadratureRule* QuadratureRuleSet::GetRule(Geometry::ReferenceTetrahedron* refGeo, int order)
    {
        for (listOfQR_3DType::iterator it = listOfQR_Tetrahedron.begin();
              it != listOfQR_Tetrahedron.end(); ++it)
          if (((GaussQuadratureRule*)*it)->order() >= order) return *it;

        return NULL;
    }
    
    GaussQuadratureRule* QuadratureRuleSet::GetRule(Geometry::ReferencePyramid* refGeo, int order)
    {
        for (listOfQR_3DType::iterator it = listOfQR_Pyramid.begin();
              it != listOfQR_Pyramid.end(); ++it)
          if (((GaussQuadratureRule*)*it)->order() >= order) return *it;

        return NULL;
    }
    
    GaussQuadratureRule* QuadratureRuleSet::GetRule(Geometry::ReferenceTriangularPrism* refGeo, int order)
    {
        for (listOfQR_3DType::iterator it = listOfQR_TriangularPrism.begin();
              it != listOfQR_TriangularPrism.end(); ++it)
          if (((GaussQuadratureRule*)*it)->order() >= order) return *it;

        return NULL;
    }
    
    GaussQuadratureRule* QuadratureRuleSet::GetRule(Geometry::ReferenceCube* refGeo, int order)
    {
        for (listOfQR_3DType::iterator it = listOfQR_Cube.begin();
              it != listOfQR_Cube.end(); ++it)
          if (((GaussQuadratureRule*)*it)->order() >= order) return *it;

        return NULL;
    }
    
    GaussQuadratureRule* QuadratureRuleSet::GetRule(Geometry::ReferenceHypercube* refGeo, int order)
    {
        for (listOfQR_4DType::iterator it = listOfQR_Hypercube.begin();
              it != listOfQR_Hypercube.end(); ++it)
          if (((GaussQuadratureRule*)*it)->order() >= order) return *it;

        return NULL;
    }

};
