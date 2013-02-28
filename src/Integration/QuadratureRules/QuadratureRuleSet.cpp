//------------------------------------------------------------------------------
// File: QuadratureRuleSet.hpp 
// Implementation for class QuadratureRuleSet to keep quadrature rules into a set.
// M.T. Julianto, Wed Feb 13 10:45:06 UTC 2013
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
    
    void QuadratureRuleSet::AddRule(GaussQuadratureRule<1>* qr)
    {
        listOfQR_1DType::iterator it = listOfQR_Line.begin();
        while (it != listOfQR_Line.end())
        {
          if (((GaussQuadratureRule<1>*)*it)->order() < qr->order()) ++it;
          else break;
        }
        listOfQR_Line.insert(it,qr);
    }
    
    void QuadratureRuleSet::AddRule(GaussQuadratureRule<2>* qr)
    {
        if (qr->forReferenceGeometry() == &Geometry::ReferenceTriangle::Instance())
        {
          listOfQR_2DType::iterator it = listOfQR_Triangle.begin();
          while (it != listOfQR_Triangle.end())
          {
            if (((GaussQuadratureRule<2>*)*it)->order() < qr->order()) ++it;
            else break;
          }
          listOfQR_Triangle.insert(it,qr);
        }
        else
        {
          listOfQR_2DType::iterator it = listOfQR_Square.begin();
          while (it != listOfQR_Square.end())
          {
            if (((GaussQuadratureRule<2>*)*it)->order() < qr->order()) ++it;
            else break;
          }
          listOfQR_Square.insert(it,qr);
        }
    }

    void QuadratureRuleSet::AddRule(GaussQuadratureRule<3>* qr)
    {
        if (qr->forReferenceGeometry() == &Geometry::ReferenceTetrahedron::Instance())
        {
          listOfQR_3DType::iterator it = listOfQR_Tetrahedron.begin();
          while (it != listOfQR_Tetrahedron.end())
          {
            if (((GaussQuadratureRule<3>*)*it)->order() < qr->order()) ++it;
            else break;
          }
          listOfQR_Tetrahedron.insert(it,qr);
        }
        else if (qr->forReferenceGeometry() == &Geometry::ReferencePyramid::Instance())
        {
          listOfQR_3DType::iterator it = listOfQR_Pyramid.begin();
          while (it != listOfQR_Pyramid.end())
          {
            if (((GaussQuadratureRule<3>*)*it)->order() < qr->order()) ++it;
            else break;
          }
          listOfQR_Pyramid.insert(it,qr);
        }
        else if (qr->forReferenceGeometry() == &Geometry::ReferenceTriangularPrism::Instance())
        {
          listOfQR_3DType::iterator it = listOfQR_TriangularPrism.begin();
          while (it != listOfQR_TriangularPrism.end())
          {
            if (((GaussQuadratureRule<3>*)*it)->order() < qr->order()) ++it;
            else break;
          }
          listOfQR_TriangularPrism.insert(it,qr);
        }
        else
        {
          listOfQR_3DType::iterator it = listOfQR_Cube.begin();
          while (it != listOfQR_Cube.end())
          {
            if (((GaussQuadratureRule<3>*)*it)->order() < qr->order()) ++it;
            else break;
          }
          listOfQR_Cube.insert(it,qr);
        }
    }
    
    void QuadratureRuleSet::AddRule(GaussQuadratureRule<4>* qr)
    {
        listOfQR_4DType::iterator it = listOfQR_Hypercube.begin();
        while (it != listOfQR_Hypercube.end())
        {
          if (((GaussQuadratureRule<4>*)*it)->order() < qr->order()) ++it;
          else break;
        }
        listOfQR_Hypercube.insert(it,qr);
    }

    GaussQuadratureRule<1>* QuadratureRuleSet::GetRule(Geometry::ReferenceLine* refGeo, int order)
    {
        for (listOfQR_1DType::iterator it = listOfQR_Line.begin(); 
              it != listOfQR_Line.end(); ++it)
          if (((GaussQuadratureRule<1>*)*it)->order() >= order) return *it;

        return NULL;
    }
    
    GaussQuadratureRule<2>* QuadratureRuleSet::GetRule(Geometry::ReferenceTriangle* refGeo, int order)
    {
        for (listOfQR_2DType::iterator it = listOfQR_Triangle.begin(); 
              it != listOfQR_Triangle.end(); ++it)
          if (((GaussQuadratureRule<1>*)*it)->order() >= order) return *it;

        return NULL;
    }

    GaussQuadratureRule<2>* QuadratureRuleSet::GetRule(Geometry::ReferenceSquare* refGeo, int order)
    {
        for (listOfQR_2DType::iterator it = listOfQR_Square.begin();
              it != listOfQR_Square.end(); ++it)
          if (((GaussQuadratureRule<1>*)*it)->order() >= order) return *it;

        return NULL;
    }


    GaussQuadratureRule<3>* QuadratureRuleSet::GetRule(Geometry::ReferenceTetrahedron* refGeo, int order)
    {
        for (listOfQR_3DType::iterator it = listOfQR_Tetrahedron.begin();
              it != listOfQR_Tetrahedron.end(); ++it)
          if (((GaussQuadratureRule<1>*)*it)->order() >= order) return *it;

        return NULL;
    }
    
    GaussQuadratureRule<3>* QuadratureRuleSet::GetRule(Geometry::ReferencePyramid* refGeo, int order)
    {
        for (listOfQR_3DType::iterator it = listOfQR_Pyramid.begin();
              it != listOfQR_Pyramid.end(); ++it)
          if (((GaussQuadratureRule<1>*)*it)->order() >= order) return *it;

        return NULL;
    }
    
    GaussQuadratureRule<3>* QuadratureRuleSet::GetRule(Geometry::ReferenceTriangularPrism* refGeo, int order)
    {
        for (listOfQR_3DType::iterator it = listOfQR_TriangularPrism.begin();
              it != listOfQR_TriangularPrism.end(); ++it)
          if (((GaussQuadratureRule<1>*)*it)->order() >= order) return *it;

        return NULL;
    }
    
    GaussQuadratureRule<3>* QuadratureRuleSet::GetRule(Geometry::ReferenceCube* refGeo, int order)
    {
        for (listOfQR_3DType::iterator it = listOfQR_Cube.begin();
              it != listOfQR_Cube.end(); ++it)
          if (((GaussQuadratureRule<1>*)*it)->order() >= order) return *it;

        return NULL;
    }
    
    GaussQuadratureRule<4>* QuadratureRuleSet::GetRule(Geometry::ReferenceHypercube* refGeo, int order)
    {
        for (listOfQR_4DType::iterator it = listOfQR_Hypercube.begin();
              it != listOfQR_Hypercube.end(); ++it)
          if (((GaussQuadratureRule<1>*)*it)->order() >= order) return *it;

        return NULL;
    }

};
