//------------------------------------------------------------------------------
// File: QuadratureRuleSet.hpp 
// Header for class QuadratureRuleSet to keep quadrature rules into a set.
// M.T. Julianto, Wed Feb 13 10:45:06 UTC 2013
//------------------------------------------------------------------------------

#ifndef QuadratureRuleSet_hpp
#define QuadratureRuleSet_hpp

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
