//------------------------------------------------------------------------------
// File: GaussQuadratureRule.hpp 
// The abstract base class for quadrature rules that calculate the integral
// as a weighted sum of function values at certain points.
// Lars Pesch, Thu Aug 18 17:44:27 CEST 2005
//----
// Modified by M.T. Julianto, Wed Feb 13 10:45:06 UTC 2013
//     -> GaussIntegrationRule and IntegrationRuleBase classes are merged into
//        this abstract base class.
//     -> Rule's id and rule criterion stuffs are removed.
//------------------------------------------------------------------------------
#ifndef GaussQuadratureRule_hpp
#define GaussQuadratureRule_hpp

#include <string>

#include "Geometry/PointReference.hpp"
#include "Geometry/ReferenceGeometry.hpp"
#include "Integration/GlobalNamespaceIntegration.hpp"

namespace Geometry
{
    // forward declaration
    class ReferenceGeometry;
}

namespace QuadratureRules
{
    class GaussQuadratureRule
    {
      public:
        //! Return the name of the quadrature.
        virtual std::string getName() const = 0;
        
        //! Return the dimension.
        virtual unsigned int dimension() const = 0;
        
        //! Return the order of the quadrature.
        virtual unsigned int order() const = 0;
        
        //! Return the number of points used in the quadrature.
        virtual unsigned int nrOfPoints() const = 0;

        //! Return the weight attached to the function value of the requested point number.
        virtual Integration::NumType weight(unsigned int) const = 0;

        //! Return the coordinates of the point with the given index.
        virtual void getPoint(unsigned int, Geometry::PointReference&) const = 0;

        //! Each rule also knows which ReferenceGeometry it is meant for.
        virtual Geometry::ReferenceGeometry* forReferenceGeometry() const = 0;
        
        virtual ~GaussQuadratureRule() {}
    };
} // close namespace Integration
#endif
