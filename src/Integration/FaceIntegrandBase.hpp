//------------------------------------------------------------------------------
// File: FaceIntegrandBase.hpp 
// Base class for integrand.
// M.T. Julianto, Wed Feb 13 17:12:33 WET 2013
//------------------------------------------------------------------------------

#ifndef FaceIntegrandBase_hpp
#define FaceIntegrandBase_hpp

#include "Base/Element.hpp"
#include "Base/Face.hpp"
#include "Geometry/PointReference.hpp"
#include "Geometry/PointPhysical.hpp"
#include "Integration/GlobalNamespaceIntegration.hpp"

namespace Integration
{

    template < typename T>
    class FaceIntegrandBase
    {
    public:
        virtual void faceIntegrand(const Base::Face* face, const Geometry::PointPhysical& normal, const Geometry::PointReference& p, T& ret) = 0;
    };

};

#endif
