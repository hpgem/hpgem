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

    template <unsigned int DIM, typename T=LinearAlgebra::NumericalVector>
    class FaceIntegrandBase
    {
    public:
        typedef T                                   ReturnType;
        typedef Geometry::PointReference<DIM-1>     PointReferenceT;
        typedef Geometry::PointPhysical<DIM>        PointPhysicalT;
       
    public:
        ~FaceIntegrandBase() {}

        virtual void operator()(const Base::Face<DIM>& fa, const PointPhysicalT& normal, 
                                const PointReferenceT& p,  ReturnType& ret) = 0;
    };

};

#endif
