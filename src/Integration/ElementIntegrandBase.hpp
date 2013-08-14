//------------------------------------------------------------------------------
// File: ElementIntegrandBase.hpp 
// Base class for integrand.
// M.T. Julianto, Wed Feb 13 17:12:33 WET 2013
//------------------------------------------------------------------------------

#ifndef ElementIntegrandBase_hpp
#define ElementIntegrandBase_hpp

#include "Base/Element.hpp"
#include "Geometry/PointReference.hpp"
#include "Integration/GlobalNamespaceIntegration.hpp"

namespace Integration
{

    template <unsigned int DIM, typename T=LinearAlgebra::NumericalVector>
    class ElementIntegrandBase
    {
    public:
        typedef T                               ReturnType;
        typedef Geometry::PointReference<DIM>   PointReferenceT;
        typedef Base::Element<DIM>              ElementT;
        
    public:
        ~ElementIntegrandBase() {}

        virtual void operator()(const Base::Element<DIM>* element, const PointReferenceT& p, LinearAlgebra::NumericalVector& ret) = 0;
    };
};

#endif
