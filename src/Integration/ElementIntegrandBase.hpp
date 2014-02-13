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

	/**
	 * If you want to integrate over elements it is likely you already have the functions
	 * elementIntegrand(const Base::Element*, const Geometry::PointReference&, LinearAlgebra::Matrix) and
	 * elementIntegrand(const Base::Element*, const Geometry::PointReference&, LinearAlgebra::NumericalVector)
	 * implemented in some class already, so that class can simply inherit from ElementIntegrandBase<LinearAlgebra::Matrix>
	 * and ElementIntegrandBase<LinearAlgebra::NumericalVector> to signal the integrators that it does so
	 */
    template <class T>
    class ElementIntegrandBase
    {
    public:
    	///compute the contribution to the returntype of this reference point
        virtual void elementIntegrand(const Base::Element* element, const Geometry::PointReference& p, T& ret) = 0;
    };
};

#endif
