
//------------------------------------------------------------------------------
// Package configuration (namespace Integration):
#ifndef RETURNRAIT1_HH
#define RETURNRAIT1_HH
#include "GlobalNamespaceIntegration.hpp"
//------------------------------------------------------------------------------
// System includes and names imported from them:
//------------------------------------------------------------------------------
// Package includes:
#include "../Geometry/PointReference.hpp"

// #include "../Geometry/PhysSpacePoint.hh"
// using Geometry::PhysSpacePoint;
// #include "../Geometry/Element.hh"
// using Geometry::Element;

// #include "../Base/ElementDataExpVector.hh"
// using Base::ElementDataExpVector;
// #include "../Base/ElementDataExpansion.hh"
// using Base::ElementDataExpansion;
// #include "../Base/PhysicalSpaceFunction.hh"
// using Base::PhysicalSpaceFunction;
// #include "../Base/PhysGradientOfBasisFunction.hh"
// using Base::PhysGradientOfBasisFunction;
//------------------------------------------------------------------------------
namespace Geometry
{
    template <unsigned int DIM>
    class Element;
}

namespace Base
{
    template <unsigned int DIM, class FType>
    class PhysicalSpaceFunction;

    template <unsigned int DIM, class EType>
    struct PhysGradientOfBasisFunction;
}


namespace Integration 
{
    using Geometry::Element;
    using Base::PhysicalSpaceFunction;
    using Base::PhysGradientOfBasisFunction;
    
    // this is now the default.... it requires the T class to define ReturnType
    template <class T>
    struct ReturnTrait1
    {
        typedef typename T::ReturnType ReturnType;
    };

        // you can provide a function as a integrand
    template <unsigned int DIM, typename T>
    struct ReturnTrait1<void (*)(const Geometry::PointReference<DIM>&, T&)>
    {
        typedef T ReturnType;
    };
        // you can provide a function as a integrand ant return calcuclated value via return value
    template <unsigned int DIM, typename T>
    struct ReturnTrait1<T (*)(const Geometry::PointReference<DIM>&)>
    {
        typedef T ReturnType;
    };

//     template <unsigned int DIM>
//     struct ReturnTrait1<PhysGradientOfBasisFunction<DIM, Element<DIM> > >
//     {
//         typedef typename PhysGradientOfBasisFunction<DIM, Element<DIM> >::ReturnType ReturnType;
//     };

    template <class T>
    struct ReturnTrait1<T*>
    {
        typedef typename ReturnTrait1<T>::ReturnType ReturnType;
    };

} // close namespace Integration
#endif
//------------------------------------------------------------------------------
// Local variables:
// mode:c++
// comment-column: 48
// End:

