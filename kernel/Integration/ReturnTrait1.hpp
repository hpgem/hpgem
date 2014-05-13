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
// Package configuration (namespace Integration):
#ifndef RETURNRAIT1_HH
#define RETURNRAIT1_HH
#include "GlobalNamespaceIntegration.hpp"
//------------------------------------------------------------------------------
// System includes and names imported from them:
//------------------------------------------------------------------------------
// Package includes:
#include "../Geometry/PointReference.hpp"
#include "../Geometry/PointPhysical.hpp"

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
namespace Base
{
    class Element;
    
    class Face;
}

namespace Base
{
    template <class FType>
    class PhysicalSpaceFunction;

        // template <unsigned int DIM, class EType>
//    struct PhysGradientOfBasisFunction;
}


namespace Integration 
{
    using Base::Element;
    using Base::PhysicalSpaceFunction;
//      using Base::PhysGradientOfBasisFunction;
    
    // this is now the default.... it requires the T class to define ReturnType
    template <class T>
    struct ReturnTrait1
    {
        typedef typename T::ReturnType ReturnType;
    };

    // you can provide a function as a integrand
    template < typename T>
    struct ReturnTrait1<void (*)(const Base::Element*, const Geometry::PointReference&, T&)>
    {
        typedef T ReturnType;
    };

    template <class B, typename T>
    struct ReturnTrait1<void (B::*)(const Base::Element*, const Geometry::PointReference&, T&)>
    {
        typedef T ReturnType;
    };
    
    template <typename B, typename T>
    struct ReturnTrait1<void (B::*)(const Base::Face*, const Geometry::PointPhysical&, const Geometry::PointReference&, T&)>
    {
        typedef T ReturnType;
    };
    
    /*template < template<unsigned int> class B, typename T>
    struct ReturnTrait1<void (B<DIM>::*)(const Base::Element*, const Geometry::PointReference&, T&)>
    {
        typedef T ReturnType;
    };*/
        
    // you can provide a function as a integrand ant return calcuclated value via return value
    template < typename T>
    struct ReturnTrait1<T (*)(const Geometry::PointReference&)>
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
