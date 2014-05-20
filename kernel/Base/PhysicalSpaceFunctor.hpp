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
// File: PhysicalSpaceFunction.hh
// A wrapper for simple functions taking Point arguments from physical space.
// Lars Pesch,  26/03/2004
// New version: the component classes IteratorSource and ElementPtrSource had
// to be abolished because of trouble with the compilers relating to the
// template template args. Now there is a version that offers both
// functionality, but has to be started with an iterator.
//
// Extension 27/7/2004, LP: It has become necessary to extend the possible
// argument functions of PhysicalSpaceFunction to also include functors
// etc. Hence one template arg is formally deleted to generate a general case
// (which expects e.g. typedefs in the argument type) and the previously only
// version is given as a specialization of the general case.
// //------------------------------------------------------------------------------
#ifndef PHYSICALSPACEFUNCTION_HH
#define PHYSICALSPACEFUNCTION_HH

namespace Geometry
{
    class ElementGeometry;
    class PointPhysical;
    class PointReference;
}
namespace Base
{
    class Element;
}

namespace Base
{

	template<class returnType, typename FType>
	class PhysicalSpaceEvaluator<returnType,FType>;

    /*! \class PhysicalSpaceFunction
     * \brief A function that can be evaluated in physical space.
     * \details
     * PhysicalSpaceFunction objects are produced e.g. by the
     * transform2RefElement function of Element. They enable to query a
     * function that is defined in physical space at reference space
     * coordinates. The necessary transformation is carried out by the mapping
     * of the element which has to be passed in through the
     * constructor. PhysicalSpaceFunction objects are integrable since they
     * offer the necessary typedefs and the evaluation operator.
     */
    template <typename FType>
    class PhysicalSpaceFunctor
    {

    public:
        typedef typename Integration::ReturnTrait1<FType>::ReturnType ReturnType;
        
        typedef Geometry::ElementGeometry                      ElementGeometryT;

        //! Ctor, type of the wrapped function is fixed by template argument.
        PhysicalSpaceFunctor(const ElementGeometryT* const element, const FType& functor) :
            element_(element),
            functor_(functor)
        {}
        
        PhysicalSpaceFunctor(const PhysicalSpaceFunctor& other) :
            element_(other.element_),
            functor_(other.functor_)
        {}

        //! Evaluation operator for _reference_ space coordinates.
        void operator()(const Element& el, const Geometry::PointReference& pRef, ReturnType& r)
        {
            Geometry::PointPhysical pPhys(pRef.size());  // Declare and...
            el.referenceToPhysical(pRef, pPhys); // ...transform the point.
            // PhysSpaceEvaluator enables us to query different types of
            // functions/functors (regarding their argument composition)
            // with one syntax:
            PhysicalSpaceEvaluator<ReturnType, FType>::eval(el, functor_, pPhys, r);
        }

    private:
        //! The element we're on.
        const ElementGeometryT* const   element_;
        //! Physical space function/functor.
        FType                           functor_;

    };
}
#endif
