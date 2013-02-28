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

#include "ConfigurationData.hpp"
#include "Geometry/PointReference.hpp"
#include "Geometry/PointPhysical.hpp"
#include "Integration/ReturnTrait1.hpp"
#include "Base/PhysicalSpaceEvaluator.hpp"

template <unsigned int DIM>
class Element;

namespace Base
{
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
    template <unsigned int DIM, class FType>
    class PhysicalSpaceFunctor
    {

    public:
        typedef typename Integration::ReturnTrait1<FType>::RetType RetType;
        
        typedef Element<DIM>                                        ElementT;

        //! Ctor, type of the wrapped function is fixed by template argument.
        PhysicalSpaceFunctor( ElementT* element, FType functor) :
            element_(element),
            functor_(functor)
        {}

        //! Evaluation operator for _reference_ space coordinates.
        void operator()(const Geometry::PointReference<DIM>& pRef, RetType& r)
        {
            Geometry::PointPhysical<DIM> pPhys;  // Declare and...
            element_->referenceToPhysical(pRef, pPhys); // ...transform the point.
            // PhysSpaceEvaluator enables us to query different types of
            // functions/functors (regarding their argument composition)
            // with one syntax:
            PhysicalSpaceEvaluator<DIM, RetType, FType>::eval(functor_, pPhys, r);
        }

    private:
        //! The element we're on.
        ElementT* element_;
        //! Physical space function/functor.
        FType functor_;

    };
}
#endif
