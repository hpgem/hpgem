//
//  PhysicalSpaceEvaluator.hpp
//  
//
//  Created by Shavarsh Nurijanyan on 2/27/13.
//
//

#ifndef _PhysicalSpaceEvaluator_hpp
#define _PhysicalSpaceEvaluator_hpp

#include "Geometry/PointPhysical.hpp"

namespace Base
{
    template <unsigned int DIM>
    class Element;
}
    //! Default case: operator() with arguments Point and result type exists.
template <unsigned int DIM, typename ResultType, class FuncType>
struct PhysicalSpaceEvaluator
{
	typedef ResultType RetType;
    
	static void eval(const Base::Element<DIM>& el, FuncType& f, const Geometry::PointPhysical<DIM>& p, RetType& r)
    {
		f(el, p, r);
    }
};

    //! Specialization for pointer (occurs when e.g. ic's are given as pointers
    //  for flexibility).
template <unsigned int DIM, typename ResultType, class FuncType>
struct PhysicalSpaceEvaluator<DIM, ResultType, FuncType*>
{
	typedef ResultType RetType;
    
	static void eval(FuncType* f, const Geometry::PointPhysical<DIM>& p, RetType& r)
    {
		(*f)(p, r);
    }
};

    //! Specialization for a function (like basis functions).
template <unsigned int DIM, typename ResultType>
struct PhysicalSpaceEvaluator<DIM, ResultType,
ResultType (*)(const Geometry::PointPhysical<DIM>&)>
{
	typedef ResultType RetType;
	static void eval(RetType (*funcPtr)(const Geometry::PointPhysical<DIM>&),
                     const Geometry::PointPhysical<DIM>& p,
                     RetType& r)
	{
	    r = funcPtr(p);
	}
};

    //! Specialization for a function (which returns void).
template <unsigned int DIM, typename ResultType>
struct PhysicalSpaceEvaluator<DIM, ResultType,
void (*)(const Geometry::PointPhysical<DIM>&, ResultType&)>
{
	typedef ResultType RetType;
	static void eval(
                     void (*funcPtr)(const Base::Element<DIM>&, const Geometry::PointPhysical<DIM>&, ResultType&),
                     const Geometry::PointPhysical<DIM>& p,
                     RetType& r)    
	{
	    funcPtr(p, r);
	}
};
#endif
