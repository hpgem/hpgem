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
    class Element;
}
    //! Default case: operator() with arguments Point and result type exists.
template <typename ResultType, class FuncType>
struct PhysicalSpaceEvaluator
{
	typedef ResultType RetType;
    
	static void eval(const Base::Element& el, FuncType& f, const Geometry::PointPhysical& p, RetType& r)
    {
		f(el, p, r);
    }
};

    //! Specialization for pointer (occurs when e.g. ic's are given as pointers
    //  for flexibility).
template <typename ResultType, class FuncType>
struct PhysicalSpaceEvaluator<ResultType, FuncType*>
{
	typedef ResultType RetType;
    
	static void eval(FuncType* f, const Geometry::PointPhysical& p, RetType& r)
    {
		(*f)(p, r);
    }
};

    //! Specialization for a function (like basis functions).
template <typename ResultType>
struct PhysicalSpaceEvaluator<ResultType,
ResultType (*)(const Geometry::PointPhysical&)>
{
	typedef ResultType RetType;
	static void eval(RetType (*funcPtr)(const Geometry::PointPhysical&),
                     const Geometry::PointPhysical& p,
                     RetType& r)
	{
	    r = funcPtr(p);
	}
};

    //! Specialization for a function (which returns void).
template <typename ResultType>
struct PhysicalSpaceEvaluator<ResultType,
void (*)(const Geometry::PointPhysical&, ResultType&)>
{
	typedef ResultType RetType;
	static void eval(
                     void (*funcPtr)(const Base::Element&, const Geometry::PointPhysical&, ResultType&),
                     const Geometry::PointPhysical& p,
                     RetType& r)    
	{
	    funcPtr(p, r);
	}
};
#endif
