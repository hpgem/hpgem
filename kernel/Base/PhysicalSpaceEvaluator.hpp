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
#ifndef _PhysicalSpaceEvaluator_hpp
#define _PhysicalSpaceEvaluator_hpp

namespace Geometry{
	class PointPhysical;
}

namespace Base
{
    class Element;
}
    //! Default case: operator() with arguments Point and result type exists.
template <typename ResultType, class FuncType>
struct PhysicalSpaceEvaluator
{
	using RetType = ResultType;
    
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
	using RetType = ResultType;
    
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
	using RetType = ResultType;
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
	using RetType = ResultType;
	static void eval(
                     void (*funcPtr)(const Base::Element&, const Geometry::PointPhysical&, ResultType&),
                     const Geometry::PointPhysical& p,
                     RetType& r)    
	{
	    funcPtr(p, r);
	}
};
#endif
