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
#ifndef AllGaussQuadratureRules_hpp
#define AllGaussQuadratureRules_hpp

#include "Integration/QuadratureRules/GaussQuadratureRulesForLine.hpp"
#include "Integration/QuadratureRules/GaussQuadratureRulesForTriangle.hpp"
#include "Integration/QuadratureRules/GaussQuadratureRulesForSquare.hpp"
#include "Integration/QuadratureRules/GaussQuadratureRulesForTetrahedron.hpp"
#include "Integration/QuadratureRules/GaussQuadratureRulesForPyramid.hpp"
#include "Integration/QuadratureRules/GaussQuadratureRulesForTriangularPrism.hpp"
#include "Integration/QuadratureRules/GaussQuadratureRulesForCube.hpp"
#include "Integration/QuadratureRules/GaussQuadratureRulesForHypercube.hpp"
#include "Integration/QuadratureRules/GaussQuadratureRulesForPoint.hpp"

namespace QuadratureRules
{
	//solves the problem the singletons are causing by adding more singletons ;)
	/**
	 * Storage class for all the quadrature rules. If you add a rule, make sure to also add it here in the constructor.
	 * If you are integrating and want a quadrature rule, this is the appropriate place to get one
	 */
	class AllGaussQuadratureRules{
	public:
		static AllGaussQuadratureRules& instance();

		//it is possible to call this from an external location, but it is nicer to list all the rules inside this class
		void addRule(const GaussQuadratureRule* rule);

		const GaussQuadratureRule* getRule(const Geometry::ReferenceGeometry* referenceGeometry, int order);

	private:
		AllGaussQuadratureRules();

		AllGaussQuadratureRules(AllGaussQuadratureRules&);//this will generate a linker error if you try to copy
		void operator=(AllGaussQuadratureRules&);

		std::map<const Geometry::ReferenceGeometry*,std::list<const GaussQuadratureRule*>> listOfRules_;
	};

}
//---------------------------------------------------------------------------
#endif

