/*
 This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
 http://hpgem.org/about-the-code/team
 
 This code is distributed using BSD 3-Clause License. A copy of which can found below.
 
 
 Copyright (c) 2014, Univesity of Twenete
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef BASISFUNCTIONS2DH1CONFORMINGTRIANGLE_HPP_
#define BASISFUNCTIONS2DH1CONFORMINGTRIANGLE_HPP_

#include "Base/BaseBasisFunction.hpp"
#include <vector>

namespace Base{
class BasisFunctionSet;
class OrientedBasisFunctionSet;
}

namespace Geometry{
class PointReference;
}

namespace Utilities{

	class BasisFunction2DVertexTriangle: public Base::BaseBasisFunction
	{
	public:
		BasisFunction2DVertexTriangle(int node):node_(node){}

	double eval(const Geometry::PointReference& p) const;

	double evalDeriv0(const Geometry::PointReference& p) const;

	double evalDeriv1(const Geometry::PointReference& p) const;

	private:
		int node_;
	};

	class BasisFunction2DFaceTriangle: public Base::BaseBasisFunction
	{
	public:
		BasisFunction2DFaceTriangle(int node0,int node1,int polynomialOrder):node0_(node0),node1_(node1),polynomialOrder_(polynomialOrder){}

	double eval(const Geometry::PointReference& p) const;

	double evalDeriv0(const Geometry::PointReference& p) const;

	double evalDeriv1(const Geometry::PointReference& p) const;

	private:
		int node0_, node1_,polynomialOrder_;
	};

	class BasisFunction2DInteriorTriangle: public Base::BaseBasisFunction
	{
	public:
		BasisFunction2DInteriorTriangle(int polynomialOrder0,int polynomialOrder1):polynomialOrder0_(polynomialOrder0),polynomialOrder1_(polynomialOrder1){}

	double eval(const Geometry::PointReference& p) const;

	double evalDeriv0(const Geometry::PointReference& p) const;

	double evalDeriv1(const Geometry::PointReference& p) const;

	private:
		int polynomialOrder0_,polynomialOrder1_;
	};

Base::BasisFunctionSet* createDGBasisFunctionSet2DH1Triangle(int order);

Base::BasisFunctionSet* createInteriorBasisFunctionSet2DH1Triangle(int order);

void createVertexBasisFunctionSet2DH1Triangle(int order, std::vector<const Base::BasisFunctionSet*>& result);

void createFaceBasisFunctionSet2DH1Triangle(int order, std::vector<const Base::OrientedBasisFunctionSet*>& result);

}




#endif /* BASISFUNCTIONS2DH1CONFORMINGTRIANGLE_HPP_ */
