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

#ifndef BASISFUNCTIONS2DH1CONFORMINGSQUARE_HPP_
#define BASISFUNCTIONS2DH1CONFORMINGSQUARE_HPP_

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

	class BasisFunction2DVertexSquare: public Base::BaseBasisFunction
	{
	public:
		BasisFunction2DVertexSquare(int node):nodePosition0_((node%2)*2-1),nodePosition1_((node/2)*2-1){}

		double eval(const Geometry::PointReference& p) const;

		double evalDeriv0(const Geometry::PointReference& p) const;

		double evalDeriv1(const Geometry::PointReference& p) const;

	private:
		int nodePosition0_;
		int nodePosition1_;
	};

	class BasisFunction2DFaceSquare_0: public Base::BaseBasisFunction
	{
	public:
		BasisFunction2DFaceSquare_0(int node0, int node1, int polynomialOrder);

		double eval(const Geometry::PointReference& p) const;

		double evalDeriv0(const Geometry::PointReference& p) const;

		double evalDeriv1(const Geometry::PointReference& p) const;

	private:
		int edgePosition_;
		int mirroring_;
		int polynomialOrder_;
	};

	class BasisFunction2DFaceSquare_1: public Base::BaseBasisFunction
	{
	public:
		BasisFunction2DFaceSquare_1(int node0, int node1, int polynomialOrder);

		double eval(const Geometry::PointReference& p) const;

		double evalDeriv0(const Geometry::PointReference& p) const;

		double evalDeriv1(const Geometry::PointReference& p) const;

	private:
		int edgePosition_;
		int mirroring_;
		int polynomialOrder_;
	};

	class BasisFunction2DInteriorSquare: public Base::BaseBasisFunction
	{
	public:
		BasisFunction2DInteriorSquare(int polynomialOrder0,int polynomialOrder1):polynomialOrder0_(polynomialOrder0),polynomialOrder1_(polynomialOrder1){}

		double eval(const Geometry::PointReference& p) const;

		double evalDeriv0(const Geometry::PointReference& p) const;

		double evalDeriv1(const Geometry::PointReference& p) const;
	private:
		int polynomialOrder0_,polynomialOrder1_;
	};

Base::BasisFunctionSet* createDGBasisFunctionSet2DH1Square(int order);

Base::BasisFunctionSet* createInteriorBasisFunctionSet2DH1Square(int order);

std::vector<const Base::BasisFunctionSet*> createVertexBasisFunctionSet2DH1Square(int order);

std::vector<const Base::OrientedBasisFunctionSet*> createFaceBasisFunctionSet2DH1Square(int order);

}



#endif /* BASISFUNCTIONS2DH1CONFORMINGSQUARE_HPP_ */
