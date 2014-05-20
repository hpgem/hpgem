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

#include "BasisFunctions1DH1ConformingLine.hpp"
#include "helperFunctions.hpp"
#include "Base/BasisFunctionSet.hpp"
#include "Geometry/PointReference.hpp"

namespace Utilities {

	double BasisFunction1DVertexLine::eval(const Geometry::PointReference& p) const {
		return (1. + nodePosition_ * p[0]) / 2.;
	}

	double BasisFunction1DVertexLine::evalDeriv0(const Geometry::PointReference& p) const {
		return nodePosition_ / 2.;
	}

	double BasisFunction1DInteriorLine::eval(const Geometry::PointReference& p) const {
		return (1 + p[0]) * (1 - p[0]) * LobattoPolynomial(polynomialOrder_, p[0]) / 4.;
	}

	double BasisFunction1DInteriorLine::evalDeriv0(const Geometry::PointReference& p) const {
		return -p[0] * LobattoPolynomial(polynomialOrder_, p[0]) / 2. + (1 + p[0]) * (1 - p[0])
						* LobattoPolynomialDerivative(polynomialOrder_, p[0]) / 4.;
	}

	Base::BasisFunctionSet* createDGBasisFunctionSet1DH1Line(int polynomialOrder) {
		Base::BasisFunctionSet* result(new Base::BasisFunctionSet(polynomialOrder));
		result->addBasisFunction(new BasisFunction1DVertexLine(0));
		result->addBasisFunction(new BasisFunction1DVertexLine(1));
		for (int i = 0; i <= polynomialOrder - 2; ++i) {
			result->addBasisFunction(new BasisFunction1DInteriorLine(i));
		}
		return result;
	}

	Base::BasisFunctionSet* createInteriorBasisFunctionSet1DH1Line(int polynomialOrder) {
		Base::BasisFunctionSet* result(new Base::BasisFunctionSet(polynomialOrder));
		for (int i = 0; i <= polynomialOrder - 2; ++i) {
			result->addBasisFunction(new BasisFunction1DInteriorLine(i));
		}
		return result;
	}

	void createVertexBasisFunctionSet1DH1Line(int polynomialOrder,std::vector<const Base::BasisFunctionSet*>& result) {
		Base::BasisFunctionSet* set(new Base::BasisFunctionSet(polynomialOrder));
		set->addBasisFunction(new BasisFunction1DVertexLine(0));
		result.push_back(set);
		set = new Base::BasisFunctionSet(polynomialOrder);
		set->addBasisFunction(new BasisFunction1DVertexLine(1));
		result.push_back(set);
	}

}
