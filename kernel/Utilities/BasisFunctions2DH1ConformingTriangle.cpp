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

#include "BasisFunctions2DH1ConformingTriangle.hpp"
#include "helperFunctions.hpp"
#include "Base/BasisFunctionSet.hpp"
#include "Base/OrientedBasisFunctionSet.hpp"
#include "Geometry/ReferenceTriangle.hpp"
#include "Geometry/PointReference.hpp"

namespace Utilities {

	double BasisFunction2DVertexTriangle::evalDeriv0(const Geometry::PointReference& p) const {
		return baricentricDeriv(node_, 0);
	}

	double BasisFunction2DInteriorTriangle::evalDeriv0(const Geometry::PointReference& p) const {
		double x0 = baricentric_2D(0, p) - baricentric_2D(1, p);
		double x1 = baricentric_2D(1, p) - baricentric_2D(2, p);
		return baricentricDeriv(0, 0) * (baricentric_2D(1, p) * baricentric_2D(2, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1)
						+ baricentric_2D(0, p) * baricentric_2D(1, p) * baricentric_2D(2, p) * LobattoPolynomialDerivative(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1))
				+ baricentricDeriv(1, 0) * (baricentric_2D(0, p) * baricentric_2D(2, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1)
						- baricentric_2D(0, p) * baricentric_2D(1, p) * baricentric_2D(2, p) * LobattoPolynomialDerivative(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1)
						+ baricentric_2D(0, p) * baricentric_2D(1, p) * baricentric_2D(2, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomialDerivative(polynomialOrder1_, x1));
	}

	double BasisFunction2DVertexTriangle::eval(const Geometry::PointReference& p) const {
		return baricentric_2D(node_, p);
	}

	double BasisFunction2DInteriorTriangle::eval(const Geometry::PointReference& p) const {
		double x0 = baricentric_2D(0, p) - baricentric_2D(1, p);
		double x1 = baricentric_2D(1, p) - baricentric_2D(2, p);
		return baricentric_2D(0, p) * baricentric_2D(1, p) * baricentric_2D(2, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1);
	}

	double BasisFunction2DVertexTriangle::evalDeriv1(const Geometry::PointReference& p) const {
		return baricentricDeriv(node_, 1);
	}

	double BasisFunction2DFaceTriangle::evalDeriv1(const Geometry::PointReference& p) const {
		return baricentricDeriv(node0_, 1) * (baricentric_2D(node1_, p) * LobattoPolynomial(polynomialOrder_, baricentric_2D(node1_, p) - baricentric_2D(node0_, p))
						- baricentric_2D(node0_, p) * baricentric_2D(node1_, p) * LobattoPolynomialDerivative(polynomialOrder_, baricentric_2D(node1_, p) - baricentric_2D(node0_, p)))
				+ baricentricDeriv(node1_, 1) * (baricentric_2D(node0_, p) * LobattoPolynomial(polynomialOrder_, baricentric_2D(node1_, p) - baricentric_2D(node0_, p))
						+ baricentric_2D(node1_, p) * baricentric_2D(node0_, p) * LobattoPolynomialDerivative(polynomialOrder_, baricentric_2D(node1_, p) - baricentric_2D(node0_, p)));
	}

	double BasisFunction2DFaceTriangle::evalDeriv0(const Geometry::PointReference& p) const {
		return baricentricDeriv(node0_, 0) * (baricentric_2D(node1_, p) * LobattoPolynomial(polynomialOrder_,baricentric_2D(node1_, p) - baricentric_2D(node0_, p))
						- baricentric_2D(node0_, p) * baricentric_2D(node1_, p) * LobattoPolynomialDerivative(polynomialOrder_, baricentric_2D(node1_, p) - baricentric_2D(node0_, p)))
				+ baricentricDeriv(node1_, 0) * (baricentric_2D(node0_, p) * LobattoPolynomial(polynomialOrder_, baricentric_2D(node1_, p) - baricentric_2D(node0_, p))
						+ baricentric_2D(node1_, p) * baricentric_2D(node0_, p) * LobattoPolynomialDerivative(polynomialOrder_,baricentric_2D(node1_, p) - baricentric_2D(node0_, p)));
	}

	double BasisFunction2DFaceTriangle::eval(const Geometry::PointReference& p) const {
		return baricentric_2D(node0_, p) * baricentric_2D(node1_, p) * LobattoPolynomial(polynomialOrder_, baricentric_2D(node1_, p) - baricentric_2D(node0_, p));
	}

	double BasisFunction2DInteriorTriangle::evalDeriv1(const Geometry::PointReference& p) const {
		double x0 = baricentric_2D(0, p) - baricentric_2D(1, p);
		double x1 = baricentric_2D(1, p) - baricentric_2D(2, p);
		return baricentricDeriv(0, 1) * (baricentric_2D(1, p) * baricentric_2D(2, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1)
						+ baricentric_2D(0, p) * baricentric_2D(1, p) * baricentric_2D(2, p) * LobattoPolynomialDerivative(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1))
				+ baricentricDeriv(2, 1) * (baricentric_2D(0, p) * baricentric_2D(1, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1)
						- baricentric_2D(0, p) * baricentric_2D(1, p) * baricentric_2D(2, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomialDerivative(polynomialOrder1_, x1));
	}

	Base::BasisFunctionSet* createDGBasisFunctionSet2DH1Triangle(std::size_t order) {
		Base::BasisFunctionSet* result(new Base::BasisFunctionSet(order));
		for (std::size_t i = 0; i < 3; ++i) {
			result->addBasisFunction(new BasisFunction2DVertexTriangle(i));
		}
		for (std::size_t k = 0; k + 2 <= order; ++k) {
			for(std::size_t i=0;i<3;++i){
				for (std::size_t j = 0; j < i; ++j) {
					result->addBasisFunction(new BasisFunction2DFaceTriangle(i, j, k));
				}
			}
			if(k + 3 <=order){
				for (std::size_t j = 0; j <= k; ++j) {
					result->addBasisFunction(new BasisFunction2DInteriorTriangle(k-j, j));
				}
			}
		}
		return result;
	}

	Base::BasisFunctionSet* createInteriorBasisFunctionSet2DH1Triangle(std::size_t order) {
		Base::BasisFunctionSet* result(new Base::BasisFunctionSet(order));
		for (std::size_t i = 0; i + 3 <= order; ++i) {
			for (std::size_t j = 0; i + j + 3 <= order; ++j) {
				result->addBasisFunction(new BasisFunction2DInteriorTriangle(i, j));
			}
		}
		return result;
	}

	std::vector<const Base::BasisFunctionSet*> createVertexBasisFunctionSet2DH1Triangle(std::size_t order) {
        std::vector<const Base::BasisFunctionSet*> result;
		Base::BasisFunctionSet* set;
		for (int i = 0; i < 3; ++i) {
			set = new Base::BasisFunctionSet(order);
			set->addBasisFunction(new BasisFunction2DVertexTriangle(i));
			result.push_back(set);
		}
        return result;
	}

	std::vector<const Base::OrientedBasisFunctionSet*> createFaceBasisFunctionSet2DH1Triangle(std::size_t order) {
        std::vector<const Base::OrientedBasisFunctionSet*> result;
		Base::OrientedBasisFunctionSet* set;
		Geometry::ReferenceTriangle& triangle = Geometry::ReferenceTriangle::Instance();
		std::vector<std::size_t> vertexindices(2);
		for (std::size_t i = 0; i < 3; ++i) {
			set = new Base::OrientedBasisFunctionSet(order, 0, i);
			vertexindices = triangle.getCodim1EntityLocalIndices(i);
			for (std::size_t j = 0; j + 2 <= order; ++j) {
				set->addBasisFunction(new BasisFunction2DFaceTriangle(vertexindices[0], vertexindices[1], j));
			}
			result.push_back(set);
			set = new Base::OrientedBasisFunctionSet(order, 1, i);
			for (std::size_t j = 0; j + 2 <= order; ++j) {
				set->addBasisFunction(new BasisFunction2DFaceTriangle(vertexindices[1], vertexindices[0], j));
			}
			result.push_back(set);
		}
        return result;
	}

}
