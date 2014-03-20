#include "BasisFunctions2DH1ConformingTriangle.hpp"
#include "helperFunctions.hpp"
#include "Base/BasisFunctionSet.hpp"
#include "Base/OrientedBasisFunctionSet.hpp"
#include "Geometry/ReferenceTriangle.hpp"

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

	Base::BasisFunctionSet* createDGBasisFunctionSet2DH1Triangle(int order) {
		Base::BasisFunctionSet* result(new Base::BasisFunctionSet(order));
		for (int i = 0; i < 3; ++i) {
			result->addBasisFunction(new BasisFunction2DVertexTriangle(i));
			for (int j = 0; j < i; ++j) {
				for (int k = 0; k <= order - 2; ++k) {
					result->addBasisFunction(new BasisFunction2DFaceTriangle(i, j, k));
				}
			}
		}
		for (int i = 0; i <= order - 3; ++i) {
			for (int j = 0; i + j <= order - 3; ++j) {
				result->addBasisFunction(new BasisFunction2DInteriorTriangle(i, j));
			}
		}
		return result;
	}

	Base::BasisFunctionSet* createInteriorBasisFunctionSet2DH1Triangle(int order) {
		Base::BasisFunctionSet* result(new Base::BasisFunctionSet(order));
		for (int i = 0; i <= order - 3; ++i) {
			for (int j = 0; i + j <= order - 3; ++j) {
				result->addBasisFunction(new BasisFunction2DInteriorTriangle(i, j));
			}
		}
		return result;
	}

	void createVertexBasisFunctionSet2DH1Triangle(int order,std::vector<const Base::BasisFunctionSet*>& result) {
		Base::BasisFunctionSet* set;
		for (int i = 0; i < 3; ++i) {
			set = new Base::BasisFunctionSet(order);
			set->addBasisFunction(new BasisFunction2DVertexTriangle(i));
			result.push_back(set);
		}
	}

	void createFaceBasisFunctionSet2DH1Triangle(int order,std::vector<const Base::OrientedBasisFunctionSet*>& result) {
		Base::OrientedBasisFunctionSet* set;
		Geometry::ReferenceTriangle triangle = Geometry::ReferenceTriangle::Instance();
		std::vector<unsigned int> vertexindices(2);
		for (int i = 0; i < 3; ++i) {
			set = new Base::OrientedBasisFunctionSet(order, 0, i);
			triangle.getCodim1EntityLocalIndices(i, vertexindices);
			for (int j = 0; j <= order - 2; ++j) {
				set->addBasisFunction(new BasisFunction2DFaceTriangle(vertexindices[0], vertexindices[1], j));
			}
			result.push_back(set);
			set = new Base::OrientedBasisFunctionSet(order, 1, i);
			for (int j = 0; j <= order - 2; ++j) {
				set->addBasisFunction(new BasisFunction2DFaceTriangle(vertexindices[1], vertexindices[0], j));
			}
			result.push_back(set);
		}
	}

}
