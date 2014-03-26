#include "BasisFunctions1DH1ConformingLine.hpp"
#include "helperFunctions.hpp"
#include "Base/BasisFunctionSet.hpp"

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
