#include "BasisFunctions2DH1ConformingSquare.hpp"
#include "helperFunctions.hpp"
#include "Base/BasisFunctionSet.hpp"
#include "Base/OrientedBasisFunctionSet.hpp"
#include "Geometry/ReferenceSquare.hpp"

namespace Utilities {

double BasisFunction2DVertexSquare::eval(const Geometry::PointReference& p) const {
	return (1 + nodePosition0_ * p[0]) * (1 + nodePosition1_ * p[1]) / 4.;
}

double BasisFunction2DVertexSquare::evalDeriv0(const Geometry::PointReference& p) const {
	return nodePosition0_ * (1 + nodePosition1_ * p[1]) / 4.;
}

double BasisFunction2DVertexSquare::evalDeriv1(const Geometry::PointReference& p) const {
	return nodePosition1_ * (1 + nodePosition0_ * p[0]) / 4.;
}

BasisFunction2DFaceSquare_0::BasisFunction2DFaceSquare_0(int node0, int node1,int polynomialOrder) :
		polynomialOrder_(polynomialOrder) {
	TestErrorDebug((node0 + node1) % 2 == 1,"please use BasisFunction2DFaceSquare_1 for edges that are aligned vertically");
	mirroring_ = (node0 > node1) ? -1 : 1;
	edgePosition_ = (node0 + node1 < 3) ? -1 : 1;
}

double BasisFunction2DFaceSquare_0::eval(const Geometry::PointReference& p) const {
	return (1 + edgePosition_ * p[1]) * (1 - p[0]) * (1 + p[0]) * LobattoPolynomial(polynomialOrder_, mirroring_ * p[0]) / 8.;
}

double BasisFunction2DFaceSquare_0::evalDeriv0(const Geometry::PointReference& p) const {
	return (1 + edgePosition_ * p[1]) * (-p[0] * LobattoPolynomial(polynomialOrder_, mirroring_ * p[0]) + (1 - p[0]) * (1 + p[0])
							* LobattoPolynomialDerivative(polynomialOrder_, mirroring_ * p[0]) * mirroring_ / 2.) / 4.;
}

double BasisFunction2DFaceSquare_0::evalDeriv1(const Geometry::PointReference& p) const {
	return edgePosition_ * (1 - p[0]) * (1 + p[0]) * LobattoPolynomial(polynomialOrder_, mirroring_ * p[0]) / 8.;
}

BasisFunction2DFaceSquare_1::BasisFunction2DFaceSquare_1(int node0, int node1, int polynomialOrder) :
		polynomialOrder_(polynomialOrder) {
	TestErrorDebug((node0 + node1) % 2 == 0,"please use BasisFunction2DFaceSquare_0 for edges that are aligned horizontally");
	mirroring_ = (node0 > node1) ? -1 : 1;
	edgePosition_ = (node0 + node1 < 3) ? -1 : 1;
}

double BasisFunction2DFaceSquare_1::eval(const Geometry::PointReference& p) const {
	return (1 + edgePosition_ * p[0]) * (1 - p[1]) * (1 + p[1]) * LobattoPolynomial(polynomialOrder_, mirroring_ * p[1]) / 8.;
}

double BasisFunction2DFaceSquare_1::evalDeriv0(const Geometry::PointReference& p) const {
	return edgePosition_ * (1 - p[1]) * (1 + p[1]) * LobattoPolynomial(polynomialOrder_, mirroring_ * p[1]) / 8.;
}

double BasisFunction2DFaceSquare_1::evalDeriv1(const Geometry::PointReference& p) const {
	return (1 + edgePosition_ * p[0]) * (-p[1] * LobattoPolynomial(polynomialOrder_, mirroring_ * p[1]) + (1 - p[1]) * (1 + p[1])
							* LobattoPolynomialDerivative(polynomialOrder_, mirroring_ * p[1]) * mirroring_ / 2.) / 4.;
}

double BasisFunction2DInteriorSquare::eval(const Geometry::PointReference& p) const {
	return (1 - p[0]) * (1 + p[0]) * (1 - p[1]) * (1 + p[1]) * LobattoPolynomial(polynomialOrder0_, p[0]) * LobattoPolynomial(polynomialOrder1_, p[1]) / 16.;
}

double BasisFunction2DInteriorSquare::evalDeriv0(const Geometry::PointReference& p) const {
	return LobattoPolynomial(polynomialOrder1_, p[1]) * (1 - p[1]) * (1 + p[1])/ 4.
			* (-p[0] * LobattoPolynomial(polynomialOrder0_, p[0]) / 2. + (1 - p[0]) * (1 + p[0]) * LobattoPolynomialDerivative(polynomialOrder0_, p[0]) / 4.);
}

double BasisFunction2DInteriorSquare::evalDeriv1(const Geometry::PointReference& p) const {
	return LobattoPolynomial(polynomialOrder0_, p[0]) * (1 - p[0]) * (1 + p[0])/ 4.
			* (-p[1] * LobattoPolynomial(polynomialOrder1_, p[1]) / 2. + (1 - p[1]) * (1 + p[1]) * LobattoPolynomialDerivative(polynomialOrder1_, p[1]) / 4.);
}

Base::BasisFunctionSet* createDGBasisFunctionSet2DH1Square(int order) {
	Base::BasisFunctionSet* result(new Base::BasisFunctionSet(order));
	for (int i = 0; i < 4; ++i) {
		result->addBasisFunction(new BasisFunction2DVertexSquare(i));
	}
	for (int i = 0; i <= order - 2; ++i) {
		result->addBasisFunction(new BasisFunction2DFaceSquare_0(0, 1, i));
		result->addBasisFunction(new BasisFunction2DFaceSquare_0(2, 3, i));
		result->addBasisFunction(new BasisFunction2DFaceSquare_1(0, 2, i));
		result->addBasisFunction(new BasisFunction2DFaceSquare_1(1, 3, i));
		for (int j = 0; j <= order - 2; ++j) {
			result->addBasisFunction(new BasisFunction2DInteriorSquare(i, j));
		}
	}
	return result;
}

Base::BasisFunctionSet* createInteriorBasisFunctionSet2DH1Square(int order) {
	Base::BasisFunctionSet* result(new Base::BasisFunctionSet(order));
	for (int i = 0; i <= order - 2; ++i) {
		for (int j = 0; j <= order - 2; ++j) {
			result->addBasisFunction(new BasisFunction2DInteriorSquare(i, j));
		}
	}
	return result;
}

void createVertexBasisFunctionSet2DH1Square(int order, std::vector<const Base::BasisFunctionSet*>& result) {
	Base::BasisFunctionSet* set;
	for (int i = 0; i < 4; ++i) {
		set = new Base::BasisFunctionSet(order);
		set->addBasisFunction(new BasisFunction2DVertexSquare(i));
		result.push_back(set);
	}
}

void createFaceBasisFunctionSet2DH1Square(int order, std::vector<const Base::OrientedBasisFunctionSet*>& result) {
	Geometry::ReferenceSquare& square = Geometry::ReferenceSquare::Instance();
	Base::OrientedBasisFunctionSet* set;
	std::vector<unsigned int> vertexindices(2);
	for (int i = 0; i < 4; ++i) {
		set = new Base::OrientedBasisFunctionSet(order, 0, i);
		square.getCodim1EntityLocalIndices(i, vertexindices);
		for (int j = 0; j <= order - 2; ++j) {
			if ((vertexindices[0] + vertexindices[1]) % 2 == 1)
				set->addBasisFunction(new BasisFunction2DFaceSquare_0(vertexindices[0], vertexindices[1], j));
			else
				set->addBasisFunction(new BasisFunction2DFaceSquare_1(vertexindices[0], vertexindices[1], j));
		}
		result.push_back(set);
		set = new Base::OrientedBasisFunctionSet(order, 1, i);
		for (int j = 0; j <= order - 2; ++j) {
			if ((vertexindices[0] + vertexindices[1]) % 2 == 1)
				set->addBasisFunction(new BasisFunction2DFaceSquare_0(vertexindices[1], vertexindices[0], j));
			else
				set->addBasisFunction(new BasisFunction2DFaceSquare_1(vertexindices[1], vertexindices[0], j));
		}
		result.push_back(set);
	}
}

}
