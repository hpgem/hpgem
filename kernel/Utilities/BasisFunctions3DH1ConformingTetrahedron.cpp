#include "BasisFunctions3DH1ConformingTetrahedron.hpp"
#include "helperFunctions.hpp"
#include "Base/BasisFunctionSet.hpp"
#include "Base/OrientedBasisFunctionSet.hpp"
#include "Geometry/ReferenceTetrahedron.hpp"

namespace Utilities {

double BasisFunction3DVertexTetrahedron::eval(const Geometry::PointReference& p) const {
	return baricentric_3D(node_, p);
}

double BasisFunction3DVertexTetrahedron::evalDeriv0(const Geometry::PointReference& p) const {
	return baricentricDeriv(node_, 0);
}

double BasisFunction3DVertexTetrahedron::evalDeriv1(const Geometry::PointReference& p) const {
	return baricentricDeriv(node_, 1);
}

double BasisFunction3DVertexTetrahedron::evalDeriv2(const Geometry::PointReference& p) const {
	return baricentricDeriv(node_, 2);
}

double BasisFunction3DEdgeTetrahedron::eval(const Geometry::PointReference& p) const {
	return baricentric_3D(node0_, p) * baricentric_3D(node1_, p) * LobattoPolynomial(polynomialOrder_, baricentric_3D(node0_, p) - baricentric_3D(node1_, p));
}

double BasisFunction3DEdgeTetrahedron::evalDeriv0(const Geometry::PointReference& p) const {
	return baricentricDeriv(node0_, 0) * (baricentric_3D(node1_, p) * LobattoPolynomial(polynomialOrder_, baricentric_3D(node0_, p) - baricentric_3D(node1_, p))
					+ baricentric_3D(node0_, p) * baricentric_3D(node1_, p) * LobattoPolynomialDerivative(polynomialOrder_, baricentric_3D(node0_, p) - baricentric_3D(node1_, p)))
			+ baricentricDeriv(node1_, 0) * (baricentric_3D(node0_, p) * LobattoPolynomial(polynomialOrder_, baricentric_3D(node0_, p) - baricentric_3D(node1_, p))
					- baricentric_3D(node0_, p) * baricentric_3D(node1_, p) * LobattoPolynomialDerivative(polynomialOrder_, baricentric_3D(node0_, p) - baricentric_3D(node1_, p)));
}

double BasisFunction3DEdgeTetrahedron::evalDeriv1(const Geometry::PointReference& p) const {
	return baricentricDeriv(node0_, 1) * (baricentric_3D(node1_, p) * LobattoPolynomial(polynomialOrder_, baricentric_3D(node0_, p) - baricentric_3D(node1_, p))
					+ baricentric_3D(node0_, p) * baricentric_3D(node1_, p)	* LobattoPolynomialDerivative(polynomialOrder_, baricentric_3D(node0_, p) - baricentric_3D(node1_, p)))
			+ baricentricDeriv(node1_, 1) * (baricentric_3D(node0_, p) * LobattoPolynomial(polynomialOrder_, baricentric_3D(node0_, p) - baricentric_3D(node1_, p))
					- baricentric_3D(node0_, p) * baricentric_3D(node1_, p) * LobattoPolynomialDerivative(polynomialOrder_, baricentric_3D(node0_, p) - baricentric_3D(node1_, p)));
}

double BasisFunction3DEdgeTetrahedron::evalDeriv2(const Geometry::PointReference& p) const {
	return baricentricDeriv(node0_, 2) * (baricentric_3D(node1_, p) * LobattoPolynomial(polynomialOrder_, baricentric_3D(node0_, p) - baricentric_3D(node1_, p))
					+ baricentric_3D(node0_, p) * baricentric_3D(node1_, p) * LobattoPolynomialDerivative(polynomialOrder_, baricentric_3D(node0_, p) - baricentric_3D(node1_, p)))
			+ baricentricDeriv(node1_, 2) * (baricentric_3D(node0_, p) * LobattoPolynomial(polynomialOrder_, baricentric_3D(node0_, p) - baricentric_3D(node1_, p))
					- baricentric_3D(node0_, p) * baricentric_3D(node1_, p) * LobattoPolynomialDerivative(polynomialOrder_, baricentric_3D(node0_, p) - baricentric_3D(node1_, p)));
}

double BasisFunction3DFaceTetrahedron::eval(const Geometry::PointReference& p) const {
	double x0(baricentric_3D(node0_, p) - baricentric_3D(node1_, p)), x1(baricentric_3D(node1_, p) - baricentric_3D(node2_, p));
	return baricentric_3D(node0_, p) * baricentric_3D(node1_, p) * baricentric_3D(node2_, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1);
}

double BasisFunction3DFaceTetrahedron::evalDeriv0(const Geometry::PointReference& p) const {
	double x0(baricentric_3D(node0_, p) - baricentric_3D(node1_, p)), x1(baricentric_3D(node1_, p) - baricentric_3D(node2_, p));
	return baricentricDeriv(node0_, 0) * (baricentric_3D(node1_, p) * baricentric_3D(node2_, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1)
					+ baricentric_3D(node0_, p) * baricentric_3D(node1_, p) * baricentric_3D(node2_, p) * LobattoPolynomialDerivative(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1))
			+ baricentricDeriv(node1_, 0) * (baricentric_3D(node0_, p) * baricentric_3D(node2_, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1)
					- baricentric_3D(node0_, p) * baricentric_3D(node1_, p) * baricentric_3D(node2_, p) * LobattoPolynomialDerivative(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1)
					+ baricentric_3D(node0_, p) * baricentric_3D(node1_, p) * baricentric_3D(node2_, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomialDerivative(polynomialOrder1_, x1))
			+ baricentricDeriv(node2_, 0) * (baricentric_3D(node0_, p) * baricentric_3D(node1_, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1)
					- baricentric_3D(node0_, p) * baricentric_3D(node1_, p) * baricentric_3D(node2_, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomialDerivative(polynomialOrder1_, x1));
}

double BasisFunction3DFaceTetrahedron::evalDeriv1(const Geometry::PointReference& p) const {
	double x0(baricentric_3D(node0_, p) - baricentric_3D(node1_, p)), x1(baricentric_3D(node1_, p) - baricentric_3D(node2_, p));
	return baricentricDeriv(node0_, 1) * (baricentric_3D(node1_, p) * baricentric_3D(node2_, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1)
					+ baricentric_3D(node0_, p) * baricentric_3D(node1_, p) * baricentric_3D(node2_, p) * LobattoPolynomialDerivative(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1))
			+ baricentricDeriv(node1_, 1) * (baricentric_3D(node0_, p) * baricentric_3D(node2_, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1)
					- baricentric_3D(node0_, p) * baricentric_3D(node1_, p) * baricentric_3D(node2_, p) * LobattoPolynomialDerivative(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1)
					+ baricentric_3D(node0_, p) * baricentric_3D(node1_, p) * baricentric_3D(node2_, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomialDerivative(polynomialOrder1_, x1))
			+ baricentricDeriv(node2_, 1) * (baricentric_3D(node0_, p) * baricentric_3D(node1_, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1)
					- baricentric_3D(node0_, p) * baricentric_3D(node1_, p) * baricentric_3D(node2_, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomialDerivative(polynomialOrder1_, x1));
}

double BasisFunction3DFaceTetrahedron::evalDeriv2(const Geometry::PointReference& p) const {
	double x0(baricentric_3D(node0_, p) - baricentric_3D(node1_, p)), x1(baricentric_3D(node1_, p) - baricentric_3D(node2_, p));
	return baricentricDeriv(node0_, 2) * (baricentric_3D(node1_, p) * baricentric_3D(node2_, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1)
					+ baricentric_3D(node0_, p) * baricentric_3D(node1_, p) * baricentric_3D(node2_, p) * LobattoPolynomialDerivative(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1))
			+ baricentricDeriv(node1_, 2) * (baricentric_3D(node0_, p) * baricentric_3D(node2_, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1)
					- baricentric_3D(node0_, p) * baricentric_3D(node1_, p) * baricentric_3D(node2_, p) * LobattoPolynomialDerivative(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1)
					+ baricentric_3D(node0_, p) * baricentric_3D(node1_, p) * baricentric_3D(node2_, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomialDerivative(polynomialOrder1_, x1))
			+ baricentricDeriv(node2_, 2) * (baricentric_3D(node0_, p) * baricentric_3D(node1_, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1)
					- baricentric_3D(node0_, p) * baricentric_3D(node1_, p) * baricentric_3D(node2_, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomialDerivative(polynomialOrder1_, x1));
}

double BasisFunction3DInteriorTetrahedron::eval(const Geometry::PointReference& p) const {
	double x0(baricentric_3D(0, p) - baricentric_3D(1, p)), x1(baricentric_3D(1, p) - baricentric_3D(2, p)), x2(baricentric_3D(2, p) - baricentric_3D(3, p));
	return baricentric_3D(0, p) * baricentric_3D(1, p) * baricentric_3D(2, p) * baricentric_3D(3, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) * LobattoPolynomial(polynomialOrder2_, x2);
}

double BasisFunction3DInteriorTetrahedron::evalDeriv0(const Geometry::PointReference& p) const {
	double x0(baricentric_3D(0, p) - baricentric_3D(1, p)), x1(baricentric_3D(1, p) - baricentric_3D(2, p)), x2(baricentric_3D(2, p) - baricentric_3D(3, p));
	return baricentricDeriv(0, 0) * (baricentric_3D(1, p) * baricentric_3D(2, p) * baricentric_3D(3, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) * LobattoPolynomial(polynomialOrder2_, x2)
					+ baricentric_3D(0, p) * baricentric_3D(1, p) * baricentric_3D(2, p) * baricentric_3D(3, p) * LobattoPolynomialDerivative(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) * LobattoPolynomial(polynomialOrder2_, x2))
			+ baricentricDeriv(1, 0) * (baricentric_3D(0, p) * baricentric_3D(2, p) * baricentric_3D(3, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) * LobattoPolynomial(polynomialOrder2_, x2)
					- baricentric_3D(0, p) * baricentric_3D(1, p) * baricentric_3D(2, p) * baricentric_3D(3, p) * LobattoPolynomialDerivative(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) * LobattoPolynomial(polynomialOrder2_, x2)
					+ baricentric_3D(0, p) * baricentric_3D(1, p) * baricentric_3D(2, p) * baricentric_3D(3, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomialDerivative(polynomialOrder1_, x1) * LobattoPolynomial(polynomialOrder2_, x2));
}

double BasisFunction3DInteriorTetrahedron::evalDeriv1(const Geometry::PointReference& p) const {
	double x0(baricentric_3D(0, p) - baricentric_3D(1, p)), x1(baricentric_3D(1, p) - baricentric_3D(2, p)), x2(baricentric_3D(2, p) - baricentric_3D(3, p));
	return baricentricDeriv(0, 1) * (baricentric_3D(1, p) * baricentric_3D(2, p) * baricentric_3D(3, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) * LobattoPolynomial(polynomialOrder2_, x2)
					+ baricentric_3D(0, p) * baricentric_3D(1, p) * baricentric_3D(2, p) * baricentric_3D(3, p) * LobattoPolynomialDerivative(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) * LobattoPolynomial(polynomialOrder2_, x2))
			+ baricentricDeriv(2, 1) * (baricentric_3D(0, p) * baricentric_3D(1, p) * baricentric_3D(3, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) * LobattoPolynomial(polynomialOrder2_, x2)
					- baricentric_3D(0, p) * baricentric_3D(1, p) * baricentric_3D(2, p) * baricentric_3D(3, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomialDerivative(polynomialOrder1_, x1) * LobattoPolynomial(polynomialOrder2_, x2)
					+ baricentric_3D(0, p) * baricentric_3D(1, p) * baricentric_3D(2, p) * baricentric_3D(3, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) * LobattoPolynomialDerivative(polynomialOrder2_, x2));
}

double BasisFunction3DInteriorTetrahedron::evalDeriv2(const Geometry::PointReference& p) const {
	double x0(baricentric_3D(0, p) - baricentric_3D(1, p)), x1(baricentric_3D(1, p) - baricentric_3D(2, p)), x2(baricentric_3D(2, p) - baricentric_3D(3, p));
	return baricentricDeriv(0, 2) * (baricentric_3D(1, p) * baricentric_3D(2, p) * baricentric_3D(3, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) * LobattoPolynomial(polynomialOrder2_, x2)
					+ baricentric_3D(0, p) * baricentric_3D(1, p) * baricentric_3D(2, p) * baricentric_3D(3, p) * LobattoPolynomialDerivative(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) * LobattoPolynomial(polynomialOrder2_, x2))
			+ baricentricDeriv(3, 2) * (baricentric_3D(0, p) * baricentric_3D(1, p) * baricentric_3D(2, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) * LobattoPolynomial(polynomialOrder2_, x2)
					- baricentric_3D(0, p) * baricentric_3D(1, p) * baricentric_3D(2, p) * baricentric_3D(3, p) * LobattoPolynomial(polynomialOrder0_, x0) * LobattoPolynomial(polynomialOrder1_, x1) * LobattoPolynomialDerivative(polynomialOrder2_, x2));
}

Base::BasisFunctionSet* createDGBasisFunctionSet3DH1Tetrahedron(int order) {
	Base::BasisFunctionSet* result = new Base::BasisFunctionSet(order);
	Geometry::ReferenceTetrahedron& tetrahedron = Geometry::ReferenceTetrahedron::Instance();
	std::vector<unsigned int> vectorOfPointIndexes(3);
	for (int i = 0; i < 4; ++i) {
		result->addBasisFunction(new BasisFunction3DVertexTetrahedron(i));
	}
	for (int i = 0; i < 6; ++i) {
		tetrahedron.getCodim2EntityLocalIndices(i, vectorOfPointIndexes);
		for (int j = 0; j <= order - 2; ++j) {
			result->addBasisFunction(new BasisFunction3DEdgeTetrahedron(vectorOfPointIndexes[0], vectorOfPointIndexes[1], j));
		}
	}
	for (int i = 0; i < 4; ++i) {
		tetrahedron.getCodim1EntityLocalIndices(i, vectorOfPointIndexes);
		for (int j = 0; j <= order - 3; ++j) {
			for (int k = 0; j + k <= order - 3; ++k) {
				result->addBasisFunction(new BasisFunction3DFaceTetrahedron(vectorOfPointIndexes[0], vectorOfPointIndexes[1], vectorOfPointIndexes[2], j, k));
			}
		}
	}
	for (int i = 0; i <= order - 4; ++i) {
		for (int j = 0; i + j <= order - 4; ++j) {
			for (int k = 0; i + j + k <= order - 4; ++k) {
				result->addBasisFunction(new BasisFunction3DInteriorTetrahedron(i, j, k));
			}
		}
	}
	return result;
}

Base::BasisFunctionSet* createInteriorBasisFunctionSet3DH1Tetrahedron(int order) {
	Base::BasisFunctionSet* result = new Base::BasisFunctionSet(order);
	for (int i = 0; i <= order - 4; ++i) {
		for (int j = 0; i + j <= order - 4; ++j) {
			for (int k = 0; i + j + k <= order - 4; ++k) {
				result->addBasisFunction(new BasisFunction3DInteriorTetrahedron(i, j, k));
			}
		}
	}
	return result;
}

void createVertexBasisFunctionSet3DH1Tetrahedron(int order, std::vector<const Base::BasisFunctionSet*>& result) {
	Base::BasisFunctionSet* set;
	for (int i = 0; i < 4; ++i) {
		set = new Base::BasisFunctionSet(order);
		set->addBasisFunction(new BasisFunction3DVertexTetrahedron(i));
		result.push_back(set);
	}
}

void createEdgeBasisFunctionSet3DH1Tetrahedron(int order, std::vector<const Base::OrientedBasisFunctionSet*>& result) {
	Base::OrientedBasisFunctionSet* set;
	Geometry::ReferenceTetrahedron& tetrahedron = Geometry::ReferenceTetrahedron::Instance();
	std::vector<unsigned int> vectorOfPointIndexes(2);
	for (int i = 0; i < 6; ++i) {
		set = new Base::OrientedBasisFunctionSet(order, 0, i);
		tetrahedron.getCodim2EntityLocalIndices(i, vectorOfPointIndexes);
		for (int j = 0; j <= order - 2; ++j) {
			set->addBasisFunction(new BasisFunction3DEdgeTetrahedron(vectorOfPointIndexes[0], vectorOfPointIndexes[1], j));
		}
		result.push_back(set);
		set = new Base::OrientedBasisFunctionSet(order, 1, i);
		for (int j = 0; j <= order - 2; ++j) {
			set->addBasisFunction(new BasisFunction3DEdgeTetrahedron(vectorOfPointIndexes[1], vectorOfPointIndexes[0], j));
		}
		result.push_back(set);
	}
}

void createFaceBasisFunctionSet3DH1Tetrahedron(int order, std::vector<const Base::OrientedBasisFunctionSet*>& result) {
	Base::OrientedBasisFunctionSet* set;
	Geometry::ReferenceTetrahedron& tetrahedron= Geometry::ReferenceTetrahedron::Instance();
	std::vector<unsigned int> vectorOfPointIndexes(3);
	for (int i = 0; i < 4; ++i) {
		tetrahedron.getCodim1EntityLocalIndices(i, vectorOfPointIndexes);
		set = new Base::OrientedBasisFunctionSet(order, 0, i);
		for (int j = 0; j <= order - 3; ++j) {
			for (int k = 0; j + k <= order - 3; ++k) {
				set->addBasisFunction(new BasisFunction3DFaceTetrahedron(vectorOfPointIndexes[0], vectorOfPointIndexes[1], vectorOfPointIndexes[2], j, k));
			}
		}
		result.push_back(set);
		set = new Base::OrientedBasisFunctionSet(order, 1, i);
		for (int j = 0; j <= order - 3; ++j) {
			for (int k = 0; j + k <= order - 3; ++k) {
				set->addBasisFunction(new BasisFunction3DFaceTetrahedron(vectorOfPointIndexes[0], vectorOfPointIndexes[2], vectorOfPointIndexes[1], j, k));
			}
		}
		result.push_back(set);
		set = new Base::OrientedBasisFunctionSet(order, 2, i);
		for (int j = 0; j <= order - 3; ++j) {
			for (int k = 0; j + k <= order - 3; ++k) {
				set->addBasisFunction(new BasisFunction3DFaceTetrahedron(vectorOfPointIndexes[1], vectorOfPointIndexes[2], vectorOfPointIndexes[0], j, k));
			}
		}
		result.push_back(set);
		set = new Base::OrientedBasisFunctionSet(order, 3, i);
		for (int j = 0; j <= order - 3; ++j) {
			for (int k = 0; j + k <= order - 3; ++k) {
				set->addBasisFunction(new BasisFunction3DFaceTetrahedron(vectorOfPointIndexes[1], vectorOfPointIndexes[0], vectorOfPointIndexes[2], j, k));
			}
		}
		result.push_back(set);
		set = new Base::OrientedBasisFunctionSet(order, 4, i);
		for (int j = 0; j <= order - 3; ++j) {
			for (int k = 0; j + k <= order - 3; ++k) {
				set->addBasisFunction(new BasisFunction3DFaceTetrahedron(vectorOfPointIndexes[2], vectorOfPointIndexes[1], vectorOfPointIndexes[0], j, k));
			}
		}
		result.push_back(set);
		set = new Base::OrientedBasisFunctionSet(order, 5, i);
		for (int j = 0; j <= order - 3; ++j) {
			for (int k = 0; j + k <= order - 3; ++k) {
				set->addBasisFunction(new BasisFunction3DFaceTetrahedron(vectorOfPointIndexes[2], vectorOfPointIndexes[0], vectorOfPointIndexes[1], j, k));
			}
		}
		result.push_back(set);
	}
}

}
