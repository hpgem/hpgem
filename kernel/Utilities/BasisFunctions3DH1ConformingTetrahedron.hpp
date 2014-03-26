/*
 * BasisFunctions3DH1ConformingCube.hpp
 *
 *  Created on: Mar 6, 2014
 *      Author: brinkf
 */

#ifndef BASISFUNCTIONS3DH1CONFORMINGTETRAHEDRON_HPP_
#define BASISFUNCTIONS3DH1CONFORMINGTETRAHEDRON_HPP_

#include "Base/BaseBasisFunction.hpp"
#include <vector>

namespace Base{
class BasisFunctionSet;
class OrientedBasisFunctionSet;
}

namespace Geometry{
class PointReference;
}

namespace Utilities
{

	class BasisFunction3DVertexTetrahedron : public Base::BaseBasisFunction
	{
	public:
		BasisFunction3DVertexTetrahedron(int node):node_(node){}

	double eval(const Geometry::PointReference& p) const;

	double evalDeriv0(const Geometry::PointReference& p) const;

	double evalDeriv1(const Geometry::PointReference& p) const;

	double evalDeriv2(const Geometry::PointReference& p) const;

	private:
		int node_;
	};

	class BasisFunction3DEdgeTetrahedron : public Base::BaseBasisFunction
	{
	public:
		BasisFunction3DEdgeTetrahedron(int node0,int node1,int polynomialOrder):node0_(node0),node1_(node1),polynomialOrder_(polynomialOrder){}

	double eval(const Geometry::PointReference& p) const;

	double evalDeriv0(const Geometry::PointReference& p) const;

	double evalDeriv1(const Geometry::PointReference& p) const;

	double evalDeriv2(const Geometry::PointReference& p) const;

	private:
		int node0_,node1_,polynomialOrder_;
	};

	class BasisFunction3DFaceTetrahedron : public Base::BaseBasisFunction
	{
	public:
		BasisFunction3DFaceTetrahedron(int node0,int node1,int node2,int polynomialOrder0,int polynomialOrder1):node0_(node0),node1_(node1),node2_(node2),polynomialOrder0_(polynomialOrder0),polynomialOrder1_(polynomialOrder1){}

	double eval(const Geometry::PointReference& p) const;

	double evalDeriv0(const Geometry::PointReference& p) const;

	double evalDeriv1(const Geometry::PointReference& p) const;

	double evalDeriv2(const Geometry::PointReference& p) const;

	private:
		int node0_,node1_,node2_,polynomialOrder0_,polynomialOrder1_;
	};

	class BasisFunction3DInteriorTetrahedron : public Base::BaseBasisFunction
	{
	public:
		BasisFunction3DInteriorTetrahedron(int polynomialOrder0,int polynomialOrder1,int polynomialOrder2):polynomialOrder0_(polynomialOrder0),polynomialOrder1_(polynomialOrder1),polynomialOrder2_(polynomialOrder2){}

	double eval(const Geometry::PointReference& p) const;

	double evalDeriv0(const Geometry::PointReference& p) const;

	double evalDeriv1(const Geometry::PointReference& p) const;

	double evalDeriv2(const Geometry::PointReference& p) const;

	private:
		int polynomialOrder0_,polynomialOrder1_,polynomialOrder2_;
	};

Base::BasisFunctionSet* createDGBasisFunctionSet3DH1Tetrahedron(int order);

Base::BasisFunctionSet* createInteriorBasisFunctionSet3DH1Tetrahedron(int order);

void createVertexBasisFunctionSet3DH1Tetrahedron(int order, std::vector<const Base::BasisFunctionSet*>& result);

void createEdgeBasisFunctionSet3DH1Tetrahedron(int order, std::vector<const Base::OrientedBasisFunctionSet*>& result);

void createFaceBasisFunctionSet3DH1Tetrahedron(int order, std::vector<const Base::OrientedBasisFunctionSet*>& result);

}

#endif
