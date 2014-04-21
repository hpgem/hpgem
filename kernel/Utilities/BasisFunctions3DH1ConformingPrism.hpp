/*
 * BasisFunctions3DH1ConformingPrism.hpp
 *
 *  Created on: Mar 7, 2014
 *      Author: brinkf
 */

#ifndef BASISFUNCTIONS3DH1CONFORMINGPRISM_HPP_
#define BASISFUNCTIONS3DH1CONFORMINGPRISM_HPP_

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

	class BasisFunction3DVertexPrism: public Base::BaseBasisFunction
	{
	public:
		BasisFunction3DVertexPrism(int node);

		double eval(const Geometry::PointReference& p) const;

		double evalDeriv0(const Geometry::PointReference& p) const;

		double evalDeriv1(const Geometry::PointReference& p) const;

		double evalDeriv2(const Geometry::PointReference& p) const;

	private:
		int nodePosition_,node_;//node is number inside triangle
	};

	class BasisFunction3DEdgePrism_0: public Base::BaseBasisFunction
	{
	public:
		BasisFunction3DEdgePrism_0(int node0,int node1,int polynomialOrder):node0_(node0%3),node1_(node1%3),polynomialOrder_(polynomialOrder),edgePosition_((node0/3)*2-1){}

		double eval(const Geometry::PointReference& p) const;

		double evalDeriv0(const Geometry::PointReference& p) const;

		double evalDeriv1(const Geometry::PointReference& p) const;

		double evalDeriv2(const Geometry::PointReference& p) const;

	private:
		int edgePosition_,node0_,node1_,polynomialOrder_;
	};

	class BasisFunction3DEdgePrism_1: public Base::BaseBasisFunction
	{
	public:
		BasisFunction3DEdgePrism_1(int node0, int node1, int polynomialOrder);

		double eval(const Geometry::PointReference& p) const;

		double evalDeriv0(const Geometry::PointReference& p) const;

		double evalDeriv1(const Geometry::PointReference& p) const;

		double evalDeriv2(const Geometry::PointReference& p) const;

	private:
		int mirroring_,node_,polynomialOrder_;
	};

	class BasisFunction3DFacePrism_0: public Base::BaseBasisFunction
	{
	public:
		BasisFunction3DFacePrism_0(int node0, int node1, int node2, int polynomialOrder0, int polynomialOrder1);

		double eval(const Geometry::PointReference& p) const;

		double evalDeriv0(const Geometry::PointReference& p) const;

		double evalDeriv1(const Geometry::PointReference& p) const;

		double evalDeriv2(const Geometry::PointReference& p) const;

	private:
		int facePosition_,polynomialOrder0_,polynomialOrder1_,node0_,node1_,node2_;
	};

	class BasisFunction3DFacePrism_1: public Base::BaseBasisFunction
	{
	public:
		BasisFunction3DFacePrism_1(int node0, int node1, int node2, int polynomialOrder0, int polynomialOrder1);

		double eval(const Geometry::PointReference& p) const;

		double evalDeriv0(const Geometry::PointReference& p) const;

		double evalDeriv1(const Geometry::PointReference& p) const;

		double evalDeriv2(const Geometry::PointReference& p) const;

	private:
		int mirroring_,node0_,node1_,polynomialOrder0_,polynomialOrder1_;
	};

	class BasisFunction3DInteriorPrism: public Base::BaseBasisFunction
	{
	public:
		BasisFunction3DInteriorPrism(int polynomialOrder0,int polynomialOrder1,int polynomialOrder2):polnomialOrder0_(polynomialOrder0),polynomialOrder1_(polynomialOrder1),polynomialOrder2_(polynomialOrder2){}

		double eval(const Geometry::PointReference& p) const;

		double evalDeriv0(const Geometry::PointReference& p) const;

		double evalDeriv1(const Geometry::PointReference& p) const;

		double evalDeriv2(const Geometry::PointReference& p) const;

	private:
		int polnomialOrder0_,polynomialOrder1_,polynomialOrder2_;
	};

	Base::BasisFunctionSet* createDGBasisFunctionSet3DH1ConformingPrism(int order);

	Base::BasisFunctionSet* createInteriorBasisFunctionSet3DH1ConformingPrism(int order);

	void createVertexBasisFunctionSet3DH1ConformingPrism(int order, std::vector<const Base::BasisFunctionSet*>& result);

	void createEdgeBasisFunctionSet3DH1ConformingPrism(int order, std::vector<const Base::OrientedBasisFunctionSet*>& result);

	void CreateFaceBasisFunctionSet3DH1ConformingPrism(int order, std::vector<const Base::OrientedBasisFunctionSet*>& result);

}


#endif /* BASISFUNCTIONS3DH1CONFORMINGPRISM_HPP_ */
