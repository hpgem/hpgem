/*
 * BasisFunctions3DH1ConformingCube.hpp
 *
 *  Created on: Mar 6, 2014
 *      Author: brinkf
 */

#ifndef BASISFUNCTIONS3DH1CONFORMINGCUBE_HPP_
#define BASISFUNCTIONS3DH1CONFORMINGCUBE_HPP_

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

	class BasisFunction3DVertexCube : public Base::BaseBasisFunction
	{
	public:
		BasisFunction3DVertexCube(int node);

		double eval(const Geometry::PointReference& p) const;

		double evalDeriv0(const Geometry::PointReference& p) const;

		double evalDeriv1(const Geometry::PointReference& p) const;

		double evalDeriv2(const Geometry::PointReference& p) const;

	private:
		int nodePosition0_,nodePosition1_,nodePosition2_;
	};

	class BasisFunction3DEdgeCube_0:public Base::BaseBasisFunction
	{
	public:
		BasisFunction3DEdgeCube_0(int node0, int node1, int polynomialOrder);

		double eval(const Geometry::PointReference& p) const;

		double evalDeriv0(const Geometry::PointReference& p) const;

		double evalDeriv1(const Geometry::PointReference& p) const;

		double evalDeriv2(const Geometry::PointReference& p) const;

	private:
		int edgePosition1_;
		int edgePosition2_;
		int mirroring_;
		int polynomialOrder_;
	};

	class BasisFunction3DEdgeCube_1:public Base::BaseBasisFunction
	{
	public:
		BasisFunction3DEdgeCube_1(int node0, int node1, int polynomialOrder);

		double eval(const Geometry::PointReference& p) const;

		double evalDeriv0(const Geometry::PointReference& p) const;

		double evalDeriv1(const Geometry::PointReference& p) const;

		double evalDeriv2(const Geometry::PointReference& p) const;

	private:
		int edgePosition0_;
		int edgePosition2_;
		int mirroring_;
		int polynomialOrder_;
	};

	class BasisFunction3DEdgeCube_2:public Base::BaseBasisFunction
	{
	public:
		BasisFunction3DEdgeCube_2(int node0, int node1, int polynomialOrder);

		double eval(const Geometry::PointReference& p) const;

		double evalDeriv0(const Geometry::PointReference& p) const;

		double evalDeriv1(const Geometry::PointReference& p) const;

		double evalDeriv2(const Geometry::PointReference& p) const;

	private:
		int edgePosition0_;
		int edgePosition1_;
		int mirroring_;
		int polynomialOrder_;
	};

	class BasisFunction3DFaceCube_0:public Base::BaseBasisFunction
	{
	public:
		BasisFunction3DFaceCube_0(int node0, int node1, int node2, int polynomialOrder1, int polynomialOrder2);

		double eval(const Geometry::PointReference& p) const;

		double evalDeriv0(const Geometry::PointReference& p) const;

		double evalDeriv1(const Geometry::PointReference& p) const;

		double evalDeriv2(const Geometry::PointReference& p) const;
	private:
		int facePosition_;
		int mirroring1_;
		int mirroring2_;
		int polynomialOrder1_;
		int polynomialOrder2_;
	};

	class BasisFunction3DFaceCube_1:public Base::BaseBasisFunction
	{
	public:
		BasisFunction3DFaceCube_1(int node0, int node1, int node2, int polynomialOrder0, int polynomialOrder2);

		double eval(const Geometry::PointReference& p) const;

		double evalDeriv0(const Geometry::PointReference& p) const;

		double evalDeriv1(const Geometry::PointReference& p) const;

		double evalDeriv2(const Geometry::PointReference& p) const;
	private:
		int facePosition_;
		int mirroring0_;
		int mirroring2_;
		int polynomialOrder0_;
		int polynomialOrder2_;
	};

	class BasisFunction3DFaceCube_2:public Base::BaseBasisFunction
	{
	public:
		BasisFunction3DFaceCube_2(int node0, int node1, int node2, int polynomialOrder0, int polynomialOrder1);

		double eval(const Geometry::PointReference& p) const;

		double evalDeriv0(const Geometry::PointReference& p) const;

		double evalDeriv1(const Geometry::PointReference& p) const;

		double evalDeriv2(const Geometry::PointReference& p) const;
	private:
		int facePosition_;
		int mirroring0_;
		int mirroring1_;
		int polynomialOrder0_;
		int polynomialOrder1_;
	};

	class BasisFunction3DInteriorCube: public Base::BaseBasisFunction
	{
	public:
		BasisFunction3DInteriorCube(int polynomialOrder0,int polynomialOrder1, int polynomialOrder2):polynomialOrder0_(polynomialOrder0),polynomialOrder1_(polynomialOrder1),polynomialOrder2_(polynomialOrder2){}

		double eval(const Geometry::PointReference& p) const;

		double evalDeriv0(const Geometry::PointReference& p) const;

		double evalDeriv1(const Geometry::PointReference& p) const;

		double evalDeriv2(const Geometry::PointReference& p) const;

	private:
		int polynomialOrder0_,polynomialOrder1_,polynomialOrder2_;
	};

	Base::BasisFunctionSet* createDGBasisFunctionSet3DH1Cube(int order);

	Base::BasisFunctionSet* createInteriorBasisFunctionSet3DH1Cube(int order);

	void createVertexBasisFunctionSet3DH1Cube(int order, std::vector<const Base::BasisFunctionSet*>& result);

	void createEdgeBasisFunctionSet3DH1Cube(int order, std::vector<const Base::OrientedBasisFunctionSet*>& result);

	void createFaceBasisFunctionSet3DH1Cube(int order, std::vector<const Base::OrientedBasisFunctionSet*>& result);
}



#endif /* BASISFUNCTIONS3DH1CONFORMINGCUBE_HPP_ */
