/*
 * BasisFunctions2DH1ConformingTriangle.hpp
 *
 *  Created on: Mar 6, 2014
 *      Author: brinkf
 */

#ifndef BASISFUNCTIONS2DH1CONFORMINGTRIANGLE_HPP_
#define BASISFUNCTIONS2DH1CONFORMINGTRIANGLE_HPP_

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

	class BasisFunction2DVertexTriangle: public Base::BaseBasisFunction
	{
	public:
		BasisFunction2DVertexTriangle(int node):node_(node){}

	double eval(const Geometry::PointReference& p) const;

	double evalDeriv0(const Geometry::PointReference& p) const;

	double evalDeriv1(const Geometry::PointReference& p) const;

	private:
		int node_;
	};

	class BasisFunction2DFaceTriangle: public Base::BaseBasisFunction
	{
	public:
		BasisFunction2DFaceTriangle(int node0,int node1,int polynomialOrder):node0_(node0),node1_(node1),polynomialOrder_(polynomialOrder){}

	double eval(const Geometry::PointReference& p) const;

	double evalDeriv0(const Geometry::PointReference& p) const;

	double evalDeriv1(const Geometry::PointReference& p) const;

	private:
		int node0_, node1_,polynomialOrder_;
	};

	class BasisFunction2DInteriorTriangle: public Base::BaseBasisFunction
	{
	public:
		BasisFunction2DInteriorTriangle(int polynomialOrder0,int polynomialOrder1):polynomialOrder0_(polynomialOrder0),polynomialOrder1_(polynomialOrder1){}

	double eval(const Geometry::PointReference& p) const;

	double evalDeriv0(const Geometry::PointReference& p) const;

	double evalDeriv1(const Geometry::PointReference& p) const;

	private:
		int polynomialOrder0_,polynomialOrder1_;
	};

Base::BasisFunctionSet* createDGBasisFunctionSet2DH1Triangle(int order);

Base::BasisFunctionSet* createInteriorBasisFunctionSet2DH1Triangle(int order);

void createVertexBasisFunctionSet2DH1Triangle(int order, std::vector<const Base::BasisFunctionSet*>& result);

void createFaceBasisFunctionSet2DH1Triangle(int order, std::vector<const Base::OrientedBasisFunctionSet*>& result);

}




#endif /* BASISFUNCTIONS2DH1CONFORMINGTRIANGLE_HPP_ */
