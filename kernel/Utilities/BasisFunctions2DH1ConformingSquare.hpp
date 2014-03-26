/*
 * BasisFunctions2DH1ConformingSquare.hpp
 *
 *  Created on: Mar 6, 2014
 *      Author: brinkf
 */

#ifndef BASISFUNCTIONS2DH1CONFORMINGSQUARE_HPP_
#define BASISFUNCTIONS2DH1CONFORMINGSQUARE_HPP_

#include "Base/TestErrorDebug.hpp"
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

	class BasisFunction2DVertexSquare: public Base::BaseBasisFunction
	{
	public:
		BasisFunction2DVertexSquare(int node):nodePosition0_((node%2)*2-1),nodePosition1_((node/2)*2-1){}

		double eval(const Geometry::PointReference& p) const;

		double evalDeriv0(const Geometry::PointReference& p) const;

		double evalDeriv1(const Geometry::PointReference& p) const;

	private:
		int nodePosition0_;
		int nodePosition1_;
	};

	class BasisFunction2DFaceSquare_0: public Base::BaseBasisFunction
	{
	public:
		BasisFunction2DFaceSquare_0(int node0, int node1, int polynomialOrder);

		double eval(const Geometry::PointReference& p) const;

		double evalDeriv0(const Geometry::PointReference& p) const;

		double evalDeriv1(const Geometry::PointReference& p) const;

	private:
		int edgePosition_;
		int mirroring_;
		int polynomialOrder_;
	};

	class BasisFunction2DFaceSquare_1: public Base::BaseBasisFunction
	{
	public:
		BasisFunction2DFaceSquare_1(int node0, int node1, int polynomialOrder);

		double eval(const Geometry::PointReference& p) const;

		double evalDeriv0(const Geometry::PointReference& p) const;

		double evalDeriv1(const Geometry::PointReference& p) const;

	private:
		int edgePosition_;
		int mirroring_;
		int polynomialOrder_;
	};

	class BasisFunction2DInteriorSquare: public Base::BaseBasisFunction
	{
	public:
		BasisFunction2DInteriorSquare(int polynomialOrder0,int polynomialOrder1):polynomialOrder0_(polynomialOrder0),polynomialOrder1_(polynomialOrder1){}

		double eval(const Geometry::PointReference& p) const;

		double evalDeriv0(const Geometry::PointReference& p) const;

		double evalDeriv1(const Geometry::PointReference& p) const;
	private:
		int polynomialOrder0_,polynomialOrder1_;
	};

Base::BasisFunctionSet* createDGBasisFunctionSet2DH1Square(int order);

Base::BasisFunctionSet* createInteriorBasisFunctionSet2DH1Square(int order);

void createVertexBasisFunctionSet2DH1Square(int order, std::vector<const Base::BasisFunctionSet*>& result);

void createFaceBasisFunctionSet2DH1Square(int order, std::vector<const Base::OrientedBasisFunctionSet*>& result);

}



#endif /* BASISFUNCTIONS2DH1CONFORMINGSQUARE_HPP_ */
