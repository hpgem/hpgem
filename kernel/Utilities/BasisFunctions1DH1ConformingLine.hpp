/*
 * BasisFunctions1DConformingLine.cpp
 *
 *  Created on: Mar 5, 2014
 *      Author: brinkf
 */

#ifndef BASISFUNCTIONS1DCONFORMINGLINE_CPP_
#define BASISFUNCTIONS1DCONFORMINGLINE_CPP_

#include "Base/BaseBasisFunction.hpp"
#include <vector>

namespace Base{
class BasisFunctionSet;
}

namespace Geometry{
class PointReference;
}

namespace Utilities{

	class BasisFunction1DVertexLine:public Base::BaseBasisFunction
	{
	public:
		BasisFunction1DVertexLine(int node):nodePosition_(2*node-1){}
	double eval(const Geometry::PointReference& p) const;

	double evalDeriv0(const Geometry::PointReference& p) const;
	private:
		int nodePosition_;
	};

	class BasisFunction1DInteriorLine:public Base::BaseBasisFunction
	{
	public:
		BasisFunction1DInteriorLine(int polynomialOrder):polynomialOrder_(polynomialOrder){}

	double eval(const Geometry::PointReference& p) const;

	double evalDeriv0(const Geometry::PointReference& p) const;
	private:
		int polynomialOrder_;
	};

Base::BasisFunctionSet* createDGBasisFunctionSet1DH1Line(int polynomialOrder);

Base::BasisFunctionSet* createInteriorBasisFunctionSet1DH1Line(int polynomialOrder);

void createVertexBasisFunctionSet1DH1Line(int polynomialOrder,std::vector<const Base::BasisFunctionSet*>& result);

}


#endif /* BASISFUNCTIONS1DCONFORMINGLINE_CPP_ */
