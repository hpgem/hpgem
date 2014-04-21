/*
 * 020GaussQuadratureRulesForLine_UnitTest.cpp
 *
 *  Created on: Apr 14, 2014
 *      Author: brinkf
 */

//naming convention: <Digit><ClassName>_UnitTest.cpp where <Digit> is a number that will make sure
//the unit tests are ordered such that the first failing unit test indicate the culprit class and
//other 'unit' tests may assume correct execution of all prior unit tests

#include "Integration/QuadratureRules/GaussQuadratureRulesForLine.hpp"
#include "cassert"

#include "Utilities/BasisFunctions1DH1ConformingLine.hpp"
#include "Geometry/ReferenceLine.hpp"
#include "Base/BasisFunctionSet.hpp"

void testRule(QuadratureRules::GaussQuadratureRule& test,int expectedOrder){
	cout<<test.getName();
	assert(("dimension",test.dimension()==1));
	assert(("order",test.order()>=expectedOrder));
	assert(("forReferenceGeometry",typeid(*test.forReferenceGeometry())==typeid(Geometry::ReferenceLine)));
	Geometry::PointReference point(1);

	Base::BasisFunctionSet* functions = Utilities::createDGBasisFunctionSet1DH1Line(expectedOrder);

	for(int i=0;i<functions->size();++i){
		double integrated=0;
		for(int j=0;j<test.nrOfPoints();++j){
			test.getPoint(j,point);
			integrated+=test.weight(j)*functions->eval(i,point);
		}
		if(i<2){
			assert(("integration",fabs(integrated-1)<1e-12));
		}else if(i==2){
			assert(("integration",fabs(integrated+sqrt(2./3.))<1e-12));
		}else{
			assert(("integration",fabs(integrated)<1e-12));
		}

	}

	delete functions;
}

int main(){

	testRule(QuadratureRules::Cn1_1_1::Instance(),1);
	testRule(QuadratureRules::Cn1_3_4::Instance(),3);
	testRule(QuadratureRules::Cn1_5_9::Instance(),5);
	testRule(QuadratureRules::C1_7_x::Instance(),7);
	testRule(QuadratureRules::C1_9_25::Instance(),9);
	testRule(QuadratureRules::C1_11_36::Instance(),10);///\BUG there are no 11th order polynomials yet...

	return 0;
}


