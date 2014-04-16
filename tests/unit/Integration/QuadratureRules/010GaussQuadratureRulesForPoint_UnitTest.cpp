/*
 * 010GaussQuadratureRulesForPoint_unitTest.cpp
 *
 *  Created on: Apr 14, 2014
 *      Author: brinkf
 */

//naming convention: <Digit><ClassName>_UnitTest.cpp where <Digit> is a number that will make sure
//the unit tests are ordered such that the first failing unit test indicate the culprit class and
//other 'unit' tests may assume correct execution of all prior unit tests

#include "Integration/QuadratureRules/GaussQuadratureRulesForPoint.hpp"

void testRule(QuadratureRules::GaussQuadratureRule& test){
	cout<<test.getName();
	assert(("dimension",test.dimension()==0));
	assert(("order",test.order()>11));
	assert(("forReferenceGeometry",typeid(*test.forReferenceGeometry())==typeid(Geometry::ReferencePoint)));
	Geometry::PointReference point(0);
	//0D Quadrature rules are special
	double integrated=0;
	for(int i=0;i<test.nrOfPoints();++i){
		integrated+=test.weight(i);
		test.getPoint(i,point);
	}
	assert(("integration",fabs(integrated-1)<1e-12));
}

int main(){

	testRule(QuadratureRules::Cn0_inf_1::Instance());

	return 0;
}

