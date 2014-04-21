/*
 * 050GaussQuadratureRuleForCube_UnitTest.cpp
 *
 *  Created on: Apr 14, 2014
 *      Author: brinkf
 */

//naming convention: <Digit><ClassName>_UnitTest.cpp where <Digit> is a number that will make sure
//the unit tests are ordered such that the first failing unit test indicate the culprit class and
//other 'unit' tests may assume correct execution of all prior unit tests

#include "Integration/QuadratureRules/GaussQuadratureRulesForCube.hpp"
#include "cassert"

#include "Utilities/BasisFunctions3DH1ConformingCube.hpp"
#include "Geometry/ReferenceCube.hpp"
#include "Base/BasisFunctionSet.hpp"

void testRule(QuadratureRules::GaussQuadratureRule& test,int expectedOrder){
	cout<<test.getName()<<endl;
	assert(("dimension",test.dimension()==3));
	assert(("order",test.order()>=expectedOrder));
	assert(("forReferenceGeometry",typeid(*test.forReferenceGeometry())==typeid(Geometry::ReferenceCube)));
	Geometry::PointReference point(3);

	Base::BasisFunctionSet* functions = Utilities::createDGBasisFunctionSet3DH1Cube(expectedOrder);
	cout.precision(14);
	for(int i=0;i<functions->size();++i){
		double integrated=0;
		for(int j=0;j<test.nrOfPoints();++j){
			test.getPoint(j,point);
			integrated+=test.weight(j)*functions->eval(i,point);
		}
		if(i<8){
			assert(("integration",fabs(integrated-1)<1e-12));
		}else if(i<20){
			assert(("integration",fabs(integrated+sqrt(2./3.))<1e-12));
		}else if(i<26){
			assert(("integration",fabs(integrated-2./3.)<1e-12));
		}else if(i==26){
			assert(("integration",fabs(integrated+sqrt(2./3.)*2./3.)<1e-12));
		}else{
			assert(("integration",fabs(integrated)<1e-12));
		}

	}

	delete functions;
}

int main(){

	testRule(QuadratureRules::Cn3_1_1::Instance(),1);
	testRule(QuadratureRules::Cn3_3_4::Instance(),3);
	testRule(QuadratureRules::Cn3_5_9::Instance(),5);
	testRule(QuadratureRules::C3_7_2::Instance(),7);
	testRule(QuadratureRules::C3_9_2::Instance(),9);
	testRule(QuadratureRules::C3_11_2::Instance(),10);///\BUG there are no 11th order polynomials yet...

	return 0;
}


