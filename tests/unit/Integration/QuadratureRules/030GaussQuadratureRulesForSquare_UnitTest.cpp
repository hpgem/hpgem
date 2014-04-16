/*
 * 030QaussQuadratureRulesForSquare_UnitTest.cpp
 *
 *  Created on: Apr 14, 2014
 *      Author: brinkf
 */

//naming convention: <Digit><ClassName>_UnitTest.cpp where <Digit> is a number that will make sure
//the unit tests are ordered such that the first failing unit test indicate the culprit class and
//other 'unit' tests may assume correct execution of all prior unit tests

#include "Integration/QuadratureRules/GaussQuadratureRulesForSquare.hpp"
#include "cassert"

#include "Utilities/BasisFunctions2DH1ConformingSquare.hpp"
#include "Geometry/ReferenceSquare.hpp"
#include "Base/BasisFunctionSet.hpp"

void testRule(QuadratureRules::GaussQuadratureRule& test,int expectedOrder){
	cout<<test.getName()<<endl;
	assert(("dimension",test.dimension()==2));
	assert(("order",test.order()>=expectedOrder));
	assert(("forReferenceGeometry",typeid(*test.forReferenceGeometry())==typeid(Geometry::ReferenceSquare)));
	Geometry::PointReference point(2);

	Base::BasisFunctionSet* functions = Utilities::createDGBasisFunctionSet2DH1Square(expectedOrder);

	for(int i=0;i<functions->size();++i){
		double integrated=0;
		for(int j=0;j<test.nrOfPoints();++j){
			test.getPoint(j,point);
			integrated+=test.weight(j)*functions->eval(i,point);
		}
		if(i<4){
			assert(("integration",fabs(integrated-1)<1e-12));
		}else if(i<8){
			assert(("integration",fabs(integrated+sqrt(2./3.))<1e-12));
		}else if(i==8){
			assert(("integration",fabs(integrated-2./3.)<1e-12));
		}else{
			assert(("integration",fabs(integrated)<1e-12));
		}

	}

	delete functions;
}

int main(){

	testRule(QuadratureRules::Cn2_1_1::Instance(),1);
	testRule(QuadratureRules::Cn2_3_4::Instance(),3);
	testRule(QuadratureRules::Cn2_5_9::Instance(),5);
	testRule(QuadratureRules::C2_7_4::Instance(),7);
	testRule(QuadratureRules::C2_9_5::Instance(),9);
	testRule(QuadratureRules::C2_11_6::Instance(),10);///\BUG there are no 11th order polynomials yet...

	return 0;
}


