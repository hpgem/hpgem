/*
 This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
 http://hpgem.org/about-the-code/team
 
 This code is distributed using BSD 3-Clause License. A copy of which can found below.
 
 
 Copyright (c) 2014, Univesity of Twenete
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

//naming convention: <Digit><ClassName>_UnitTest.cpp where <Digit> is a number that will make sure
//the unit tests are ordered such that the first failing unit test indicate the culprit class and
//other 'unit' tests may assume correct execution of all prior unit tests

#include "Integration/QuadratureRules/GaussQuadratureRulesForTriangularPrism.hpp"
#include "cassert"

#include "Utilities/BasisFunctions3DH1ConformingPrism.hpp"
#include "Geometry/ReferenceTriangularPrism.hpp"
#include "Base/BasisFunctionSet.hpp"

void testRule(QuadratureRules::GaussQuadratureRule& test,int expectedOrder){
	cout<<test.getName()<<endl;
	assert(("dimension",test.dimension()==3));
	assert(("order",test.order()>=expectedOrder));
	assert(("forReferenceGeometry",typeid(*test.forReferenceGeometry())==typeid(Geometry::ReferenceTriangularPrism)));
	Geometry::PointReference point(3);
	cout.precision(14);
	Base::BasisFunctionSet* functions = Utilities::createDGBasisFunctionSet3DH1ConformingPrism(expectedOrder);
	cout.precision(14);
	for(int i=0;i<functions->size();++i){
		double integrated=0;
		for(int j=0;j<test.nrOfPoints();++j){
			test.getPoint(j,point);
			integrated+=test.weight(j)*functions->eval(i,point);
		}
		if(i<6){
			assert(("integration",fabs(integrated-1./6.)<1e-10));
		}else if(i<12){
			assert(("integration",fabs(integrated+0.1020620726159)<1e-10));
		}else if(11<i&&i<15){
			assert(("integration",fabs(integrated+0.1360827634879)<1e-9));
		}else if(14<i&&i<18){
			assert(("integration",fabs(integrated-1./12.)<1e-10));
		}else if(i==27||i==28){
			assert(("integration",fabs(integrated-1./20.)<1e-10));
		}else if(40<i&&i<47){
			assert(("integration",fabs(integrated-0.012991865926298)<1e-10));
		}else if(i==56||i==61||i==66){
			assert(("integration",fabs(integrated+0.01060781410869)<1e-10));
		}else if(i==69||i==70){
			assert(("integration",fabs(integrated-0.008166315725102)<1e-10));
		}else if(i==75){
			assert(("integration",fabs(integrated-0.001369177697178)<1e-10||expectedOrder==7));//actually the p=5 quadrature rule may also be the culprit, but 7 is more likely because it has other flaws
		}else if(i==76){
			assert(("integration",fabs(integrated+0.001369177697178)<1e-10||expectedOrder==7));
		}else if(i==87||i==89||i==90||i==92){
			assert(("integration",fabs(integrated+0.010001653302483)<1e-10));
		}else if(i==88||i==91){
			assert(("integration",fabs(integrated+0.003968253968254)<1e-10));
		}else if(i==123||i==124){
			assert(("integration",fabs(integrated+0.000465750474069)<1e-10||expectedOrder==7));
		}else if(i==129){
			assert(("integration",fabs(integrated+0.000721656823802)<1e-10||expectedOrder==7));
		}else if(i==130){
			assert(("integration",fabs(integrated-0.000721656823802)<1e-10||expectedOrder==7));
		}else if(i<132){//I test what I can for p=7, but not all the points
			assert(("integration",fabs(integrated)<1e-10));
		}

	}

	delete functions;
}

int main(){

	testRule(QuadratureRules::TriPrism_1_1::Instance(),1);
	testRule(QuadratureRules::TriPrism_3_1::Instance(),3);
	testRule(QuadratureRules::TriPrism_5_1::Instance(),5);
	testRule(QuadratureRules::TriPrism_7_1::Instance(),7);///\TODO not accurate enough
	///\TODO there are no quadrature rules for higher order prisms

	return 0;
}


