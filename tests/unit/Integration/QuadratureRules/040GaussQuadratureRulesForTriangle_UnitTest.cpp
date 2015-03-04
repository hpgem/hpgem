/*
 This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
 http://hpgem.org/about-the-code/team
 
 This code is distributed using BSD 3-Clause License. A copy of which can found below.
 
 
 Copyright (c) 2014, University of Twente
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

#include "Integration/QuadratureRules/GaussQuadratureRulesForTriangle.hpp"
#include <cassert>
#include <typeinfo>

#include "Utilities/BasisFunctions2DH1ConformingTriangle.hpp"
#include "Geometry/ReferenceTriangle.hpp"
#include "Base/BasisFunctionSet.hpp"
#include "Geometry/PointReference.hpp"
#include "LinearAlgebra/NumericalVector.hpp"
#include <cmath>

void testRule(QuadratureRules::GaussQuadratureRule& test,std::size_t expectedOrder){
	std::cout<<test.getName()<<std::endl;
	assert(("dimension",test.dimension()==2));
	assert(("order",test.order()>=expectedOrder));
	assert(("forReferenceGeometry",typeid(*test.forReferenceGeometry())==typeid(Geometry::ReferenceTriangle)));
	Geometry::PointReference point(2);
	std::cout.precision(14);
	Base::BasisFunctionSet* functions = Utilities::createDGBasisFunctionSet2DH1Triangle(expectedOrder);
	for(std::size_t i=0;i<functions->size();++i){
		double integrated=0;
		for(std::size_t j=0;j<test.nrOfPoints();++j){
			point = test.getPoint(j);
			integrated+=test.weight(j)*functions->eval(i,point);
		}
		if(i<3){
			assert(("integration",std::abs(integrated-1./6.)<1e-12));
		}else if(i<6){
			assert(("integration",std::abs(integrated+0.102062072616)<1e-12));
		}else if(i==6){
			assert(("integration",std::abs(integrated-1./20.)<1e-12));
		}else if(11<i&&i<15){
			assert(("integration",std::abs(integrated-0.0129918659263)<1e-12));
		}else if(i==15||i==17){
			assert(("integration",std::abs(integrated+0.010001653302)<1e-12));
		}else if(i==16){
			assert(("integration",std::abs(integrated+0.00396825396825)<1e-12));
		}else if(i==22){
			assert(("integration",std::abs(integrated-0.001467281692237)<1e-12));
		}else if(i==23){
			assert(("integration",std::abs(integrated+0.001467281692237)<1e-12));
		}else if(24<i&&i<28){
			assert(("integration",std::abs(integrated+0.000814308291636)<1e-12));
		}else if(i==28||i==32){
			assert(("integration",std::abs(integrated-0.001994639807826)<1e-12));
		}else if(i==29||i==31){
			assert(("integration",std::abs(integrated-0.001663741054687)<1e-12));
		}else if(i==30){
			assert(("integration",std::abs(integrated-0.0025173611111111)<1e-12));
		}else if(i==37){
			assert(("integration",std::abs(integrated+0.0010300275676522)<1e-12));
		}else if(i==38){
			assert(("integration",std::abs(integrated-0.0000984282481784)<1e-12));
		}else if(i==39){
			assert(("integration",std::abs(integrated+0.0000984282481784)<1e-12));
		}else if(i==40){
			assert(("integration",std::abs(integrated-0.0010300275676522)<1e-12));
		}else if(41<i&&i<45){
			assert(("integration",std::abs(integrated-0.0001528243743039)<1e-12));
		}else if(i==45||i==51){
			assert(("integration",std::abs(integrated+0.0003743417373047)<1e-12));
		}else if(i==46||i==50){
			assert(("integration",std::abs(integrated+0.0045990061560235)<1e-12));
		}else if(i==47||i==49){
			assert(("integration",std::abs(integrated+0.0005242977910035)<1e-12));
		}else if(i==48){
			assert(("integration",std::abs(integrated+0.001020698051948)<1e-12));
		}else if(i==56){
			assert(("integration",std::abs(integrated-0.0004027275873253)<1e-12));
		}else if(i==57){
			assert(("integration",std::abs(integrated+0.001833588494787)<1e-12));
		}else if(i==58){
			assert(("integration",std::abs(integrated-0.0002669631696782)<1e-12));
		}else if(i==59){
			assert(("integration",std::abs(integrated+0.0002669631696782)<1e-12));
		}else if(i==60){
			assert(("integration",std::abs(integrated-0.001833588494787)<1e-12));
		}else if(i==61){
			assert(("integration",std::abs(integrated+0.0004027275873253)<1e-12));
		}else if(62<i&&i<66){
			assert(("integration",std::abs(integrated+0.0000445921151835)<1e-12));
		}else{
			assert(("integration",std::abs(integrated)<1e-12));
		}

	}

	delete functions;
}

int main(){

	testRule(QuadratureRules::Tn2_1_1::Instance(),1);
	testRule(QuadratureRules::Tn2_2_1::Instance(),2);
	testRule(QuadratureRules::Tn2_3_1::Instance(),3);
	testRule(QuadratureRules::Tn2_4_1::Instance(),4);
	testRule(QuadratureRules::T2_5_1::Instance(),5);
	testRule(QuadratureRules::T2_6_1::Instance(),6);
	testRule(QuadratureRules::T2_7_1::Instance(),7);
	testRule(QuadratureRules::T2_8_1::Instance(),8);
	testRule(QuadratureRules::T2_9_1::Instance(),9);
	testRule(QuadratureRules::T2_10_1::Instance(),10);
	testRule(QuadratureRules::T2_11_1::Instance(),10);///\BUG there are no 11th order polynomials yet...

	return 0;
}


