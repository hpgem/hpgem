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

#include "Geometry/PointReference.hpp"
#include <iostream>
#include "cassert"
#include <cmath>
using Geometry::PointReference;

int main(){

	double coord0[]={};
	double coord1[]={1.1};
	double coord2[]={1.2,2.2};
	double coord3[]={1.3,2.3,3.3};
	double coord4[]={1.4,2.4,3.4,4.4};

	LinearAlgebra::NumericalVector vec0(coord0,0),vec1(coord1,1),vec2(coord2,2),vec3(coord3,3),vec4(coord4,4);

	Geometry::Point pb0(0),pb1(1),pb2(2),pb3(3),pb4(4);

	//testing constructors up to DIM=4
	PointReference p0(0),p1(1),p2(2),p3(3),p4(4),pp0(pb0),pp1(pb1),pp2(pb2),pp3(pb3),pp4(pb4),
			pc0(coord0,0),pc1(coord1,1),pc2(coord2,2),pc3(coord3,3),pc4(coord4,4),
			pv0(vec0),pv1(vec1),pv2(vec2),pv3(vec3),pv4(vec4);

	assert(("1D default constructor",p1[0]==0.));
		for(std::size_t i=0;i<2;++i){
			assert(("2D default constructor",p2[i]==0.));
		}
		for(std::size_t i=0;i<3;++i){
			assert(("3D default constructor",p3[i]==0.));
		}
		for(std::size_t i=0;i<4;++i){
			assert(("4D default constructor",p4[i]==0.));
		}

		assert(("1D copy constructor",pp1[0]==0.));
		for(std::size_t i=0;i<2;++i){
			assert(("2D copy constructor",pp2[i]==0.));
		}
		for(std::size_t i=0;i<3;++i){
			assert(("3D copy constructor",pp3[i]==0.));
		}
		for(std::size_t i=0;i<4;++i){
			assert(("4D copy constructor",pp4[i]==0.));
		}

		assert(("1D from array constructor",std::abs(pc1[0]-1.1)<1e-12));
		for(std::size_t i=0;i<2;++i){
			assert(("2D from array constructor",std::abs(pc2[i]-1.2-i)<1e-12));
		}
		for(std::size_t i=0;i<3;++i){
			assert(("3D from array constructor",std::abs(pc3[i]-1.3-i)<1e-12));
		}
		for(std::size_t i=0;i<4;++i){
			assert(("4D from array constructor",std::abs(pc4[i]-1.4-i)<1e-12));
		}

		assert(("1D from NumericalVector constructor",std::abs(pv1[0]-1.1)<1e-12));
		for(std::size_t i=0;i<2;++i){
			assert(("2D from NumericalVector constructor",std::abs(pv2[i]-1.2-i)<1e-12));
		}
		for(std::size_t i=0;i<3;++i){
			assert(("3D from NumericalVector constructor",std::abs(pv3[i]-1.3-i)<1e-12));
		}
		for(std::size_t i=0;i<4;++i){
			assert(("4D from NumericalVector constructor",std::abs(pv4[i]-1.4-i)<1e-12));
		}

		//testing operators

		const PointReference pr0 = pp0 = pv0;
		const PointReference pr1 = pp1 = pv1;
		const PointReference pr2 = pp2 = pv2;
		const PointReference pr3 = pp3 = pv3;
		const PointReference pr4 = pp4 = pv4;


		pv0*6.;
		assert(("1D multiplication",std::abs((pv1*5.)[0]-5.5)<1e-12));
		for(std::size_t i=0;i<2;++i){
			assert(("2D multiplication",std::abs((pv2*4.)[i]-4.8-4*i)<1e-12));
		}
		for(std::size_t i=0;i<3;++i){
			assert(("3D multiplication",std::abs((pv3*3.)[i]-3.9-3*i)<1e-12));
		}
		for(std::size_t i=0;i<4;++i){
			assert(("4D multiplication",std::abs((pv4*2.)[i]-2.8-2*i)<1e-12));
		}

		assert(("0D multiplication",(pr0*0.)==p0));
		assert(("1D multiplication",(pr1*0.)==p1));
		assert(("2D multiplication",(pr2*0.)==p2));
		assert(("3D multiplication",(pr3*0.)==p3));
		assert(("4D multiplication",(pr4*0.)==p4));

		pc0+pv0;
		assert(("1D addition",std::abs((pc1+pv1)[0]-2.2)<1e-12));
		for(std::size_t i=0;i<2;++i){
			assert(("2D addition",std::abs((pc2+pv2)[i]-2.4-2*i)<1e-12));
		}
		for(std::size_t i=0;i<3;++i){
			assert(("3D addition",std::abs((pc3+pv3)[i]-2.6-2*i)<1e-12));
		}
		for(std::size_t i=0;i<4;++i){
			assert(("4D addition",std::abs((pc4+pv4)[i]-2.8-2*i)<1e-12));
		}

		pr0+pv0;
		assert(("1D addition",std::abs((pr1+pv1)[0]-2.2)<1e-12));
		for(std::size_t i=0;i<2;++i){
			assert(("2D addition",std::abs((pr2+pv2)[i]-2.4-2*i)<1e-12));
		}
		for(std::size_t i=0;i<3;++i){
			assert(("3D addition",std::abs((pr3+pv3)[i]-2.6-2*i)<1e-12));
		}
		for(std::size_t i=0;i<4;++i){
			assert(("4D addition",std::abs((pr4+pv4)[i]-2.8-2*i)<1e-12));
		}

		pc0-pv0;
		assert(("1D subtraction",std::abs((pc1-pv1)[0])<1e-12));
		for(std::size_t i=0;i<2;++i){
			assert(("2D subtraction",std::abs((pc2-pv2)[i])<1e-12));
		}
		for(std::size_t i=0;i<3;++i){
			assert(("3D subtraction",std::abs((pc3-pv3)[i])<1e-12));
		}
		for(std::size_t i=0;i<4;++i){
			assert(("4D subtraction",std::abs((pc4-pv4)[i])<1e-12));
		}

		pr0-pv0;
		assert(("1D subtraction",std::abs((pr1-pv1)[0])<1e-12));
		for(std::size_t i=0;i<2;++i){
			assert(("2D subtraction",std::abs((pr2-pv2)[i])<1e-12));
		}
		for(std::size_t i=0;i<3;++i){
			assert(("3D subtraction",std::abs((pr3-pv3)[i])<1e-12));
		}
		for(std::size_t i=0;i<4;++i){
			assert(("4D subtraction",std::abs((pr4-pv4)[i])<1e-12));
		}

		//and point already works so done

}


