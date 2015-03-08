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

#include "Geometry/Point.hpp"
#include <iostream>
#include "cassert"
#include <cmath>
using std::fabs;
using Geometry::Point;

int main(){

	double coord0[]={};
	double coord1[]={1.1};
	double coord2[]={1.2,2.2};
	double coord3[]={1.3,2.3,3.3};
	double coord4[]={1.4,2.4,3.4,4.4};

	LinearAlgebra::NumericalVector vec0(coord0,0),vec1(coord1,1),vec2(coord2,2),vec3(coord3,3),vec4(coord4,4);

	//testing constructors up to DIM=4
	Point p0(0),p1(1),p2(2),p3(3),p4(4),pc0(coord0,0),pc1(coord1,1),pc2(coord2,2),pc3(coord3,3),pc4(coord4,4),pp0(p0),pp1(p1),pp2(p2),pp3(p3),pp4(p4),pv0(vec0),pv1(vec1),pv2(vec2),pv3(vec3),pv4(vec4);

	//testing operator[]

	assert(("1D default constructor or operator[] of Point",p1[0]==0.));
	for(int i=0;i<2;++i){
		assert(("2D default constructor or operator[] of Point",p2[i]==0.));
	}
	for(int i=0;i<3;++i){
		assert(("3D default constructor or operator[] of Point",p3[i]==0.));
	}
	for(int i=0;i<4;++i){
		assert(("4D default constructor or operator[] of Point",p4[i]==0.));
	}

	assert(("1D copy constructor",pp1[0]==0.));
	for(int i=0;i<2;++i){
		assert(("2D copy constructor",pp2[i]==0.));
	}
	for(int i=0;i<3;++i){
		assert(("3D copy constructor",pp3[i]==0.));
	}
	for(int i=0;i<4;++i){
		assert(("4D copy constructor",pp4[i]==0.));
	}

	assert(("1D from array constructor or operator[] of Point",fabs(pc1[0]-1.1)<1e-12));
	for(int i=0;i<2;++i){
		assert(("2D from array constructor or operator[] of Point",fabs(pc2[i]-1.2-i)<1e-12));
	}
	for(int i=0;i<3;++i){
		assert(("3D from array constructor or operator[] of Point",fabs(pc3[i]-1.3-i)<1e-12));
	}
	for(int i=0;i<4;++i){
		assert(("4D from array constructor or operator[] of Point",fabs(pc4[i]-1.4-i)<1e-12));
	}

	assert(("1D from NumericalVector constructor",fabs(pv1[0]-1.1)<1e-12));
	for(int i=0;i<2;++i){
		assert(("2D from NumericalVector constructor",fabs(pv2[i]-1.2-i)<1e-12));
	}
	for(int i=0;i<3;++i){
		assert(("3D from NumericalVector constructor",fabs(pv3[i]-1.3-i)<1e-12));
	}
	for(int i=0;i<4;++i){
		assert(("4D from NumericalVector constructor",fabs(pv4[i]-1.4-i)<1e-12));
	}

	//testing setCoordinates and setCoordinate

	p1.setCoordinates(vec1);
	assert(("1D setCoordinates",fabs(p1[0]-1.1)<1e-12));
	p2.setCoordinates(vec2);
	for(int i=0;i<2;++i){
		assert(("2D setCoordinates",fabs(p2[i]-1.2-i)<1e-12));
	}
	p3.setCoordinates(vec3);
	for(int i=0;i<3;++i){
		assert(("3D setCoordinates",fabs(p3[i]-1.3-i)<1e-12));
	}
	p4.setCoordinates(vec4);
	for(int i=0;i<4;++i){
		assert(("4D setCoordinates",fabs(p4[i]-1.4-i)<1e-12));
	}

	p1.setCoordinate(0,0.9);
	assert(("1D setCoordinate",fabs(p1[0]-0.9)<1e-12));
	for(int i=0;i<2;++i){
		p2.setCoordinate(i,0.8+double(i));
		for(int j=0;j<=i;++j){
			assert(("2D setCoordinate",fabs(p2[j]-0.8-j)<1e-12));
		}
		for(int j=i+1;j<2;++j){
			assert(("2D setCoordinate",fabs(p2[j]-1.2-j)<1e-12));
		}
	}
	for(int i=0;i<3;++i){
		p3.setCoordinate(i,0.7+i);
		for(int j=0;j<=i;++j){
			assert(("3D setCoordinate",fabs(p3[j]-0.7-j)<1e-12));
		}
		for(int j=i+1;j<3;++j){
			assert(("3D setCoordinate",fabs(p3[j]-1.3-j)<1e-12));
		}
	}
	for(int i=0;i<4;++i){
		p4.setCoordinate(i,0.6+i);
		for(int j=0;j<=i;++j){
			assert(("4D setCoordinate",fabs(p4[j]-0.6-j)<1e-12));
		}
		for(int j=i+1;j<4;++j){
			assert(("4D setCoordinate",fabs(p4[j]-1.4-j)<1e-12));
		}
	}

	//testing operators

	const Point pr0 = pc0 = p0;
	const Point pr1 = pc1 = p1;
	assert(("1D assignment operator",fabs(pc1[0]-0.9)<1e-12));
	assert(("1D assignment operator",fabs(pr1[0]-0.9)<1e-12));
	const Point pr2 = pc2 = p2;
	for(int i=0;i<2;++i){
		assert(("2D assignment operator",fabs(pc2[i]-0.8-i)<1e-12));
		assert(("2D assignment operator",fabs(pr2[i]-0.8-i)<1e-12));
	}
	const Point pr3 = pc3 = p3;
	for(int i=0;i<3;++i){
		assert(("3D assignment operator",fabs(pc3[i]-0.7-i)<1e-12));
		assert(("3D assignment operator",fabs(pr3[i]-0.7-i)<1e-12));
	}
	const Point pr4 = pc4 = p4;
	for(int i=0;i<4;++i){
		assert(("4D assignment operator",fabs(pc4[i]-0.6-i)<1e-12));
		assert(("4D assignment operator",fabs(pr4[i]-0.6-i)<1e-12));
	}

	assert(("0D equality operator",pr0==pc0&&pc0==pr0&&pc0==p0));
	assert(("1D equality operator",pr1==pc1&&pc1==pr1&&pc1==p1&&!(pr1==pv1||pv1==pr1||p1==pv1)));
	assert(("2D equality operator",pr2==pc2&&pc2==pr2&&pc2==p2&&!(pr2==pv2||pv2==pr2||p2==pv2)));
	assert(("3D equality operator",pr3==pc3&&pc3==pr3&&pc3==p3&&!(pr3==pv3||pv3==pr3||p3==pv3)));
	assert(("4D equality operator",pr4==pc4&&pc4==pr4&&pc4==p4&&!(pr4==pv4||pv4==pr4||p4==pv4)));
	assert(("equality operator - different dimensions",
			!(pr0==pv1||pv1==pr0||p0==p1||pr0==pv2||pv2==pr0||p0==p2||
					pr0==pv3||pv3==pr0||p0==p3||pr0==pv4||pv4==pr0||p0==p4||
					pr1==pv2||pv2==pr1||p1==p2||
					pr1==pv3||pv3==pr1||p1==p3||pr1==pv4||pv4==pr1||p1==p4||
					pr2==pv3||pv3==pr2||p2==p3||pr2==pv4||pv4==pr2||p2==p4||
					pr3==pv4||pv4==pr3||p3==p4)));

	pc0+=p0;
	pc1+=p1;
	assert(("1D increment operator",fabs(pc1[0]-1.8)<1e-12));
	pc2+=p2;
	for(int i=0;i<2;++i){
		assert(("2D increment operator",fabs(pc2[i]-1.6-2*i)<1e-12));
	}
	pc3+=p3;
	for(int i=0;i<3;++i){
		assert(("3D increment operator",fabs(pc3[i]-1.4-2*i)<1e-12));
	}
	pc4+=p4;
	for(int i=0;i<4;++i){
		assert(("4D increment operator",fabs(pc4[i]-1.2-2*i)<1e-12));
	}

	pc0-=pv0;
	pc1-=pv1;
	assert(("1D decrement operator",fabs(pc1[0]-0.7)<1e-12));
	pc2-=pv2;
	for(int i=0;i<2;++i){
		assert(("2D decrement operator",fabs(pc2[i]-0.4-i)<1e-12));
	}
	pc3-=pv3;
	for(int i=0;i<3;++i){
		assert(("3D decrement operator",fabs(pc3[i]-0.1-i)<1e-12));
	}
	pc4-=pv4;
	for(int i=0;i<4;++i){
		assert(("4D decrement operator or negative numbers",fabs(pc4[i]+0.2-i)<1e-12));
	}

	pc0*=2.;
	pc1*=3.;
	assert(("1D multiply operator",fabs(pc1[0]-2.1)<1e-12));
	pc2*=4.;
	for(int i=0;i<2;++i){
		assert(("2D multiply operator",fabs(pc2[i]-1.6-4*i)<1e-12));
	}
	pc3*=5.;
	for(int i=0;i<3;++i){
		assert(("3D multiply operator",fabs(pc3[i]-0.5-5*i)<1e-12));
	}
	pc4*=6.;
	for(int i=0;i<4;++i){
		assert(("4D multiply operator",fabs(pc4[i]+1.2-6*i)<1e-12));
	}

	pv0*6.;
	assert(("1D multiplication",fabs((pv1*5.)[0]-5.5)<1e-12));
	for(int i=0;i<2;++i){
		assert(("2D multiplication",fabs((pv2*4.)[i]-4.8-4*i)<1e-12));
	}
	for(int i=0;i<3;++i){
		assert(("3D multiplication",fabs((pv3*3.)[i]-3.9-3*i)<1e-12));
	}
	for(int i=0;i<4;++i){
		assert(("4D multiplication",fabs((pv4*2.)[i]-2.8-2*i)<1e-12));
	}

	assert(("0D multiplication",(pr0*0.)==pp0));
	assert(("1D multiplication",(pr1*0.)==pp1));
	assert(("2D multiplication",(pr2*0.)==pp2));
	assert(("3D multiplication",(pr3*0.)==pp3));
	assert(("4D multiplication",(pr4*0.)==pp4));

	pc0+pv0;
	assert(("1D addition",fabs((pc1+pv1)[0]-3.2)<1e-12));
	for(int i=0;i<2;++i){
		assert(("2D addition",fabs((pc2+pv2)[i]-2.8-5*i)<1e-12));
	}
	for(int i=0;i<3;++i){
		assert(("3D addition",fabs((pc3+pv3)[i]-1.8-6*i)<1e-12));
	}
	for(int i=0;i<4;++i){
		assert(("4D addition",fabs((pc4+pv4)[i]-0.2-7*i)<1e-12));
	}

	pr0+pv0;
	assert(("1D addition",fabs((pr1+pv1)[0]-2.)<1e-12));
	for(int i=0;i<2;++i){
		assert(("2D addition",fabs((pr2+pv2)[i]-2.-2*i)<1e-12));
	}
	for(int i=0;i<3;++i){
		assert(("3D addition",fabs((pr3+pv3)[i]-2.-2*i)<1e-12));
	}
	for(int i=0;i<4;++i){
		assert(("4D addition",fabs((pr4+pv4)[i]-2.-2*i)<1e-12));
	}

	pc0-pv0;
	assert(("1D subtraction",fabs((pc1-pv1)[0]-1.)<1e-12));
	for(int i=0;i<2;++i){
		assert(("2D subtraction",fabs((pc2-pv2)[i]-0.4-3*i)<1e-12));
	}
	for(int i=0;i<3;++i){
		assert(("3D subtraction",fabs((pc3-pv3)[i]+0.8-4*i)<1e-12));
	}
	for(int i=0;i<4;++i){
		assert(("4D subtraction",fabs((pc4-pv4)[i]+2.6-5*i)<1e-12));
	}

	pr0-pv0;
	assert(("1D subtraction",fabs((pr1-pv1)[0]+0.2)<1e-12));
	for(int i=0;i<2;++i){
		assert(("2D subtraction",fabs((pr2-pv2)[i]+0.4)<1e-12));
	}
	for(int i=0;i<3;++i){
		assert(("3D subtraction",fabs((pr3-pv3)[i]+0.6)<1e-12));
	}
	for(int i=0;i<4;++i){
		assert(("4D subtraction",fabs((pr4-pv4)[i]+0.8)<1e-12));
	}

	//testing size

	assert(("size of a 0D point",p0.size()==0&&pp0.size()==0&&pr0.size()==0&&pc0.size()==0&&pv0.size()==0));
	assert(("size of a 1D point",p1.size()==1&&pp1.size()==1&&pr1.size()==1&&pc1.size()==1&&pv1.size()==1));
	assert(("size of a 2D point",p2.size()==2&&pp2.size()==2&&pr2.size()==2&&pc2.size()==2&&pv2.size()==2));
	assert(("size of a 3D point",p3.size()==3&&pp3.size()==3&&pr3.size()==3&&pc3.size()==3&&pv3.size()==3));
	assert(("size of a 4D point",p4.size()==4&&pp4.size()==4&&pr4.size()==4&&pc4.size()==4&&pv4.size()==4));

	//testing getCoordinate and getCoordinates

	assert(("1D getCoordinate",p1.getCoordinate(0)==p1[0]));
	for(int i=0;i<2;++i){
		assert(("2D getCoordinate",p2.getCoordinate(i)==p2[i]));
	}
	for(int i=0;i<3;++i){
		assert(("3D getCoordinate",p3.getCoordinate(i)==p3[i]));
	}
	for(int i=0;i<4;++i){
		assert(("4D getCoordinate",p4.getCoordinate(i)==p4[i]));
	}

	assert(("0D getCoordinates",pv0.getCoordinates()==vec0));
	assert(("1D getCoordinates",pv1.getCoordinates()==vec1));
	assert(("2D getCoordinates",pv2.getCoordinates()==vec2));
	assert(("3D getCoordinates",pv3.getCoordinates()==vec3));
	assert(("4D getCoordinates",pv4.getCoordinates()==vec4));

	//testing friends

	-pc0;
	assert(("1D unary -",fabs((-pc1)[0]+2.1)<1e-12));
	for(int i=0;i<2;++i){
		assert(("2D unary -",fabs((-pc2)[i]+1.6+4*i)<1e-12));
	}
	for(int i=0;i<3;++i){
		assert(("3D unary -",fabs((-pc3)[i]+0.5+5*i)<1e-12));
	}
	for(int i=0;i<4;++i){
		assert(("4D unary -",fabs((-pc4)[i]-1.2+6*i)<1e-12));
	}

	6.*pv0;
	assert(("1D left multiplication",fabs((5.*pv1)[0]-5.5)<1e-12));
	for(int i=0;i<2;++i){
		assert(("2D left multiplication",fabs((4.*pv2)[i]-4.8-4*i)<1e-12));
	}
	for(int i=0;i<3;++i){
		assert(("3D left multiplication",fabs((3.*pv3)[i]-3.9-3*i)<1e-12));
	}
	for(int i=0;i<4;++i){
		assert(("4D left multiplication",fabs((2.*pv4)[i]-2.8-2*i)<1e-12));
	}

	std::cout<<p0<<p1<<p2<<p3<<p4<<pr0<<pr1<<pr2<<pr3<<pr4<<std::endl;

	return 0;
}
