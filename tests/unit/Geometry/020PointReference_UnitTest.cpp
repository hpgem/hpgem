/*
 * 020PointReference_UnitTest.cpp
 *
 *  Created on: Mar 26, 2014
 *      Author: brinkf
 */

//naming convention: <Digit><ClassName>_UnitTest.cpp where <Digit> is a number that will make sure
//the unit tests are ordered such that the first failing unit test indicate the culprit class and
//other 'unit' tests may assume correct execution of all prior unit tests

#include "Geometry/PointReference.hpp"
#include <iostream>
#include "cassert"
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
		for(int i=0;i<2;++i){
			assert(("2D default constructor",p2[i]==0.));
		}
		for(int i=0;i<3;++i){
			assert(("3D default constructor",p3[i]==0.));
		}
		for(int i=0;i<4;++i){
			assert(("4D default constructor",p4[i]==0.));
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

		assert(("1D from array constructor",fabs(pc1[0]-1.1)<1e-12));
		for(int i=0;i<2;++i){
			assert(("2D from array constructor",fabs(pc2[i]-1.2-i)<1e-12));
		}
		for(int i=0;i<3;++i){
			assert(("3D from array constructor",fabs(pc3[i]-1.3-i)<1e-12));
		}
		for(int i=0;i<4;++i){
			assert(("4D from array constructor",fabs(pc4[i]-1.4-i)<1e-12));
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

		//testing operators

		const PointReference pr0 = pp0 = pv0;
		const PointReference pr1 = pp1 = pv1;
		const PointReference pr2 = pp2 = pv2;
		const PointReference pr3 = pp3 = pv3;
		const PointReference pr4 = pp4 = pv4;


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

		assert(("0D multiplication",(pr0*0.)==p0));
		assert(("1D multiplication",(pr1*0.)==p1));
		assert(("2D multiplication",(pr2*0.)==p2));
		assert(("3D multiplication",(pr3*0.)==p3));
		assert(("4D multiplication",(pr4*0.)==p4));

		pc0+pv0;
		assert(("1D addition",fabs((pc1+pv1)[0]-2.2)<1e-12));
		for(int i=0;i<2;++i){
			assert(("2D addition",fabs((pc2+pv2)[i]-2.4-2*i)<1e-12));
		}
		for(int i=0;i<3;++i){
			assert(("3D addition",fabs((pc3+pv3)[i]-2.6-2*i)<1e-12));
		}
		for(int i=0;i<4;++i){
			assert(("4D addition",fabs((pc4+pv4)[i]-2.8-2*i)<1e-12));
		}

		pr0+pv0;
		assert(("1D addition",fabs((pr1+pv1)[0]-2.2)<1e-12));
		for(int i=0;i<2;++i){
			assert(("2D addition",fabs((pr2+pv2)[i]-2.4-2*i)<1e-12));
		}
		for(int i=0;i<3;++i){
			assert(("3D addition",fabs((pr3+pv3)[i]-2.6-2*i)<1e-12));
		}
		for(int i=0;i<4;++i){
			assert(("4D addition",fabs((pr4+pv4)[i]-2.8-2*i)<1e-12));
		}

		pc0-pv0;
		assert(("1D subtraction",fabs((pc1-pv1)[0])<1e-12));
		for(int i=0;i<2;++i){
			assert(("2D subtraction",fabs((pc2-pv2)[i])<1e-12));
		}
		for(int i=0;i<3;++i){
			assert(("3D subtraction",fabs((pc3-pv3)[i])<1e-12));
		}
		for(int i=0;i<4;++i){
			assert(("4D subtraction",fabs((pc4-pv4)[i])<1e-12));
		}

		pr0-pv0;
		assert(("1D subtraction",fabs((pr1-pv1)[0])<1e-12));
		for(int i=0;i<2;++i){
			assert(("2D subtraction",fabs((pr2-pv2)[i])<1e-12));
		}
		for(int i=0;i<3;++i){
			assert(("3D subtraction",fabs((pr3-pv3)[i])<1e-12));
		}
		for(int i=0;i<4;++i){
			assert(("4D subtraction",fabs((pr4-pv4)[i])<1e-12));
		}

		//and point already works so done

}


