/*
 * 010Norms_UnitTest.cpp
 *
 *  Created on: Apr 15, 2014
 *      Author: brinkf
 */

//naming convention: <Digit><ClassName>_UnitTest.cpp where <Digit> is a number that will make sure
//the unit tests are ordered such that the first failing unit test indicate the culprit class and
//other 'unit' tests may assume correct execution of all prior unit tests

#include "Base/L2Norm.hpp"
#include "Base/Norm2.hpp"
#include "cassert"

int main(){

	double test0[0],test1[1],test2[2],test3[3],test4[4];

	LinearAlgebra::NumericalVector vec0D(test0,0);
	Geometry::PointPhysical point0D(vec0D);

	assert(("0D case",Base::L2Norm(vec0D)==0));
	assert(("0D case",Base::L2Norm(point0D)==0));
	assert(("0D case",Utilities::norm2(point0D)==0));

	test1[0]=1;

	LinearAlgebra::NumericalVector vec1D(test1,1);
	Geometry::PointPhysical point1D(vec1D);

	assert(("1D case, positive",fabs(Base::L2Norm(vec1D)-1)<1e-12));
	assert(("1D case, positive",fabs(Base::L2Norm(point1D)-1)<1e-12));
	assert(("1D case, positive",fabs(Utilities::norm2(point1D)-1)<1e-12));

	vec1D[0]=-1;
	point1D[0]=-1;

	assert(("1D case, negative",fabs(Base::L2Norm(vec1D)-1)<1e-12));
	assert(("1D case, negative",fabs(Base::L2Norm(point1D)-1)<1e-12));
	assert(("1D case, negative",fabs(Utilities::norm2(point1D)-1)<1e-12));

	vec1D[0]=4.38573895783677438;
	point1D[0]=4.38573895783677438;

	assert(("non-unit data",fabs(Base::L2Norm(vec1D)-4.38573895783677438)<1e-12));
	assert(("non-unit data",fabs(Base::L2Norm(point1D)-4.38573895783677438)<1e-12));
	assert(("non-unit data",fabs(Utilities::norm2(point1D)-4.38573895783677438)<1e-12));

	test2[0]=1;
	test2[1]=1;

	LinearAlgebra::NumericalVector vec2D(test2,2);
	Geometry::PointPhysical point2D(vec2D);

	assert(("2D case, positive",fabs(Base::L2Norm(vec2D)-sqrt(2.))<1e-12));
	assert(("2D case, positive",fabs(Base::L2Norm(point2D)-sqrt(2.))<1e-12));
	assert(("2D case, positive",fabs(Utilities::norm2(point2D)-sqrt(2.))<1e-12));

	vec2D[0]=-1;
	point2D[0]=-1;

	assert(("2D case, mix",fabs(Base::L2Norm(vec2D)-sqrt(2.))<1e-12));
	assert(("2D case, mix",fabs(Base::L2Norm(point2D)-sqrt(2.))<1e-12));
	assert(("2D case, mix",fabs(Utilities::norm2(point2D)-sqrt(2.))<1e-12));

	vec2D[1]=-1;
	point2D[1]=-1;

	assert(("2D case, negative",fabs(Base::L2Norm(vec2D)-sqrt(2.))<1e-12));
	assert(("2D case, negative",fabs(Base::L2Norm(point2D)-sqrt(2.))<1e-12));
	assert(("2D case, negative",fabs(Utilities::norm2(point2D)-sqrt(2.))<1e-12));

	test3[0]=1;
	test3[1]=1;
	test3[2]=2;

	LinearAlgebra::NumericalVector vec3D(test3,3);
	Geometry::PointPhysical point3D(vec3D);

	assert(("3D case, positive",fabs(Base::L2Norm(vec3D)-sqrt(6.))<1e-12));
	assert(("3D case, positive",fabs(Base::L2Norm(point3D)-sqrt(6.))<1e-12));
	assert(("3D case, positive",fabs(Utilities::norm2(point3D)-sqrt(6.))<1e-12));

	vec3D[0]=-1;
	point3D[0]=-1;

	assert(("3D case, mix",fabs(Base::L2Norm(vec3D)-sqrt(6.))<1e-12));
	assert(("3D case, mix",fabs(Base::L2Norm(point3D)-sqrt(6.))<1e-12));
	assert(("3D case, mix",fabs(Utilities::norm2(point3D)-sqrt(6.))<1e-12));

	vec3D[1]=-1;
	point3D[1]=-1;

	assert(("3D case, mix",fabs(Base::L2Norm(vec3D)-sqrt(6.))<1e-12));
	assert(("3D case, mix",fabs(Base::L2Norm(point3D)-sqrt(6.))<1e-12));
	assert(("3D case, mix",fabs(Utilities::norm2(point3D)-sqrt(6.))<1e-12));

	vec3D[2]=-2;
	point3D[2]=-2;

	assert(("3D case, negative",fabs(Base::L2Norm(vec3D)-sqrt(6.))<1e-12));
	assert(("3D case, negative",fabs(Base::L2Norm(point3D)-sqrt(6.))<1e-12));
	assert(("3D case, negative",fabs(Utilities::norm2(point3D)-sqrt(6.))<1e-12));

	return 0;
}

