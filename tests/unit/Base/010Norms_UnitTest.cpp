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

#include "Base/L2Norm.hpp"
#include "Base/Norm2.hpp"

#include <cmath>
#include <iostream>
#include "LinearAlgebra/NumericalVector.hpp"
#include "Geometry/PointPhysical.hpp"
#include "Base/Logger.h"
//#include "cassert"

#ifndef LOG_TEST_LEVEL
#define LOG_TEST_LEVEL Log::DEBUG
#endif
Logger<LOG_TEST_LEVEL> logger(std::string(__FILE__).append(" (test suite)"));

void assert(std::string message,bool comparison,int lineNumber){
	if(!comparison){
		logger.log(Log::FATAL,"\n fundamental problem detected for % at line %",message,lineNumber);
	}
}

int main(){

	loggerOutput->onFatal=[](std::string module, std::string message){
        std::cerr <<  "Module: " << module << message << std::endl;
        std::exit(-1);
	};

	double test0[0],test1[1],test2[2],test3[3],test4[4];

	LinearAlgebra::NumericalVector vec0D(test0,0);
	Geometry::PointPhysical point0D(vec0D);

	assert("0D case",Base::L2Norm(vec0D)==0,__LINE__);
	assert("0D case",Base::L2Norm(point0D)==0,__LINE__);
	assert("0D case",Utilities::norm2(point0D)==0,__LINE__);

	test1[0]=1;

	LinearAlgebra::NumericalVector vec1D(test1,1);
	Geometry::PointPhysical point1D(vec1D);

	assert("1D case, positive",fabs(Base::L2Norm(vec1D)-1)<1e-12,__LINE__);
	assert("1D case, positive",fabs(Base::L2Norm(point1D)-1)<1e-12,__LINE__);
	assert("1D case, positive",fabs(Utilities::norm2(point1D)-1)<1e-12,__LINE__);

	vec1D[0]=-1;
	point1D[0]=-1;

	assert("1D case, negative",fabs(Base::L2Norm(vec1D)-1)<1e-12,__LINE__);
	assert("1D case, negative",fabs(Base::L2Norm(point1D)-1)<1e-12,__LINE__);
	assert("1D case, negative",fabs(Utilities::norm2(point1D)-1)<1e-12,__LINE__);

	vec1D[0]=4.38573895783677438;
	point1D[0]=4.38573895783677438;

	assert("non-unit data",fabs(Base::L2Norm(vec1D)-4.38573895783677438)<1e-12,__LINE__);
	assert("non-unit data",fabs(Base::L2Norm(point1D)-4.38573895783677438)<1e-12,__LINE__);
	assert("non-unit data",fabs(Utilities::norm2(point1D)-4.38573895783677438)<1e-12,__LINE__);

	test2[0]=1;
	test2[1]=1;

	LinearAlgebra::NumericalVector vec2D(test2,2);
	Geometry::PointPhysical point2D(vec2D);

	assert("2D case, positive",fabs(Base::L2Norm(vec2D)-sqrt(2.))<1e-12,__LINE__);
	assert("2D case, positive",fabs(Base::L2Norm(point2D)-sqrt(2.))<1e-12,__LINE__);
	assert("2D case, positive",fabs(Utilities::norm2(point2D)-sqrt(2.))<1e-12,__LINE__);

	vec2D[0]=-1;
	point2D[0]=-1;

	assert("2D case, mix",fabs(Base::L2Norm(vec2D)-sqrt(2.))<1e-12,__LINE__);
	assert("2D case, mix",fabs(Base::L2Norm(point2D)-sqrt(2.))<1e-12,__LINE__);
	assert("2D case, mix",fabs(Utilities::norm2(point2D)-sqrt(2.))<1e-12,__LINE__);

	vec2D[1]=-1;
	point2D[1]=-1;

	assert("2D case, negative",fabs(Base::L2Norm(vec2D)-sqrt(2.))<1e-12,__LINE__);
	assert("2D case, negative",fabs(Base::L2Norm(point2D)-sqrt(2.))<1e-12,__LINE__);
	assert("2D case, negative",fabs(Utilities::norm2(point2D)-sqrt(2.))<1e-12,__LINE__);

	test3[0]=1;
	test3[1]=1;
	test3[2]=2;

	LinearAlgebra::NumericalVector vec3D(test3,3);
	Geometry::PointPhysical point3D(vec3D);

	assert("3D case, positive",fabs(Base::L2Norm(vec3D)-sqrt(6.))<1e-12,__LINE__);
	assert("3D case, positive",fabs(Base::L2Norm(point3D)-sqrt(6.))<1e-12,__LINE__);
	assert("3D case, positive",fabs(Utilities::norm2(point3D)-sqrt(6.))<1e-12,__LINE__);

	vec3D[0]=-1;
	point3D[0]=-1;

	assert("3D case, mix",fabs(Base::L2Norm(vec3D)-sqrt(6.))<1e-12,__LINE__);
	assert("3D case, mix",fabs(Base::L2Norm(point3D)-sqrt(6.))<1e-12,__LINE__);
	assert("3D case, mix",fabs(Utilities::norm2(point3D)-sqrt(6.))<1e-12,__LINE__);

	vec3D[1]=-1;
	point3D[1]=-1;

	assert("3D case, mix",fabs(Base::L2Norm(vec3D)-sqrt(6.))<1e-12,__LINE__);
	assert("3D case, mix",fabs(Base::L2Norm(point3D)-sqrt(6.))<1e-12,__LINE__);
	assert("3D case, mix",fabs(Utilities::norm2(point3D)-sqrt(6.))<1e-12,__LINE__);

	vec3D[2]=-2;
	point3D[2]=-2;

	assert("3D case, negative",fabs(Base::L2Norm(vec3D)-sqrt(6.))<1e-12,__LINE__);
	assert("3D case, negative",fabs(Base::L2Norm(point3D)-sqrt(6.))<1e-12,__LINE__);
	assert("3D case, negative",fabs(Utilities::norm2(point3D)-sqrt(6.))<1e-12,__LINE__);

	return 0;
}

