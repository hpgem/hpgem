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

#include "Geometry/ReferenceTriangle.hpp"
#include "Geometry/ReferenceLine.hpp"
#include "Geometry/ReferencePoint.hpp"
#include <iostream>
#include "cassert"

#include "Geometry/PointReference.hpp"
#include "Geometry/Mappings/MappingToRefLineToTriangle.hpp"
#include "Geometry/Mappings/MappingToRefTriangleToTriangle.hpp"
#include "Integration/QuadratureRules/GaussQuadratureRule.hpp"
#include <cmath>
using std::fabs;
using Geometry::ReferenceTriangle;

int main(){
	ReferenceTriangle& test=ReferenceTriangle::Instance();

	Geometry::PointReference pTest(2);

	//testing basic functionality

	for(pTest[0]=-3.141;pTest[0]<0;pTest[0]+=0.1){
		for(pTest[1]=-3.1416;pTest[1]<3.1416;pTest[1]+=0.1){
			assert(("isInternalPoint",!test.isInternalPoint(pTest)));
		}
	}
	for(;pTest[0]<1;pTest[0]+=0.1){
		for(pTest[1]=-3.1417;pTest[1]<0;pTest[1]+=0.1){
			assert(("isInternalPoint",!test.isInternalPoint(pTest)));
		}
		for(;pTest[1]<1-pTest[0];pTest[1]+=0.1){
			assert(("isInternalPoint",test.isInternalPoint(pTest)));
		}
		for(;pTest[1]<3.141;pTest[1]+=0.1){
			assert(("isInternalPoint",!test.isInternalPoint(pTest)));
		}
	}
	for(;pTest[0]<3.141;pTest[0]+=0.1){
		for(pTest[1]=-3.1416;pTest[1]<3.1416;pTest[1]+=0.1){
			assert(("isInternalPoint",!test.isInternalPoint(pTest)));
		}
	}

	test.getCenter(pTest);
	assert(("getCenter",test.isInternalPoint(pTest)&&fabs(pTest[0]-1./3.)<1e-12&&fabs(pTest[1]-1./3.)<1e-12));
	test.getNode(0,pTest);
	assert(("getNode 0",fabs(pTest[0])<1e-12&&fabs(pTest[1])<1e-12));
	test.getNode(1,pTest);
	assert(("getNode 1",fabs(pTest[0]-1)<1e-12&&fabs(pTest[1])<1e-12));
	test.getNode(2,pTest);
	assert(("getNode 2",fabs(pTest[0])<1e-12&&fabs(pTest[1]-1)<1e-12));
	std::cout<<test.getName();

	assert(("getLocalNodeIndex 0",test.getLocalNodeIndex(0,0)==0));//the nodes of the face must always be specified IN THIS SPECIFIC ORDER
	assert(("getLocalNodeIndex 0",test.getLocalNodeIndex(0,1)==1));//im not sure if I like this myself, but this should at least verify
	assert(("getLocalNodeIndex 1",test.getLocalNodeIndex(1,0)==0));//that all face nodes are specified, none are specified twice
	assert(("getLocalNodeIndex 1",test.getLocalNodeIndex(1,1)==2));//and only face nodes are specified and the ordering of the nodes is consistent
	assert(("getLocalNodeIndex 2",test.getLocalNodeIndex(2,0)==1));//across function calls
	assert(("getLocalNodeIndex 2",test.getLocalNodeIndex(2,1)==2));

	std::cout<<test;

	//testing mappings and quadrature rules

	std::vector<unsigned int> base(3),transformed(3),faceIndices(2);
	for(int i=0;i<3;++i){//doesnt test against reordering of the nodes in the first vector
		base[i]=transformed[i]=i;
	}
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(test.getCodim0MappingIndex(base,transformed))==&Geometry::MappingToRefTriangleToTriangle0::Instance()));
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(base,transformed)==&Geometry::MappingToRefTriangleToTriangle0::Instance()));
	transformed[0]=0;
	transformed[1]=2;
	transformed[2]=1;
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(test.getCodim0MappingIndex(base,transformed))==&Geometry::MappingToRefTriangleToTriangle1::Instance()));
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(base,transformed)==&Geometry::MappingToRefTriangleToTriangle1::Instance()));
	transformed[0]=1;
	transformed[1]=2;
	transformed[2]=0;
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(test.getCodim0MappingIndex(base,transformed))==&Geometry::MappingToRefTriangleToTriangle2::Instance()));
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(base,transformed)==&Geometry::MappingToRefTriangleToTriangle2::Instance()));
	transformed[0]=1;
	transformed[1]=0;
	transformed[2]=2;
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(test.getCodim0MappingIndex(base,transformed))==&Geometry::MappingToRefTriangleToTriangle3::Instance()));
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(base,transformed)==&Geometry::MappingToRefTriangleToTriangle3::Instance()));
	transformed[0]=2;
	transformed[1]=1;
	transformed[2]=0;
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(test.getCodim0MappingIndex(base,transformed))==&Geometry::MappingToRefTriangleToTriangle4::Instance()));
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(base,transformed)==&Geometry::MappingToRefTriangleToTriangle4::Instance()));
	transformed[0]=2;
	transformed[1]=0;
	transformed[2]=1;
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(test.getCodim0MappingIndex(base,transformed))==&Geometry::MappingToRefTriangleToTriangle5::Instance()));
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(base,transformed)==&Geometry::MappingToRefTriangleToTriangle5::Instance()));


	assert(("higher codimensional entities",test.getNrOfCodim1Entities()==3&&test.getNrOfCodim2Entities()==3)&&test.getNrOfCodim3Entities()==0);
	assert(("getCodim1ReferenceGeometry",test.getCodim1ReferenceGeometry(0)==&Geometry::ReferenceLine::Instance()&&
										 test.getCodim1ReferenceGeometry(1)==&Geometry::ReferenceLine::Instance()&&
										 test.getCodim1ReferenceGeometry(2)==&Geometry::ReferenceLine::Instance()));
	assert(("getCodim2ReferenceGeometry",test.getCodim2ReferenceGeometry(0)==&Geometry::ReferencePoint::Instance()&&
										 test.getCodim2ReferenceGeometry(1)==&Geometry::ReferencePoint::Instance()&&
										 test.getCodim2ReferenceGeometry(2)==&Geometry::ReferencePoint::Instance()));
	assert(("getCodim1MappingPtr",test.getCodim1MappingPtr(0)==&Geometry::MappingToRefLineToTriangle0::Instance()));
	assert(("getCodim1MappingPtr",test.getCodim1MappingPtr(1)==&Geometry::MappingToRefLineToTriangle1::Instance()));
	assert(("getCodim1MappingPtr",test.getCodim1MappingPtr(2)==&Geometry::MappingToRefLineToTriangle2::Instance()));
	test.getCodim1EntityLocalIndices(0,faceIndices);
	assert(("getCodim1EntityLocalIndices",faceIndices[0]==test.getLocalNodeIndex(0,0)));
	assert(("getCodim1EntityLocalIndices",faceIndices[1]==test.getLocalNodeIndex(0,1)));
	test.getCodim1EntityLocalIndices(1,faceIndices);
	assert(("getCodim1EntityLocalIndices",faceIndices[0]==test.getLocalNodeIndex(1,0)));
	assert(("getCodim1EntityLocalIndices",faceIndices[1]==test.getLocalNodeIndex(1,1)));
	test.getCodim1EntityLocalIndices(2,faceIndices);
	assert(("getCodim1EntityLocalIndices",faceIndices[0]==test.getLocalNodeIndex(2,0)));
	assert(("getCodim1EntityLocalIndices",faceIndices[1]==test.getLocalNodeIndex(2,1)));
	faceIndices.resize(1);
	test.getCodim2EntityLocalIndices(0,faceIndices);
	assert(("getCodim2EntityLocalIndices",faceIndices[0]==0));
	test.getCodim2EntityLocalIndices(1,faceIndices);
	assert(("getCodim2EntityLocalIndices",faceIndices[0]==1));
	test.getCodim2EntityLocalIndices(2,faceIndices);
	assert(("getCodim2EntityLocalIndices",faceIndices[0]==2));

	assert(("quadrature rules",test.getGaussQuadratureRule(3)->order()>=3));
	assert(("quadrature rules",test.getGaussQuadratureRule(5)->order()>=5));
	assert(("quadrature rules",test.getGaussQuadratureRule(7)->order()>=7));
	assert(("quadrature rules",test.getGaussQuadratureRule(9)->order()>=9));
	assert(("quadrature rules",test.getGaussQuadratureRule(11)->order()>=11));

	//testing functionality of abstract parent classes

	assert(("number of nodes",test.getNumberOfNodes()==3));
	assert(("type of geometry",test.getGeometryType()==Geometry::TRIANGLE));

	///\TODO if it is decided that getBasisFunctionValue and getBasisFucntionDerivative remain here, test them

	///\TODO testing that the refinement maps behave exactly like the forwarded calls of this class
}


