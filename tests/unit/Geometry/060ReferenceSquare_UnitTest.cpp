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

#include "Geometry/ReferenceSquare.hpp"
#include "Geometry/ReferenceLine.hpp"
#include "Geometry/ReferencePoint.hpp"
#include "Geometry/Mappings/MappingToRefSquareToSquare.hpp"
#include "Geometry/Mappings/MappingToRefLineToSquare.hpp"
#include <iostream>
#include "cassert"

#include "Geometry/PointReference.hpp"
#include "Integration/QuadratureRules/GaussQuadratureRule.hpp"
#include <cmath>
using Geometry::ReferenceSquare;

int main(){
	ReferenceSquare& test=ReferenceSquare::Instance();

	Geometry::PointReference pTest(2);

	//testing basic functionality

	for(pTest[0]=-3.141;pTest[0]<-1;pTest[0]+=0.1){
		for(pTest[1]=-3.1416;pTest[1]<3.1416;pTest[1]+=0.1){
			assert(("isInternalPoint",!test.isInternalPoint(pTest)));
		}
	}
	for(;pTest[0]<1;pTest[0]+=0.1){
		for(pTest[1]=-3.1416;pTest[1]<-1;pTest[1]+=0.1){
			assert(("isInternalPoint",!test.isInternalPoint(pTest)));
		}
		for(;pTest[1]<1.;pTest[1]+=0.1){
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
	assert(("getCenter",test.isInternalPoint(pTest)&&std::abs(pTest[0])<1e-12&&std::abs(pTest[1])<1e-12));
	test.getNode(0,pTest);
	assert(("getNode 0",std::abs(pTest[0]+1)<1e-12&&std::abs(pTest[1]+1)<1e-12));
	test.getNode(1,pTest);
	assert(("getNode 1",std::abs(pTest[0]-1)<1e-12&&std::abs(pTest[1]+1)<1e-12));
	test.getNode(2,pTest);
	assert(("getNode 2",std::abs(pTest[0]+1)<1e-12&&std::abs(pTest[1]-1)<1e-12));
	test.getNode(3,pTest);
	assert(("getNode 3",std::abs(pTest[0]-1)<1e-12&&std::abs(pTest[1]-1)<1e-12));
	std::cout<<test.getName();

	assert(("getLocalNodeIndex 0",test.getLocalNodeIndex(0,0)==0));//the nodes of the face must always be specified IN THIS SPECIFIC ORDER
	assert(("getLocalNodeIndex 0",test.getLocalNodeIndex(0,1)==1));//this is needed because the outward pointing normal vector
	assert(("getLocalNodeIndex 1",test.getLocalNodeIndex(1,0)==0));//will automatically point outward when compute it using this node ordering
	assert(("getLocalNodeIndex 1",test.getLocalNodeIndex(1,1)==2));///\BUG some orderings are wrong
	assert(("getLocalNodeIndex 2",test.getLocalNodeIndex(2,0)==1));
	assert(("getLocalNodeIndex 2",test.getLocalNodeIndex(2,1)==3));
	assert(("getLocalNodeIndex 3",test.getLocalNodeIndex(3,0)==2));
	assert(("getLocalNodeIndex 3",test.getLocalNodeIndex(3,1)==3));

	std::cout<<test;

	//testing mappings and quadrature rules

	std::vector<std::size_t> base(4),transformed(4),faceIndices(2);
	for(int i=0;i<4;++i){//doesnt test against reordering of the nodes in the first vector
		base[i]=transformed[i]=i;
	}
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(test.getCodim0MappingIndex(base,transformed))==&Geometry::MappingToRefSquareToSquare0::Instance()));
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(base,transformed)==&Geometry::MappingToRefSquareToSquare0::Instance()));
	transformed[0]=1;
	transformed[1]=3;
	transformed[2]=0;
	transformed[3]=2;
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(test.getCodim0MappingIndex(base,transformed))==&Geometry::MappingToRefSquareToSquare1::Instance()));
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(base,transformed)==&Geometry::MappingToRefSquareToSquare1::Instance()));
	transformed[0]=3;
	transformed[1]=2;
	transformed[2]=1;
	transformed[3]=0;
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(test.getCodim0MappingIndex(base,transformed))==&Geometry::MappingToRefSquareToSquare2::Instance()));
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(base,transformed)==&Geometry::MappingToRefSquareToSquare2::Instance()));
	transformed[0]=2;
	transformed[1]=0;
	transformed[2]=3;
	transformed[3]=1;
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(test.getCodim0MappingIndex(base,transformed))==&Geometry::MappingToRefSquareToSquare3::Instance()));
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(base,transformed)==&Geometry::MappingToRefSquareToSquare3::Instance()));
	transformed[0]=2;
	transformed[1]=3;
	transformed[2]=0;
	transformed[3]=1;
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(test.getCodim0MappingIndex(base,transformed))==&Geometry::MappingToRefSquareToSquare4::Instance()));
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(base,transformed)==&Geometry::MappingToRefSquareToSquare4::Instance()));
	transformed[0]=1;
	transformed[1]=0;
	transformed[2]=3;
	transformed[3]=2;
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(test.getCodim0MappingIndex(base,transformed))==&Geometry::MappingToRefSquareToSquare5::Instance()));
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(base,transformed)==&Geometry::MappingToRefSquareToSquare5::Instance()));
	transformed[0]=3;
	transformed[1]=1;
	transformed[2]=2;
	transformed[3]=0;
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(test.getCodim0MappingIndex(base,transformed))==&Geometry::MappingToRefSquareToSquare6::Instance()));
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(base,transformed)==&Geometry::MappingToRefSquareToSquare6::Instance()));
	transformed[0]=0;
	transformed[1]=2;
	transformed[2]=1;
	transformed[3]=3;
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(test.getCodim0MappingIndex(base,transformed))==&Geometry::MappingToRefSquareToSquare7::Instance()));
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(base,transformed)==&Geometry::MappingToRefSquareToSquare7::Instance()));

	assert(("higher codimensional entities",test.getNrOfCodim1Entities()==4&&test.getNrOfCodim2Entities()==4)&&test.getNrOfCodim3Entities()==0);
	assert(("getCodim1ReferenceGeometry",test.getCodim1ReferenceGeometry(0)==&Geometry::ReferenceLine::Instance()&&
										 test.getCodim1ReferenceGeometry(1)==&Geometry::ReferenceLine::Instance()&&
										 test.getCodim1ReferenceGeometry(2)==&Geometry::ReferenceLine::Instance()&&
										 test.getCodim1ReferenceGeometry(3)==&Geometry::ReferenceLine::Instance()));
	assert(("getCodim2ReferenceGeometry",test.getCodim2ReferenceGeometry(0)==&Geometry::ReferencePoint::Instance()&&
										 test.getCodim2ReferenceGeometry(1)==&Geometry::ReferencePoint::Instance()&&
										 test.getCodim2ReferenceGeometry(2)==&Geometry::ReferencePoint::Instance()&&
										 test.getCodim2ReferenceGeometry(3)==&Geometry::ReferencePoint::Instance()));
	assert(("getCodim1MappingPtr",test.getCodim1MappingPtr(0)==&Geometry::MappingToRefLineToSquare0::Instance()));
	assert(("getCodim1MappingPtr",test.getCodim1MappingPtr(1)==&Geometry::MappingToRefLineToSquare1::Instance()));
	assert(("getCodim1MappingPtr",test.getCodim1MappingPtr(2)==&Geometry::MappingToRefLineToSquare2::Instance()));
	assert(("getCodim1MappingPtr",test.getCodim1MappingPtr(3)==&Geometry::MappingToRefLineToSquare3::Instance()));
	test.getCodim1EntityLocalIndices(0,faceIndices);
	assert(("getCodim1EntityLocalIndices",faceIndices[0]==test.getLocalNodeIndex(0,0)));
	assert(("getCodim1EntityLocalIndices",faceIndices[1]==test.getLocalNodeIndex(0,1)));
	test.getCodim1EntityLocalIndices(1,faceIndices);
	assert(("getCodim1EntityLocalIndices",faceIndices[0]==test.getLocalNodeIndex(1,0)));
	assert(("getCodim1EntityLocalIndices",faceIndices[1]==test.getLocalNodeIndex(1,1)));
	test.getCodim1EntityLocalIndices(2,faceIndices);
	assert(("getCodim1EntityLocalIndices",faceIndices[0]==test.getLocalNodeIndex(2,0)));
	assert(("getCodim1EntityLocalIndices",faceIndices[1]==test.getLocalNodeIndex(2,1)));
	test.getCodim1EntityLocalIndices(3,faceIndices);
	assert(("getCodim1EntityLocalIndices",faceIndices[0]==test.getLocalNodeIndex(3,0)));
	assert(("getCodim1EntityLocalIndices",faceIndices[1]==test.getLocalNodeIndex(3,1)));
	faceIndices.resize(1);
	test.getCodim2EntityLocalIndices(0,faceIndices);
	assert(("getCodim2EntityLocalIndices",faceIndices[0]==0));
	test.getCodim2EntityLocalIndices(1,faceIndices);
	assert(("getCodim2EntityLocalIndices",faceIndices[0]==1));
	test.getCodim2EntityLocalIndices(2,faceIndices);
	assert(("getCodim2EntityLocalIndices",faceIndices[0]==2));
	test.getCodim2EntityLocalIndices(3,faceIndices);
	assert(("getCodim2EntityLocalIndices",faceIndices[0]==3));

	assert(("quadrature rules",test.getGaussQuadratureRule(3)->order()>=3));
	assert(("quadrature rules",test.getGaussQuadratureRule(5)->order()>=5));
	assert(("quadrature rules",test.getGaussQuadratureRule(7)->order()>=7));
	assert(("quadrature rules",test.getGaussQuadratureRule(9)->order()>=9));
	assert(("quadrature rules",test.getGaussQuadratureRule(11)->order()>=11));

	//testing functionality of abstract parent classes

	assert(("number of nodes",test.getNumberOfNodes()==4));
	assert(("type of geometry",test.getGeometryType()==Geometry::SQUARE));

	///\TODO if it is decided that getBasisFunctionValue and getBasisFucntionDerivative remain here, test them

	///\TODO testing that the refinement maps behave exactly like the forwarded calls of this class
}


