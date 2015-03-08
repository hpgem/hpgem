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

#include "Geometry/ReferenceCube.hpp"
#include "Geometry/ReferenceSquare.hpp"
#include "Geometry/ReferenceLine.hpp"
#include "Geometry/ReferencePoint.hpp"
#include <iostream>
#include "cassert"

#include "Geometry/PointReference.hpp"
#include "Geometry/Mappings/MappingToRefSquareToCube.hpp"
#include "Geometry/Mappings/MappingToRefCubeToCube.hpp"
#include "Integration/QuadratureRules/GaussQuadratureRule.hpp"
#include <cmath>
using std::fabs;
using Geometry::ReferenceCube;

int main(){
	ReferenceCube& test=ReferenceCube::Instance();

	Geometry::PointReference pTest(3);

	//testing basic functionality

	for(pTest[0]=-3.141;pTest[0]<-1.;pTest[0]+=0.1){
		for(pTest[1]=-3.1416;pTest[1]<3.1416;pTest[1]+=0.1){
			for(pTest[2]=-3.1416;pTest[2]<3.1416;pTest[2]+=0.1){
				assert(("isInternalPoint",!test.isInternalPoint(pTest)));
			}
		}
	}
	for(;pTest[0]<1;pTest[0]+=0.1){
		for(pTest[1]=-3.1417;pTest[1]<-1.;pTest[1]+=0.1){
			for(pTest[2]=-3.1416;pTest[2]<3.1416;pTest[2]+=0.1){
				assert(("isInternalPoint",!test.isInternalPoint(pTest)));
			}
		}
		for(;pTest[1]<1.;pTest[1]+=0.1){
			for(pTest[2]=-3.1416;pTest[2]<-1.;pTest[2]+=0.1){
				assert(("isInternalPoint",!test.isInternalPoint(pTest)));
			}
			for(;pTest[2]<1.;pTest[2]+=0.1){
				assert(("isInternalPoint",test.isInternalPoint(pTest)));
			}
			for(;pTest[2]<3.141;pTest[2]+=0.1){
				assert(("isInternalPoint",!test.isInternalPoint(pTest)));
			}
		}
		for(;pTest[1]<3.141;pTest[1]+=0.1){
			for(pTest[2]=-3.1416;pTest[2]<3.1416;pTest[2]+=0.1){
				assert(("isInternalPoint",!test.isInternalPoint(pTest)));
			}
		}
	}
	for(;pTest[0]<3.141;pTest[0]+=0.1){
		for(pTest[1]=-3.1416;pTest[1]<3.1416;pTest[1]+=0.1){
			for(pTest[2]=-3.1416;pTest[2]<3.1416;pTest[2]+=0.1){
				assert(("isInternalPoint",!test.isInternalPoint(pTest)));
			}
		}
	}

	test.getCenter(pTest);
	assert(("getCenter",test.isInternalPoint(pTest)&&fabs(pTest[0])<1e-12&&fabs(pTest[1])<1e-12)&&fabs(pTest[2])<1e-12);
	test.getNode(0,pTest);
	assert(("getNode 0",fabs(pTest[0]+1)<1e-12&&fabs(pTest[1]+1)<1e-12&&fabs(pTest[2]+1)<1e-12));
	test.getNode(1,pTest);
	assert(("getNode 1",fabs(pTest[0]-1)<1e-12&&fabs(pTest[1]+1)<1e-12&&fabs(pTest[2]+1)<1e-12));
	test.getNode(2,pTest);
	assert(("getNode 2",fabs(pTest[0]+1)<1e-12&&fabs(pTest[1]-1)<1e-12&&fabs(pTest[2]+1)<1e-12));
	test.getNode(3,pTest);
	assert(("getNode 3",fabs(pTest[0]-1)<1e-12&&fabs(pTest[1]-1)<1e-12&&fabs(pTest[2]+1)<1e-12));
	test.getNode(4,pTest);
	assert(("getNode 4",fabs(pTest[0]+1)<1e-12&&fabs(pTest[1]+1)<1e-12&&fabs(pTest[2]-1)<1e-12));
	test.getNode(5,pTest);
	assert(("getNode 5",fabs(pTest[0]-1)<1e-12&&fabs(pTest[1]+1)<1e-12&&fabs(pTest[2]-1)<1e-12));
	test.getNode(6,pTest);
	assert(("getNode 6",fabs(pTest[0]+1)<1e-12&&fabs(pTest[1]-1)<1e-12&&fabs(pTest[2]-1)<1e-12));
	test.getNode(7,pTest);
	assert(("getNode 7",fabs(pTest[0]-1)<1e-12&&fabs(pTest[1]-1)<1e-12&&fabs(pTest[2]-1)<1e-12));
	std::cout<<test.getName();

	assert(("getLocalNodeIndex 0",test.getLocalNodeIndex(0,0)==0));//the nodes of the face must always be specified IN THIS SPECIFIC ORDER
	assert(("getLocalNodeIndex 0",test.getLocalNodeIndex(0,1)==1));//im not sure if I like this myself, but this should at least verify
	assert(("getLocalNodeIndex 0",test.getLocalNodeIndex(0,2)==2));//that all face nodes are specified, none are specified twice
	assert(("getLocalNodeIndex 0",test.getLocalNodeIndex(0,3)==3));//and only face nodes are specified and the ordering of the nodes is consistent
	assert(("getLocalNodeIndex 1",test.getLocalNodeIndex(1,0)==0));//across function calls
	assert(("getLocalNodeIndex 1",test.getLocalNodeIndex(1,1)==1));
	assert(("getLocalNodeIndex 1",test.getLocalNodeIndex(1,2)==4));
	assert(("getLocalNodeIndex 1",test.getLocalNodeIndex(1,3)==5));
	assert(("getLocalNodeIndex 2",test.getLocalNodeIndex(2,0)==0));
	assert(("getLocalNodeIndex 2",test.getLocalNodeIndex(2,1)==2));
	assert(("getLocalNodeIndex 2",test.getLocalNodeIndex(2,2)==4));
	assert(("getLocalNodeIndex 2",test.getLocalNodeIndex(2,3)==6));
	assert(("getLocalNodeIndex 3",test.getLocalNodeIndex(3,0)==1));
	assert(("getLocalNodeIndex 3",test.getLocalNodeIndex(3,1)==3));
	assert(("getLocalNodeIndex 3",test.getLocalNodeIndex(3,2)==5));
	assert(("getLocalNodeIndex 3",test.getLocalNodeIndex(3,3)==7));
	assert(("getLocalNodeIndex 4",test.getLocalNodeIndex(4,0)==2));
	assert(("getLocalNodeIndex 4",test.getLocalNodeIndex(4,1)==3));
	assert(("getLocalNodeIndex 4",test.getLocalNodeIndex(4,2)==6));
	assert(("getLocalNodeIndex 4",test.getLocalNodeIndex(4,3)==7));
	assert(("getLocalNodeIndex 5",test.getLocalNodeIndex(5,0)==4));
	assert(("getLocalNodeIndex 5",test.getLocalNodeIndex(5,1)==5));
	assert(("getLocalNodeIndex 5",test.getLocalNodeIndex(5,2)==6));
	assert(("getLocalNodeIndex 5",test.getLocalNodeIndex(5,3)==7));

	std::cout<<test;

	//testing mappings and quadrature rules

	std::vector<unsigned int> base(8),transformed(8),faceIndices(4);
	for(int i=0;i<8;++i){//doesnt test against reordering of the nodes in the first vector
		base[i]=transformed[i]=i;
	}
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(test.getCodim0MappingIndex(base,transformed))==&Geometry::MappingToRefCubeToCube0::Instance()));
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(base,transformed)==&Geometry::MappingToRefCubeToCube0::Instance()));
	transformed[0]=1;
	transformed[1]=3;
	transformed[2]=0;
	transformed[3]=2;
	transformed[4]=5;
	transformed[5]=7;
	transformed[6]=4;
	transformed[7]=6;
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(test.getCodim0MappingIndex(base,transformed))==&Geometry::MappingToRefCubeToCube1::Instance()));
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(base,transformed)==&Geometry::MappingToRefCubeToCube1::Instance()));
	transformed[0]=3;
	transformed[1]=2;
	transformed[2]=1;
	transformed[3]=0;
	transformed[4]=7;
	transformed[5]=6;
	transformed[6]=5;
	transformed[7]=4;
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(test.getCodim0MappingIndex(base,transformed))==&Geometry::MappingToRefCubeToCube2::Instance()));
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(base,transformed)==&Geometry::MappingToRefCubeToCube2::Instance()));
	transformed[0]=2;
	transformed[1]=0;
	transformed[2]=3;
	transformed[3]=1;
	transformed[4]=6;
	transformed[5]=4;
	transformed[6]=7;
	transformed[7]=5;
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(test.getCodim0MappingIndex(base,transformed))==&Geometry::MappingToRefCubeToCube3::Instance()));
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(base,transformed)==&Geometry::MappingToRefCubeToCube3::Instance()));
	transformed[0]=2;
	transformed[1]=3;
	transformed[2]=0;
	transformed[3]=1;
	transformed[4]=6;
	transformed[5]=7;
	transformed[6]=4;
	transformed[7]=5;
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(test.getCodim0MappingIndex(base,transformed))==&Geometry::MappingToRefCubeToCube4::Instance()));
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(base,transformed)==&Geometry::MappingToRefCubeToCube4::Instance()));
	transformed[0]=1;
	transformed[1]=0;
	transformed[2]=3;
	transformed[3]=2;
	transformed[4]=5;
	transformed[5]=4;
	transformed[6]=7;
	transformed[7]=6;
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(test.getCodim0MappingIndex(base,transformed))==&Geometry::MappingToRefCubeToCube5::Instance()));
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(base,transformed)==&Geometry::MappingToRefCubeToCube5::Instance()));
	transformed[0]=3;
	transformed[1]=1;
	transformed[2]=2;
	transformed[3]=0;
	transformed[4]=7;
	transformed[5]=5;
	transformed[6]=6;
	transformed[7]=4;
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(test.getCodim0MappingIndex(base,transformed))==&Geometry::MappingToRefCubeToCube6::Instance()));
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(base,transformed)==&Geometry::MappingToRefCubeToCube6::Instance()));
	transformed[0]=0;
	transformed[1]=2;
	transformed[2]=1;
	transformed[3]=3;
	transformed[4]=4;
	transformed[5]=6;
	transformed[6]=5;
	transformed[7]=7;
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(test.getCodim0MappingIndex(base,transformed))==&Geometry::MappingToRefCubeToCube7::Instance()));
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(base,transformed)==&Geometry::MappingToRefCubeToCube7::Instance()));


	assert(("higher codimensional entities",test.getNrOfCodim1Entities()==6&&test.getNrOfCodim2Entities()==12)&&test.getNrOfCodim3Entities()==8);
	assert(("getCodim1ReferenceGeometry",test.getCodim1ReferenceGeometry(0)==&Geometry::ReferenceSquare::Instance()&&
										 test.getCodim1ReferenceGeometry(1)==&Geometry::ReferenceSquare::Instance()&&
										 test.getCodim1ReferenceGeometry(2)==&Geometry::ReferenceSquare::Instance()&&
										 test.getCodim1ReferenceGeometry(3)==&Geometry::ReferenceSquare::Instance()&&
										 test.getCodim1ReferenceGeometry(4)==&Geometry::ReferenceSquare::Instance()&&
										 test.getCodim1ReferenceGeometry(5)==&Geometry::ReferenceSquare::Instance()));
	assert(("getCodim2ReferenceGeometry",test.getCodim2ReferenceGeometry(0)==&Geometry::ReferenceLine::Instance()&&
										 test.getCodim2ReferenceGeometry(1)==&Geometry::ReferenceLine::Instance()&&
										 test.getCodim2ReferenceGeometry(2)==&Geometry::ReferenceLine::Instance()&&
										 test.getCodim2ReferenceGeometry(3)==&Geometry::ReferenceLine::Instance()&&
										 test.getCodim2ReferenceGeometry(4)==&Geometry::ReferenceLine::Instance()&&
										 test.getCodim2ReferenceGeometry(5)==&Geometry::ReferenceLine::Instance()&&
										 test.getCodim2ReferenceGeometry(6)==&Geometry::ReferenceLine::Instance()&&
										 test.getCodim2ReferenceGeometry(7)==&Geometry::ReferenceLine::Instance()&&
										 test.getCodim2ReferenceGeometry(8)==&Geometry::ReferenceLine::Instance()&&
										 test.getCodim2ReferenceGeometry(9)==&Geometry::ReferenceLine::Instance()&&
										 test.getCodim2ReferenceGeometry(10)==&Geometry::ReferenceLine::Instance()&&
										 test.getCodim2ReferenceGeometry(11)==&Geometry::ReferenceLine::Instance()));
	assert(("getCodim1MappingPtr",test.getCodim1MappingPtr(0)==&Geometry::MappingToRefSquareToCube0::Instance()));
	assert(("getCodim1MappingPtr",test.getCodim1MappingPtr(1)==&Geometry::MappingToRefSquareToCube1::Instance()));
	assert(("getCodim1MappingPtr",test.getCodim1MappingPtr(2)==&Geometry::MappingToRefSquareToCube2::Instance()));
	assert(("getCodim1MappingPtr",test.getCodim1MappingPtr(3)==&Geometry::MappingToRefSquareToCube3::Instance()));
	assert(("getCodim1MappingPtr",test.getCodim1MappingPtr(4)==&Geometry::MappingToRefSquareToCube4::Instance()));
	assert(("getCodim1MappingPtr",test.getCodim1MappingPtr(5)==&Geometry::MappingToRefSquareToCube5::Instance()));
	test.getCodim1EntityLocalIndices(0,faceIndices);
	assert(("getCodim1EntityLocalIndices",faceIndices[0]==test.getLocalNodeIndex(0,0)));
	assert(("getCodim1EntityLocalIndices",faceIndices[1]==test.getLocalNodeIndex(0,1)));
	assert(("getCodim1EntityLocalIndices",faceIndices[2]==test.getLocalNodeIndex(0,2)));
	assert(("getCodim1EntityLocalIndices",faceIndices[3]==test.getLocalNodeIndex(0,3)));
	test.getCodim1EntityLocalIndices(1,faceIndices);
	assert(("getCodim1EntityLocalIndices",faceIndices[0]==test.getLocalNodeIndex(1,0)));
	assert(("getCodim1EntityLocalIndices",faceIndices[1]==test.getLocalNodeIndex(1,1)));
	assert(("getCodim1EntityLocalIndices",faceIndices[2]==test.getLocalNodeIndex(1,2)));
	assert(("getCodim1EntityLocalIndices",faceIndices[3]==test.getLocalNodeIndex(1,3)));
	test.getCodim1EntityLocalIndices(2,faceIndices);
	assert(("getCodim1EntityLocalIndices",faceIndices[0]==test.getLocalNodeIndex(2,0)));
	assert(("getCodim1EntityLocalIndices",faceIndices[1]==test.getLocalNodeIndex(2,1)));
	assert(("getCodim1EntityLocalIndices",faceIndices[2]==test.getLocalNodeIndex(2,2)));
	assert(("getCodim1EntityLocalIndices",faceIndices[3]==test.getLocalNodeIndex(2,3)));
	test.getCodim1EntityLocalIndices(3,faceIndices);
	assert(("getCodim1EntityLocalIndices",faceIndices[0]==test.getLocalNodeIndex(3,0)));
	assert(("getCodim1EntityLocalIndices",faceIndices[1]==test.getLocalNodeIndex(3,1)));
	assert(("getCodim1EntityLocalIndices",faceIndices[2]==test.getLocalNodeIndex(3,2)));
	assert(("getCodim1EntityLocalIndices",faceIndices[3]==test.getLocalNodeIndex(3,3)));
	test.getCodim1EntityLocalIndices(4,faceIndices);
	assert(("getCodim1EntityLocalIndices",faceIndices[0]==test.getLocalNodeIndex(4,0)));
	assert(("getCodim1EntityLocalIndices",faceIndices[1]==test.getLocalNodeIndex(4,1)));
	assert(("getCodim1EntityLocalIndices",faceIndices[2]==test.getLocalNodeIndex(4,2)));
	assert(("getCodim1EntityLocalIndices",faceIndices[3]==test.getLocalNodeIndex(4,3)));
	test.getCodim1EntityLocalIndices(5,faceIndices);
	assert(("getCodim1EntityLocalIndices",faceIndices[0]==test.getLocalNodeIndex(5,0)));
	assert(("getCodim1EntityLocalIndices",faceIndices[1]==test.getLocalNodeIndex(5,1)));
	assert(("getCodim1EntityLocalIndices",faceIndices[2]==test.getLocalNodeIndex(5,2)));
	assert(("getCodim1EntityLocalIndices",faceIndices[3]==test.getLocalNodeIndex(5,3)));
	faceIndices.resize(2);
	test.getCodim2EntityLocalIndices(0,faceIndices);
	assert(("getCodim2EntityLocalIndices",faceIndices[0]==0));
	assert(("getCodim2EntityLocalIndices",faceIndices[1]==1));
	test.getCodim2EntityLocalIndices(1,faceIndices);
	assert(("getCodim2EntityLocalIndices",faceIndices[0]==2));
	assert(("getCodim2EntityLocalIndices",faceIndices[1]==3));
	test.getCodim2EntityLocalIndices(2,faceIndices);
	assert(("getCodim2EntityLocalIndices",faceIndices[0]==4));
	assert(("getCodim2EntityLocalIndices",faceIndices[1]==5));
	test.getCodim2EntityLocalIndices(3,faceIndices);
	assert(("getCodim2EntityLocalIndices",faceIndices[0]==6));
	assert(("getCodim2EntityLocalIndices",faceIndices[1]==7));
	test.getCodim2EntityLocalIndices(4,faceIndices);
	assert(("getCodim2EntityLocalIndices",faceIndices[0]==0));
	assert(("getCodim2EntityLocalIndices",faceIndices[1]==2));
	test.getCodim2EntityLocalIndices(5,faceIndices);
	assert(("getCodim2EntityLocalIndices",faceIndices[0]==1));
	assert(("getCodim2EntityLocalIndices",faceIndices[1]==3));
	test.getCodim2EntityLocalIndices(6,faceIndices);
	assert(("getCodim2EntityLocalIndices",faceIndices[0]==4));
	assert(("getCodim2EntityLocalIndices",faceIndices[1]==6));
	test.getCodim2EntityLocalIndices(7,faceIndices);
	assert(("getCodim2EntityLocalIndices",faceIndices[0]==5));
	assert(("getCodim2EntityLocalIndices",faceIndices[1]==7));
	test.getCodim2EntityLocalIndices(8,faceIndices);
	assert(("getCodim2EntityLocalIndices",faceIndices[0]==0));
	assert(("getCodim2EntityLocalIndices",faceIndices[1]==4));
	test.getCodim2EntityLocalIndices(9,faceIndices);
	assert(("getCodim2EntityLocalIndices",faceIndices[0]==1));
	assert(("getCodim2EntityLocalIndices",faceIndices[1]==5));
	test.getCodim2EntityLocalIndices(10,faceIndices);
	assert(("getCodim2EntityLocalIndices",faceIndices[0]==2));
	assert(("getCodim2EntityLocalIndices",faceIndices[1]==6));
	test.getCodim2EntityLocalIndices(11,faceIndices);
	assert(("getCodim2EntityLocalIndices",faceIndices[0]==3));
	assert(("getCodim2EntityLocalIndices",faceIndices[1]==7));
	faceIndices.resize(1);
	test.getCodim3EntityLocalIndices(0,faceIndices);
	assert(("getCodim3EntityLocalIndices",faceIndices[0]==0));
	test.getCodim3EntityLocalIndices(1,faceIndices);
	assert(("getCodim3EntityLocalIndices",faceIndices[0]==1));
	test.getCodim3EntityLocalIndices(2,faceIndices);
	assert(("getCodim3EntityLocalIndices",faceIndices[0]==2));
	test.getCodim3EntityLocalIndices(3,faceIndices);
	assert(("getCodim3EntityLocalIndices",faceIndices[0]==3));
	test.getCodim3EntityLocalIndices(4,faceIndices);
	assert(("getCodim3EntityLocalIndices",faceIndices[0]==4));
	test.getCodim3EntityLocalIndices(5,faceIndices);
	assert(("getCodim3EntityLocalIndices",faceIndices[0]==5));
	test.getCodim3EntityLocalIndices(6,faceIndices);
	assert(("getCodim3EntityLocalIndices",faceIndices[0]==6));
	test.getCodim3EntityLocalIndices(7,faceIndices);
	assert(("getCodim3EntityLocalIndices",faceIndices[0]==7));


	assert(("quadrature rules",test.getGaussQuadratureRule(3)->order()>=3));
	assert(("quadrature rules",test.getGaussQuadratureRule(5)->order()>=5));
	assert(("quadrature rules",test.getGaussQuadratureRule(7)->order()>=7));
	assert(("quadrature rules",test.getGaussQuadratureRule(9)->order()>=9));
	assert(("quadrature rules",test.getGaussQuadratureRule(11)->order()>=11));

	//testing functionality of abstract parent classes

	assert(("number of nodes",test.getNumberOfNodes()==8));
	assert(("type of geometry",test.getGeometryType()==Geometry::CUBE));

	///\TODO if it is decided that getBasisFunctionValue and getBasisFucntionDerivative remain here, test them

	///\TODO testing that the refinement maps behave exactly like the forwarded calls of this class
}


