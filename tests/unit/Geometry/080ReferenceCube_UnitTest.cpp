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
#include "Logger.h"

#include "Geometry/PointReference.hpp"
#include "Geometry/Mappings/MappingToRefSquareToCube.hpp"
#include "Geometry/Mappings/MappingToRefCubeToCube.hpp"
#include "Integration/QuadratureRules/GaussQuadratureRule.hpp"
#include <cmath>
using Geometry::ReferenceCube;

int main(){
	ReferenceCube& test=ReferenceCube::Instance();

	Geometry::PointReference pTest(3);

	//testing basic functionality

	for(pTest[0]=-3.141;pTest[0]<-1.;pTest[0]+=0.1){
		for(pTest[1]=-3.1416;pTest[1]<3.1416;pTest[1]+=0.1){
			for(pTest[2]=-3.1416;pTest[2]<3.1416;pTest[2]+=0.1){
				logger.assert_always((!test.isInternalPoint(pTest)),"isInternalPoint");
			}
		}
	}
	for(;pTest[0]<1;pTest[0]+=0.1){
		for(pTest[1]=-3.1417;pTest[1]<-1.;pTest[1]+=0.1){
			for(pTest[2]=-3.1416;pTest[2]<3.1416;pTest[2]+=0.1){
				logger.assert_always((!test.isInternalPoint(pTest)),"isInternalPoint");
			}
		}
		for(;pTest[1]<1.;pTest[1]+=0.1){
			for(pTest[2]=-3.1416;pTest[2]<-1.;pTest[2]+=0.1){
				logger.assert_always((!test.isInternalPoint(pTest)),"isInternalPoint");
			}
			for(;pTest[2]<1.;pTest[2]+=0.1){
				logger.assert_always((test.isInternalPoint(pTest)),"isInternalPoint");
			}
			for(;pTest[2]<3.141;pTest[2]+=0.1){
				logger.assert_always((!test.isInternalPoint(pTest)),"isInternalPoint");
			}
		}
		for(;pTest[1]<3.141;pTest[1]+=0.1){
			for(pTest[2]=-3.1416;pTest[2]<3.1416;pTest[2]+=0.1){
				logger.assert_always((!test.isInternalPoint(pTest)),"isInternalPoint");
			}
		}
	}
	for(;pTest[0]<3.141;pTest[0]+=0.1){
		for(pTest[1]=-3.1416;pTest[1]<3.1416;pTest[1]+=0.1){
			for(pTest[2]=-3.1416;pTest[2]<3.1416;pTest[2]+=0.1){
				logger.assert_always((!test.isInternalPoint(pTest)),"isInternalPoint");
			}
		}
	}

	test.getCenter(pTest);
	logger.assert_always((test.isInternalPoint(pTest)&&std::abs(pTest[0])<1e-12&&std::abs(pTest[1])<1e-12)&&std::abs(pTest[2])<1e-12,"getCenter");
	test.getNode(0,pTest);
	logger.assert_always((std::abs(pTest[0]+1)<1e-12&&std::abs(pTest[1]+1)<1e-12&&std::abs(pTest[2]+1)<1e-12),"getNode 0");
	test.getNode(1,pTest);
	logger.assert_always((std::abs(pTest[0]-1)<1e-12&&std::abs(pTest[1]+1)<1e-12&&std::abs(pTest[2]+1)<1e-12),"getNode 1");
	test.getNode(2,pTest);
	logger.assert_always((std::abs(pTest[0]+1)<1e-12&&std::abs(pTest[1]-1)<1e-12&&std::abs(pTest[2]+1)<1e-12),"getNode 2");
	test.getNode(3,pTest);
	logger.assert_always((std::abs(pTest[0]-1)<1e-12&&std::abs(pTest[1]-1)<1e-12&&std::abs(pTest[2]+1)<1e-12),"getNode 3");
	test.getNode(4,pTest);
	logger.assert_always((std::abs(pTest[0]+1)<1e-12&&std::abs(pTest[1]+1)<1e-12&&std::abs(pTest[2]-1)<1e-12),"getNode 4");
	test.getNode(5,pTest);
	logger.assert_always((std::abs(pTest[0]-1)<1e-12&&std::abs(pTest[1]+1)<1e-12&&std::abs(pTest[2]-1)<1e-12),"getNode 5");
	test.getNode(6,pTest);
	logger.assert_always((std::abs(pTest[0]+1)<1e-12&&std::abs(pTest[1]-1)<1e-12&&std::abs(pTest[2]-1)<1e-12),"getNode 6");
	test.getNode(7,pTest);
	logger.assert_always((std::abs(pTest[0]-1)<1e-12&&std::abs(pTest[1]-1)<1e-12&&std::abs(pTest[2]-1)<1e-12),"getNode 7");
	std::cout<<test.getName();

	logger.assert_always((test.getLocalNodeIndex(0,0)==0),"getLocalNodeIndex 0");//the nodes of the face must always be specified IN THIS SPECIFIC ORDER
	logger.assert_always((test.getLocalNodeIndex(0,1)==1),"getLocalNodeIndex 0");//im not sure if I like this myself, but this should at least verify
	logger.assert_always((test.getLocalNodeIndex(0,2)==2),"getLocalNodeIndex 0");//that all face nodes are specified, none are specified twice
	logger.assert_always((test.getLocalNodeIndex(0,3)==3),"getLocalNodeIndex 0");//and only face nodes are specified and the ordering of the nodes is consistent
	logger.assert_always((test.getLocalNodeIndex(1,0)==0),"getLocalNodeIndex 1");//across function calls
	logger.assert_always((test.getLocalNodeIndex(1,1)==1),"getLocalNodeIndex 1");
	logger.assert_always((test.getLocalNodeIndex(1,2)==4),"getLocalNodeIndex 1");
	logger.assert_always((test.getLocalNodeIndex(1,3)==5),"getLocalNodeIndex 1");
	logger.assert_always((test.getLocalNodeIndex(2,0)==0),"getLocalNodeIndex 2");
	logger.assert_always((test.getLocalNodeIndex(2,1)==2),"getLocalNodeIndex 2");
	logger.assert_always((test.getLocalNodeIndex(2,2)==4),"getLocalNodeIndex 2");
	logger.assert_always((test.getLocalNodeIndex(2,3)==6),"getLocalNodeIndex 2");
	logger.assert_always((test.getLocalNodeIndex(3,0)==1),"getLocalNodeIndex 3");
	logger.assert_always((test.getLocalNodeIndex(3,1)==3),"getLocalNodeIndex 3");
	logger.assert_always((test.getLocalNodeIndex(3,2)==5),"getLocalNodeIndex 3");
	logger.assert_always((test.getLocalNodeIndex(3,3)==7),"getLocalNodeIndex 3");
	logger.assert_always((test.getLocalNodeIndex(4,0)==2),"getLocalNodeIndex 4");
	logger.assert_always((test.getLocalNodeIndex(4,1)==3),"getLocalNodeIndex 4");
	logger.assert_always((test.getLocalNodeIndex(4,2)==6),"getLocalNodeIndex 4");
	logger.assert_always((test.getLocalNodeIndex(4,3)==7),"getLocalNodeIndex 4");
	logger.assert_always((test.getLocalNodeIndex(5,0)==4),"getLocalNodeIndex 5");
	logger.assert_always((test.getLocalNodeIndex(5,1)==5),"getLocalNodeIndex 5");
	logger.assert_always((test.getLocalNodeIndex(5,2)==6),"getLocalNodeIndex 5");
	logger.assert_always((test.getLocalNodeIndex(5,3)==7),"getLocalNodeIndex 5");

	std::cout<<test;

	//testing mappings and quadrature rules

	std::vector<std::size_t> base(8),transformed(8),faceIndices(4);
	for(std::size_t i=0;i<8;++i){//doesnt test against reordering of the nodes in the first vector
		base[i]=transformed[i]=i;
	}
	logger.assert_always((test.getCodim0MappingPtr(test.getCodim0MappingIndex(base,transformed))==&Geometry::MappingToRefCubeToCube0::Instance()),"getCodim0MappingIndex&Ptr");
	logger.assert_always((test.getCodim0MappingPtr(base,transformed)==&Geometry::MappingToRefCubeToCube0::Instance()),"getCodim0MappingIndex&Ptr");
	transformed[0]=1;
	transformed[1]=3;
	transformed[2]=0;
	transformed[3]=2;
	transformed[4]=5;
	transformed[5]=7;
	transformed[6]=4;
	transformed[7]=6;
	logger.assert_always((test.getCodim0MappingPtr(test.getCodim0MappingIndex(base,transformed))==&Geometry::MappingToRefCubeToCube1::Instance()),"getCodim0MappingIndex&Ptr");
	logger.assert_always((test.getCodim0MappingPtr(base,transformed)==&Geometry::MappingToRefCubeToCube1::Instance()),"getCodim0MappingIndex&Ptr");
	transformed[0]=3;
	transformed[1]=2;
	transformed[2]=1;
	transformed[3]=0;
	transformed[4]=7;
	transformed[5]=6;
	transformed[6]=5;
	transformed[7]=4;
	logger.assert_always((test.getCodim0MappingPtr(test.getCodim0MappingIndex(base,transformed))==&Geometry::MappingToRefCubeToCube2::Instance()),"getCodim0MappingIndex&Ptr");
	logger.assert_always((test.getCodim0MappingPtr(base,transformed)==&Geometry::MappingToRefCubeToCube2::Instance()),"getCodim0MappingIndex&Ptr");
	transformed[0]=2;
	transformed[1]=0;
	transformed[2]=3;
	transformed[3]=1;
	transformed[4]=6;
	transformed[5]=4;
	transformed[6]=7;
	transformed[7]=5;
	logger.assert_always((test.getCodim0MappingPtr(test.getCodim0MappingIndex(base,transformed))==&Geometry::MappingToRefCubeToCube3::Instance()),"getCodim0MappingIndex&Ptr");
	logger.assert_always((test.getCodim0MappingPtr(base,transformed)==&Geometry::MappingToRefCubeToCube3::Instance()),"getCodim0MappingIndex&Ptr");
	transformed[0]=2;
	transformed[1]=3;
	transformed[2]=0;
	transformed[3]=1;
	transformed[4]=6;
	transformed[5]=7;
	transformed[6]=4;
	transformed[7]=5;
	logger.assert_always((test.getCodim0MappingPtr(test.getCodim0MappingIndex(base,transformed))==&Geometry::MappingToRefCubeToCube4::Instance()),"getCodim0MappingIndex&Ptr");
	logger.assert_always((test.getCodim0MappingPtr(base,transformed)==&Geometry::MappingToRefCubeToCube4::Instance()),"getCodim0MappingIndex&Ptr");
	transformed[0]=1;
	transformed[1]=0;
	transformed[2]=3;
	transformed[3]=2;
	transformed[4]=5;
	transformed[5]=4;
	transformed[6]=7;
	transformed[7]=6;
	logger.assert_always((test.getCodim0MappingPtr(test.getCodim0MappingIndex(base,transformed))==&Geometry::MappingToRefCubeToCube5::Instance()),"getCodim0MappingIndex&Ptr");
	logger.assert_always((test.getCodim0MappingPtr(base,transformed)==&Geometry::MappingToRefCubeToCube5::Instance()),"getCodim0MappingIndex&Ptr");
	transformed[0]=3;
	transformed[1]=1;
	transformed[2]=2;
	transformed[3]=0;
	transformed[4]=7;
	transformed[5]=5;
	transformed[6]=6;
	transformed[7]=4;
	logger.assert_always((test.getCodim0MappingPtr(test.getCodim0MappingIndex(base,transformed))==&Geometry::MappingToRefCubeToCube6::Instance()),"getCodim0MappingIndex&Ptr");
	logger.assert_always((test.getCodim0MappingPtr(base,transformed)==&Geometry::MappingToRefCubeToCube6::Instance()),"getCodim0MappingIndex&Ptr");
	transformed[0]=0;
	transformed[1]=2;
	transformed[2]=1;
	transformed[3]=3;
	transformed[4]=4;
	transformed[5]=6;
	transformed[6]=5;
	transformed[7]=7;
	logger.assert_always((test.getCodim0MappingPtr(test.getCodim0MappingIndex(base,transformed))==&Geometry::MappingToRefCubeToCube7::Instance()),"getCodim0MappingIndex&Ptr");
	logger.assert_always((test.getCodim0MappingPtr(base,transformed)==&Geometry::MappingToRefCubeToCube7::Instance()),"getCodim0MappingIndex&Ptr");


	logger.assert_always((test.getNrOfCodim1Entities()==6&&test.getNrOfCodim2Entities()==12)&&test.getNrOfCodim3Entities()==8,"higher codimensional entities");
	logger.assert_always((test.getCodim1ReferenceGeometry(0)==&Geometry::ReferenceSquare::Instance()&&
						test.getCodim1ReferenceGeometry(1)==&Geometry::ReferenceSquare::Instance()&&
						test.getCodim1ReferenceGeometry(2)==&Geometry::ReferenceSquare::Instance()&&
						test.getCodim1ReferenceGeometry(3)==&Geometry::ReferenceSquare::Instance()&&
						test.getCodim1ReferenceGeometry(4)==&Geometry::ReferenceSquare::Instance()&&
						test.getCodim1ReferenceGeometry(5)==&Geometry::ReferenceSquare::Instance()),"getCodim1ReferenceGeometry");
	logger.assert_always((test.getCodim2ReferenceGeometry(0)==&Geometry::ReferenceLine::Instance()&&
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
						test.getCodim2ReferenceGeometry(11)==&Geometry::ReferenceLine::Instance()),"getCodim2ReferenceGeometry");
	logger.assert_always((test.getCodim1MappingPtr(0)==&Geometry::MappingToRefSquareToCube0::Instance()),"getCodim1MappingPtr");
	logger.assert_always((test.getCodim1MappingPtr(1)==&Geometry::MappingToRefSquareToCube1::Instance()),"getCodim1MappingPtr");
	logger.assert_always((test.getCodim1MappingPtr(2)==&Geometry::MappingToRefSquareToCube2::Instance()),"getCodim1MappingPtr");
	logger.assert_always((test.getCodim1MappingPtr(3)==&Geometry::MappingToRefSquareToCube3::Instance()),"getCodim1MappingPtr");
	logger.assert_always((test.getCodim1MappingPtr(4)==&Geometry::MappingToRefSquareToCube4::Instance()),"getCodim1MappingPtr");
	logger.assert_always((test.getCodim1MappingPtr(5)==&Geometry::MappingToRefSquareToCube5::Instance()),"getCodim1MappingPtr");
	test.getCodim1EntityLocalIndices(0,faceIndices);
	logger.assert_always((faceIndices[0]==test.getLocalNodeIndex(0,0)),"getCodim1EntityLocalIndices");
	logger.assert_always((faceIndices[1]==test.getLocalNodeIndex(0,1)),"getCodim1EntityLocalIndices");
	logger.assert_always((faceIndices[2]==test.getLocalNodeIndex(0,2)),"getCodim1EntityLocalIndices");
	logger.assert_always((faceIndices[3]==test.getLocalNodeIndex(0,3)),"getCodim1EntityLocalIndices");
	test.getCodim1EntityLocalIndices(1,faceIndices);
	logger.assert_always((faceIndices[0]==test.getLocalNodeIndex(1,0)),"getCodim1EntityLocalIndices");
	logger.assert_always((faceIndices[1]==test.getLocalNodeIndex(1,1)),"getCodim1EntityLocalIndices");
	logger.assert_always((faceIndices[2]==test.getLocalNodeIndex(1,2)),"getCodim1EntityLocalIndices");
	logger.assert_always((faceIndices[3]==test.getLocalNodeIndex(1,3)),"getCodim1EntityLocalIndices");
	test.getCodim1EntityLocalIndices(2,faceIndices);
	logger.assert_always((faceIndices[0]==test.getLocalNodeIndex(2,0)),"getCodim1EntityLocalIndices");
	logger.assert_always((faceIndices[1]==test.getLocalNodeIndex(2,1)),"getCodim1EntityLocalIndices");
	logger.assert_always((faceIndices[2]==test.getLocalNodeIndex(2,2)),"getCodim1EntityLocalIndices");
	logger.assert_always((faceIndices[3]==test.getLocalNodeIndex(2,3)),"getCodim1EntityLocalIndices");
	test.getCodim1EntityLocalIndices(3,faceIndices);
	logger.assert_always((faceIndices[0]==test.getLocalNodeIndex(3,0)),"getCodim1EntityLocalIndices");
	logger.assert_always((faceIndices[1]==test.getLocalNodeIndex(3,1)),"getCodim1EntityLocalIndices");
	logger.assert_always((faceIndices[2]==test.getLocalNodeIndex(3,2)),"getCodim1EntityLocalIndices");
	logger.assert_always((faceIndices[3]==test.getLocalNodeIndex(3,3)),"getCodim1EntityLocalIndices");
	test.getCodim1EntityLocalIndices(4,faceIndices);
	logger.assert_always((faceIndices[0]==test.getLocalNodeIndex(4,0)),"getCodim1EntityLocalIndices");
	logger.assert_always((faceIndices[1]==test.getLocalNodeIndex(4,1)),"getCodim1EntityLocalIndices");
	logger.assert_always((faceIndices[2]==test.getLocalNodeIndex(4,2)),"getCodim1EntityLocalIndices");
	logger.assert_always((faceIndices[3]==test.getLocalNodeIndex(4,3)),"getCodim1EntityLocalIndices");
	test.getCodim1EntityLocalIndices(5,faceIndices);
	logger.assert_always((faceIndices[0]==test.getLocalNodeIndex(5,0)),"getCodim1EntityLocalIndices");
	logger.assert_always((faceIndices[1]==test.getLocalNodeIndex(5,1)),"getCodim1EntityLocalIndices");
	logger.assert_always((faceIndices[2]==test.getLocalNodeIndex(5,2)),"getCodim1EntityLocalIndices");
	logger.assert_always((faceIndices[3]==test.getLocalNodeIndex(5,3)),"getCodim1EntityLocalIndices");
	faceIndices.resize(2);
	test.getCodim2EntityLocalIndices(0,faceIndices);
	logger.assert_always((faceIndices[0]==0),"getCodim2EntityLocalIndices");
	logger.assert_always((faceIndices[1]==1),"getCodim2EntityLocalIndices");
	test.getCodim2EntityLocalIndices(1,faceIndices);
	logger.assert_always((faceIndices[0]==2),"getCodim2EntityLocalIndices");
	logger.assert_always((faceIndices[1]==3),"getCodim2EntityLocalIndices");
	test.getCodim2EntityLocalIndices(2,faceIndices);
	logger.assert_always((faceIndices[0]==4),"getCodim2EntityLocalIndices");
	logger.assert_always((faceIndices[1]==5),"getCodim2EntityLocalIndices");
	test.getCodim2EntityLocalIndices(3,faceIndices);
	logger.assert_always((faceIndices[0]==6),"getCodim2EntityLocalIndices");
	logger.assert_always((faceIndices[1]==7),"getCodim2EntityLocalIndices");
	test.getCodim2EntityLocalIndices(4,faceIndices);
	logger.assert_always((faceIndices[0]==0),"getCodim2EntityLocalIndices");
	logger.assert_always((faceIndices[1]==2),"getCodim2EntityLocalIndices");
	test.getCodim2EntityLocalIndices(5,faceIndices);
	logger.assert_always((faceIndices[0]==1),"getCodim2EntityLocalIndices");
	logger.assert_always((faceIndices[1]==3),"getCodim2EntityLocalIndices");
	test.getCodim2EntityLocalIndices(6,faceIndices);
	logger.assert_always((faceIndices[0]==4),"getCodim2EntityLocalIndices");
	logger.assert_always((faceIndices[1]==6),"getCodim2EntityLocalIndices");
	test.getCodim2EntityLocalIndices(7,faceIndices);
	logger.assert_always((faceIndices[0]==5),"getCodim2EntityLocalIndices");
	logger.assert_always((faceIndices[1]==7),"getCodim2EntityLocalIndices");
	test.getCodim2EntityLocalIndices(8,faceIndices);
	logger.assert_always((faceIndices[0]==0),"getCodim2EntityLocalIndices");
	logger.assert_always((faceIndices[1]==4),"getCodim2EntityLocalIndices");
	test.getCodim2EntityLocalIndices(9,faceIndices);
	logger.assert_always((faceIndices[0]==1),"getCodim2EntityLocalIndices");
	logger.assert_always((faceIndices[1]==5),"getCodim2EntityLocalIndices");
	test.getCodim2EntityLocalIndices(10,faceIndices);
	logger.assert_always((faceIndices[0]==2),"getCodim2EntityLocalIndices");
	logger.assert_always((faceIndices[1]==6),"getCodim2EntityLocalIndices");
	test.getCodim2EntityLocalIndices(11,faceIndices);
	logger.assert_always((faceIndices[0]==3),"getCodim2EntityLocalIndices");
	logger.assert_always((faceIndices[1]==7),"getCodim2EntityLocalIndices");
	faceIndices.resize(1);
	test.getCodim3EntityLocalIndices(0,faceIndices);
	logger.assert_always((faceIndices[0]==0),"getCodim3EntityLocalIndices");
	test.getCodim3EntityLocalIndices(1,faceIndices);
	logger.assert_always((faceIndices[0]==1),"getCodim3EntityLocalIndices");
	test.getCodim3EntityLocalIndices(2,faceIndices);
	logger.assert_always((faceIndices[0]==2),"getCodim3EntityLocalIndices");
	test.getCodim3EntityLocalIndices(3,faceIndices);
	logger.assert_always((faceIndices[0]==3),"getCodim3EntityLocalIndices");
	test.getCodim3EntityLocalIndices(4,faceIndices);
	logger.assert_always((faceIndices[0]==4),"getCodim3EntityLocalIndices");
	test.getCodim3EntityLocalIndices(5,faceIndices);
	logger.assert_always((faceIndices[0]==5),"getCodim3EntityLocalIndices");
	test.getCodim3EntityLocalIndices(6,faceIndices);
	logger.assert_always((faceIndices[0]==6),"getCodim3EntityLocalIndices");
	test.getCodim3EntityLocalIndices(7,faceIndices);
	logger.assert_always((faceIndices[0]==7),"getCodim3EntityLocalIndices");


	logger.assert_always((test.getGaussQuadratureRule(3)->order()>=3),"quadrature rules");
	logger.assert_always((test.getGaussQuadratureRule(5)->order()>=5),"quadrature rules");
	logger.assert_always((test.getGaussQuadratureRule(7)->order()>=7),"quadrature rules");
	logger.assert_always((test.getGaussQuadratureRule(9)->order()>=9),"quadrature rules");
	logger.assert_always((test.getGaussQuadratureRule(11)->order()>=11),"quadrature rules");

	//testing functionality of abstract parent classes

	logger.assert_always((test.getNumberOfNodes()==8),"number of nodes");
	logger.assert_always((test.getGeometryType()==Geometry::CUBE),"type of geometry");

	///\TODO if it is decided that getBasisFunctionValue and getBasisFucntionDerivative remain here, test them

	///\TODO testing that the refinement maps behave exactly like the forwarded calls of this class
}


