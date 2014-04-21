/*
 * 120ReferenceHyperCube_UnitTest.cpp
 *
 *  Created on: Mar 28, 2014
 *      Author: brinkf
 */

//naming convention: <Digit><ClassName>_UnitTest.cpp where <Digit> is a number that will make sure
//the unit tests are ordered such that the first failing unit test indicate the culprit class and
//other 'unit' tests may assume correct execution of all prior unit tests

#include "Geometry/ReferenceHypercube.hpp"
#include <iostream>
#include "cassert"

#include "Geometry/PointReference.hpp"
using Geometry::ReferenceHypercube;

int main(){
	ReferenceHypercube& test=ReferenceHypercube::Instance();

	Geometry::PointReference pTest(4);

	//testing basic functionality

	for(pTest[0]=-3.141;pTest[0]<-1.;pTest[0]+=0.1){
		for(pTest[1]=-3.1416;pTest[1]<3.1416;pTest[1]+=0.1){
			for(pTest[2]=-3.1416;pTest[2]<3.1416;pTest[2]+=0.1){
				for(pTest[3]=-3.1416;pTest[3]<3.1416;pTest[3]+=0.1){
					assert(("isInternalPoint",!test.isInternalPoint(pTest)));
				}
			}
		}
	}
	for(;pTest[0]<1;pTest[0]+=0.1){
		for(pTest[1]=-3.1417;pTest[1]<-1.;pTest[1]+=0.1){
			for(pTest[2]=-3.1416;pTest[2]<3.1416;pTest[2]+=0.1){
				for(pTest[3]=-3.1416;pTest[3]<3.1416;pTest[3]+=0.1){
					assert(("isInternalPoint",!test.isInternalPoint(pTest)));
				}
			}
		}
		for(;pTest[1]<1.;pTest[1]+=0.1){
			for(pTest[2]=-3.1416;pTest[2]<-1.;pTest[2]+=0.1){
				for(pTest[3]=-3.1416;pTest[3]<3.1416;pTest[3]+=0.1){
					assert(("isInternalPoint",!test.isInternalPoint(pTest)));
				}
			}
			for(;pTest[2]<1.;pTest[2]+=0.1){
				for(pTest[3]=-3.1416;pTest[3]<-1.;pTest[3]+=0.1){
					assert(("isInternalPoint",!test.isInternalPoint(pTest)));
				}
				for(;pTest[3]<1.;pTest[3]+=0.1){
					assert(("isInternalPoint",test.isInternalPoint(pTest)));
				}
				for(;pTest[3]<3.1416;pTest[3]+=0.1){
					assert(("isInternalPoint",!test.isInternalPoint(pTest)));
				}
			}
			for(;pTest[2]<3.141;pTest[2]+=0.1){
				for(pTest[3]=-3.1416;pTest[3]<3.1416;pTest[3]+=0.1){
					assert(("isInternalPoint",!test.isInternalPoint(pTest)));
				}
			}
		}
		for(;pTest[1]<3.141;pTest[1]+=0.1){
			for(pTest[2]=-3.1416;pTest[2]<3.1416;pTest[2]+=0.1){
				for(pTest[3]=-3.1416;pTest[3]<3.1416;pTest[3]+=0.1){
					assert(("isInternalPoint",!test.isInternalPoint(pTest)));
				}
			}
		}
	}
	for(;pTest[0]<3.141;pTest[0]+=0.1){
		for(pTest[1]=-3.1416;pTest[1]<3.1416;pTest[1]+=0.1){
			for(pTest[2]=-3.1416;pTest[2]<3.1416;pTest[2]+=0.1){
				for(pTest[3]=-3.1416;pTest[3]<3.1416;pTest[3]+=0.1){
					assert(("isInternalPoint",!test.isInternalPoint(pTest)));
				}
			}
		}
	}

	test.getCenter(pTest);
	assert(("getCenter",test.isInternalPoint(pTest)&&fabs(pTest[0])<1e-12&&fabs(pTest[1])<1e-12)&&fabs(pTest[2])<1e-12&&fabs(pTest[3])<1e-12);
	test.getNode(0,pTest);
	assert(("getNode 0",fabs(pTest[0]+1)<1e-12&&fabs(pTest[1]+1)<1e-12&&fabs(pTest[2]+1)<1e-12&&fabs(pTest[3]+1)<1e-12));
	test.getNode(1,pTest);
	assert(("getNode 1",fabs(pTest[0]-1)<1e-12&&fabs(pTest[1]+1)<1e-12&&fabs(pTest[2]+1)<1e-12&&fabs(pTest[3]+1)<1e-12));
	test.getNode(2,pTest);
	assert(("getNode 2",fabs(pTest[0]+1)<1e-12&&fabs(pTest[1]-1)<1e-12&&fabs(pTest[2]+1)<1e-12&&fabs(pTest[3]+1)<1e-12));
	test.getNode(3,pTest);
	assert(("getNode 3",fabs(pTest[0]-1)<1e-12&&fabs(pTest[1]-1)<1e-12&&fabs(pTest[2]+1)<1e-12&&fabs(pTest[3]+1)<1e-12));
	test.getNode(4,pTest);
	assert(("getNode 4",fabs(pTest[0]+1)<1e-12&&fabs(pTest[1]+1)<1e-12&&fabs(pTest[2]-1)<1e-12&&fabs(pTest[3]+1)<1e-12));
	test.getNode(5,pTest);
	assert(("getNode 5",fabs(pTest[0]-1)<1e-12&&fabs(pTest[1]+1)<1e-12&&fabs(pTest[2]-1)<1e-12&&fabs(pTest[3]+1)<1e-12));
	test.getNode(6,pTest);
	assert(("getNode 6",fabs(pTest[0]+1)<1e-12&&fabs(pTest[1]-1)<1e-12&&fabs(pTest[2]-1)<1e-12&&fabs(pTest[3]+1)<1e-12));
	test.getNode(7,pTest);
	assert(("getNode 7",fabs(pTest[0]-1)<1e-12&&fabs(pTest[1]-1)<1e-12&&fabs(pTest[2]-1)<1e-12&&fabs(pTest[3]+1)<1e-12));
	test.getNode(8,pTest);
	assert(("getNode 8",fabs(pTest[0]+1)<1e-12&&fabs(pTest[1]+1)<1e-12&&fabs(pTest[2]+1)<1e-12&&fabs(pTest[3]-1)<1e-12));
	test.getNode(9,pTest);
	assert(("getNode 9",fabs(pTest[0]-1)<1e-12&&fabs(pTest[1]+1)<1e-12&&fabs(pTest[2]+1)<1e-12&&fabs(pTest[3]-1)<1e-12));
	test.getNode(10,pTest);
	assert(("getNode 10",fabs(pTest[0]+1)<1e-12&&fabs(pTest[1]-1)<1e-12&&fabs(pTest[2]+1)<1e-12&&fabs(pTest[3]-1)<1e-12));
	test.getNode(11,pTest);
	assert(("getNode 11",fabs(pTest[0]-1)<1e-12&&fabs(pTest[1]-1)<1e-12&&fabs(pTest[2]+1)<1e-12&&fabs(pTest[3]-1)<1e-12));
	test.getNode(12,pTest);
	assert(("getNode 12",fabs(pTest[0]+1)<1e-12&&fabs(pTest[1]+1)<1e-12&&fabs(pTest[2]-1)<1e-12&&fabs(pTest[3]-1)<1e-12));
	test.getNode(13,pTest);
	assert(("getNode 13",fabs(pTest[0]-1)<1e-12&&fabs(pTest[1]+1)<1e-12&&fabs(pTest[2]-1)<1e-12&&fabs(pTest[3]-1)<1e-12));
	test.getNode(14,pTest);
	assert(("getNode 14",fabs(pTest[0]+1)<1e-12&&fabs(pTest[1]-1)<1e-12&&fabs(pTest[2]-1)<1e-12&&fabs(pTest[3]-1)<1e-12));
	test.getNode(15,pTest);
	assert(("getNode 15",fabs(pTest[0]-1)<1e-12&&fabs(pTest[1]-1)<1e-12&&fabs(pTest[2]-1)<1e-12&&fabs(pTest[3]-1)<1e-12));
	cout<<test.getName();

	assert(("getLocalNodeIndex 0",test.getLocalNodeIndex(0,0)==0));//the nodes of the face must always be specified IN THIS SPECIFIC ORDER
	assert(("getLocalNodeIndex 0",test.getLocalNodeIndex(0,1)==1));//im not sure if I like this myself, but this should at least verify
	assert(("getLocalNodeIndex 0",test.getLocalNodeIndex(0,2)==2));//that all face nodes are specified, none are specified twice
	assert(("getLocalNodeIndex 0",test.getLocalNodeIndex(0,3)==3));//and only face nodes are specified and the ordering of the nodes is consistent
	assert(("getLocalNodeIndex 0",test.getLocalNodeIndex(0,4)==4));//across function calls
	assert(("getLocalNodeIndex 0",test.getLocalNodeIndex(0,5)==5));
	assert(("getLocalNodeIndex 0",test.getLocalNodeIndex(0,6)==6));
	assert(("getLocalNodeIndex 0",test.getLocalNodeIndex(0,7)==7));
	assert(("getLocalNodeIndex 1",test.getLocalNodeIndex(1,0)==0));
	assert(("getLocalNodeIndex 1",test.getLocalNodeIndex(1,1)==1));
	assert(("getLocalNodeIndex 1",test.getLocalNodeIndex(1,2)==2));
	assert(("getLocalNodeIndex 1",test.getLocalNodeIndex(1,3)==3));
	assert(("getLocalNodeIndex 1",test.getLocalNodeIndex(1,4)==8));
	assert(("getLocalNodeIndex 1",test.getLocalNodeIndex(1,5)==9));
	assert(("getLocalNodeIndex 1",test.getLocalNodeIndex(1,6)==10));
	assert(("getLocalNodeIndex 1",test.getLocalNodeIndex(1,7)==11));
	assert(("getLocalNodeIndex 2",test.getLocalNodeIndex(2,0)==0));
	assert(("getLocalNodeIndex 2",test.getLocalNodeIndex(2,1)==1));
	assert(("getLocalNodeIndex 2",test.getLocalNodeIndex(2,2)==4));
	assert(("getLocalNodeIndex 2",test.getLocalNodeIndex(2,3)==5));
	assert(("getLocalNodeIndex 2",test.getLocalNodeIndex(2,4)==8));
	assert(("getLocalNodeIndex 2",test.getLocalNodeIndex(2,5)==9));
	assert(("getLocalNodeIndex 2",test.getLocalNodeIndex(2,6)==12));
	assert(("getLocalNodeIndex 2",test.getLocalNodeIndex(2,7)==13));
	assert(("getLocalNodeIndex 3",test.getLocalNodeIndex(3,0)==0));
	assert(("getLocalNodeIndex 3",test.getLocalNodeIndex(3,1)==2));
	assert(("getLocalNodeIndex 3",test.getLocalNodeIndex(3,2)==4));
	assert(("getLocalNodeIndex 3",test.getLocalNodeIndex(3,3)==6));
	assert(("getLocalNodeIndex 3",test.getLocalNodeIndex(3,4)==8));
	assert(("getLocalNodeIndex 3",test.getLocalNodeIndex(3,5)==10));
	assert(("getLocalNodeIndex 3",test.getLocalNodeIndex(3,6)==12));
	assert(("getLocalNodeIndex 3",test.getLocalNodeIndex(3,7)==14));
	assert(("getLocalNodeIndex 4",test.getLocalNodeIndex(4,0)==1));
	assert(("getLocalNodeIndex 4",test.getLocalNodeIndex(4,1)==3));
	assert(("getLocalNodeIndex 4",test.getLocalNodeIndex(4,2)==5));
	assert(("getLocalNodeIndex 4",test.getLocalNodeIndex(4,3)==7));
	assert(("getLocalNodeIndex 4",test.getLocalNodeIndex(4,4)==9));
	assert(("getLocalNodeIndex 4",test.getLocalNodeIndex(4,5)==11));
	assert(("getLocalNodeIndex 4",test.getLocalNodeIndex(4,6)==13));
	assert(("getLocalNodeIndex 4",test.getLocalNodeIndex(4,7)==15));
	assert(("getLocalNodeIndex 5",test.getLocalNodeIndex(5,0)==2));
	assert(("getLocalNodeIndex 5",test.getLocalNodeIndex(5,1)==3));
	assert(("getLocalNodeIndex 5",test.getLocalNodeIndex(5,2)==6));
	assert(("getLocalNodeIndex 5",test.getLocalNodeIndex(5,3)==7));
	assert(("getLocalNodeIndex 5",test.getLocalNodeIndex(5,4)==10));
	assert(("getLocalNodeIndex 5",test.getLocalNodeIndex(5,5)==11));
	assert(("getLocalNodeIndex 5",test.getLocalNodeIndex(5,6)==14));
	assert(("getLocalNodeIndex 5",test.getLocalNodeIndex(5,7)==15));
	assert(("getLocalNodeIndex 6",test.getLocalNodeIndex(6,0)==4));
	assert(("getLocalNodeIndex 6",test.getLocalNodeIndex(6,1)==5));
	assert(("getLocalNodeIndex 6",test.getLocalNodeIndex(6,2)==6));
	assert(("getLocalNodeIndex 6",test.getLocalNodeIndex(6,3)==7));
	assert(("getLocalNodeIndex 6",test.getLocalNodeIndex(6,4)==12));
	assert(("getLocalNodeIndex 6",test.getLocalNodeIndex(6,5)==13));
	assert(("getLocalNodeIndex 6",test.getLocalNodeIndex(6,6)==14));
	assert(("getLocalNodeIndex 6",test.getLocalNodeIndex(6,7)==15));
	assert(("getLocalNodeIndex 7",test.getLocalNodeIndex(7,0)==8));
	assert(("getLocalNodeIndex 7",test.getLocalNodeIndex(7,1)==9));
	assert(("getLocalNodeIndex 7",test.getLocalNodeIndex(7,2)==10));
	assert(("getLocalNodeIndex 7",test.getLocalNodeIndex(7,3)==11));
	assert(("getLocalNodeIndex 7",test.getLocalNodeIndex(7,4)==12));
	assert(("getLocalNodeIndex 7",test.getLocalNodeIndex(7,5)==13));
	assert(("getLocalNodeIndex 7",test.getLocalNodeIndex(7,6)==14));
	assert(("getLocalNodeIndex 7",test.getLocalNodeIndex(7,7)==15));

	cout<<test;

	//testing mappings and quadrature rules

	std::vector<unsigned int> faceIndices(8);
	//there is no 5D element so codim0mappings are not needed


	assert(("higher codimensional entities",test.getNrOfCodim1Entities()==8&&test.getNrOfCodim2Entities()==24)&&test.getNrOfCodim3Entities()==32);
	assert(("getCodim1ReferenceGeometry",test.getCodim1ReferenceGeometry(0)==&Geometry::ReferenceCube::Instance()&&
										 test.getCodim1ReferenceGeometry(1)==&Geometry::ReferenceCube::Instance()&&
										 test.getCodim1ReferenceGeometry(2)==&Geometry::ReferenceCube::Instance()&&
										 test.getCodim1ReferenceGeometry(3)==&Geometry::ReferenceCube::Instance()&&
										 test.getCodim1ReferenceGeometry(4)==&Geometry::ReferenceCube::Instance()&&
										 test.getCodim1ReferenceGeometry(5)==&Geometry::ReferenceCube::Instance()&&
										 test.getCodim1ReferenceGeometry(6)==&Geometry::ReferenceCube::Instance()&&
										 test.getCodim1ReferenceGeometry(7)==&Geometry::ReferenceCube::Instance()));
	assert(("getCodim2ReferenceGeometry",test.getCodim2ReferenceGeometry(0)==&Geometry::ReferenceSquare::Instance()&&
										 test.getCodim2ReferenceGeometry(1)==&Geometry::ReferenceSquare::Instance()&&
										 test.getCodim2ReferenceGeometry(2)==&Geometry::ReferenceSquare::Instance()&&
										 test.getCodim2ReferenceGeometry(3)==&Geometry::ReferenceSquare::Instance()&&
										 test.getCodim2ReferenceGeometry(4)==&Geometry::ReferenceSquare::Instance()&&
										 test.getCodim2ReferenceGeometry(5)==&Geometry::ReferenceSquare::Instance()&&
										 test.getCodim2ReferenceGeometry(6)==&Geometry::ReferenceSquare::Instance()&&
										 test.getCodim2ReferenceGeometry(7)==&Geometry::ReferenceSquare::Instance()&&
										 test.getCodim2ReferenceGeometry(8)==&Geometry::ReferenceSquare::Instance()&&
										 test.getCodim2ReferenceGeometry(9)==&Geometry::ReferenceSquare::Instance()&&
										 test.getCodim2ReferenceGeometry(10)==&Geometry::ReferenceSquare::Instance()&&
										 test.getCodim2ReferenceGeometry(11)==&Geometry::ReferenceSquare::Instance()&&
										 test.getCodim2ReferenceGeometry(12)==&Geometry::ReferenceSquare::Instance()&&
										 test.getCodim2ReferenceGeometry(13)==&Geometry::ReferenceSquare::Instance()&&
										 test.getCodim2ReferenceGeometry(14)==&Geometry::ReferenceSquare::Instance()&&
										 test.getCodim2ReferenceGeometry(15)==&Geometry::ReferenceSquare::Instance()&&
										 test.getCodim2ReferenceGeometry(16)==&Geometry::ReferenceSquare::Instance()&&
										 test.getCodim2ReferenceGeometry(17)==&Geometry::ReferenceSquare::Instance()&&
										 test.getCodim2ReferenceGeometry(18)==&Geometry::ReferenceSquare::Instance()&&
										 test.getCodim2ReferenceGeometry(19)==&Geometry::ReferenceSquare::Instance()&&
										 test.getCodim2ReferenceGeometry(20)==&Geometry::ReferenceSquare::Instance()&&
										 test.getCodim2ReferenceGeometry(21)==&Geometry::ReferenceSquare::Instance()&&
										 test.getCodim2ReferenceGeometry(22)==&Geometry::ReferenceSquare::Instance()&&
										 test.getCodim2ReferenceGeometry(23)==&Geometry::ReferenceSquare::Instance()));
	assert(("getCodim1MappingPtr"&&test.getCodim1MappingPtr(0)==&Geometry::MappingToRefCubeToHypercube0::Instance()));
	assert(("getCodim1MappingPtr",test.getCodim1MappingPtr(1)==&Geometry::MappingToRefCubeToHypercube1::Instance()));
	assert(("getCodim1MappingPtr",test.getCodim1MappingPtr(2)==&Geometry::MappingToRefCubeToHypercube2::Instance()));
	assert(("getCodim1MappingPtr",test.getCodim1MappingPtr(3)==&Geometry::MappingToRefCubeToHypercube3::Instance()));
	assert(("getCodim1MappingPtr",test.getCodim1MappingPtr(4)==&Geometry::MappingToRefCubeToHypercube4::Instance()));
	assert(("getCodim1MappingPtr",test.getCodim1MappingPtr(5)==&Geometry::MappingToRefCubeToHypercube5::Instance()));
	assert(("getCodim1MappingPtr",test.getCodim1MappingPtr(6)==&Geometry::MappingToRefCubeToHypercube6::Instance()));
	assert(("getCodim1MappingPtr",test.getCodim1MappingPtr(7)==&Geometry::MappingToRefCubeToHypercube7::Instance()));
	test.getCodim1EntityLocalIndices(0,faceIndices);
	assert(("getCodim1EntityLocalIndices",faceIndices[0]==test.getLocalNodeIndex(0,0)));
	assert(("getCodim1EntityLocalIndices",faceIndices[1]==test.getLocalNodeIndex(0,1)));
	assert(("getCodim1EntityLocalIndices",faceIndices[2]==test.getLocalNodeIndex(0,2)));
	assert(("getCodim1EntityLocalIndices",faceIndices[3]==test.getLocalNodeIndex(0,3)));
	assert(("getCodim1EntityLocalIndices",faceIndices[4]==test.getLocalNodeIndex(0,4)));
	assert(("getCodim1EntityLocalIndices",faceIndices[5]==test.getLocalNodeIndex(0,5)));
	assert(("getCodim1EntityLocalIndices",faceIndices[6]==test.getLocalNodeIndex(0,6)));
	assert(("getCodim1EntityLocalIndices",faceIndices[7]==test.getLocalNodeIndex(0,7)));
	test.getCodim1EntityLocalIndices(1,faceIndices);
	assert(("getCodim1EntityLocalIndices",faceIndices[0]==test.getLocalNodeIndex(1,0)));
	assert(("getCodim1EntityLocalIndices",faceIndices[1]==test.getLocalNodeIndex(1,1)));
	assert(("getCodim1EntityLocalIndices",faceIndices[2]==test.getLocalNodeIndex(1,2)));
	assert(("getCodim1EntityLocalIndices",faceIndices[3]==test.getLocalNodeIndex(1,3)));
	assert(("getCodim1EntityLocalIndices",faceIndices[4]==test.getLocalNodeIndex(1,4)));
	assert(("getCodim1EntityLocalIndices",faceIndices[5]==test.getLocalNodeIndex(1,5)));
	assert(("getCodim1EntityLocalIndices",faceIndices[6]==test.getLocalNodeIndex(1,6)));
	assert(("getCodim1EntityLocalIndices",faceIndices[7]==test.getLocalNodeIndex(1,7)));
	test.getCodim1EntityLocalIndices(2,faceIndices);
	assert(("getCodim1EntityLocalIndices",faceIndices[0]==test.getLocalNodeIndex(2,0)));
	assert(("getCodim1EntityLocalIndices",faceIndices[1]==test.getLocalNodeIndex(2,1)));
	assert(("getCodim1EntityLocalIndices",faceIndices[2]==test.getLocalNodeIndex(2,2)));
	assert(("getCodim1EntityLocalIndices",faceIndices[3]==test.getLocalNodeIndex(2,3)));
	assert(("getCodim1EntityLocalIndices",faceIndices[4]==test.getLocalNodeIndex(2,4)));
	assert(("getCodim1EntityLocalIndices",faceIndices[5]==test.getLocalNodeIndex(2,5)));
	assert(("getCodim1EntityLocalIndices",faceIndices[6]==test.getLocalNodeIndex(2,6)));
	assert(("getCodim1EntityLocalIndices",faceIndices[7]==test.getLocalNodeIndex(2,7)));
	test.getCodim1EntityLocalIndices(3,faceIndices);
	assert(("getCodim1EntityLocalIndices",faceIndices[0]==test.getLocalNodeIndex(3,0)));
	assert(("getCodim1EntityLocalIndices",faceIndices[1]==test.getLocalNodeIndex(3,1)));
	assert(("getCodim1EntityLocalIndices",faceIndices[2]==test.getLocalNodeIndex(3,2)));
	assert(("getCodim1EntityLocalIndices",faceIndices[3]==test.getLocalNodeIndex(3,3)));
	assert(("getCodim1EntityLocalIndices",faceIndices[4]==test.getLocalNodeIndex(3,4)));
	assert(("getCodim1EntityLocalIndices",faceIndices[5]==test.getLocalNodeIndex(3,5)));
	assert(("getCodim1EntityLocalIndices",faceIndices[6]==test.getLocalNodeIndex(3,6)));
	assert(("getCodim1EntityLocalIndices",faceIndices[7]==test.getLocalNodeIndex(3,7)));
	test.getCodim1EntityLocalIndices(4,faceIndices);
	assert(("getCodim1EntityLocalIndices",faceIndices[0]==test.getLocalNodeIndex(4,0)));
	assert(("getCodim1EntityLocalIndices",faceIndices[1]==test.getLocalNodeIndex(4,1)));
	assert(("getCodim1EntityLocalIndices",faceIndices[2]==test.getLocalNodeIndex(4,2)));
	assert(("getCodim1EntityLocalIndices",faceIndices[3]==test.getLocalNodeIndex(4,3)));
	assert(("getCodim1EntityLocalIndices",faceIndices[4]==test.getLocalNodeIndex(4,4)));
	assert(("getCodim1EntityLocalIndices",faceIndices[5]==test.getLocalNodeIndex(4,5)));
	assert(("getCodim1EntityLocalIndices",faceIndices[6]==test.getLocalNodeIndex(4,6)));
	assert(("getCodim1EntityLocalIndices",faceIndices[7]==test.getLocalNodeIndex(4,7)));
	test.getCodim1EntityLocalIndices(5,faceIndices);
	assert(("getCodim1EntityLocalIndices",faceIndices[0]==test.getLocalNodeIndex(5,0)));
	assert(("getCodim1EntityLocalIndices",faceIndices[1]==test.getLocalNodeIndex(5,1)));
	assert(("getCodim1EntityLocalIndices",faceIndices[2]==test.getLocalNodeIndex(5,2)));
	assert(("getCodim1EntityLocalIndices",faceIndices[3]==test.getLocalNodeIndex(5,3)));
	assert(("getCodim1EntityLocalIndices",faceIndices[4]==test.getLocalNodeIndex(5,4)));
	assert(("getCodim1EntityLocalIndices",faceIndices[5]==test.getLocalNodeIndex(5,5)));
	assert(("getCodim1EntityLocalIndices",faceIndices[6]==test.getLocalNodeIndex(5,6)));
	assert(("getCodim1EntityLocalIndices",faceIndices[7]==test.getLocalNodeIndex(5,7)));
	test.getCodim1EntityLocalIndices(6,faceIndices);
	assert(("getCodim1EntityLocalIndices",faceIndices[0]==test.getLocalNodeIndex(6,0)));
	assert(("getCodim1EntityLocalIndices",faceIndices[1]==test.getLocalNodeIndex(6,1)));
	assert(("getCodim1EntityLocalIndices",faceIndices[2]==test.getLocalNodeIndex(6,2)));
	assert(("getCodim1EntityLocalIndices",faceIndices[3]==test.getLocalNodeIndex(6,3)));
	assert(("getCodim1EntityLocalIndices",faceIndices[4]==test.getLocalNodeIndex(6,4)));
	assert(("getCodim1EntityLocalIndices",faceIndices[5]==test.getLocalNodeIndex(6,5)));
	assert(("getCodim1EntityLocalIndices",faceIndices[6]==test.getLocalNodeIndex(6,6)));
	assert(("getCodim1EntityLocalIndices",faceIndices[7]==test.getLocalNodeIndex(6,7)));
	test.getCodim1EntityLocalIndices(7,faceIndices);
	assert(("getCodim1EntityLocalIndices",faceIndices[0]==test.getLocalNodeIndex(7,0)));
	assert(("getCodim1EntityLocalIndices",faceIndices[1]==test.getLocalNodeIndex(7,1)));
	assert(("getCodim1EntityLocalIndices",faceIndices[2]==test.getLocalNodeIndex(7,2)));
	assert(("getCodim1EntityLocalIndices",faceIndices[3]==test.getLocalNodeIndex(7,3)));
	assert(("getCodim1EntityLocalIndices",faceIndices[4]==test.getLocalNodeIndex(7,4)));
	assert(("getCodim1EntityLocalIndices",faceIndices[5]==test.getLocalNodeIndex(7,5)));
	assert(("getCodim1EntityLocalIndices",faceIndices[6]==test.getLocalNodeIndex(7,6)));
	assert(("getCodim1EntityLocalIndices",faceIndices[7]==test.getLocalNodeIndex(7,7)));

	//other codimensions are not implemented


	assert(("quadrature rules",test.getGaussQuadratureRule(3)->order()>=3));
	//assert(("quadrature rules",test.getGaussQuadratureRule(5)->order()>=5));///\TODO implement more quadrature rules
	//assert(("quadrature rules",test.getGaussQuadratureRule(7)->order()>=7));
	//assert(("quadrature rules",test.getGaussQuadratureRule(9)->order()>=9));
	//assert(("quadrature rules",test.getGaussQuadratureRule(11)->order()>=11));

	//testing functionality of abstract parent classes

	assert(("number of nodes",test.getNumberOfNodes()==16));
	assert(("type of geometry",test.getGeometryType()==Geometry::HYPERCUBE));

	///\TODO if it is decided that getBasisFunctionValue and getBasisFucntionDerivative remain here, test them

	///\TODO testing that the refinement maps behave exactly like the forwarded calls of this class
}





