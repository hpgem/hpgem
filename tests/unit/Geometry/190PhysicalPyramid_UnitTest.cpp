/*
 * 190PhysicalPyramid_UnitTest.cpp
 *
 *  Created on: Mar 31, 2014
 *      Author: brinkf
 */

//naming convention: <Digit><ClassName>_UnitTest.cpp where <Digit> is a number that will make sure
//the unit tests are ordered such that the first failing unit test indicate the culprit class and
//other 'unit' tests may assume correct execution of all prior unit tests

#include "Geometry/PhysicalPyramid.hpp"
#include "cassert"

using Geometry::PhysicalPyramid;

int main(){

	std::vector<unsigned int> pointIndexes;
	std::vector<Geometry::PointPhysical> nodes;

	Geometry::PointPhysical point(3);

	pointIndexes.push_back(4);
	pointIndexes.push_back(7);
	pointIndexes.push_back(10);
	pointIndexes.push_back(11);
	pointIndexes.push_back(12);

	for(double i=0.;i<10;++i){
		point[0]=1.+i/10.;
		point[1]=2.+i/10.;
		point[2]=3.+i/10.;
		nodes.push_back(point);
	}

	point[0]=3.5;
	point[1]=4.6;
	point[2]=5.4;
	nodes.push_back(point);
	point[0]=6.7;
	point[1]=2.8;
	point[2]=5.7;
	nodes.push_back(point);
	point[0]=1.4;
	point[1]=2.4;
	point[2]=5.4;
	nodes.push_back(point);

	PhysicalPyramid test(pointIndexes,nodes,&Geometry::ReferencePyramid::Instance());

	cout<<test;

	pointIndexes = test.getNodeIndexes();

	assert(("getNodeIndexes",pointIndexes[0]==4&&pointIndexes[1]==7&&pointIndexes[2]==10&&pointIndexes[3]==11&&pointIndexes[4]==12));
	assert(("getNodes",nodes==test.getNodes()));
	assert(("getNodeIndex",test.getNodeIndex(0)==4&&test.getNodeIndex(1)==7&&test.getNodeIndex(2)==10&&test.getNodeIndex(3)==11&&test.getNodeIndex(4)==12));

	cout<<test.getName();

	point = *test.getNodePtr(test.getNodeIndex(0));
	assert(("getNodePtr",fabs(point[0]-1.4)<1e-12));
	assert(("getNodePtr",fabs(point[1]-2.4)<1e-12));
	assert(("getNodePtr",fabs(point[2]-3.4)<1e-12));
	point = *test.getNodePtr(test.getNodeIndex(1));
	assert(("getNodePtr",fabs(point[0]-1.7)<1e-12));
	assert(("getNodePtr",fabs(point[1]-2.7)<1e-12));
	assert(("getNodePtr",fabs(point[2]-3.7)<1e-12));

	assert(("getNumberOfNodes",test.getNumberOfNodes()==5));

	test.getNodeCoordinates(0,point);
	assert(("getNodeCoordinates",fabs(point[0]-1.4)<1e-12));
	assert(("getNodeCoordinates",fabs(point[1]-2.4)<1e-12));
	assert(("getNodeCoordinates",fabs(point[2]-3.4)<1e-12));
	test.getNodeCoordinates(1,point);
	assert(("getNodeCoordinates",fabs(point[0]-1.7)<1e-12));
	assert(("getNodeCoordinates",fabs(point[1]-2.7)<1e-12));
	assert(("getNodeCoordinates",fabs(point[2]-3.7)<1e-12));
	test.getNodeCoordinates(2,point);
	assert(("getNodeCoordinates",fabs(point[0]-3.5)<1e-12));
	assert(("getNodeCoordinates",fabs(point[1]-4.6)<1e-12));
	assert(("getNodeCoordinates",fabs(point[2]-5.4)<1e-12));
	test.getNodeCoordinates(3,point);
	assert(("getNodeCoordinates",fabs(point[0]-6.7)<1e-12));
	assert(("getNodeCoordinates",fabs(point[1]-2.8)<1e-12));
	assert(("getNodeCoordinates",fabs(point[2]-5.7)<1e-12));
	test.getNodeCoordinates(4,point);
	assert(("getNodeCoordinates",fabs(point[0]-1.4)<1e-12));
	assert(("getNodeCoordinates",fabs(point[1]-2.4)<1e-12));
	assert(("getNodeCoordinates",fabs(point[2]-5.4)<1e-12));

	test.getLocalNodeCoordinates(0,point);
	assert(("getLocalNodeCoordinates",fabs(point[0]-1.4)<1e-12));
	assert(("getLocalNodeCoordinates",fabs(point[1]-2.4)<1e-12));
	assert(("getLocalNodeCoordinates",fabs(point[2]-3.4)<1e-12));
	test.getLocalNodeCoordinates(1,point);
	assert(("getLocalNodeCoordinates",fabs(point[0]-1.7)<1e-12));
	assert(("getLocalNodeCoordinates",fabs(point[1]-2.7)<1e-12));
	assert(("getLocalNodeCoordinates",fabs(point[2]-3.7)<1e-12));
	test.getLocalNodeCoordinates(2,point);
	assert(("getLocalNodeCoordinates",fabs(point[0]-3.5)<1e-12));
	assert(("getLocalNodeCoordinates",fabs(point[1]-4.6)<1e-12));
	assert(("getLocalNodeCoordinates",fabs(point[2]-5.4)<1e-12));
	test.getLocalNodeCoordinates(3,point);
	assert(("getLocalNodeCoordinates",fabs(point[0]-6.7)<1e-12));
	assert(("getLocalNodeCoordinates",fabs(point[1]-2.8)<1e-12));
	assert(("getLocalNodeCoordinates",fabs(point[2]-5.7)<1e-12));
	test.getLocalNodeCoordinates(4,point);
	assert(("getLocalNodeCoordinates",fabs(point[0]-1.4)<1e-12));
	assert(("getLocalNodeCoordinates",fabs(point[1]-2.4)<1e-12));
	assert(("getLocalNodeCoordinates",fabs(point[2]-5.4)<1e-12));

	for(double i=0;i<10;++i){
		test.getGlobalNodeCoordinates(i,point);
		assert(("getGlobalNodeCoordinates",fabs(point[0]-1.-i/10.)<1e-12));
		assert(("getGlobalNodeCoordinates",fabs(point[1]-2.-i/10.)<1e-12));
		assert(("getGlobalNodeCoordinates",fabs(point[2]-3.-i/10.)<1e-12));
	}

	pointIndexes.resize(4);
	for(int i=0;i<1;++i){
		test.getGlobalFaceNodeIndices(i,pointIndexes);
		for(int j=0;j<4;++j){
			assert(("getGlobalFaceNodeIndices",pointIndexes[j]==test.getNodeIndex(test.getRefGeometry()->getLocalNodeIndex(i,j))));
		}
	}
	pointIndexes.resize(3);
	for(int i=1;i<5;++i){
		test.getGlobalFaceNodeIndices(i,pointIndexes);
		for(int j=0;j<3;++j){
			assert(("getGlobalFaceNodeIndices",pointIndexes[j]==test.getNodeIndex(test.getRefGeometry()->getLocalNodeIndex(i,j))));
		}
	}

	pointIndexes.resize(4);
	for(int i=0;i<1;++i){
		test.getLocalFaceNodeIndices(i,pointIndexes);
		for(int j=0;j<4;++j){
			assert(("getLocalFaceNodeIndices",pointIndexes[j]==test.getRefGeometry()->getLocalNodeIndex(i,j)));
		}
	}
	pointIndexes.resize(3);
	for(int i=1;i<5;++i){
		test.getLocalFaceNodeIndices(i,pointIndexes);
		for(int j=0;j<3;++j){
			assert(("getLocalFaceNodeIndices",pointIndexes[j]==test.getRefGeometry()->getLocalNodeIndex(i,j)));
		}
	}

	assert(("getNrOfFaces",test.getNrOfFaces()==5));

	assert(("getRefGeometry",test.getRefGeometry()==&Geometry::ReferencePyramid::Instance()));


	return 0;
}


