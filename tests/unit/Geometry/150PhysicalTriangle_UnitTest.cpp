/*
 * 150PhysicalTriangle_UnitTest.cpp
 *
 *  Created on: Mar 31, 2014
 *      Author: brinkf
 */

//naming convention: <Digit><ClassName>_UnitTest.cpp where <Digit> is a number that will make sure
//the unit tests are ordered such that the first failing unit test indicate the culprit class and
//other 'unit' tests may assume correct execution of all prior unit tests

#include "Geometry/PhysicalTriangle.hpp"
#include "cassert"

using Geometry::PhysicalTriangle;

int main(){

	std::vector<unsigned int> pointIndexes;
	std::vector<Geometry::PointPhysical> nodes;

	Geometry::PointPhysical point(2);

	pointIndexes.push_back(4);
	pointIndexes.push_back(7);
	pointIndexes.push_back(10);

	for(double i=0.;i<10;++i){
		point[0]=1.+i/10.;
		point[1]=2.+i/10.;
		nodes.push_back(point);
	}

	point[0]=3.5;
	point[1]=4.6;
	nodes.push_back(point);

	PhysicalTriangle test(pointIndexes,nodes,&Geometry::ReferenceTriangle::Instance());

	cout<<test;

	pointIndexes = test.getNodeIndexes();

	assert(("getNodeIndexes",pointIndexes[0]==4&&pointIndexes[1]==7&&pointIndexes[2]==10));
	assert(("getNodes",nodes==test.getNodes()));
	assert(("getNodeIndex",test.getNodeIndex(0)==4&&test.getNodeIndex(1)==7&&test.getNodeIndex(2)==10));

	cout<<test.getName();

	point = *test.getNodePtr(test.getNodeIndex(0));
	assert(("getNodePtr",fabs(point[0]-1.4)<1e-12));
	assert(("getNodePtr",fabs(point[1]-2.4)<1e-12));
	point = *test.getNodePtr(test.getNodeIndex(1));
	assert(("getNodePtr",fabs(point[0]-1.7)<1e-12));
	assert(("getNodePtr",fabs(point[1]-2.7)<1e-12));

	assert(("getNumberOfNodes",test.getNumberOfNodes()==3));

	test.getNodeCoordinates(0,point);
	assert(("getNodeCoordinates",fabs(point[0]-1.4)<1e-12));
	assert(("getNodeCoordinates",fabs(point[1]-2.4)<1e-12));
	test.getNodeCoordinates(1,point);
	assert(("getNodeCoordinates",fabs(point[0]-1.7)<1e-12));
	assert(("getNodeCoordinates",fabs(point[1]-2.7)<1e-12));
	test.getNodeCoordinates(2,point);
	assert(("getNodeCoordinates",fabs(point[0]-3.5)<1e-12));
	assert(("getNodeCoordinates",fabs(point[1]-4.6)<1e-12));

	test.getLocalNodeCoordinates(0,point);
	assert(("getLocalNodeCoordinates",fabs(point[0]-1.4)<1e-12));
	assert(("getLocalNodeCoordinates",fabs(point[1]-2.4)<1e-12));
	test.getLocalNodeCoordinates(1,point);
	assert(("getLocalNodeCoordinates",fabs(point[0]-1.7)<1e-12));
	assert(("getLocalNodeCoordinates",fabs(point[1]-2.7)<1e-12));
	test.getLocalNodeCoordinates(2,point);
	assert(("getLocalNodeCoordinates",fabs(point[0]-3.5)<1e-12));
	assert(("getLocalNodeCoordinates",fabs(point[1]-4.6)<1e-12));

	for(double i=0;i<10;++i){
		test.getGlobalNodeCoordinates(i,point);
		assert(("getGlobalNodeCoordinates",fabs(point[0]-1.-i/10.)<1e-12));
		assert(("getGlobalNodeCoordinates",fabs(point[1]-2.-i/10.)<1e-12));
	}

	pointIndexes.resize(2);

	test.getGlobalFaceNodeIndices(0,pointIndexes);
	assert(("getGlobalFaceNodeIndices",pointIndexes[0]==test.getNodeIndex(test.getRefGeometry()->getLocalNodeIndex(0,0))));
	assert(("getGlobalFaceNodeIndices",pointIndexes[1]==test.getNodeIndex(test.getRefGeometry()->getLocalNodeIndex(0,1))));
	test.getGlobalFaceNodeIndices(1,pointIndexes);
	assert(("getGlobalFaceNodeIndices",pointIndexes[0]==test.getNodeIndex(test.getRefGeometry()->getLocalNodeIndex(1,0))));
	assert(("getGlobalFaceNodeIndices",pointIndexes[1]==test.getNodeIndex(test.getRefGeometry()->getLocalNodeIndex(1,1))));
	test.getGlobalFaceNodeIndices(2,pointIndexes);
	assert(("getGlobalFaceNodeIndices",pointIndexes[0]==test.getNodeIndex(test.getRefGeometry()->getLocalNodeIndex(2,0))));
	assert(("getGlobalFaceNodeIndices",pointIndexes[1]==test.getNodeIndex(test.getRefGeometry()->getLocalNodeIndex(2,1))));

	test.getLocalFaceNodeIndices(0,pointIndexes);
	assert(("getLocalFaceNodeIndices",pointIndexes[0]==test.getRefGeometry()->getLocalNodeIndex(0,0)));
	assert(("getLocalFaceNodeIndices",pointIndexes[1]==test.getRefGeometry()->getLocalNodeIndex(0,1)));
	test.getLocalFaceNodeIndices(1,pointIndexes);
	assert(("getLocalFaceNodeIndices",pointIndexes[0]==test.getRefGeometry()->getLocalNodeIndex(1,0)));
	assert(("getLocalFaceNodeIndices",pointIndexes[1]==test.getRefGeometry()->getLocalNodeIndex(1,1)));
	test.getLocalFaceNodeIndices(2,pointIndexes);
	assert(("getLocalFaceNodeIndices",pointIndexes[0]==test.getRefGeometry()->getLocalNodeIndex(2,0)));
	assert(("getLocalFaceNodeIndices",pointIndexes[1]==test.getRefGeometry()->getLocalNodeIndex(2,1)));

	assert(("getNrOfFaces",test.getNrOfFaces()==3));

	assert(("getRefGeometry",test.getRefGeometry()==&Geometry::ReferenceTriangle::Instance()));


	return 0;
}
