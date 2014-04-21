/*
 * 130PhysicalLine_UnitTest.cpp
 *
 *  Created on: Mar 28, 2014
 *      Author: brinkf
 */

//naming convention: <Digit><ClassName>_UnitTest.cpp where <Digit> is a number that will make sure
//the unit tests are ordered such that the first failing unit test indicate the culprit class and
//other 'unit' tests may assume correct execution of all prior unit tests

#include "Geometry/PhysicalLine.hpp"
#include "cassert"

using Geometry::PhysicalLine;

int main(){

	std::vector<unsigned int> pointIndexes;
	std::vector<Geometry::PointPhysical> nodes;

	Geometry::PointPhysical point(1);

	pointIndexes.push_back(4);
	pointIndexes.push_back(7);

	for(double i=0.;i<10;++i){
		point[0]=1.+i/10.;
		nodes.push_back(point);
	}

	PhysicalLine test(pointIndexes,nodes,&Geometry::ReferenceLine::Instance());

	cout<<test;

	pointIndexes = test.getNodeIndexes();

	assert(("getNodeIndexes",pointIndexes[0]==4&&pointIndexes[1]==7));
	assert(("getNodes",nodes==test.getNodes()));
	assert(("getNodeIndex",test.getNodeIndex(0)==4&&test.getNodeIndex(1)==7));

	cout<<test.getName();

	point = *test.getNodePtr(test.getNodeIndex(0));
	assert(("getNodePtr",fabs(point[0]-1.4)<1e-12));
	point = *test.getNodePtr(test.getNodeIndex(1));
	assert(("getNodePtr",fabs(point[0]-1.7)<1e-12));

	assert(("getNumberOfNodes",test.getNumberOfNodes()==2));

	test.getNodeCoordinates(0,point);
	assert(("getNodeCoordinates",fabs(point[0]-1.4)<1e-12));
	test.getNodeCoordinates(1,point);
	assert(("getNodeCoordinates",fabs(point[0]-1.7)<1e-12));

	test.getLocalNodeCoordinates(0,point);
	assert(("getLocalNodeCoordinates",fabs(point[0]-1.4)<1e-12));
	test.getLocalNodeCoordinates(1,point);
	assert(("getLocalNodeCoordinates",fabs(point[0]-1.7)<1e-12));

	for(double i=0;i<10;++i){
		test.getGlobalNodeCoordinates(i,point);
		assert(("getGlobalNodeCoordinates",fabs(point[0]-1.-i/10.)<1e-12));
	}

	pointIndexes.resize(1);

	test.getGlobalFaceNodeIndices(0,pointIndexes);
	assert(("getGlobalFaceNodeIndices",pointIndexes[0]==4));
	test.getGlobalFaceNodeIndices(1,pointIndexes);
	assert(("getGlobalFaceNodeIndices",pointIndexes[0]==7));

	test.getLocalFaceNodeIndices(0,pointIndexes);
	assert(("getLocalFaceNodeIndices",pointIndexes[0]==0));
	test.getLocalFaceNodeIndices(1,pointIndexes);
	assert(("getLocalFaceNodeIndices",pointIndexes[0]==1));

	assert(("getNrOfFaces",test.getNrOfFaces()==2));

	assert(("getRefGeometry",test.getRefGeometry()==&Geometry::ReferenceLine::Instance()));


	return 0;
}

