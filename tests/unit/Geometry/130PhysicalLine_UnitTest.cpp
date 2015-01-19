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

#include "Geometry/PhysicalLine.hpp"
#include "cassert"
#include "Geometry/PointPhysical.hpp"
#include "Geometry/PointReference.hpp"
#include "Geometry/ReferenceLine.hpp"

#include <cmath>
using Geometry::PhysicalLine;

int main(){

	std::vector<std::size_t> pointIndexes;
	std::vector<Geometry::PointPhysical> nodes;

	Geometry::PointPhysical point(1);

	pointIndexes.push_back(4);
	pointIndexes.push_back(7);

	for(double i=0.;i<10;++i){
		point[0]=1.+i/10.;
		nodes.push_back(point);
	}

	PhysicalLine test(pointIndexes,nodes,&Geometry::ReferenceLine::Instance());

	std::cout<<test;

	pointIndexes = test.getNodeIndexes();

	assert(("getNodeIndexes",pointIndexes[0]==4&&pointIndexes[1]==7));
	assert(("getNodes",nodes==test.getNodes()));
	assert(("getNodeIndex",test.getNodeIndex(0)==4&&test.getNodeIndex(1)==7));

	std::cout<<test.getName();

	point = *test.getNodePtr(test.getNodeIndex(0));
	assert(("getNodePtr",std::abs(point[0]-1.4)<1e-12));
	point = *test.getNodePtr(test.getNodeIndex(1));
	assert(("getNodePtr",std::abs(point[0]-1.7)<1e-12));

	assert(("getNumberOfNodes",test.getNumberOfNodes()==2));

	test.getNodeCoordinates(0,point);
	assert(("getNodeCoordinates",std::abs(point[0]-1.4)<1e-12));
	test.getNodeCoordinates(1,point);
	assert(("getNodeCoordinates",std::abs(point[0]-1.7)<1e-12));

	test.getLocalNodeCoordinates(0,point);
	assert(("getLocalNodeCoordinates",std::abs(point[0]-1.4)<1e-12));
	test.getLocalNodeCoordinates(1,point);
	assert(("getLocalNodeCoordinates",std::abs(point[0]-1.7)<1e-12));

	for(double i=0;i<10;++i){
		test.getGlobalNodeCoordinates(i,point);
		assert(("getGlobalNodeCoordinates",std::abs(point[0]-1.-i/10.)<1e-12));
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

