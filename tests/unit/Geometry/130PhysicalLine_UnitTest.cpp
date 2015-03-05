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
#include "Logger.h"
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

	for(double i = 0.0; i < 1; i += 0.1){
		point[0]=1.+i;
		nodes.push_back(point);
	}

	PhysicalLine test(pointIndexes,nodes);

	std::cout<<test;

	pointIndexes = test.getNodeIndexes();

	logger.assert_always((pointIndexes[0]==4&&pointIndexes[1]==7),"getNodeIndexes");
	logger.assert_always((nodes==test.getNodes()),"getNodes");
	logger.assert_always((test.getNodeIndex(0)==4&&test.getNodeIndex(1)==7),"getNodeIndex");

	std::cout<<test.getName();

	point = *test.getNodePtr(test.getNodeIndex(0));
	logger.assert_always((std::abs(point[0]-1.4)<1e-12),"getNodePtr");
	point = *test.getNodePtr(test.getNodeIndex(1));
	logger.assert_always((std::abs(point[0]-1.7)<1e-12),"getNodePtr");

	logger.assert_always((test.getNumberOfNodes()==2),"getNumberOfNodes");

	point = test.getNodeCoordinates(0);
	logger.assert_always((std::abs(point[0]-1.4)<1e-12),"getNodeCoordinates");
	point = test.getNodeCoordinates(1);
	logger.assert_always((std::abs(point[0]-1.7)<1e-12),"getNodeCoordinates");

	point = test.getLocalNodeCoordinates(0);
	logger.assert_always((std::abs(point[0]-1.4)<1e-12),"getLocalNodeCoordinates");
	point = test.getLocalNodeCoordinates(1);
	logger.assert_always((std::abs(point[0]-1.7)<1e-12),"getLocalNodeCoordinates");

	for(std::size_t i=0;i<10;++i){
		point = test.getGlobalNodeCoordinates(i);
		logger.assert_always((std::abs(point[0]-1.-i/10.)<1e-12),"getGlobalNodeCoordinates");
	}

	pointIndexes.resize(1);

	pointIndexes = test.getGlobalFaceNodeIndices(0);
	logger.assert_always((pointIndexes[0]==4),"getGlobalFaceNodeIndices");
	pointIndexes = test.getGlobalFaceNodeIndices(1);
	logger.assert_always((pointIndexes[0]==7),"getGlobalFaceNodeIndices");

	pointIndexes = test.getLocalFaceNodeIndices(0);
	logger.assert_always((pointIndexes[0]==0),"getLocalFaceNodeIndices");
	pointIndexes = test.getLocalFaceNodeIndices(1);
	logger.assert_always((pointIndexes[0]==1),"getLocalFaceNodeIndices");

	logger.assert_always((test.getNrOfFaces()==2),"getNrOfFaces");

	logger.assert_always((test.getRefGeometry()==&Geometry::ReferenceLine::Instance()),"getRefGeometry");


	return 0;
}

