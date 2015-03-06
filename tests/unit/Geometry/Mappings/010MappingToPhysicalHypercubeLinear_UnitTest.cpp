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

#include "Geometry/Mappings/MappingToPhysHypercubeLinear.hpp"
#include "Logger.h"

#include "Geometry/ReferenceLine.hpp"
#include "Geometry/PhysicalLine.hpp"
#include "Geometry/ReferenceSquare.hpp"
#include "Geometry/PhysicalQuadrilateral.hpp"
#include "Geometry/ReferenceCube.hpp"
#include "Geometry/PhysicalHexahedron.hpp"
#include "Geometry/PointReference.hpp"
#include "Geometry/Jacobian.hpp"
#include <cmath>
//transformations should map internal points to internal points, external points to external points
//and nodes to nodes so construct the physical geometries such that this can be checked :(

bool isInternal1D(const Geometry::PointPhysical& p) {
	return p[0]>1.4&&p[0]<1.7;
}

bool isInternal2D(const Geometry::PointPhysical& p) {
	return (p[1]-p[0])>1.&&(p[1]+p[0])<7.&&(2.4*p[0]-0.8*p[1])>1.44&&(-1.5*p[0]+1.1*p[1])>.42;
}

bool isInternal3D(const Geometry::PointPhysical& p) {
	return (p[1]-p[0])>1.&&(p[1]+p[0])<7.&&(2.4*p[0]-0.8*p[1])>1.44&&(-1.5*p[0]+1.1*p[1])>.42&&p[2]>0.&&p[2]<4.;
}

int main() {

	//dim 1

	std::vector<std::size_t> pointIndexes;
	std::vector<Geometry::PointPhysical> nodes1D;

	Geometry::PointPhysical point1D(1),compare1D(1);
	Geometry::PointReference refPoint1D(1);

	pointIndexes.push_back(4);
	pointIndexes.push_back(6);

	for(double i=0.;i<10;++i){
		point1D[0]=1.+i/10.;
		nodes1D.push_back(point1D);
	}

	Geometry::ReferenceLine& rGeom = Geometry::ReferenceLine::Instance();

    Geometry::PhysicalLine oops(pointIndexes, nodes1D);
	pointIndexes[1]=7;
	Geometry::PhysicalLine pGeom(pointIndexes,nodes1D);

	Geometry::MappingToPhysHypercubeLinear<1> mapping1D(&pGeom),reinit1D(&oops);
	reinit1D.reinit(&pGeom);

	Geometry::Jacobian jac1D(1,1);

	for(refPoint1D[0]=-2.8189;refPoint1D[0]<3.141;refPoint1D[0]+=0.1) {
		point1D = mapping1D.transform(refPoint1D);
		logger.assert_always((rGeom.isInternalPoint(refPoint1D)==isInternal1D(point1D)),"transform");
		point1D = reinit1D.transform(refPoint1D);
		logger.assert_always((rGeom.isInternalPoint(refPoint1D)==isInternal1D(point1D)),"reinit");

		refPoint1D[0]+=-1.e-8;
		compare1D = mapping1D.transform(refPoint1D);
		refPoint1D[0]+=2.e-8;
		point1D = mapping1D.transform(refPoint1D);

		refPoint1D[0]+=-1e-8;
		jac1D = mapping1D.calcJacobian(refPoint1D);
		logger.assert_always((std::abs(jac1D[0]-5.e7*(point1D[0]-compare1D[0]))<1e-5),"jacobian");//estimate is a bit rough, but should work for most mappings
	}

	for(int i=0;i<rGeom.getNumberOfNodes();++i){
		refPoint1D = rGeom.getNode(i);
		compare1D = pGeom.getNodeCoordinates(i);
		point1D = mapping1D.transform(refPoint1D);
		logger.assert_always((std::abs(point1D[0]-compare1D[0])<1e-12),"transform");
	}

	logger.assert_always((mapping1D.getTargetDimension()==1),"getTargetDimension");

	for(int i=0;i<10;++i){
		compare1D = mapping1D.getNodeCoordinates(i);
		point1D = pGeom.getGlobalNodeCoordinates(i);
		logger.assert_always((compare1D==point1D),"getNodeCoordinates");
	}

	//dim2

	std::vector<Geometry::PointPhysical> nodes2D;

	Geometry::PointPhysical point2D(2),compare2D(2);
	Geometry::PointReference refPoint2D(2);

	pointIndexes.push_back(12);
	pointIndexes.push_back(13);

	for(double i=0.;i<10;++i){
		point2D[0]=1.+i/10.;
		point2D[1]=2.+i/10.;
		nodes2D.push_back(point2D);
	}
	for(double i=0.;i<10;++i){
		point2D[0]=2.+i/10.;
		point2D[1]=5.-i/10.;
		nodes2D.push_back(point2D);
	}

	Geometry::ReferenceSquare& rGeom2D = Geometry::ReferenceSquare::Instance();

	Geometry::PhysicalQuadrilateral oops2D(pointIndexes,nodes2D);
	pointIndexes[3]=18;
	Geometry::PhysicalQuadrilateral pGeom2D(pointIndexes,nodes2D);

	Geometry::MappingToPhysHypercubeLinear<2> mapping2D(&pGeom2D),reinit2D(&oops2D);
	reinit2D.reinit(&pGeom2D);

	Geometry::Jacobian jac2D(2,2);

	for(refPoint2D[0]=-1.5189;refPoint2D[0]<1.541;refPoint2D[0]+=0.1) {
		for(refPoint2D[1]=-1.5188;refPoint2D[1]<1.541;refPoint2D[1]+=0.1){
			point2D = mapping2D.transform(refPoint2D);
			logger.assert_always(((rGeom2D.isInternalPoint(refPoint2D)&&isInternal2D(point2D))||(!rGeom2D.isInternalPoint(refPoint2D)&&!isInternal2D(point2D))),"transform");
			point2D = reinit2D.transform(refPoint2D);
			logger.assert_always((rGeom2D.isInternalPoint(refPoint2D)==isInternal2D(point2D)),"reinit");

			refPoint2D[0]+=-1.e-8;
			compare2D = mapping2D.transform(refPoint2D);
			refPoint2D[0]+=2.e-8;
			point2D = mapping2D.transform(refPoint2D);

			refPoint2D[0]+=-1e-8;
			jac2D = mapping2D.calcJacobian(refPoint2D);
			logger.assert_always((std::abs(jac2D[0]-5.e7*(point2D[0]-compare2D[0]))<1e-5),"jacobian");//estimate is a bit rough, but should work for most mappings
			logger.assert_always((std::abs(jac2D[1]-5.e7*(point2D[1]-compare2D[1]))<1e-5),"jacobian");//implementations are strongly recommended to be more accurate

			refPoint2D[1]+=-1.e-8;
			compare2D = mapping2D.transform(refPoint2D);
			refPoint2D[1]+=2.e-8;
			point2D = mapping2D.transform(refPoint2D);

			refPoint2D[1]+=-1e-8;
			jac2D = mapping2D.calcJacobian(refPoint2D);
			logger.assert_always((std::abs(jac2D[2]-5.e7*(point2D[0]-compare2D[0]))<1e-5),"jacobian");
			logger.assert_always((std::abs(jac2D[3]-5.e7*(point2D[1]-compare2D[1]))<1e-5),"jacobian");
		}
	}

	for(int i=0;i<rGeom2D.getNumberOfNodes();++i){
		refPoint2D = rGeom2D.getNode(i);
		compare2D = pGeom2D.getNodeCoordinates(i);
		point2D = mapping2D.transform(refPoint2D);
		logger.assert_always((std::abs(point2D[0]-compare2D[0])<1e-12)&&std::abs(point2D[1]-compare2D[1])<1e-12,"transform");
	}

	logger.assert_always((mapping2D.getTargetDimension()==2),"getTargetDimension");

	for(std::size_t i=0;i<20;++i){
		compare2D = mapping2D.getNodeCoordinates(i);
		point2D = pGeom2D.getGlobalNodeCoordinates(i);
		logger.assert_always((compare2D==point2D),"getNodeCoordinates");
	}

	//dim3

	std::vector<Geometry::PointPhysical> nodes3D;

	Geometry::PointPhysical point3D(3),compare3D(3);
	Geometry::PointReference refPoint3D(3);

	pointIndexes.push_back(24);
	pointIndexes.push_back(27);
	pointIndexes.push_back(32);
	pointIndexes.push_back(13);

	for(double i=0.;i<10;++i){
		point3D[0]=1.+i/10.;
		point3D[1]=2.+i/10.;
		point3D[2]=0.;
		nodes3D.push_back(point3D);
	}
	for(double i=0.;i<10;++i){
		point3D[0]=2.+i/10.;
		point3D[1]=5.-i/10.;
		point3D[2]=0.;
		nodes3D.push_back(point3D);
	}
	for(double i=0.;i<10;++i){
		point3D[0]=1.+i/10.;
		point3D[1]=2.+i/10.;
		point3D[2]=4.;
		nodes3D.push_back(point3D);
	}
	for(double i=0.;i<10;++i){
		point3D[0]=2.+i/10.;
		point3D[1]=5.-i/10.;
		point3D[2]=4.;
		nodes3D.push_back(point3D);
	}

	Geometry::ReferenceCube& rGeom3D = Geometry::ReferenceCube::Instance();

	Geometry::PhysicalHexahedron oops3D(pointIndexes,nodes3D);
	pointIndexes[7]=38;
	Geometry::PhysicalHexahedron pGeom3D(pointIndexes,nodes3D);

	Geometry::MappingToPhysHypercubeLinear<3> mapping3D(&pGeom3D),reinit3D(&oops3D);
	reinit3D.reinit(&pGeom3D);

	Geometry::Jacobian jac3D(3,3);

	for(refPoint3D[0]=-1.5189;refPoint3D[0]<1.541;refPoint3D[0]+=0.2) {
		for(refPoint3D[1]=-1.5188;refPoint3D[1]<1.541;refPoint3D[1]+=0.2){
			for(refPoint3D[2]=-1.5188;refPoint3D[2]<1.541;refPoint3D[2]+=0.2){
				point3D = mapping3D.transform(refPoint3D);
				logger.assert_always((rGeom3D.isInternalPoint(refPoint3D)==isInternal3D(point3D)),"transform");
				point3D = reinit3D.transform(refPoint3D);
				logger.assert_always((rGeom3D.isInternalPoint(refPoint3D)==isInternal3D(point3D)),"reinit");

				refPoint3D[0]+=-1.e-8;
				compare3D = mapping3D.transform(refPoint3D);
				refPoint3D[0]+=2.e-8;
				point3D = mapping3D.transform(refPoint3D);

				refPoint3D[0]+=-1e-8;
				jac3D = mapping3D.calcJacobian(refPoint3D);
				logger.assert_always((std::abs(jac3D[0]-5.e7*(point3D[0]-compare3D[0]))<1e-5),"jacobian");//estimate is a bit rough, but should work for most mappings
				logger.assert_always((std::abs(jac3D[1]-5.e7*(point3D[1]-compare3D[1]))<1e-5),"jacobian");//implementations are strongly recommended to be more accurate
				logger.assert_always((std::abs(jac3D[2]-5.e7*(point3D[2]-compare3D[2]))<1e-5),"jacobian");

				refPoint3D[1]+=-1.e-8;
				compare3D = mapping3D.transform(refPoint3D);
				refPoint3D[1]+=2.e-8;
				point3D = mapping3D.transform(refPoint3D);

				refPoint3D[1]+=-1e-8;
				jac3D = mapping3D.calcJacobian(refPoint3D);
				logger.assert_always((std::abs(jac3D[3]-5.e7*(point3D[0]-compare3D[0]))<1e-5),"jacobian");
				logger.assert_always((std::abs(jac3D[4]-5.e7*(point3D[1]-compare3D[1]))<1e-5),"jacobian");
				logger.assert_always((std::abs(jac3D[5]-5.e7*(point3D[2]-compare3D[2]))<1e-5),"jacobian");

				refPoint3D[2]+=-1.e-8;
				compare3D = mapping3D.transform(refPoint3D);
				refPoint3D[2]+=2.e-8;
				point3D = mapping3D.transform(refPoint3D);

				refPoint3D[2]+=-1e-8;
				jac3D = mapping3D.calcJacobian(refPoint3D);
				logger.assert_always((std::abs(jac3D[6]-5.e7*(point3D[0]-compare3D[0]))<1e-5),"jacobian");
				logger.assert_always((std::abs(jac3D[7]-5.e7*(point3D[1]-compare3D[1]))<1e-5),"jacobian");
				logger.assert_always((std::abs(jac3D[8]-5.e7*(point3D[2]-compare3D[2]))<1e-5),"jacobian");
			}
		}
	}

	for(int i=0;i<rGeom3D.getNumberOfNodes();++i){
		refPoint3D = rGeom3D.getNode(i);
		compare3D = pGeom3D.getNodeCoordinates(i);
		point3D = mapping3D.transform(refPoint3D);
		logger.assert_always((std::abs(point3D[0]-compare3D[0])<1e-12)&&std::abs(point3D[1]-compare3D[1])<1e-12&&std::abs(point3D[2]-compare3D[2])<1e-12,"transform");
	}

	logger.assert_always((mapping3D.getTargetDimension()==3),"getTargetDimension");

	for(int i=0;i<40;++i){
		compare3D = mapping3D.getNodeCoordinates(i);
		point3D = pGeom3D.getGlobalNodeCoordinates(i);
		logger.assert_always((compare3D==point3D),"getNodeCoordinates");
	}

	return 0;
}
