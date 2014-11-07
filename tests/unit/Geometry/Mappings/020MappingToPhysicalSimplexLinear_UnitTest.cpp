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

#include "Geometry/Mappings/MappingToPhysSimplexLinear.hpp"
#include "cassert"

#include "Geometry/ReferenceLine.hpp"
#include "Geometry/PhysicalLine.hpp"
#include "Geometry/ReferenceTriangle.hpp"
#include "Geometry/PhysicalTriangle.hpp"
#include "Geometry/ReferenceTetrahedron.hpp"
#include "Geometry/PhysicalTetrahedron.hpp"
#include "Geometry/PointPhysical.hpp"
#include "Geometry/PhysicalGeometry.hpp"
#include <cmath>
//transformations should map internal points to internal points, external points to external points
//and nodes to nodes so construct the physical geometries such that this can be checked :(

//bool isInternal1D(const Geometry::PointPhysical& p) {
//	return p[0]>1.4&&p[1]<1.7;
//}

bool isInternal2D(const Geometry::PointPhysical& p) {
	return (p[1]-p[0])>1.&&(1.8*p[0]-1.4*p[1])>-0.84&&(1.5*p[0]-1.1*p[1])<-0.42;
}

bool isInternal3D(const Geometry::PointPhysical& p) {
	return p[2]>0&&(p[1]-p[0])>(1+p[2]/10.)&&(1.8*p[0]-1.4*p[1])>-0.84&&(1.5*p[0]-1.1*p[1])<-0.42;
}

int main() {

	//dim 1

	std::vector<unsigned int> pointIndexes;
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

	Geometry::PhysicalLine oops(pointIndexes,nodes1D,&rGeom);
	pointIndexes[1]=7;
	Geometry::PhysicalLine pGeom(pointIndexes,nodes1D,&rGeom);

	Geometry::MappingToPhysSimplexLinear<1> mapping1D(&pGeom),reinit1D(&oops);
	reinit1D.reinit(&pGeom);

	Geometry::Jacobian jac1D(1,1);

	/*for(refPoint1D[0]=-2.8189;refPoint1D[0]<3.141;refPoint1D[0]+=0.1) {
		mapping1D.transform(refPoint1D,point1D);
		assert(("transform",rGeom.isInternalPoint(refPoint1D)==isInternal1D(point1D)));
		reinit1D.transform(refPoint1D,point1D);
		assert(("reinit",rGeom.isInternalPoint(refPoint1D)==isInternal1D(point1D)));

		refPoint1D[0]+=-1.e-8;
		mapping1D.transform(refPoint1D,compare1D);
		refPoint1D[0]+=2.e-8;
		mapping1D.transform(refPoint1D,point1D);

		refPoint1D[0]+=-1e-8;
		mapping1D.calcJacobian(refPoint1D,jac1D);
		assert(("jacobian",std::abs(jac1D[0]-5.e7*(point1D[0]-compare1D[0]))<1e-5));//estimate is a bit rough, but should work for most mappings
	}

	for(int i=0;i<rGeom.getNumberOfNodes();++i){
		rGeom.getNode(i,refPoint1D);
		pGeom.getNodeCoordinates(i,compare1D);
		mapping1D.transform(refPoint1D,point1D);
		assert(("transform",std::abs(point1D[0]-compare1D[0])<1e-12));
	}

	assert(("getTargetDimension",mapping1D.getTargetDimension()==1));

	for(int i=0;i<10;++i){
		mapping1D.getNodeCoordinates(i,compare1D);
		pGeom.getGlobalNodeCoordinates(i,point1D);
		assert(("getNodeCoordinates",compare1D==point1D));
	}*/

	//dim2

	std::vector<Geometry::PointPhysical> nodes2D;

	Geometry::PointPhysical point2D(2),compare2D(2);
	Geometry::PointReference refPoint2D(2);

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

	Geometry::ReferenceTriangle& rGeom2D = Geometry::ReferenceTriangle::Instance();

	Geometry::PhysicalTriangle oops2D(pointIndexes,nodes2D,&rGeom2D);
	pointIndexes[2]=18;
	Geometry::PhysicalTriangle pGeom2D(pointIndexes,nodes2D,&rGeom2D);

	Geometry::MappingToPhysSimplexLinear<2> mapping2D(&pGeom2D),reinit2D(&oops2D);
	reinit2D.reinit(&pGeom2D);

	Geometry::Jacobian jac2D(2,2);

	for(refPoint2D[0]=-1.5189;refPoint2D[0]<1.541;refPoint2D[0]+=0.1) {
		for(refPoint2D[1]=-1.5188;refPoint2D[1]<1.541;refPoint2D[1]+=0.1){
			mapping2D.transform(refPoint2D,point2D);
			assert(("transform"&&(rGeom2D.isInternalPoint(refPoint2D)&&isInternal2D(point2D))||(!rGeom2D.isInternalPoint(refPoint2D)&&!isInternal2D(point2D))));
			reinit2D.transform(refPoint2D,point2D);
			assert(("reinit",rGeom2D.isInternalPoint(refPoint2D)==isInternal2D(point2D)));

			refPoint2D[0]+=-1.e-8;
			mapping2D.transform(refPoint2D,compare2D);
			refPoint2D[0]+=2.e-8;
			mapping2D.transform(refPoint2D,point2D);

			refPoint2D[0]+=-1e-8;
			mapping2D.calcJacobian(refPoint2D,jac2D);
			assert(("jacobian",std::abs(jac2D[0]-5.e7*(point2D[0]-compare2D[0]))<1e-5));//estimate is a bit rough, but should work for most mappings
			assert(("jacobian",std::abs(jac2D[1]-5.e7*(point2D[1]-compare2D[1]))<1e-5));//implementations are strongly recommended to be more accurate

			refPoint2D[1]+=-1.e-8;
			mapping2D.transform(refPoint2D,compare2D);
			refPoint2D[1]+=2.e-8;
			mapping2D.transform(refPoint2D,point2D);

			refPoint2D[1]+=-1e-8;
			mapping2D.calcJacobian(refPoint2D,jac2D);
			assert(("jacobian",std::abs(jac2D[2]-5.e7*(point2D[0]-compare2D[0]))<1e-5));
			assert(("jacobian",std::abs(jac2D[3]-5.e7*(point2D[1]-compare2D[1]))<1e-5));
		}
	}

	for(int i=0;i<rGeom2D.getNumberOfNodes();++i){
		rGeom2D.getNode(i,refPoint2D);
		pGeom2D.getNodeCoordinates(i,compare2D);
		mapping2D.transform(refPoint2D,point2D);
		assert(("transform",std::abs(point2D[0]-compare2D[0])<1e-12)&&std::abs(point2D[1]-compare2D[1])<1e-12);
	}

	assert(("getTargetDimension",mapping2D.getTargetDimension()==2));

	for(int i=0;i<20;++i){
		mapping2D.getNodeCoordinates(i,compare2D);
		pGeom2D.getGlobalNodeCoordinates(i,point2D);
		assert(("getNodeCoordinates",compare2D==point2D));
	}

	//dim3

	std::vector<Geometry::PointPhysical> nodes3D;

	Geometry::PointPhysical point3D(3),compare3D(3);
	Geometry::PointReference refPoint3D(3);

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

	Geometry::ReferenceTetrahedron& rGeom3D = Geometry::ReferenceTetrahedron::Instance();

	Geometry::PhysicalTetrahedron oops3D(pointIndexes,nodes3D,&rGeom3D);
	pointIndexes[3]=38;
	Geometry::PhysicalTetrahedron pGeom3D(pointIndexes,nodes3D,&rGeom3D);

	Geometry::MappingToPhysSimplexLinear<3> mapping3D(&pGeom3D),reinit3D(&oops3D);
	reinit3D.reinit(&pGeom3D);

	Geometry::Jacobian jac3D(3,3);

	for(refPoint3D[0]=-1.5189;refPoint3D[0]<1.541;refPoint3D[0]+=0.2) {
		for(refPoint3D[1]=-1.5188;refPoint3D[1]<1.541;refPoint3D[1]+=0.2){
			for(refPoint3D[2]=-1.5188;refPoint3D[2]<1.541;refPoint3D[2]+=0.2){
				mapping3D.transform(refPoint3D,point3D);
				assert(("transform",rGeom3D.isInternalPoint(refPoint3D)==isInternal3D(point3D)));
				reinit3D.transform(refPoint3D,point3D);
				assert(("reinit",rGeom3D.isInternalPoint(refPoint3D)==isInternal3D(point3D)));

				refPoint3D[0]+=-1.e-8;
				mapping3D.transform(refPoint3D,compare3D);
				refPoint3D[0]+=2.e-8;
				mapping3D.transform(refPoint3D,point3D);

				refPoint3D[0]+=-1e-8;
				mapping3D.calcJacobian(refPoint3D,jac3D);
				assert(("jacobian",std::abs(jac3D[0]-5.e7*(point3D[0]-compare3D[0]))<1e-5));//estimate is a bit rough, but should work for most mappings
				assert(("jacobian",std::abs(jac3D[1]-5.e7*(point3D[1]-compare3D[1]))<1e-5));//implementations are strongly recommended to be more accurate
				assert(("jacobian",std::abs(jac3D[2]-5.e7*(point3D[2]-compare3D[2]))<1e-5));

				refPoint3D[1]+=-1.e-8;
				mapping3D.transform(refPoint3D,compare3D);
				refPoint3D[1]+=2.e-8;
				mapping3D.transform(refPoint3D,point3D);

				refPoint3D[1]+=-1e-8;
				mapping3D.calcJacobian(refPoint3D,jac3D);
				assert(("jacobian",std::abs(jac3D[3]-5.e7*(point3D[0]-compare3D[0]))<1e-5));
				assert(("jacobian",std::abs(jac3D[4]-5.e7*(point3D[1]-compare3D[1]))<1e-5));
				assert(("jacobian",std::abs(jac3D[5]-5.e7*(point3D[2]-compare3D[2]))<1e-5));

				refPoint3D[2]+=-1.e-8;
				mapping3D.transform(refPoint3D,compare3D);
				refPoint3D[2]+=2.e-8;
				mapping3D.transform(refPoint3D,point3D);

				refPoint3D[2]+=-1e-8;
				mapping3D.calcJacobian(refPoint3D,jac3D);
				assert(("jacobian",std::abs(jac3D[6]-5.e7*(point3D[0]-compare3D[0]))<1e-5));
				assert(("jacobian",std::abs(jac3D[7]-5.e7*(point3D[1]-compare3D[1]))<1e-5));
				assert(("jacobian",std::abs(jac3D[8]-5.e7*(point3D[2]-compare3D[2]))<1e-5));
			}
		}
	}

	for(int i=0;i<rGeom3D.getNumberOfNodes();++i){
		rGeom3D.getNode(i,refPoint3D);
		pGeom3D.getNodeCoordinates(i,compare3D);
		mapping3D.transform(refPoint3D,point3D);
		assert(("transform",std::abs(point3D[0]-compare3D[0])<1e-12)&&std::abs(point3D[1]-compare3D[1])<1e-12&&std::abs(point3D[2]-compare3D[2])<1e-12);
	}

	assert(("getTargetDimension",mapping3D.getTargetDimension()==3));

	for(int i=0;i<40;++i){
		mapping3D.getNodeCoordinates(i,compare3D);
		pGeom3D.getGlobalNodeCoordinates(i,point3D);
		assert(("getNodeCoordinates",compare3D==point3D));
	}

	return 0;
}


