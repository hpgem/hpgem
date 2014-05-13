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

#include "Geometry/Mappings/MappingToPhysPyramid.hpp"
#include "cassert"

#include "Geometry/ReferencePyramid.hpp"
#include "Geometry/PhysicalPyramid.hpp"
//transformations should map internal points to internal points, external points to external points
//and nodes to nodes so construct the physical geometries such that this can be checked :(

bool isInternal3D(const Geometry::PointPhysical& p) {
	return (p[1]-p[0])>(1.)&&(p[1]+p[0])<(7.-p[2]*0.65)&&(2.4*p[0]-0.8*p[1])>(1.44+0.12*p[2])&&(-1.5*p[0]+1.1*p[1])>(.42)&&p[2]>0;
}

int main() {

	std::vector<unsigned int> pointIndexes;

	pointIndexes.push_back(27);
	pointIndexes.push_back(4);
	pointIndexes.push_back(7);
	pointIndexes.push_back(12);

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

	Geometry::ReferencePyramid& rGeom3D = Geometry::ReferencePyramid::Instance();

	Geometry::PhysicalPyramid oops3D(pointIndexes,nodes3D,&rGeom3D);
	pointIndexes[4]=18;
	Geometry::PhysicalPyramid pGeom3D(pointIndexes,nodes3D,&rGeom3D);

	Geometry::MappingToPhysPyramid mapping3D(&pGeom3D),reinit3D(&oops3D);
	reinit3D.reinit(&pGeom3D);

	Geometry::Jacobian jac3D(3,3);

	for(refPoint3D[0]=-1.5189;refPoint3D[0]<1.541;refPoint3D[0]+=0.2) {
		for(refPoint3D[1]=-1.5188;refPoint3D[1]<1.541;refPoint3D[1]+=0.2){
			for(refPoint3D[2]=-1.5188;refPoint3D[2]<1.541;refPoint3D[2]+=0.2){
				mapping3D.transform(refPoint3D,point3D);
				if(rGeom3D.isInternalPoint(refPoint3D))//not perfect, but the degenerate cube face makes the mapping less linear than desired
					assert(("transform",isInternal3D(point3D)));
				reinit3D.transform(refPoint3D,point3D);
				if(rGeom3D.isInternalPoint(refPoint3D))
					assert(("reinit",isInternal3D(point3D)));

				refPoint3D[0]+=-1.e-8;
				mapping3D.transform(refPoint3D,compare3D);
				refPoint3D[0]+=2.e-8;
				mapping3D.transform(refPoint3D,point3D);

				refPoint3D[0]+=-1e-8;
				mapping3D.calcJacobian(refPoint3D,jac3D);
				assert(("jacobian",fabs(jac3D[0]-5.e7*(point3D[0]-compare3D[0]))<1e-5));//estimate is a bit rough, but should work for most mappings
				assert(("jacobian",fabs(jac3D[1]-5.e7*(point3D[1]-compare3D[1]))<1e-5));//implementations are strongly recommended to be more accurate
				assert(("jacobian",fabs(jac3D[2]-5.e7*(point3D[2]-compare3D[2]))<1e-5));

				refPoint3D[1]+=-1.e-8;
				mapping3D.transform(refPoint3D,compare3D);
				refPoint3D[1]+=2.e-8;
				mapping3D.transform(refPoint3D,point3D);

				refPoint3D[1]+=-1e-8;
				mapping3D.calcJacobian(refPoint3D,jac3D);
				assert(("jacobian",fabs(jac3D[3]-5.e7*(point3D[0]-compare3D[0]))<1e-5));
				assert(("jacobian",fabs(jac3D[4]-5.e7*(point3D[1]-compare3D[1]))<1e-5));
				assert(("jacobian",fabs(jac3D[5]-5.e7*(point3D[2]-compare3D[2]))<1e-5));

				refPoint3D[2]+=-1.e-8;
				mapping3D.transform(refPoint3D,compare3D);
				refPoint3D[2]+=2.e-8;
				mapping3D.transform(refPoint3D,point3D);

				refPoint3D[2]+=-1e-8;
				mapping3D.calcJacobian(refPoint3D,jac3D);
				assert(("jacobian",fabs(jac3D[6]-5.e7*(point3D[0]-compare3D[0]))<1e-5));
				assert(("jacobian",fabs(jac3D[7]-5.e7*(point3D[1]-compare3D[1]))<1e-5));
				assert(("jacobian",fabs(jac3D[8]-5.e7*(point3D[2]-compare3D[2]))<1e-5));
			}
		}
	}

	for(int i=0;i<rGeom3D.getNumberOfNodes();++i){
		rGeom3D.getNode(i,refPoint3D);
		pGeom3D.getNodeCoordinates(i,compare3D);
		mapping3D.transform(refPoint3D,point3D);
		assert(("transform",fabs(point3D[0]-compare3D[0])<1e-12)&&fabs(point3D[1]-compare3D[1])<1e-12&&fabs(point3D[2]-compare3D[2])<1e-12);
	}

	assert(("getTargetDimension",mapping3D.getTargetDimension()==3));

	for(int i=0;i<40;++i){
		mapping3D.getNodeCoordinates(i,compare3D);
		pGeom3D.getGlobalNodeCoordinates(i,point3D);
		assert(("getNodeCoordinates",compare3D==point3D));
	}

	return 0;
}



