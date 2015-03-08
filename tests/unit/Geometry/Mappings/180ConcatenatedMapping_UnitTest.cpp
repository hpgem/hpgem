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

#include "Geometry/Mappings/ConcatenatedMapping.hpp"
#include "Geometry/Mappings/MappingToRefPointToPoint.hpp"
#include "Geometry/Mappings/MappingToRefPointToLine.hpp"
#include "Geometry/Mappings/MappingToRefLineToLine.hpp"
#include "Geometry/Mappings/MappingToRefLineToSquare.hpp"
#include "Geometry/Mappings/MappingToRefSquareToSquare.hpp"
#include "Geometry/Mappings/MappingToRefSquareToCube.hpp"
#include "Geometry/Mappings/MappingToRefCubeToCube.hpp"
#include "Geometry/Mappings/MappingToRefCubeToHypercube.hpp"
#include "cassert"

#include "Geometry/ReferencePoint.hpp"
#include "Geometry/ReferenceLine.hpp"
#include "Geometry/ReferenceSquare.hpp"
#include "Geometry/ReferenceCube.hpp"
#include "Geometry/ReferenceHypercube.hpp"
#include "Geometry/PointReference.hpp"
#include "Geometry/Jacobian.hpp"
#include "LinearAlgebra/NumericalVector.hpp"

int main(){

	Geometry::PointReference point0D(0),point1D(1),point2D(2),point3D(3),point4D(4),
			compare0D(0),compare1D(1),compare2D(2),compare3D(3),compare4D(4),
			orig0D(0),orig1D(1),orig2D(2),orig3D(3),orig4D(4);

	//0->0->0

	Geometry::ReferenceGeometry* source = &Geometry::ReferencePoint::Instance();
	Geometry::ReferenceGeometry* target = &Geometry::ReferencePoint::Instance();

	Geometry::Jacobian jac(0,0);

	std::vector<int> nodesAfterTransformation(1);

	Geometry::ConcatenatedMapping* test = new Geometry::ConcatenatedMapping(Geometry::MappingToRefPointToPoint::Instance(),Geometry::MappingToRefPointToPoint::Instance());
	nodesAfterTransformation[0]=0;

	test->transform(orig0D,point0D);
	assert(("transform",source->isInternalPoint(orig0D)==target->isInternalPoint(point0D)));

	test->calcJacobian(orig0D,jac);

	for(int i=0;i<source->getNumberOfNodes();++i){
		source->getNode(i,orig0D);
		target->getNode(nodesAfterTransformation[i],compare0D);
		test->transform(orig0D,point0D);
	}

	assert(("getTargetDimension",test->getTargetDimension()==0));

	//0->0->1

	target = &Geometry::ReferenceLine::Instance();
	jac.resize(0,1);
	delete test;
	test=new Geometry::ConcatenatedMapping(Geometry::MappingToRefPointToPoint::Instance(),Geometry::MappingToRefPointToLine0::Instance());

	test->transform(orig0D,point1D);
	assert(("transform",source->isInternalPoint(orig0D)==target->isInternalPoint(point1D)));

	test->calcJacobian(orig0D,jac);

	for(int i=0;i<source->getNumberOfNodes();++i){
		source->getNode(i,orig0D);
		target->getNode(nodesAfterTransformation[i],compare1D);
		test->transform(orig0D,point1D);
	}

	assert(("getTargetDimension",test->getTargetDimension()==1));

	//0->1->1

	delete test;
	test=new Geometry::ConcatenatedMapping(Geometry::MappingToRefPointToLine1::Instance(),Geometry::MappingToRefLineToLine1::Instance());

	test->transform(orig0D,point1D);
	assert(("transform",source->isInternalPoint(orig0D)==target->isInternalPoint(point1D)));

	test->calcJacobian(orig0D,jac);

	for(int i=0;i<source->getNumberOfNodes();++i){
		source->getNode(i,orig0D);
		target->getNode(nodesAfterTransformation[i],compare1D);
		test->transform(orig0D,point1D);
	}

	assert(("getTargetDimension",test->getTargetDimension()==1));

	//0->1->2

	target = &Geometry::ReferenceSquare::Instance();
	jac.resize(0,2);
	delete test;
	test=new Geometry::ConcatenatedMapping(Geometry::MappingToRefPointToLine1::Instance(),Geometry::MappingToRefLineToSquare3::Instance());
	nodesAfterTransformation[0]=3;

	test->transform(orig0D,point2D);
	assert(("transform",source->isInternalPoint(orig0D)==target->isInternalPoint(point2D)));

	test->calcJacobian(orig0D,jac);

	for(int i=0;i<source->getNumberOfNodes();++i){
		source->getNode(i,orig0D);
		target->getNode(nodesAfterTransformation[i],compare2D);
		test->transform(orig0D,point2D);
	}

	assert(("getTargetDimension",test->getTargetDimension()==2));

	//1->1->1

	source = &Geometry::ReferenceLine::Instance();
	target = &Geometry::ReferenceLine::Instance();
	jac.resize(1,1);
	nodesAfterTransformation.resize(2);
	delete test;
	test=new Geometry::ConcatenatedMapping(Geometry::MappingToRefLineToLine1::Instance(),Geometry::MappingToRefLineToLine0::Instance());
	nodesAfterTransformation[0]=1;
	nodesAfterTransformation[1]=0;

	for(orig1D[0]=-2.8189;orig1D[0]<3.141;orig1D[0]+=0.1) {
		test->transform(orig1D,point1D);
		assert(("transform",source->isInternalPoint(orig1D)==target->isInternalPoint(point1D)));

		orig1D[0]+=-1.e-8;
		test->transform(orig1D,compare1D);
		orig1D[0]+=2.e-8;
		test->transform(orig1D,point1D);

		orig1D[0]+=-1e-8;
		test->calcJacobian(orig1D,jac);
		assert(("jacobian",fabs(jac[0]-5.e7*(point1D[0]-compare1D[0]))<1e-5));//estimate is a bit rough, but should work for most mappings
	}

	for(int i=0;i<source->getNumberOfNodes();++i){
		source->getNode(i,orig1D);
		target->getNode(nodesAfterTransformation[i],compare1D);
		test->transform(orig1D,point1D);
		assert(("transform",fabs(point1D[0]-compare1D[0])<1e-12));
	}

	assert(("getTargetDimension",test->getTargetDimension()==1));

	//1->1->2

	target = &Geometry::ReferenceSquare::Instance();
	jac.resize(1,2);
	delete test;
	test=new Geometry::ConcatenatedMapping(Geometry::MappingToRefLineToLine1::Instance(),Geometry::MappingToRefLineToSquare3::Instance());
	nodesAfterTransformation[0]=3;
	nodesAfterTransformation[1]=2;

	for(orig1D[0]=-2.8189;orig1D[0]<3.141;orig1D[0]+=0.1) {
		test->transform(orig1D,point2D);
		assert(("transform",source->isInternalPoint(orig1D)==target->isInternalPoint(point2D)));

		orig1D[0]+=-1.e-8;
		test->transform(orig1D,compare2D);
		orig1D[0]+=2.e-8;
		test->transform(orig1D,point2D);

		orig1D[0]+=-1e-8;
		test->calcJacobian(orig1D,jac);
		assert(("jacobian",fabs(jac[0]-5.e7*(point2D[0]-compare2D[0]))<1e-5));//estimate is a bit rough, but should work for most mappings
		assert(("jacobian",fabs(jac[1]-5.e7*(point2D[1]-compare2D[1]))<1e-5));//implementations are very strongly recommended to be more accurate
	}

	for(int i=0;i<source->getNumberOfNodes();++i){
		source->getNode(i,orig1D);
		target->getNode(nodesAfterTransformation[i],compare2D);
		test->transform(orig1D,point2D);
		assert(("transform",fabs(point2D[0]-compare2D[0])<1e-12));
		assert(("transform",fabs(point2D[1]-compare2D[1])<1e-12));
	}

	assert(("getTargetDimension",test->getTargetDimension()==2));

	//1->2->2

	delete test;
	test=new Geometry::ConcatenatedMapping(Geometry::MappingToRefLineToSquare3::Instance(),Geometry::MappingToRefSquareToSquare7::Instance());
	nodesAfterTransformation[0]=1;
	nodesAfterTransformation[1]=3;

	for(orig1D[0]=-2.8189;orig1D[0]<3.141;orig1D[0]+=0.1) {
		test->transform(orig1D,point2D);
		assert(("transform",source->isInternalPoint(orig1D)==target->isInternalPoint(point2D)));

		orig1D[0]+=-1.e-8;
		test->transform(orig1D,compare2D);
		orig1D[0]+=2.e-8;
		test->transform(orig1D,point2D);

		orig1D[0]+=-1e-8;
		test->calcJacobian(orig1D,jac);
		assert(("jacobian",fabs(jac[0]-5.e7*(point2D[0]-compare2D[0]))<1e-5));//estimate is a bit rough, but should work for most mappings
		assert(("jacobian",fabs(jac[1]-5.e7*(point2D[1]-compare2D[1]))<1e-5));//implementations are very strongly recommended to be more accurate
	}

	for(int i=0;i<source->getNumberOfNodes();++i){
		source->getNode(i,orig1D);
		target->getNode(nodesAfterTransformation[i],compare2D);
		test->transform(orig1D,point2D);
		assert(("transform",fabs(point2D[0]-compare2D[0])<1e-12));
		assert(("transform",fabs(point2D[1]-compare2D[1])<1e-12));
	}

	assert(("getTargetDimension",test->getTargetDimension()==2));

	//1->2->3

	target = &Geometry::ReferenceCube::Instance();
	jac.resize(1,3);
	delete test;
	test=new Geometry::ConcatenatedMapping(Geometry::MappingToRefLineToSquare1::Instance(),Geometry::MappingToRefSquareToCube2::Instance());
	nodesAfterTransformation[0]=0;
	nodesAfterTransformation[1]=4;

	for(orig1D[0]=-2.8189;orig1D[0]<3.141;orig1D[0]+=0.1) {
		test->transform(orig1D,point3D);
		assert(("transform",source->isInternalPoint(orig1D)==target->isInternalPoint(point3D)));

		orig1D[0]+=-1.e-8;
		test->transform(orig1D,compare3D);
		orig1D[0]+=2.e-8;
		test->transform(orig1D,point3D);

		orig1D[0]+=-1e-8;
		test->calcJacobian(orig1D,jac);
		assert(("jacobian",fabs(jac[0]-5.e7*(point3D[0]-compare3D[0]))<1e-5));//estimate is a bit rough, but should work for most mappings
		assert(("jacobian",fabs(jac[1]-5.e7*(point3D[1]-compare3D[1]))<1e-5));//implementations are very strongly recommended to be more accurate
	}

	for(int i=0;i<source->getNumberOfNodes();++i){
		source->getNode(i,orig1D);
		target->getNode(nodesAfterTransformation[i],compare3D);
		test->transform(orig1D,point3D);
		assert(("transform",fabs(point3D[0]-compare3D[0])<1e-12));
		assert(("transform",fabs(point3D[1]-compare3D[1])<1e-12));
	}

	assert(("getTargetDimension",test->getTargetDimension()==3));

	//2->2->2

	source = &Geometry::ReferenceSquare::Instance();
	target = &Geometry::ReferenceSquare::Instance();
	nodesAfterTransformation.resize(4);
	jac.resize(2,2);
	delete test;
	test=new Geometry::ConcatenatedMapping(Geometry::MappingToRefSquareToSquare3::Instance(),Geometry::MappingToRefSquareToSquare2::Instance());
	nodesAfterTransformation[0]=1;
	nodesAfterTransformation[1]=3;
	nodesAfterTransformation[2]=0;
	nodesAfterTransformation[3]=2;

	for(orig2D[0]=-2.8189;orig2D[0]<3.141;orig2D[0]+=0.1) {
		for(orig2D[1]=-2.8189;orig2D[1]<3.141;orig2D[1]+=0.1) {
			test->transform(orig2D,point2D);
			assert(("transform",source->isInternalPoint(orig2D)==target->isInternalPoint(point2D)));

			orig2D[0]+=-1.e-8;
			test->transform(orig2D,compare2D);
			orig2D[0]+=2.e-8;
			test->transform(orig2D,point2D);

			orig2D[0]+=-1e-8;
			test->calcJacobian(orig2D,jac);
			assert(("jacobian",fabs(jac[0]-5.e7*(point2D[0]-compare2D[0]))<1e-5));//estimate is a bit rough, but should work for most mappings
			assert(("jacobian",fabs(jac[1]-5.e7*(point2D[1]-compare2D[1]))<1e-5));//implementations are very strongly recommended to be more accurate

			orig2D[1]+=-1.e-8;
			test->transform(orig2D,compare2D);
			orig2D[1]+=2.e-8;
			test->transform(orig2D,point2D);

			orig2D[1]+=-1e-8;
			test->calcJacobian(orig2D,jac);
			assert(("jacobian",fabs(jac[2]-5.e7*(point2D[0]-compare2D[0]))<1e-5));
			assert(("jacobian",fabs(jac[3]-5.e7*(point2D[1]-compare2D[1]))<1e-5));
		}
	}

	for(int i=0;i<source->getNumberOfNodes();++i){
		source->getNode(i,orig2D);
		target->getNode(nodesAfterTransformation[i],compare2D);
		test->transform(orig2D,point2D);
		assert(("transform",fabs(point2D[0]-compare2D[0])<1e-12));
		assert(("transform",fabs(point2D[1]-compare2D[1])<1e-12));
	}

	assert(("getTargetDimension",test->getTargetDimension()==2));

	//2->2->3

	target = &Geometry::ReferenceCube::Instance();
	jac.resize(2,3);
	delete test;
	test=new Geometry::ConcatenatedMapping(Geometry::MappingToRefSquareToSquare3::Instance(),Geometry::MappingToRefSquareToCube2::Instance());
	nodesAfterTransformation[0]=4;
	nodesAfterTransformation[1]=0;
	nodesAfterTransformation[2]=6;
	nodesAfterTransformation[3]=2;

	for(orig2D[0]=-2.8189;orig2D[0]<3.141;orig2D[0]+=0.1) {
		for(orig2D[1]=-2.8189;orig2D[1]<3.141;orig2D[1]+=0.1) {
			test->transform(orig2D,point3D);
			assert(("transform",source->isInternalPoint(orig2D)==target->isInternalPoint(point3D)));

			orig2D[0]+=-1.e-8;
			test->transform(orig2D,compare3D);
			orig2D[0]+=2.e-8;
			test->transform(orig2D,point3D);

			orig2D[0]+=-1e-8;
			test->calcJacobian(orig2D,jac);
			assert(("jacobian",fabs(jac[0]-5.e7*(point3D[0]-compare3D[0]))<1e-5));//estimate is a bit rough, but should work for most mappings
			assert(("jacobian",fabs(jac[1]-5.e7*(point3D[1]-compare3D[1]))<1e-5));//implementations are very strongly recommended to be more accurate
			assert(("jacobian",fabs(jac[2]-5.e7*(point3D[2]-compare3D[2]))<1e-5));

			orig2D[1]+=-1.e-8;
			test->transform(orig2D,compare3D);
			orig2D[1]+=2.e-8;
			test->transform(orig2D,point3D);

			orig2D[1]+=-1e-8;
			test->calcJacobian(orig2D,jac);
			assert(("jacobian",fabs(jac[3]-5.e7*(point3D[0]-compare3D[0]))<1e-5));
			assert(("jacobian",fabs(jac[4]-5.e7*(point3D[1]-compare3D[1]))<1e-5));
			assert(("jacobian",fabs(jac[5]-5.e7*(point3D[2]-compare3D[2]))<1e-5));
		}
	}

	for(int i=0;i<source->getNumberOfNodes();++i){
		source->getNode(i,orig2D);
		target->getNode(nodesAfterTransformation[i],compare3D);
		test->transform(orig2D,point3D);
		assert(("transform",fabs(point3D[0]-compare3D[0])<1e-12));
		assert(("transform",fabs(point3D[1]-compare3D[1])<1e-12));
		assert(("transform",fabs(point3D[2]-compare3D[2])<1e-12));
	}

	assert(("getTargetDimension",test->getTargetDimension()==3));

	//2->3->3

	delete test;
	test=new Geometry::ConcatenatedMapping(Geometry::MappingToRefSquareToCube3::Instance(),Geometry::MappingToRefCubeToCube2::Instance());
	nodesAfterTransformation[0]=4;
	nodesAfterTransformation[1]=0;
	nodesAfterTransformation[2]=6;
	nodesAfterTransformation[3]=2;

	for(orig2D[0]=-2.8189;orig2D[0]<3.141;orig2D[0]+=0.1) {
		for(orig2D[1]=-2.8189;orig2D[1]<3.141;orig2D[1]+=0.1) {
			test->transform(orig2D,point3D);
			assert(("transform",source->isInternalPoint(orig2D)==target->isInternalPoint(point3D)));

			orig2D[0]+=-1.e-8;
			test->transform(orig2D,compare3D);
			orig2D[0]+=2.e-8;
			test->transform(orig2D,point3D);

			orig2D[0]+=-1e-8;
			test->calcJacobian(orig2D,jac);
			assert(("jacobian",fabs(jac[0]-5.e7*(point3D[0]-compare3D[0]))<1e-5));//estimate is a bit rough, but should work for most mappings
			assert(("jacobian",fabs(jac[1]-5.e7*(point3D[1]-compare3D[1]))<1e-5));//implementations are very strongly recommended to be more accurate
			assert(("jacobian",fabs(jac[2]-5.e7*(point3D[2]-compare3D[2]))<1e-5));

			orig2D[1]+=-1.e-8;
			test->transform(orig2D,compare3D);
			orig2D[1]+=2.e-8;
			test->transform(orig2D,point3D);

			orig2D[1]+=-1e-8;
			test->calcJacobian(orig2D,jac);
			assert(("jacobian",fabs(jac[3]-5.e7*(point3D[0]-compare3D[0]))<1e-5));
			assert(("jacobian",fabs(jac[4]-5.e7*(point3D[1]-compare3D[1]))<1e-5));
			assert(("jacobian",fabs(jac[5]-5.e7*(point3D[2]-compare3D[2]))<1e-5));
		}
	}

	/*for(int i=0;i<source->getNumberOfNodes();++i){///\TODO figure out how the cube->cube maps actually map their nodes
		source->getNode(i,orig2D);
		target->getNode(nodesAfterTransformation[i],compare3D);
		test->transform(orig2D,point3D);
		assert(("transform",fabs(point3D[0]-compare3D[0])<1e-12));
		assert(("transform",fabs(point3D[1]-compare3D[1])<1e-12));
		assert(("transform",fabs(point3D[2]-compare3D[2])<1e-12));
	}*/

	assert(("getTargetDimension",test->getTargetDimension()==3));

	//2->3->4
	jac.resize(2,4);

	target = &Geometry::ReferenceHypercube::Instance();
	delete test;
	test=new Geometry::ConcatenatedMapping(Geometry::MappingToRefSquareToCube3::Instance(),Geometry::MappingToRefCubeToHypercube4::Instance());
	nodesAfterTransformation[0]=3;
	nodesAfterTransformation[1]=7;
	nodesAfterTransformation[2]=11;
	nodesAfterTransformation[3]=15;

	for(orig2D[0]=-2.8189;orig2D[0]<3.141;orig2D[0]+=0.1) {
		for(orig2D[1]=-2.8189;orig2D[1]<3.141;orig2D[1]+=0.1) {
			test->transform(orig2D,point4D);
			assert(("transform",source->isInternalPoint(orig2D)==target->isInternalPoint(point4D)));

			orig2D[0]+=-1.e-8;
			test->transform(orig2D,compare4D);
			orig2D[0]+=2.e-8;
			test->transform(orig2D,point4D);

			orig2D[0]+=-1e-8;
			test->calcJacobian(orig2D,jac);
			assert(("jacobian",fabs(jac[0]-5.e7*(point4D[0]-compare4D[0]))<1e-5));//estimate is a bit rough, but should work for most mappings
			assert(("jacobian",fabs(jac[1]-5.e7*(point4D[1]-compare4D[1]))<1e-5));//implementations are very strongly recommended to be more accurate
			assert(("jacobian",fabs(jac[2]-5.e7*(point4D[2]-compare4D[2]))<1e-5));
			assert(("jacobian",fabs(jac[3]-5.e7*(point4D[3]-compare4D[3]))<1e-5));

			orig2D[1]+=-1.e-8;
			test->transform(orig2D,compare4D);
			orig2D[1]+=2.e-8;
			test->transform(orig2D,point4D);

			orig2D[1]+=-1e-8;
			test->calcJacobian(orig2D,jac);
			assert(("jacobian",fabs(jac[4]-5.e7*(point4D[0]-compare4D[0]))<1e-5));
			assert(("jacobian",fabs(jac[5]-5.e7*(point4D[1]-compare4D[1]))<1e-5));
			assert(("jacobian",fabs(jac[6]-5.e7*(point4D[2]-compare4D[2]))<1e-5));
			assert(("jacobian",fabs(jac[7]-5.e7*(point4D[3]-compare4D[3]))<1e-5));
		}
	}

	for(int i=0;i<source->getNumberOfNodes();++i){
		source->getNode(i,orig2D);
		target->getNode(nodesAfterTransformation[i],compare4D);
		test->transform(orig2D,point4D);
		assert(("transform",fabs(point4D[0]-compare4D[0])<1e-12));
		assert(("transform",fabs(point4D[1]-compare4D[1])<1e-12));
		assert(("transform",fabs(point4D[2]-compare4D[2])<1e-12));
		assert(("transform",fabs(point4D[3]-compare4D[3])<1e-12));
	}

	assert(("getTargetDimension",test->getTargetDimension()==4));

	//3->3->3

	source = &Geometry::ReferenceCube::Instance();
	target = &Geometry::ReferenceCube::Instance();
	nodesAfterTransformation.resize(8);
	jac.resize(3,3);
	delete test;
	test=new Geometry::ConcatenatedMapping(Geometry::MappingToRefCubeToCube3::Instance(),Geometry::MappingToRefCubeToCube7::Instance());
	nodesAfterTransformation[0]=4;
	nodesAfterTransformation[1]=0;
	nodesAfterTransformation[2]=6;
	nodesAfterTransformation[3]=2;

	for(orig3D[0]=-2.8189;orig3D[0]<3.141;orig3D[0]+=0.2) {
		for(orig3D[1]=-2.8189;orig3D[1]<3.141;orig3D[1]+=0.2) {
			for(orig3D[2]=-2.8189;orig3D[2]<3.141;orig3D[2]+=0.2) {
				test->transform(orig3D,point3D);
				assert(("transform",source->isInternalPoint(orig3D)==target->isInternalPoint(point3D)));

				orig3D[0]+=-1.e-8;
				test->transform(orig3D,compare3D);
				orig3D[0]+=2.e-8;
				test->transform(orig3D,point3D);

				orig3D[0]+=-1e-8;
				test->calcJacobian(orig3D,jac);
				assert(("jacobian",fabs(jac[0]-5.e7*(point3D[0]-compare3D[0]))<1e-5));//estimate is a bit rough, but should work for most mappings
				assert(("jacobian",fabs(jac[1]-5.e7*(point3D[1]-compare3D[1]))<1e-5));//implementations are very strongly recommended to be more accurate
				assert(("jacobian",fabs(jac[2]-5.e7*(point3D[2]-compare3D[2]))<1e-5));

				orig3D[1]+=-1.e-8;
				test->transform(orig3D,compare3D);
				orig3D[1]+=2.e-8;
				test->transform(orig3D,point3D);

				orig3D[1]+=-1e-8;
				test->calcJacobian(orig3D,jac);
				assert(("jacobian",fabs(jac[3]-5.e7*(point3D[0]-compare3D[0]))<1e-5));
				assert(("jacobian",fabs(jac[4]-5.e7*(point3D[1]-compare3D[1]))<1e-5));
				assert(("jacobian",fabs(jac[5]-5.e7*(point3D[2]-compare3D[2]))<1e-5));

				orig3D[2]+=-1.e-8;
				test->transform(orig3D,compare3D);
				orig3D[2]+=2.e-8;
				test->transform(orig3D,point3D);

				orig3D[2]+=-1e-8;
				test->calcJacobian(orig3D,jac);
				assert(("jacobian",fabs(jac[6]-5.e7*(point3D[0]-compare3D[0]))<1e-5));
				assert(("jacobian",fabs(jac[7]-5.e7*(point3D[1]-compare3D[1]))<1e-5));
				assert(("jacobian",fabs(jac[8]-5.e7*(point3D[2]-compare3D[2]))<1e-5));
			}
		}
	}

	/*for(int i=0;i<source->getNumberOfNodes();++i){///\TODO figure out how the cube->cube maps actually map their nodes
		source->getNode(i,orig3D);
		target->getNode(nodesAfterTransformation[i],compare3D);
		test->transform(orig3D,point3D);
		assert(("transform",fabs(point3D[0]-compare3D[0])<1e-12));
		assert(("transform",fabs(point3D[1]-compare3D[1])<1e-12));
		assert(("transform",fabs(point3D[2]-compare3D[2])<1e-12));
	}*/

	assert(("getTargetDimension",test->getTargetDimension()==3));

	//3->3->4

	target = &Geometry::ReferenceHypercube::Instance();
	jac.resize(3,4);
	delete test;
	test=new Geometry::ConcatenatedMapping(Geometry::MappingToRefCubeToCube3::Instance(),Geometry::MappingToRefCubeToHypercube7::Instance());
	nodesAfterTransformation[0]=4;
	nodesAfterTransformation[1]=0;
	nodesAfterTransformation[2]=6;
	nodesAfterTransformation[3]=2;

	for(orig3D[0]=-2.8189;orig3D[0]<3.141;orig3D[0]+=0.2) {
		for(orig3D[1]=-2.8189;orig3D[1]<3.141;orig3D[1]+=0.2) {
			for(orig3D[2]=-2.8189;orig3D[2]<3.141;orig3D[2]+=0.2) {
				test->transform(orig3D,point4D);
				assert(("transform",source->isInternalPoint(orig3D)==target->isInternalPoint(point4D)));

				orig3D[0]+=-1.e-8;
				test->transform(orig3D,compare4D);
				orig3D[0]+=2.e-8;
				test->transform(orig3D,point4D);

				orig3D[0]+=-1e-8;
				test->calcJacobian(orig3D,jac);
				assert(("jacobian",fabs(jac[0]-5.e7*(point4D[0]-compare4D[0]))<1e-5));//estimate is a bit rough, but should work for most mappings
				assert(("jacobian",fabs(jac[1]-5.e7*(point4D[1]-compare4D[1]))<1e-5));//implementations are very strongly recommended to be more accurate
				assert(("jacobian",fabs(jac[2]-5.e7*(point4D[2]-compare4D[2]))<1e-5));
				assert(("jacobian",fabs(jac[3]-5.e7*(point4D[3]-compare4D[3]))<1e-5));

				orig3D[1]+=-1.e-8;
				test->transform(orig3D,compare4D);
				orig3D[1]+=2.e-8;
				test->transform(orig3D,point4D);

				orig3D[1]+=-1e-8;
				test->calcJacobian(orig3D,jac);
				assert(("jacobian",fabs(jac[4]-5.e7*(point4D[0]-compare4D[0]))<1e-5));
				assert(("jacobian",fabs(jac[5]-5.e7*(point4D[1]-compare4D[1]))<1e-5));
				assert(("jacobian",fabs(jac[6]-5.e7*(point4D[2]-compare4D[2]))<1e-5));
				assert(("jacobian",fabs(jac[7]-5.e7*(point4D[3]-compare4D[3]))<1e-5));

				orig3D[2]+=-1.e-8;
				test->transform(orig3D,compare4D);
				orig3D[2]+=2.e-8;
				test->transform(orig3D,point4D);

				orig3D[2]+=-1e-8;
				test->calcJacobian(orig3D,jac);
				assert(("jacobian",fabs(jac[8]-5.e7*(point4D[0]-compare4D[0]))<1e-5));
				assert(("jacobian",fabs(jac[9]-5.e7*(point4D[1]-compare4D[1]))<1e-5));
				assert(("jacobian",fabs(jac[10]-5.e7*(point4D[2]-compare4D[2]))<1e-5));
				assert(("jacobian",fabs(jac[11]-5.e7*(point4D[3]-compare4D[3]))<1e-5));
			}
		}
	}

	/*for(int i=0;i<source->getNumberOfNodes();++i){///\TODO figure out how the cube->cube maps actually map their nodes
		source->getNode(i,orig3D);
		target->getNode(nodesAfterTransformation[i],compare4D);
		test->transform(orig3D,point4D);
		assert(("transform",fabs(point4D[0]-compare4D[0])<1e-12));
		assert(("transform",fabs(point4D[1]-compare4D[1])<1e-12));
		assert(("transform",fabs(point4D[2]-compare4D[2])<1e-12));
		assert(("transform",fabs(point4D[3]-compare4D[3])<1e-12));
	}*/

	assert(("getTargetDimension",test->getTargetDimension()==4));

	//0(->1->)2(->3->)4 (chaining)

	source = &Geometry::ReferencePoint::Instance();
	jac.resize(0,4);
	nodesAfterTransformation.resize(1);
	delete test;

	Geometry::ConcatenatedMapping map1(Geometry::MappingToRefPointToLine1::Instance(),Geometry::MappingToRefLineToSquare1::Instance());
	Geometry::ConcatenatedMapping map2(Geometry::MappingToRefSquareToCube1::Instance(),Geometry::MappingToRefCubeToHypercube1::Instance());

	
	test=new Geometry::ConcatenatedMapping(map1,map2);
	nodesAfterTransformation[0]=8;

	test->transform(orig0D,point4D);
	assert(("transform",source->isInternalPoint(orig0D)==target->isInternalPoint(point4D)));

	test->calcJacobian(orig3D,jac);


	for(int i=0;i<source->getNumberOfNodes();++i){///\TODO figure out how the cube->cube maps actually map their nodes
		source->getNode(i,orig0D);
		target->getNode(nodesAfterTransformation[i],compare4D);
		test->transform(orig0D,point4D);
		assert(("transform",fabs(point4D[0]-compare4D[0])<1e-12));
		assert(("transform",fabs(point4D[1]-compare4D[1])<1e-12));
		assert(("transform",fabs(point4D[2]-compare4D[2])<1e-12));
		assert(("transform",fabs(point4D[3]-compare4D[3])<1e-12));
	}

	assert(("getTargetDimension",test->getTargetDimension()==4));
	return 0;
}

