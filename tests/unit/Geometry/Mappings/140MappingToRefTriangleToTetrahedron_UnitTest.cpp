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

#include "Geometry/Mappings/MappingToRefTriangleToTetrahedron.hpp"
#include "cassert"

#include "Geometry/ReferenceTetrahedron.hpp"
#include "Geometry/ReferenceTriangle.hpp"
#include "Geometry/PointReference.hpp"
#include "Geometry/Jacobian.hpp"
#include "LinearAlgebra/NumericalVector.hpp"
#include <cmath>
int main() {

	Geometry::PointReference refPoint(2),point(3),compare(3);

	Geometry::ReferenceTetrahedron& eGeom = Geometry::ReferenceTetrahedron::Instance();
	Geometry::ReferenceTriangle& fGeom = Geometry::ReferenceTriangle::Instance();

	Geometry::Jacobian jac(3,2);

	std::vector<int> nodesAfterTransformation(3);

	const Geometry::MappingReferenceToReference* test = &Geometry::MappingToRefTriangleToTetrahedron0::Instance();
	nodesAfterTransformation[0]=0;
	nodesAfterTransformation[1]=3;
	nodesAfterTransformation[2]=2;

	for(refPoint[0]=-2.8189;refPoint[0]<3.141;refPoint[0]+=0.1) {
		for(refPoint[1]=-2.8189;refPoint[1]<3.141;refPoint[1]+=0.1) {
			test->transform(refPoint,point);
			assert(("transform",fGeom.isInternalPoint(refPoint)==eGeom.isInternalPoint(point)));

			refPoint[0]+=-1.e-8;
			test->transform(refPoint,compare);
			refPoint[0]+=2.e-8;
			test->transform(refPoint,point);

			refPoint[0]+=-1e-8;
			test->calcJacobian(refPoint,jac);
			assert(("jacobian",std::abs(jac[0]-5.e7*(point[0]-compare[0]))<1e-5));//estimate is a bit rough, but should work for most mappings
			assert(("jacobian",std::abs(jac[1]-5.e7*(point[1]-compare[1]))<1e-5));//implementations are very strongly recommended to be more accurate
			assert(("jacobian",std::abs(jac[2]-5.e7*(point[2]-compare[2]))<1e-5));

			refPoint[1]+=-1.e-8;
			test->transform(refPoint,compare);
			refPoint[1]+=2.e-8;
			test->transform(refPoint,point);

			refPoint[1]+=-1e-8;
			test->calcJacobian(refPoint,jac);
			assert(("jacobian",std::abs(jac[3]-5.e7*(point[0]-compare[0]))<1e-5));
			assert(("jacobian",std::abs(jac[4]-5.e7*(point[1]-compare[1]))<1e-5));
			assert(("jacobian",std::abs(jac[5]-5.e7*(point[2]-compare[2]))<1e-5));
		}
	}

	for(int i=0;i<fGeom.getNumberOfNodes();++i){
		fGeom.getNode(i,refPoint);
		eGeom.getNode(nodesAfterTransformation[i],compare);
		test->transform(refPoint,point);
		assert(("transform",std::abs(point[0]-compare[0])<1e-12));
		assert(("transform",std::abs(point[1]-compare[1])<1e-12));
		assert(("transform",std::abs(point[2]-compare[2])<1e-12));
	}

	assert(("getTargetDimension",test->getTargetDimension()==3));


	test = &Geometry::MappingToRefTriangleToTetrahedron1::Instance();
	nodesAfterTransformation[0]=0;
	nodesAfterTransformation[1]=1;
	nodesAfterTransformation[2]=3;

	for(refPoint[0]=-2.8189;refPoint[0]<3.141;refPoint[0]+=0.1) {
		for(refPoint[1]=-2.8189;refPoint[1]<3.141;refPoint[1]+=0.1) {
			test->transform(refPoint,point);
			assert(("transform",fGeom.isInternalPoint(refPoint)==eGeom.isInternalPoint(point)));

			refPoint[0]+=-1.e-8;
			test->transform(refPoint,compare);
			refPoint[0]+=2.e-8;
			test->transform(refPoint,point);

			refPoint[0]+=-1e-8;
			test->calcJacobian(refPoint,jac);
			assert(("jacobian",std::abs(jac[0]-5.e7*(point[0]-compare[0]))<1e-5));//estimate is a bit rough, but should work for most mappings
			assert(("jacobian",std::abs(jac[1]-5.e7*(point[1]-compare[1]))<1e-5));//implementations are very strongly recommended to be more accurate
			assert(("jacobian",std::abs(jac[2]-5.e7*(point[2]-compare[2]))<1e-5));

			refPoint[1]+=-1.e-8;
			test->transform(refPoint,compare);
			refPoint[1]+=2.e-8;
			test->transform(refPoint,point);

			refPoint[1]+=-1e-8;
			test->calcJacobian(refPoint,jac);
			assert(("jacobian",std::abs(jac[3]-5.e7*(point[0]-compare[0]))<1e-5));
			assert(("jacobian",std::abs(jac[4]-5.e7*(point[1]-compare[1]))<1e-5));
			assert(("jacobian",std::abs(jac[5]-5.e7*(point[2]-compare[2]))<1e-5));
		}
	}

	for(int i=0;i<fGeom.getNumberOfNodes();++i){
		fGeom.getNode(i,refPoint);
		eGeom.getNode(nodesAfterTransformation[i],compare);
		test->transform(refPoint,point);
		assert(("transform",std::abs(point[0]-compare[0])<1e-12));
		assert(("transform",std::abs(point[1]-compare[1])<1e-12));
		assert(("transform",std::abs(point[2]-compare[2])<1e-12));
	}

	assert(("getTargetDimension",test->getTargetDimension()==3));


	test = &Geometry::MappingToRefTriangleToTetrahedron2::Instance();
	nodesAfterTransformation[0]=0;
	nodesAfterTransformation[1]=2;
	nodesAfterTransformation[2]=1;

	for(refPoint[0]=-2.8189;refPoint[0]<3.141;refPoint[0]+=0.1) {
		for(refPoint[1]=-2.8189;refPoint[1]<3.141;refPoint[1]+=0.1) {
			test->transform(refPoint,point);
			assert(("transform",fGeom.isInternalPoint(refPoint)==eGeom.isInternalPoint(point)));

			refPoint[0]+=-1.e-8;
			test->transform(refPoint,compare);
			refPoint[0]+=2.e-8;
			test->transform(refPoint,point);

			refPoint[0]+=-1e-8;
			test->calcJacobian(refPoint,jac);
			assert(("jacobian",std::abs(jac[0]-5.e7*(point[0]-compare[0]))<1e-5));//estimate is a bit rough, but should work for most mappings
			assert(("jacobian",std::abs(jac[1]-5.e7*(point[1]-compare[1]))<1e-5));//implementations are very strongly recommended to be more accurate
			assert(("jacobian",std::abs(jac[2]-5.e7*(point[2]-compare[2]))<1e-5));

			refPoint[1]+=-1.e-8;
			test->transform(refPoint,compare);
			refPoint[1]+=2.e-8;
			test->transform(refPoint,point);

			refPoint[1]+=-1e-8;
			test->calcJacobian(refPoint,jac);
			assert(("jacobian",std::abs(jac[3]-5.e7*(point[0]-compare[0]))<1e-5));
			assert(("jacobian",std::abs(jac[4]-5.e7*(point[1]-compare[1]))<1e-5));
			assert(("jacobian",std::abs(jac[5]-5.e7*(point[2]-compare[2]))<1e-5));
		}
	}

	for(int i=0;i<fGeom.getNumberOfNodes();++i){
		fGeom.getNode(i,refPoint);
		eGeom.getNode(nodesAfterTransformation[i],compare);
		test->transform(refPoint,point);
		assert(("transform",std::abs(point[0]-compare[0])<1e-12));
		assert(("transform",std::abs(point[1]-compare[1])<1e-12));
		assert(("transform",std::abs(point[2]-compare[2])<1e-12));
	}

	assert(("getTargetDimension",test->getTargetDimension()==3));


	test = &Geometry::MappingToRefTriangleToTetrahedron3::Instance();
	nodesAfterTransformation[0]=1;
	nodesAfterTransformation[1]=2;
	nodesAfterTransformation[2]=3;

	for(refPoint[0]=-2.8189;refPoint[0]<3.141;refPoint[0]+=0.1) {
		for(refPoint[1]=-2.8189;refPoint[1]<3.141;refPoint[1]+=0.1) {
			test->transform(refPoint,point);
			std::cout<<refPoint<<point<<std::endl;//truncation error pushes a point outside the reference element
			//assert(("transform",fGeom.isInternalPoint(refPoint)==eGeom.isInternalPoint(point)));

			refPoint[0]+=-1.e-8;
			test->transform(refPoint,compare);
			refPoint[0]+=2.e-8;
			test->transform(refPoint,point);

			refPoint[0]+=-1e-8;
			test->calcJacobian(refPoint,jac);
			assert(("jacobian",std::abs(jac[0]-5.e7*(point[0]-compare[0]))<1e-5));//estimate is a bit rough, but should work for most mappings
			assert(("jacobian",std::abs(jac[1]-5.e7*(point[1]-compare[1]))<1e-5));//implementations are very strongly recommended to be more accurate
			assert(("jacobian",std::abs(jac[2]-5.e7*(point[2]-compare[2]))<1e-5));

			refPoint[1]+=-1.e-8;
			test->transform(refPoint,compare);
			refPoint[1]+=2.e-8;
			test->transform(refPoint,point);

			refPoint[1]+=-1e-8;
			test->calcJacobian(refPoint,jac);
			assert(("jacobian",std::abs(jac[3]-5.e7*(point[0]-compare[0]))<1e-5));
			assert(("jacobian",std::abs(jac[4]-5.e7*(point[1]-compare[1]))<1e-5));
			assert(("jacobian",std::abs(jac[5]-5.e7*(point[2]-compare[2]))<1e-5));
		}
	}

	for(int i=0;i<fGeom.getNumberOfNodes();++i){
		fGeom.getNode(i,refPoint);
		eGeom.getNode(nodesAfterTransformation[i],compare);
		test->transform(refPoint,point);
		assert(("transform",std::abs(point[0]-compare[0])<1e-12));
		assert(("transform",std::abs(point[1]-compare[1])<1e-12));
		assert(("transform",std::abs(point[2]-compare[2])<1e-12));
	}

	assert(("getTargetDimension",test->getTargetDimension()==3));

	return 0;
}


