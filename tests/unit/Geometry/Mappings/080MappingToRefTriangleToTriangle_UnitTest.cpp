/*
 * 080MappingToRefTriangleToTriangle_UnitTest.cpp
 *
 *  Created on: Apr 10, 2014
 *      Author: brinkf
 */

//naming convention: <Digit><ClassName>_UnitTest.cpp where <Digit> is a number that will make sure
//the unit tests are ordered such that the first failing unit test indicate the culprit class and
//other 'unit' tests may assume correct execution of all prior unit tests

#include "Geometry/Mappings/MappingToRefTriangleToTriangle.hpp"
#include "cassert"

#include "Geometry/ReferenceTriangle.hpp"
int main() {

	Geometry::PointReference refPoint(2),point(2),compare(2);

	Geometry::ReferenceTriangle& geom = Geometry::ReferenceTriangle::Instance();

	Geometry::Jacobian jac(2,2);

	std::vector<int> nodesAfterTransformation(3);

	const Geometry::MappingReferenceToReference* test = &Geometry::MappingToRefTriangleToTriangle0::Instance();
	nodesAfterTransformation[0]=0;
	nodesAfterTransformation[1]=1;
	nodesAfterTransformation[2]=2;

	for(refPoint[0]=-2.8189;refPoint[0]<3.141;refPoint[0]+=0.1) {
		for(refPoint[1]=-2.8189;refPoint[1]<3.141;refPoint[1]+=0.1) {
			test->transform(refPoint,point);
			assert(("transform",geom.isInternalPoint(refPoint)==geom.isInternalPoint(point)));

			refPoint[0]+=-1.e-8;
			test->transform(refPoint,compare);
			refPoint[0]+=2.e-8;
			test->transform(refPoint,point);

			refPoint[0]+=-1e-8;
			test->calcJacobian(refPoint,jac);
			assert(("jacobian",fabs(jac[0]-5.e7*(point[0]-compare[0]))<1e-5));//estimate is a bit rough, but should work for most mappings
			assert(("jacobian",fabs(jac[1]-5.e7*(point[1]-compare[1]))<1e-5));//implementations are very strongly recommended to be more accurate

			refPoint[1]+=-1.e-8;
			test->transform(refPoint,compare);
			refPoint[1]+=2.e-8;
			test->transform(refPoint,point);

			refPoint[1]+=-1e-8;
			test->calcJacobian(refPoint,jac);
			assert(("jacobian",fabs(jac[2]-5.e7*(point[0]-compare[0]))<1e-5));
			assert(("jacobian",fabs(jac[3]-5.e7*(point[1]-compare[1]))<1e-5));
		}
	}

	for(int i=0;i<geom.getNumberOfNodes();++i){
		geom.getNode(i,refPoint);
		geom.getNode(nodesAfterTransformation[i],compare);
		test->transform(refPoint,point);
		assert(("transform",fabs(point[0]-compare[0])<1e-12));
		assert(("transform",fabs(point[1]-compare[1])<1e-12));
	}

	assert(("getTargetDimension",test->getTargetDimension()==2));


	test = &Geometry::MappingToRefTriangleToTriangle1::Instance();
	nodesAfterTransformation[0]=0;
	nodesAfterTransformation[1]=2;
	nodesAfterTransformation[2]=1;

	for(refPoint[0]=-2.8189;refPoint[0]<3.141;refPoint[0]+=0.1) {
		for(refPoint[1]=-2.8189;refPoint[1]<3.141;refPoint[1]+=0.1) {
			test->transform(refPoint,point);
			assert(("transform",geom.isInternalPoint(refPoint)==geom.isInternalPoint(point)));

			refPoint[0]+=-1.e-8;
			test->transform(refPoint,compare);
			refPoint[0]+=2.e-8;
			test->transform(refPoint,point);

			refPoint[0]+=-1e-8;
			test->calcJacobian(refPoint,jac);
			assert(("jacobian",fabs(jac[0]-5.e7*(point[0]-compare[0]))<1e-5));//estimate is a bit rough, but should work for most mappings
			assert(("jacobian",fabs(jac[1]-5.e7*(point[1]-compare[1]))<1e-5));//implementations are very strongly recommended to be more accurate

			refPoint[1]+=-1.e-8;
			test->transform(refPoint,compare);
			refPoint[1]+=2.e-8;
			test->transform(refPoint,point);

			refPoint[1]+=-1e-8;
			test->calcJacobian(refPoint,jac);
			assert(("jacobian",fabs(jac[2]-5.e7*(point[0]-compare[0]))<1e-5));
			assert(("jacobian",fabs(jac[3]-5.e7*(point[1]-compare[1]))<1e-5));
		}
	}

	for(int i=0;i<geom.getNumberOfNodes();++i){
		geom.getNode(i,refPoint);
		geom.getNode(nodesAfterTransformation[i],compare);
		test->transform(refPoint,point);
		assert(("transform",fabs(point[0]-compare[0])<1e-12));
		assert(("transform",fabs(point[1]-compare[1])<1e-12));
	}

	assert(("getTargetDimension",test->getTargetDimension()==2));


	test = &Geometry::MappingToRefTriangleToTriangle2::Instance();
	nodesAfterTransformation[0]=1;
	nodesAfterTransformation[1]=2;
	nodesAfterTransformation[2]=0;

	for(refPoint[0]=-2.8189;refPoint[0]<3.141;refPoint[0]+=0.1) {
		for(refPoint[1]=-2.8189;refPoint[1]<3.141;refPoint[1]+=0.1) {
			test->transform(refPoint,point);
			assert(("transform",geom.isInternalPoint(refPoint)==geom.isInternalPoint(point)));

			refPoint[0]+=-1.e-8;
			test->transform(refPoint,compare);
			refPoint[0]+=2.e-8;
			test->transform(refPoint,point);

			refPoint[0]+=-1e-8;
			test->calcJacobian(refPoint,jac);
			assert(("jacobian",fabs(jac[0]-5.e7*(point[0]-compare[0]))<1e-5));//estimate is a bit rough, but should work for most mappings
			assert(("jacobian",fabs(jac[1]-5.e7*(point[1]-compare[1]))<1e-5));//implementations are very strongly recommended to be more accurate

			refPoint[1]+=-1.e-8;
			test->transform(refPoint,compare);
			refPoint[1]+=2.e-8;
			test->transform(refPoint,point);

			refPoint[1]+=-1e-8;
			test->calcJacobian(refPoint,jac);
			assert(("jacobian",fabs(jac[2]-5.e7*(point[0]-compare[0]))<1e-5));
			assert(("jacobian",fabs(jac[3]-5.e7*(point[1]-compare[1]))<1e-5));
		}
	}

	for(int i=0;i<geom.getNumberOfNodes();++i){
		geom.getNode(i,refPoint);
		geom.getNode(nodesAfterTransformation[i],compare);
		test->transform(refPoint,point);
		assert(("transform",fabs(point[0]-compare[0])<1e-12));
		assert(("transform",fabs(point[1]-compare[1])<1e-12));
	}

	assert(("getTargetDimension",test->getTargetDimension()==2));


	test = &Geometry::MappingToRefTriangleToTriangle3::Instance();
	nodesAfterTransformation[0]=1;
	nodesAfterTransformation[1]=0;
	nodesAfterTransformation[2]=2;

	for(refPoint[0]=-2.8189;refPoint[0]<3.141;refPoint[0]+=0.1) {
		for(refPoint[1]=-2.8189;refPoint[1]<3.141;refPoint[1]+=0.1) {
			test->transform(refPoint,point);
			assert(("transform",geom.isInternalPoint(refPoint)==geom.isInternalPoint(point)));

			refPoint[0]+=-1.e-8;
			test->transform(refPoint,compare);
			refPoint[0]+=2.e-8;
			test->transform(refPoint,point);

			refPoint[0]+=-1e-8;
			test->calcJacobian(refPoint,jac);
			assert(("jacobian",fabs(jac[0]-5.e7*(point[0]-compare[0]))<1e-5));//estimate is a bit rough, but should work for most mappings
			assert(("jacobian",fabs(jac[1]-5.e7*(point[1]-compare[1]))<1e-5));//implementations are very strongly recommended to be more accurate

			refPoint[1]+=-1.e-8;
			test->transform(refPoint,compare);
			refPoint[1]+=2.e-8;
			test->transform(refPoint,point);

			refPoint[1]+=-1e-8;
			test->calcJacobian(refPoint,jac);
			assert(("jacobian",fabs(jac[2]-5.e7*(point[0]-compare[0]))<1e-5));
			assert(("jacobian",fabs(jac[3]-5.e7*(point[1]-compare[1]))<1e-5));
		}
	}

	for(int i=0;i<geom.getNumberOfNodes();++i){
		geom.getNode(i,refPoint);
		geom.getNode(nodesAfterTransformation[i],compare);
		test->transform(refPoint,point);
		assert(("transform",fabs(point[0]-compare[0])<1e-12));
		assert(("transform",fabs(point[1]-compare[1])<1e-12));
	}

	assert(("getTargetDimension",test->getTargetDimension()==2));


	test = &Geometry::MappingToRefTriangleToTriangle4::Instance();
	nodesAfterTransformation[0]=2;
	nodesAfterTransformation[1]=1;
	nodesAfterTransformation[2]=0;

	for(refPoint[0]=-2.8189;refPoint[0]<3.141;refPoint[0]+=0.1) {
		for(refPoint[1]=-2.8189;refPoint[1]<3.141;refPoint[1]+=0.1) {
			test->transform(refPoint,point);
			assert(("transform",geom.isInternalPoint(refPoint)==geom.isInternalPoint(point)));

			refPoint[0]+=-1.e-8;
			test->transform(refPoint,compare);
			refPoint[0]+=2.e-8;
			test->transform(refPoint,point);

			refPoint[0]+=-1e-8;
			test->calcJacobian(refPoint,jac);
			assert(("jacobian",fabs(jac[0]-5.e7*(point[0]-compare[0]))<1e-5));//estimate is a bit rough, but should work for most mappings
			assert(("jacobian",fabs(jac[1]-5.e7*(point[1]-compare[1]))<1e-5));//implementations are very strongly recommended to be more accurate

			refPoint[1]+=-1.e-8;
			test->transform(refPoint,compare);
			refPoint[1]+=2.e-8;
			test->transform(refPoint,point);

			refPoint[1]+=-1e-8;
			test->calcJacobian(refPoint,jac);
			assert(("jacobian",fabs(jac[2]-5.e7*(point[0]-compare[0]))<1e-5));
			assert(("jacobian",fabs(jac[3]-5.e7*(point[1]-compare[1]))<1e-5));
		}
	}

	for(int i=0;i<geom.getNumberOfNodes();++i){
		geom.getNode(i,refPoint);
		geom.getNode(nodesAfterTransformation[i],compare);
		test->transform(refPoint,point);
		assert(("transform",fabs(point[0]-compare[0])<1e-12));
		assert(("transform",fabs(point[1]-compare[1])<1e-12));
	}

	assert(("getTargetDimension",test->getTargetDimension()==2));


	test = &Geometry::MappingToRefTriangleToTriangle5::Instance();
	nodesAfterTransformation[0]=2;
	nodesAfterTransformation[1]=0;
	nodesAfterTransformation[2]=1;

	for(refPoint[0]=-2.8189;refPoint[0]<3.141;refPoint[0]+=0.1) {
		for(refPoint[1]=-2.8189;refPoint[1]<3.141;refPoint[1]+=0.1) {
			test->transform(refPoint,point);
			assert(("transform",geom.isInternalPoint(refPoint)==geom.isInternalPoint(point)));

			refPoint[0]+=-1.e-8;
			test->transform(refPoint,compare);
			refPoint[0]+=2.e-8;
			test->transform(refPoint,point);

			refPoint[0]+=-1e-8;
			test->calcJacobian(refPoint,jac);
			assert(("jacobian",fabs(jac[0]-5.e7*(point[0]-compare[0]))<1e-5));//estimate is a bit rough, but should work for most mappings
			assert(("jacobian",fabs(jac[1]-5.e7*(point[1]-compare[1]))<1e-5));//implementations are very strongly recommended to be more accurate

			refPoint[1]+=-1.e-8;
			test->transform(refPoint,compare);
			refPoint[1]+=2.e-8;
			test->transform(refPoint,point);

			refPoint[1]+=-1e-8;
			test->calcJacobian(refPoint,jac);
			assert(("jacobian",fabs(jac[2]-5.e7*(point[0]-compare[0]))<1e-5));
			assert(("jacobian",fabs(jac[3]-5.e7*(point[1]-compare[1]))<1e-5));
		}
	}

	for(int i=0;i<geom.getNumberOfNodes();++i){
		geom.getNode(i,refPoint);
		geom.getNode(nodesAfterTransformation[i],compare);
		test->transform(refPoint,point);
		assert(("transform",fabs(point[0]-compare[0])<1e-12));
		assert(("transform",fabs(point[1]-compare[1])<1e-12));
	}

	assert(("getTargetDimension",test->getTargetDimension()==2));

	return 0;
}


