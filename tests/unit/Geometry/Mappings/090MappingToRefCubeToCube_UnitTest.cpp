/*
 * 090MappingToRefCubeToCube_UnitTest.cpp
 *
 *  Created on: Apr 10, 2014
 *      Author: brinkf
 */

//naming convention: <Digit><ClassName>_UnitTest.cpp where <Digit> is a number that will make sure
//the unit tests are ordered such that the first failing unit test indicate the culprit class and
//other 'unit' tests may assume correct execution of all prior unit tests

#include "Geometry/Mappings/MappingToRefCubeToCube.hpp"
#include "cassert"

#include "Geometry/ReferenceCube.hpp"
int main() {

	Geometry::PointReference refPoint(3),point(3),compare(3);

	Geometry::ReferenceCube& geom = Geometry::ReferenceCube::Instance();

	Geometry::Jacobian jac(3,3);

	std::vector<int> nodesAfterTransformation(4);

	const Geometry::MappingReferenceToReference* test = &Geometry::MappingToRefCubeToCube0::Instance();
	nodesAfterTransformation[0]=0;
	nodesAfterTransformation[1]=1;
	nodesAfterTransformation[2]=2;
	nodesAfterTransformation[3]=3;
	nodesAfterTransformation[4]=4;
	nodesAfterTransformation[5]=5;
	nodesAfterTransformation[6]=6;
	nodesAfterTransformation[7]=7;

	for(refPoint[0]=-2.8189;refPoint[0]<3.141;refPoint[0]+=0.1) {
		for(refPoint[1]=-2.8189;refPoint[1]<3.141;refPoint[1]+=0.1) {
			for(refPoint[2]=-2.8189;refPoint[2]<3.141;refPoint[2]+=0.1) {
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
				assert(("jacobian",fabs(jac[2]-5.e7*(point[2]-compare[2]))<1e-5));

				refPoint[1]+=-1.e-8;
				test->transform(refPoint,compare);
				refPoint[1]+=2.e-8;
				test->transform(refPoint,point);

				refPoint[1]+=-1e-8;
				test->calcJacobian(refPoint,jac);
				assert(("jacobian",fabs(jac[3]-5.e7*(point[0]-compare[0]))<1e-5));
				assert(("jacobian",fabs(jac[4]-5.e7*(point[1]-compare[1]))<1e-5));
				assert(("jacobian",fabs(jac[5]-5.e7*(point[2]-compare[2]))<1e-5));

				refPoint[2]+=-1.e-8;
				test->transform(refPoint,compare);
				refPoint[2]+=2.e-8;
				test->transform(refPoint,point);

				refPoint[2]+=-1e-8;
				test->calcJacobian(refPoint,jac);
				assert(("jacobian",fabs(jac[6]-5.e7*(point[0]-compare[0]))<1e-5));
				assert(("jacobian",fabs(jac[7]-5.e7*(point[1]-compare[1]))<1e-5));
				assert(("jacobian",fabs(jac[8]-5.e7*(point[2]-compare[2]))<1e-5));
			}
		}
	}

	for(int i=0;i<geom.getNumberOfNodes();++i){
		geom.getNode(i,refPoint);
		geom.getNode(nodesAfterTransformation[i],compare);
		test->transform(refPoint,point);
		assert(("transform",fabs(point[0]-compare[0])<1e-10 ));
		assert(("transform",fabs(point[1]-compare[1])<1e-10));
		assert(("transform",fabs(point[2]-compare[2])<1e-10));
	}

	return (0);

	assert(("getTargetDimension",test->getTargetDimension()==3));


	test = &Geometry::MappingToRefCubeToCube1::Instance();
	nodesAfterTransformation[0]=1;
	nodesAfterTransformation[1]=3;
	nodesAfterTransformation[2]=0;
	nodesAfterTransformation[3]=2;
	nodesAfterTransformation[4]=5;
	nodesAfterTransformation[5]=7;
	nodesAfterTransformation[6]=4;
	nodesAfterTransformation[7]=6;

	for(refPoint[0]=-2.8189;refPoint[0]<3.141;refPoint[0]+=0.1) {
		for(refPoint[1]=-2.8189;refPoint[1]<3.141;refPoint[1]+=0.1) {
			for(refPoint[2]=-2.8189;refPoint[2]<3.141;refPoint[2]+=0.1) {
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
				assert(("jacobian",fabs(jac[2]-5.e7*(point[2]-compare[2]))<1e-5));

				refPoint[1]+=-1.e-8;
				test->transform(refPoint,compare);
				refPoint[1]+=2.e-8;
				test->transform(refPoint,point);

				refPoint[1]+=-1e-8;
				test->calcJacobian(refPoint,jac);
				assert(("jacobian",fabs(jac[3]-5.e7*(point[0]-compare[0]))<1e-5));
				assert(("jacobian",fabs(jac[4]-5.e7*(point[1]-compare[1]))<1e-5));
				assert(("jacobian",fabs(jac[5]-5.e7*(point[2]-compare[2]))<1e-5));

				refPoint[2]+=-1.e-8;
				test->transform(refPoint,compare);
				refPoint[2]+=2.e-8;
				test->transform(refPoint,point);

				refPoint[2]+=-1e-8;
				test->calcJacobian(refPoint,jac);
				assert(("jacobian",fabs(jac[6]-5.e7*(point[0]-compare[0]))<1e-5));
				assert(("jacobian",fabs(jac[7]-5.e7*(point[1]-compare[1]))<1e-5));
				assert(("jacobian",fabs(jac[8]-5.e7*(point[2]-compare[2]))<1e-5));
			}
		}
	}

	/*for(int i=0;i<geom.getNumberOfNodes();++i){
		geom.getNode(i,refPoint);
		geom.getNode(nodesAfterTransformation[i],compare);
		test->transform(refPoint,point);
		assert(("transform",fabs(point[0]-compare[0])<1e-12));
		assert(("transform",fabs(point[1]-compare[1])<1e-12));
		assert(("transform",fabs(point[2]-compare[2])<1e-12));
	}*/

	assert(("getTargetDimension",test->getTargetDimension()==3));


	test = &Geometry::MappingToRefCubeToCube2::Instance();
	nodesAfterTransformation[0]=3;
	nodesAfterTransformation[1]=2;
	nodesAfterTransformation[2]=1;
	nodesAfterTransformation[3]=0;
	nodesAfterTransformation[4]=7;
	nodesAfterTransformation[5]=6;
	nodesAfterTransformation[6]=5;
	nodesAfterTransformation[7]=4;

	for(refPoint[0]=-2.8189;refPoint[0]<3.141;refPoint[0]+=0.1) {
		for(refPoint[1]=-2.8189;refPoint[1]<3.141;refPoint[1]+=0.1) {
			for(refPoint[2]=-2.8189;refPoint[2]<3.141;refPoint[2]+=0.1) {
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
				assert(("jacobian",fabs(jac[2]-5.e7*(point[2]-compare[2]))<1e-5));

				refPoint[1]+=-1.e-8;
				test->transform(refPoint,compare);
				refPoint[1]+=2.e-8;
				test->transform(refPoint,point);

				refPoint[1]+=-1e-8;
				test->calcJacobian(refPoint,jac);
				assert(("jacobian",fabs(jac[3]-5.e7*(point[0]-compare[0]))<1e-5));
				assert(("jacobian",fabs(jac[4]-5.e7*(point[1]-compare[1]))<1e-5));
				assert(("jacobian",fabs(jac[5]-5.e7*(point[2]-compare[2]))<1e-5));

				refPoint[2]+=-1.e-8;
				test->transform(refPoint,compare);
				refPoint[2]+=2.e-8;
				test->transform(refPoint,point);

				refPoint[2]+=-1e-8;
				test->calcJacobian(refPoint,jac);
				assert(("jacobian",fabs(jac[6]-5.e7*(point[0]-compare[0]))<1e-5));
				assert(("jacobian",fabs(jac[7]-5.e7*(point[1]-compare[1]))<1e-5));
				assert(("jacobian",fabs(jac[8]-5.e7*(point[2]-compare[2]))<1e-5));
			}
		}
	}

	/*for(int i=0;i<geom.getNumberOfNodes();++i){
		geom.getNode(i,refPoint);
		geom.getNode(nodesAfterTransformation[i],compare);
		test->transform(refPoint,point);
		assert(("transform",fabs(point[0]-compare[0])<1e-12));
		assert(("transform",fabs(point[1]-compare[1])<1e-12));
		assert(("transform",fabs(point[2]-compare[2])<1e-12));
	}*/

	assert(("getTargetDimension",test->getTargetDimension()==3));


	test = &Geometry::MappingToRefCubeToCube3::Instance();
	nodesAfterTransformation[0]=2;
	nodesAfterTransformation[1]=0;
	nodesAfterTransformation[2]=3;
	nodesAfterTransformation[3]=1;
	nodesAfterTransformation[4]=6;
	nodesAfterTransformation[5]=4;
	nodesAfterTransformation[6]=7;
	nodesAfterTransformation[7]=5;

	for(refPoint[0]=-2.8189;refPoint[0]<3.141;refPoint[0]+=0.1) {
		for(refPoint[1]=-2.8189;refPoint[1]<3.141;refPoint[1]+=0.1) {
			for(refPoint[2]=-2.8189;refPoint[2]<3.141;refPoint[2]+=0.1) {
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
				assert(("jacobian",fabs(jac[2]-5.e7*(point[2]-compare[2]))<1e-5));

				refPoint[1]+=-1.e-8;
				test->transform(refPoint,compare);
				refPoint[1]+=2.e-8;
				test->transform(refPoint,point);

				refPoint[1]+=-1e-8;
				test->calcJacobian(refPoint,jac);
				assert(("jacobian",fabs(jac[3]-5.e7*(point[0]-compare[0]))<1e-5));
				assert(("jacobian",fabs(jac[4]-5.e7*(point[1]-compare[1]))<1e-5));
				assert(("jacobian",fabs(jac[5]-5.e7*(point[2]-compare[2]))<1e-5));

				refPoint[2]+=-1.e-8;
				test->transform(refPoint,compare);
				refPoint[2]+=2.e-8;
				test->transform(refPoint,point);

				refPoint[2]+=-1e-8;
				test->calcJacobian(refPoint,jac);
				assert(("jacobian",fabs(jac[6]-5.e7*(point[0]-compare[0]))<1e-5));
				assert(("jacobian",fabs(jac[7]-5.e7*(point[1]-compare[1]))<1e-5));
				assert(("jacobian",fabs(jac[8]-5.e7*(point[2]-compare[2]))<1e-5));
			}
		}
	}

	/*for(int i=0;i<geom.getNumberOfNodes();++i){
		geom.getNode(i,refPoint);
		geom.getNode(nodesAfterTransformation[i],compare);
		test->transform(refPoint,point);
		assert(("transform",fabs(point[0]-compare[0])<1e-12));
		assert(("transform",fabs(point[1]-compare[1])<1e-12));
		assert(("transform",fabs(point[2]-compare[2])<1e-12));
	}*/

	assert(("getTargetDimension",test->getTargetDimension()==3));


	test = &Geometry::MappingToRefCubeToCube4::Instance();
	nodesAfterTransformation[0]=2;
	nodesAfterTransformation[1]=3;
	nodesAfterTransformation[2]=0;
	nodesAfterTransformation[3]=1;
	nodesAfterTransformation[4]=6;
	nodesAfterTransformation[5]=7;
	nodesAfterTransformation[6]=4;
	nodesAfterTransformation[7]=5;

	for(refPoint[0]=-2.8189;refPoint[0]<3.141;refPoint[0]+=0.1) {
		for(refPoint[1]=-2.8189;refPoint[1]<3.141;refPoint[1]+=0.1) {
			for(refPoint[2]=-2.8189;refPoint[2]<3.141;refPoint[2]+=0.1) {
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
				assert(("jacobian",fabs(jac[2]-5.e7*(point[2]-compare[2]))<1e-5));

				refPoint[1]+=-1.e-8;
				test->transform(refPoint,compare);
				refPoint[1]+=2.e-8;
				test->transform(refPoint,point);

				refPoint[1]+=-1e-8;
				test->calcJacobian(refPoint,jac);
				assert(("jacobian",fabs(jac[3]-5.e7*(point[0]-compare[0]))<1e-5));
				assert(("jacobian",fabs(jac[4]-5.e7*(point[1]-compare[1]))<1e-5));
				assert(("jacobian",fabs(jac[5]-5.e7*(point[2]-compare[2]))<1e-5));

				refPoint[2]+=-1.e-8;
				test->transform(refPoint,compare);
				refPoint[2]+=2.e-8;
				test->transform(refPoint,point);

				refPoint[2]+=-1e-8;
				test->calcJacobian(refPoint,jac);
				assert(("jacobian",fabs(jac[6]-5.e7*(point[0]-compare[0]))<1e-5));
				assert(("jacobian",fabs(jac[7]-5.e7*(point[1]-compare[1]))<1e-5));
				assert(("jacobian",fabs(jac[8]-5.e7*(point[2]-compare[2]))<1e-5));
			}
		}
	}

	/*for(int i=0;i<geom.getNumberOfNodes();++i){
		geom.getNode(i,refPoint);
		geom.getNode(nodesAfterTransformation[i],compare);
		test->transform(refPoint,point);
		assert(("transform",fabs(point[0]-compare[0])<1e-12));
		assert(("transform",fabs(point[1]-compare[1])<1e-12));
		assert(("transform",fabs(point[2]-compare[2])<1e-12));
	}*/

	assert(("getTargetDimension",test->getTargetDimension()==3));


	test = &Geometry::MappingToRefCubeToCube5::Instance();
	nodesAfterTransformation[0]=1;
	nodesAfterTransformation[1]=0;
	nodesAfterTransformation[2]=3;
	nodesAfterTransformation[3]=2;
	nodesAfterTransformation[4]=5;
	nodesAfterTransformation[5]=4;
	nodesAfterTransformation[6]=7;
	nodesAfterTransformation[7]=6;

	for(refPoint[0]=-2.8189;refPoint[0]<3.141;refPoint[0]+=0.1) {
		for(refPoint[1]=-2.8189;refPoint[1]<3.141;refPoint[1]+=0.1) {
			for(refPoint[2]=-2.8189;refPoint[2]<3.141;refPoint[2]+=0.1) {
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
				assert(("jacobian",fabs(jac[2]-5.e7*(point[2]-compare[2]))<1e-5));

				refPoint[1]+=-1.e-8;
				test->transform(refPoint,compare);
				refPoint[1]+=2.e-8;
				test->transform(refPoint,point);

				refPoint[1]+=-1e-8;
				test->calcJacobian(refPoint,jac);
				assert(("jacobian",fabs(jac[3]-5.e7*(point[0]-compare[0]))<1e-5));
				assert(("jacobian",fabs(jac[4]-5.e7*(point[1]-compare[1]))<1e-5));
				assert(("jacobian",fabs(jac[5]-5.e7*(point[2]-compare[2]))<1e-5));

				refPoint[2]+=-1.e-8;
				test->transform(refPoint,compare);
				refPoint[2]+=2.e-8;
				test->transform(refPoint,point);

				refPoint[2]+=-1e-8;
				test->calcJacobian(refPoint,jac);
				assert(("jacobian",fabs(jac[6]-5.e7*(point[0]-compare[0]))<1e-5));
				assert(("jacobian",fabs(jac[7]-5.e7*(point[1]-compare[1]))<1e-5));
				assert(("jacobian",fabs(jac[8]-5.e7*(point[2]-compare[2]))<1e-5));
			}
		}
	}

	/*for(int i=0;i<geom.getNumberOfNodes();++i){
		geom.getNode(i,refPoint);
		geom.getNode(nodesAfterTransformation[i],compare);
		test->transform(refPoint,point);
		assert(("transform",fabs(point[0]-compare[0])<1e-12));
		assert(("transform",fabs(point[1]-compare[1])<1e-12));
		assert(("transform",fabs(point[2]-compare[2])<1e-12));
	}*/

	assert(("getTargetDimension",test->getTargetDimension()==3));


	test = &Geometry::MappingToRefCubeToCube6::Instance();
	nodesAfterTransformation[0]=3;
	nodesAfterTransformation[1]=1;
	nodesAfterTransformation[2]=2;
	nodesAfterTransformation[3]=0;
	nodesAfterTransformation[4]=7;
	nodesAfterTransformation[5]=5;
	nodesAfterTransformation[6]=6;
	nodesAfterTransformation[7]=4;

	for(refPoint[0]=-2.8189;refPoint[0]<3.141;refPoint[0]+=0.1) {
		for(refPoint[1]=-2.8189;refPoint[1]<3.141;refPoint[1]+=0.1) {
			for(refPoint[2]=-2.8189;refPoint[2]<3.141;refPoint[2]+=0.1) {
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
				assert(("jacobian",fabs(jac[2]-5.e7*(point[2]-compare[2]))<1e-5));

				refPoint[1]+=-1.e-8;
				test->transform(refPoint,compare);
				refPoint[1]+=2.e-8;
				test->transform(refPoint,point);

				refPoint[1]+=-1e-8;
				test->calcJacobian(refPoint,jac);
				assert(("jacobian",fabs(jac[3]-5.e7*(point[0]-compare[0]))<1e-5));
				assert(("jacobian",fabs(jac[4]-5.e7*(point[1]-compare[1]))<1e-5));
				assert(("jacobian",fabs(jac[5]-5.e7*(point[2]-compare[2]))<1e-5));

				refPoint[2]+=-1.e-8;
				test->transform(refPoint,compare);
				refPoint[2]+=2.e-8;
				test->transform(refPoint,point);

				refPoint[2]+=-1e-8;
				test->calcJacobian(refPoint,jac);
				assert(("jacobian",fabs(jac[6]-5.e7*(point[0]-compare[0]))<1e-5));
				assert(("jacobian",fabs(jac[7]-5.e7*(point[1]-compare[1]))<1e-5));
				assert(("jacobian",fabs(jac[8]-5.e7*(point[2]-compare[2]))<1e-5));
			}
		}
	}

	/*for(int i=0;i<geom.getNumberOfNodes();++i){
		geom.getNode(i,refPoint);
		geom.getNode(nodesAfterTransformation[i],compare);
		test->transform(refPoint,point);
		assert(("transform",fabs(point[0]-compare[0])<1e-12));
		assert(("transform",fabs(point[1]-compare[1])<1e-12));
		assert(("transform",fabs(point[2]-compare[2])<1e-12));
	}*/

	assert(("getTargetDimension",test->getTargetDimension()==3));


	test = &Geometry::MappingToRefCubeToCube7::Instance();
	nodesAfterTransformation[0]=0;
	nodesAfterTransformation[1]=2;
	nodesAfterTransformation[2]=1;
	nodesAfterTransformation[3]=3;
	nodesAfterTransformation[4]=4;
	nodesAfterTransformation[5]=6;
	nodesAfterTransformation[6]=5;
	nodesAfterTransformation[7]=7;

	for(refPoint[0]=-2.8189;refPoint[0]<3.141;refPoint[0]+=0.1) {
		for(refPoint[1]=-2.8189;refPoint[1]<3.141;refPoint[1]+=0.1) {
			for(refPoint[2]=-2.8189;refPoint[2]<3.141;refPoint[2]+=0.1) {
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
				assert(("jacobian",fabs(jac[2]-5.e7*(point[2]-compare[2]))<1e-5));

				refPoint[1]+=-1.e-8;
				test->transform(refPoint,compare);
				refPoint[1]+=2.e-8;
				test->transform(refPoint,point);

				refPoint[1]+=-1e-8;
				test->calcJacobian(refPoint,jac);
				assert(("jacobian",fabs(jac[3]-5.e7*(point[0]-compare[0]))<1e-5));
				assert(("jacobian",fabs(jac[4]-5.e7*(point[1]-compare[1]))<1e-5));
				assert(("jacobian",fabs(jac[5]-5.e7*(point[2]-compare[2]))<1e-5));

				refPoint[2]+=-1.e-8;
				test->transform(refPoint,compare);
				refPoint[2]+=2.e-8;
				test->transform(refPoint,point);

				refPoint[2]+=-1e-8;
				test->calcJacobian(refPoint,jac);
				assert(("jacobian",fabs(jac[6]-5.e7*(point[0]-compare[0]))<1e-5));
				assert(("jacobian",fabs(jac[7]-5.e7*(point[1]-compare[1]))<1e-5));
				assert(("jacobian",fabs(jac[8]-5.e7*(point[2]-compare[2]))<1e-5));
			}
		}
	}

	/*for(int i=0;i<geom.getNumberOfNodes();++i){
		geom.getNode(i,refPoint);
		geom.getNode(nodesAfterTransformation[i],compare);
		test->transform(refPoint,point);
		assert(("transform",fabs(point[0]-compare[0])<1e-12));
		assert(("transform",fabs(point[1]-compare[1])<1e-12));
		assert(("transform",fabs(point[2]-compare[2])<1e-12));
	}*/

	assert(("getTargetDimension",test->getTargetDimension()==3));

	return 0;
}


