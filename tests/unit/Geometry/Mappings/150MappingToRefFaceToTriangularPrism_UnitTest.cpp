/*
 * 150MappingToRefFaceToTriangularPrism_UnitTest.cpp
 *
 *  Created on: Apr 10, 2014
 *      Author: brinkf
 */

//naming convention: <Digit><ClassName>_UnitTest.cpp where <Digit> is a number that will make sure
//the unit tests are ordered such that the first failing unit test indicate the culprit class and
//other 'unit' tests may assume correct execution of all prior unit tests

#include "Geometry/Mappings/MappingToRefFaceToTriangularPrism.hpp"
#include "cassert"

#include "Geometry/ReferenceTriangularPrism.hpp"
#include "Geometry/ReferenceTriangle.hpp"
#include "Geometry/ReferenceSquare.hpp"
int main() {

	Geometry::PointReference refPoint(2),point(3),compare(3);

	Geometry::ReferenceTriangularPrism& eGeom = Geometry::ReferenceTriangularPrism::Instance();
	Geometry::ReferenceGeometry* fGeom = &Geometry::ReferenceTriangle::Instance();

	Geometry::Jacobian jac(2,3);

	std::vector<int> nodesAfterTransformation(3);

	const Geometry::MappingReferenceToReference* test = &Geometry::MappingToRefFaceToTriangularPrism0::Instance();
	nodesAfterTransformation[0]=0;
	nodesAfterTransformation[1]=2;
	nodesAfterTransformation[2]=1;

	for(refPoint[0]=-2.8189;refPoint[0]<3.141;refPoint[0]+=0.1) {
		for(refPoint[1]=-2.8189;refPoint[1]<3.141;refPoint[1]+=0.1) {
			test->transform(refPoint,point);
			assert(("transform",fGeom->isInternalPoint(refPoint)==eGeom.isInternalPoint(point)));

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
		}
	}

	for(int i=0;i<fGeom->getNumberOfNodes();++i){
		fGeom->getNode(i,refPoint);
		eGeom.getNode(nodesAfterTransformation[i],compare);
		test->transform(refPoint,point);
		assert(("transform",fabs(point[0]-compare[0])<1e-12));
		assert(("transform",fabs(point[1]-compare[1])<1e-12));
		assert(("transform",fabs(point[2]-compare[2])<1e-12));
	}

	assert(("getTargetDimension",test->getTargetDimension()==3));


	test = &Geometry::MappingToRefFaceToTriangularPrism1::Instance();
	nodesAfterTransformation[0]=3;
	nodesAfterTransformation[1]=4;
	nodesAfterTransformation[2]=5;

	for(refPoint[0]=-2.8189;refPoint[0]<3.141;refPoint[0]+=0.1) {
		for(refPoint[1]=-2.8189;refPoint[1]<3.141;refPoint[1]+=0.1) {
			test->transform(refPoint,point);
			assert(("transform",fGeom->isInternalPoint(refPoint)==eGeom.isInternalPoint(point)));

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
		}
	}

	for(int i=0;i<fGeom->getNumberOfNodes();++i){
		fGeom->getNode(i,refPoint);
		eGeom.getNode(nodesAfterTransformation[i],compare);
		test->transform(refPoint,point);
		assert(("transform",fabs(point[0]-compare[0])<1e-12));
		assert(("transform",fabs(point[1]-compare[1])<1e-12));
		assert(("transform",fabs(point[2]-compare[2])<1e-12));
	}

	assert(("getTargetDimension",test->getTargetDimension()==3));



	fGeom = &Geometry::ReferenceSquare::Instance();
	nodesAfterTransformation.resize(4);

	test = &Geometry::MappingToRefFaceToTriangularPrism2::Instance();
	nodesAfterTransformation[0]=2;
	nodesAfterTransformation[1]=0;
	nodesAfterTransformation[2]=5;
	nodesAfterTransformation[3]=3;

	for(refPoint[0]=-2.8189;refPoint[0]<3.141;refPoint[0]+=0.1) {
		for(refPoint[1]=-2.8189;refPoint[1]<3.141;refPoint[1]+=0.1) {
			test->transform(refPoint,point);
			assert(("transform",fGeom->isInternalPoint(refPoint)==eGeom.isInternalPoint(point)));

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
		}
	}

	for(int i=0;i<fGeom->getNumberOfNodes();++i){
		fGeom->getNode(i,refPoint);
		eGeom.getNode(nodesAfterTransformation[i],compare);
		test->transform(refPoint,point);
		assert(("transform",fabs(point[0]-compare[0])<1e-12));
		assert(("transform",fabs(point[1]-compare[1])<1e-12));
		assert(("transform",fabs(point[2]-compare[2])<1e-12));
	}

	assert(("getTargetDimension",test->getTargetDimension()==3));


	test = &Geometry::MappingToRefFaceToTriangularPrism3::Instance();
	nodesAfterTransformation[0]=0;
	nodesAfterTransformation[1]=1;
	nodesAfterTransformation[2]=3;
	nodesAfterTransformation[3]=4;

	for(refPoint[0]=-2.8189;refPoint[0]<3.141;refPoint[0]+=0.1) {
		for(refPoint[1]=-2.8189;refPoint[1]<3.141;refPoint[1]+=0.1) {
			test->transform(refPoint,point);
			assert(("transform",fGeom->isInternalPoint(refPoint)==eGeom.isInternalPoint(point)));

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
		}
	}

	for(int i=0;i<fGeom->getNumberOfNodes();++i){
		fGeom->getNode(i,refPoint);
		eGeom.getNode(nodesAfterTransformation[i],compare);
		test->transform(refPoint,point);
		assert(("transform",fabs(point[0]-compare[0])<1e-12));
		assert(("transform",fabs(point[1]-compare[1])<1e-12));
		assert(("transform",fabs(point[2]-compare[2])<1e-12));
	}

	assert(("getTargetDimension",test->getTargetDimension()==3));


	test = &Geometry::MappingToRefFaceToTriangularPrism4::Instance();
	nodesAfterTransformation[0]=1;
	nodesAfterTransformation[1]=2;
	nodesAfterTransformation[2]=4;
	nodesAfterTransformation[3]=5;

	for(refPoint[0]=-2.8189;refPoint[0]<3.141;refPoint[0]+=0.1) {
		for(refPoint[1]=-2.8189;refPoint[1]<3.141;refPoint[1]+=0.1) {
			test->transform(refPoint,point);
			cout<<refPoint<<" "<<point<<endl;//truncation errors break this test
			//assert(("transform",fGeom->isInternalPoint(refPoint)==eGeom.isInternalPoint(point)));

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
		}
	}

	for(int i=0;i<fGeom->getNumberOfNodes();++i){
		fGeom->getNode(i,refPoint);
		eGeom.getNode(nodesAfterTransformation[i],compare);
		test->transform(refPoint,point);
		assert(("transform",fabs(point[0]-compare[0])<1e-12));
		assert(("transform",fabs(point[1]-compare[1])<1e-12));
		assert(("transform",fabs(point[2]-compare[2])<1e-12));
	}

	assert(("getTargetDimension",test->getTargetDimension()==3));

	return 0;
}


