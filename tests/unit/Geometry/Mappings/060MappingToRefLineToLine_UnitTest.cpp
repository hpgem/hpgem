/*
 * 060MappingToRefLineToLine_UnitTest.cpp
 *
 *  Created on: Apr 10, 2014
 *      Author: brinkf
 */

//naming convention: <Digit><ClassName>_UnitTest.cpp where <Digit> is a number that will make sure
//the unit tests are ordered such that the first failing unit test indicate the culprit class and
//other 'unit' tests may assume correct execution of all prior unit tests

#include "Geometry/Mappings/MappingToRefLineToLine.hpp"
#include "cassert"

#include "Geometry/ReferenceLine.hpp"
int main() {

	Geometry::PointReference refPoint(1),point(1),compare(1);

	const Geometry::ReferenceLine& geom = Geometry::ReferenceLine::Instance();

	Geometry::Jacobian jac(1,1);

	std::vector<int> nodesAfterTransformation(2);

	const Geometry::MappingReferenceToReference* test = &Geometry::MappingToRefLineToLine0::Instance();
	nodesAfterTransformation[0]=0;
	nodesAfterTransformation[1]=1;

	for(refPoint[0]=-2.8189;refPoint[0]<3.141;refPoint[0]+=0.1) {
		test->transform(refPoint,point);
		assert(("transform",geom.isInternalPoint(refPoint)==geom.isInternalPoint(point)));

		refPoint[0]+=-1.e-8;
		test->transform(refPoint,compare);
		refPoint[0]+=2.e-8;
		test->transform(refPoint,point);

		refPoint[0]+=-1e-8;
		test->calcJacobian(refPoint,jac);
		assert(("jacobian",fabs(jac[0]-5.e7*(point[0]-compare[0]))<1e-5));//estimate is a bit rough, but should work for most mappings
	}

	for(int i=0;i<geom.getNumberOfNodes();++i){
		geom.getNode(i,refPoint);
		geom.getNode(nodesAfterTransformation[i],compare);
		test->transform(refPoint,point);
		assert(("transform",fabs(point[0]-compare[0])<1e-12));
	}

	assert(("getTargetDimension",test->getTargetDimension()==1));


	test = &Geometry::MappingToRefLineToLine1::Instance();
	nodesAfterTransformation[0]=1;
	nodesAfterTransformation[1]=0;

	for(refPoint[0]=-2.8189;refPoint[0]<3.141;refPoint[0]+=0.1) {
		test->transform(refPoint,point);
		assert(("transform",geom.isInternalPoint(refPoint)==geom.isInternalPoint(point)));

		refPoint[0]+=-1.e-8;
		test->transform(refPoint,compare);
		refPoint[0]+=2.e-8;
		test->transform(refPoint,point);

		refPoint[0]+=-1e-8;
		test->calcJacobian(refPoint,jac);
		assert(("jacobian",fabs(jac[0]-5.e7*(point[0]-compare[0]))<1e-5));//estimate is a bit rough, but should work for most mappings
	}

	for(int i=0;i<geom.getNumberOfNodes();++i){
		geom.getNode(i,refPoint);
		geom.getNode(nodesAfterTransformation[i],compare);
		test->transform(refPoint,point);
		assert(("transform",fabs(point[0]-compare[0])<1e-12));
	}

	assert(("getTargetDimension",test->getTargetDimension()==1));

	return 0;
}


