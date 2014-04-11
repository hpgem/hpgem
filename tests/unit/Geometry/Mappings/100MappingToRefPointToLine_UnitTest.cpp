/*
 * 100MappingToRefPointToLine_UnitTest.cpp
 *
 *  Created on: Apr 10, 2014
 *      Author: brinkf
 */

//naming convention: <Digit><ClassName>_UnitTest.cpp where <Digit> is a number that will make sure
//the unit tests are ordered such that the first failing unit test indicate the culprit class and
//other 'unit' tests may assume correct execution of all prior unit tests

#include "Geometry/Mappings/MappingToRefPointToLine.hpp"
#include "cassert"

#include "Geometry/ReferenceLine.hpp"
#include "Geometry/ReferencePoint.hpp"
int main() {

	Geometry::PointReference refPoint(0),point(1),compare(1);

	Geometry::ReferenceLine& eGeom = Geometry::ReferenceLine::Instance();
	Geometry::ReferencePoint& fGeom = Geometry::ReferencePoint::Instance();

	Geometry::Jacobian jac(0,1);

	std::vector<int> nodesAfterTransformation(1);

	const Geometry::MappingReferenceToReference* test = &Geometry::MappingToRefPointToLine0::Instance();
	nodesAfterTransformation[0]=0;

	test->transform(refPoint,point);
	assert(("transform",fGeom.isInternalPoint(refPoint)==eGeom.isInternalPoint(point)));

	test->calcJacobian(refPoint,jac);

	for(int i=0;i<fGeom.getNumberOfNodes();++i){
		fGeom.getNode(i,refPoint);
		eGeom.getNode(nodesAfterTransformation[i],compare);
		test->transform(refPoint,point);
		assert(("transform",fabs(point[0]-compare[0])<1e-12));
	}

	assert(("getTargetDimension",test->getTargetDimension()==1));


	test = &Geometry::MappingToRefPointToLine1::Instance();
	nodesAfterTransformation[0]=1;

	test->transform(refPoint,point);
	assert(("transform",fGeom.isInternalPoint(refPoint)==eGeom.isInternalPoint(point)));

	test->calcJacobian(refPoint,jac);

	for(int i=0;i<fGeom.getNumberOfNodes();++i){
		fGeom.getNode(i,refPoint);
		eGeom.getNode(nodesAfterTransformation[i],compare);
		test->transform(refPoint,point);
	}

	assert(("getTargetDimension",test->getTargetDimension()==1));

	return 0;
}

