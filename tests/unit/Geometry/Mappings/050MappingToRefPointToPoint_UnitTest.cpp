/*
 * 050MappingToRefPointToPoint_UnitTest.cpp
 *
 *  Created on: Apr 10, 2014
 *      Author: brinkf
 */

//naming convention: <Digit><ClassName>_UnitTest.cpp where <Digit> is a number that will make sure
//the unit tests are ordered such that the first failing unit test indicate the culprit class and
//other 'unit' tests may assume correct execution of all prior unit tests

#include "Geometry/Mappings/MappingToRefPointToPoint.hpp"
#include "cassert"

#include "Geometry/ReferencePoint.hpp"

int main() {//The 0D case is mostly testing if there are any crashing functions

	Geometry::PointReference refPoint(0),point(0),compare(0);

	Geometry::ReferencePoint& geom = Geometry::ReferencePoint::Instance();

	const Geometry::MappingReferenceToReference* test = &Geometry::MappingToRefPointToPoint::Instance();

	Geometry::Jacobian jac(0,0);

	test->transform(refPoint,point);
	assert(("transform",geom.isInternalPoint(refPoint)==geom.isInternalPoint(point)));

	test->calcJacobian(refPoint,jac);

	for(int i=0;i<geom.getNumberOfNodes();++i){
		geom.getNode(i,refPoint);
		geom.getNode(i,compare);
		test->transform(refPoint,point);
	}

	assert(("getTargetDimension",test->getTargetDimension()==0));

	return 0;
}



