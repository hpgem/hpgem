/*
 * 040ReferencePoint_UnitTest.cpp
 *
 *  Created on: Mar 27, 2014
 *      Author: brinkf
 */

//naming convention: <Digit><ClassName>_UnitTest.cpp where <Digit> is a number that will make sure
//the unit tests are ordered such that the first failing unit test indicate the culprit class and
//other 'unit' tests may assume correct execution of all prior unit tests

#include "Geometry/ReferencePoint.hpp"
#include <iostream>
#include "cassert"

#include "Geometry/PointReference.hpp"
#include "Geometry/Mappings/MappingToRefPointToPoint.hpp"
using Geometry::ReferencePoint;

int main(){
	ReferencePoint& test=ReferencePoint::Instance();

	Geometry::PointReference pTest(0);

	//testing basic functionality

	assert(("isInternalPoint",test.isInternalPoint(pTest)));
	test.getCenter(pTest);
	test.getNode(0,pTest);
	cout<<test.getName();

	//getLocalNodeIndex should always break since dimension -1 entities dont exist

	std::vector<unsigned int> base(1),transformed(1);

	//testing mappings and quadrature rules

	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(test.getCodim0MappingIndex(base,transformed))==&Geometry::MappingToRefPointToPoint::Instance()));
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(base,transformed)==&Geometry::MappingToRefPointToPoint::Instance()));
	assert(("higher codimensional entities",test.getNrOfCodim1Entities()==0&&test.getNrOfCodim2Entities()==0)&&test.getNrOfCodim3Entities()==0);

	assert(("quadrature rules",test.getGaussQuadratureRule(3)->order()>=3));
	assert(("quadrature rules",test.getGaussQuadratureRule(5)->order()>=5));
	assert(("quadrature rules",test.getGaussQuadratureRule(7)->order()>=7));
	assert(("quadrature rules",test.getGaussQuadratureRule(9)->order()>=9));
	assert(("quadrature rules",test.getGaussQuadratureRule(11)->order()>=11));

	//testing functionality of abstract parent classes

	assert(("number of nodes",test.getNumberOfNodes()==1));
	assert(("type of geometry",test.getGeometryType()==Geometry::POINT));

	//getBasisFunctionValue and getBasisFunctionDerivative require 0D basisfunctions to be implemented

	///\TODO testing that the refinement maps behave exactly like the forwarded calls of this class
}


