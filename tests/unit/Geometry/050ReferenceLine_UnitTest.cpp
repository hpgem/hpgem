/*
 * 050ReferenceLine_UnitTest.cpp
 *
 *  Created on: Mar 27, 2014
 *      Author: brinkf
 */

//naming convention: <Digit><ClassName>_UnitTest.cpp where <Digit> is a number that will make sure
//the unit tests are ordered such that the first failing unit test indicate the culprit class and
//other 'unit' tests may assume correct execution of all prior unit tests

#include "Geometry/ReferenceLine.hpp"
#include <iostream>
#include "cassert"

#include "Geometry/PointReference.hpp"
#include "Geometry/Mappings/MappingToRefLineToLine.hpp"
#include "Geometry/Mappings/MappingToRefPointToLine.hpp"
using Geometry::ReferenceLine;

int main(){
	ReferenceLine& test=ReferenceLine::Instance();

	Geometry::PointReference pTest(1);

	//testing basic functionality

	for(pTest[0]=-3.141;pTest[0]<-1.;pTest[0]+=0.1){
		assert(("isInternalPoint",!test.isInternalPoint(pTest)));
	}
	for(;pTest[0]<1;pTest[0]+=0.1){
		assert(("isInternalPoint",test.isInternalPoint(pTest)));
	}
	for(;pTest[0]<3.141;pTest[0]+=0.1){
		assert(("isInternalPoint",!test.isInternalPoint(pTest)));
	}

	test.getCenter(pTest);
	assert(("getCenter",test.isInternalPoint(pTest)&&fabs(pTest[0])<1e-12));
	test.getNode(0,pTest);
	assert(("getNode 0",fabs(pTest[0]+1)<1e-12));
	test.getNode(1,pTest);
	assert(("getNode 1",fabs(pTest[0]-1)<1e-12));
	cout<<test.getName();

	assert(("getLocalNodeIndex 0",test.getLocalNodeIndex(0,0)==0));
	assert(("getLocalNodeIndex 1",test.getLocalNodeIndex(1,0)==1));

	cout<<test;

	//testing mappings and quadrature rules

	std::vector<unsigned int> base(2),transformed(2),faceIndices(1);
	for(int i=0;i<2;++i){
		base[i]=transformed[i]=i;
	}
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(test.getCodim0MappingIndex(base,transformed))==&Geometry::MappingToRefLineToLine0::Instance()));
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(base,transformed)==&Geometry::MappingToRefLineToLine0::Instance()));
	transformed[0]=1;
	transformed[1]=0;
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(test.getCodim0MappingIndex(base,transformed))==&Geometry::MappingToRefLineToLine1::Instance()));
	assert(("getCodim0MappingIndex&Ptr",test.getCodim0MappingPtr(base,transformed)==&Geometry::MappingToRefLineToLine1::Instance()));

	assert(("higher codimensional entities",test.getNrOfCodim1Entities()==2&&test.getNrOfCodim2Entities()==0)&&test.getNrOfCodim3Entities()==0);
	assert(("getCodim1ReferenceGeometry",test.getCodim1ReferenceGeometry(0)==&Geometry::ReferencePoint::Instance()&&
										 test.getCodim1ReferenceGeometry(1)==&Geometry::ReferencePoint::Instance()));
	assert(("getCodim1MappingPtr",test.getCodim1MappingPtr(0)==&Geometry::MappingToRefPointToLine0::Instance()));
	assert(("getCodim1MappingPtr",test.getCodim1MappingPtr(1)==&Geometry::MappingToRefPointToLine1::Instance()));
	test.getCodim1EntityLocalIndices(0,faceIndices);
	assert(("getCodim1EntityLocalIndices",faceIndices[0]==test.getLocalNodeIndex(0,0)));
	test.getCodim1EntityLocalIndices(1,faceIndices);
	assert(("getCodim1EntityLocalIndices",faceIndices[0]==test.getLocalNodeIndex(1,0)));

	assert(("quadrature rules",test.getGaussQuadratureRule(3)->order()>=3));
	assert(("quadrature rules",test.getGaussQuadratureRule(5)->order()>=5));
	assert(("quadrature rules",test.getGaussQuadratureRule(7)->order()>=7));
	assert(("quadrature rules",test.getGaussQuadratureRule(9)->order()>=9));
	assert(("quadrature rules",test.getGaussQuadratureRule(11)->order()>=11));

	//testing functionality of abstract parent classes

	assert(("number of nodes",test.getNumberOfNodes()==2));
	assert(("type of geometry",test.getGeometryType()==Geometry::LINE));

	///\TODO if it is decided that getBasisFunctionValue and getBasisFucntionDerivative remain here, test them

	///\TODO testing that the refinement maps behave exactly like the forwarded calls of this class
}



