/*
 * 030MappingToPhysicalTriangularPrism_UnitTest.cpp
 *
 *  Created on: Apr 9, 2014
 *      Author: brinkf
 */

//naming convention: <Digit><ClassName>_UnitTest.cpp where <Digit> is a number that will make sure
//the unit tests are ordered such that the first failing unit test indicate the culprit class and
//other 'unit' tests may assume correct execution of all prior unit tests

#include "Geometry/Mappings/MappingToPhysTriangularPrism.hpp"
#include "cassert"

#include "Geometry/ReferenceTriangularPrism.hpp"
#include "Geometry/PhysicalTriangularPrism.hpp"
//transformations should map internal points to internal points, external points to external points
//and nodes to nodes so construct the physical geometries such that this can be checked :(

bool isInternal3D(const Geometry::PointPhysical& p) {
	return p[2]>0&&(p[1]-p[0])>(1)&&(1.8*p[0]-1.4*p[1])>-0.84&&(1.5*p[0]-1.1*p[1])<-0.42&&p[2]<4.;
}

int main() {

	std::vector<unsigned int> pointIndexes;

	pointIndexes.push_back(4);
	pointIndexes.push_back(7);
	pointIndexes.push_back(18);

	std::vector<Geometry::PointPhysical> nodes3D;

	Geometry::PointPhysical point3D(3),compare3D(3);
	Geometry::PointReference refPoint3D(3);

	pointIndexes.push_back(13);
	pointIndexes.push_back(27);
	pointIndexes.push_back(38);

	for(double i=0.;i<10;++i){
		point3D[0]=1.+i/10.;
		point3D[1]=2.+i/10.;
		point3D[2]=0.;
		nodes3D.push_back(point3D);
	}
	for(double i=0.;i<10;++i){
		point3D[0]=2.+i/10.;
		point3D[1]=5.-i/10.;
		point3D[2]=0.;
		nodes3D.push_back(point3D);
	}
	for(double i=0.;i<10;++i){
		point3D[0]=1.+i/10.;
		point3D[1]=2.+i/10.;
		point3D[2]=4.;
		nodes3D.push_back(point3D);
	}
	for(double i=0.;i<10;++i){
		point3D[0]=2.+i/10.;
		point3D[1]=5.-i/10.;
		point3D[2]=4.;
		nodes3D.push_back(point3D);
	}

	Geometry::ReferenceTriangularPrism& rGeom3D = Geometry::ReferenceTriangularPrism::Instance();

	Geometry::PhysicalTriangularPrism oops3D(pointIndexes,nodes3D,&rGeom3D);
	pointIndexes[3]=24;
	Geometry::PhysicalTriangularPrism pGeom3D(pointIndexes,nodes3D,&rGeom3D);

	Geometry::MappingToPhysTriangularPrism mapping3D(&pGeom3D),reinit3D(&oops3D);
	reinit3D.reinit(&pGeom3D);

	Geometry::Jacobian jac3D(3,3);

	for(refPoint3D[0]=-1.5189;refPoint3D[0]<1.541;refPoint3D[0]+=0.2) {
		for(refPoint3D[1]=-1.5188;refPoint3D[1]<1.541;refPoint3D[1]+=0.2){
			for(refPoint3D[2]=-1.5188;refPoint3D[2]<1.541;refPoint3D[2]+=0.2){
				mapping3D.transform(refPoint3D,point3D);
				assert(("transform",rGeom3D.isInternalPoint(refPoint3D)==isInternal3D(point3D)));
				reinit3D.transform(refPoint3D,point3D);
				assert(("reinit",rGeom3D.isInternalPoint(refPoint3D)==isInternal3D(point3D)));

				refPoint3D[0]+=-1.e-8;
				mapping3D.transform(refPoint3D,compare3D);
				refPoint3D[0]+=2.e-8;
				mapping3D.transform(refPoint3D,point3D);

				refPoint3D[0]+=-1e-8;
				mapping3D.calcJacobian(refPoint3D,jac3D);
				assert(("jacobian",fabs(jac3D[0]-5.e7*(point3D[0]-compare3D[0]))<1e-5));//estimate is a bit rough, but should work for most mappings
				assert(("jacobian",fabs(jac3D[1]-5.e7*(point3D[1]-compare3D[1]))<1e-5));//implementations are strongly recommended to be more accurate
				assert(("jacobian",fabs(jac3D[2]-5.e7*(point3D[2]-compare3D[2]))<1e-5));

				refPoint3D[1]+=-1.e-8;
				mapping3D.transform(refPoint3D,compare3D);
				refPoint3D[1]+=2.e-8;
				mapping3D.transform(refPoint3D,point3D);

				refPoint3D[1]+=-1e-8;
				mapping3D.calcJacobian(refPoint3D,jac3D);
				assert(("jacobian",fabs(jac3D[3]-5.e7*(point3D[0]-compare3D[0]))<1e-5));
				assert(("jacobian",fabs(jac3D[4]-5.e7*(point3D[1]-compare3D[1]))<1e-5));
				assert(("jacobian",fabs(jac3D[5]-5.e7*(point3D[2]-compare3D[2]))<1e-5));

				refPoint3D[2]+=-1.e-8;
				mapping3D.transform(refPoint3D,compare3D);
				refPoint3D[2]+=2.e-8;
				mapping3D.transform(refPoint3D,point3D);

				refPoint3D[2]+=-1e-8;
				mapping3D.calcJacobian(refPoint3D,jac3D);
				assert(("jacobian",fabs(jac3D[6]-5.e7*(point3D[0]-compare3D[0]))<1e-5));
				assert(("jacobian",fabs(jac3D[7]-5.e7*(point3D[1]-compare3D[1]))<1e-5));
				assert(("jacobian",fabs(jac3D[8]-5.e7*(point3D[2]-compare3D[2]))<1e-5));
			}
		}
	}

	for(int i=0;i<rGeom3D.getNumberOfNodes();++i){
		rGeom3D.getNode(i,refPoint3D);
		pGeom3D.getNodeCoordinates(i,compare3D);
		mapping3D.transform(refPoint3D,point3D);
		assert(("transform",fabs(point3D[0]-compare3D[0])<1e-12)&&fabs(point3D[1]-compare3D[1])<1e-12&&fabs(point3D[2]-compare3D[2])<1e-12);
	}

	assert(("getTargetDimension",mapping3D.getTargetDimension()==3));

	for(int i=0;i<40;++i){
		mapping3D.getNodeCoordinates(i,compare3D);
		pGeom3D.getGlobalNodeCoordinates(i,point3D);
		assert(("getNodeCoordinates",compare3D==point3D));
	}

	return 0;
}
