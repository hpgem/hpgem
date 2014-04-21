/*
 * 010MappingToPhysicalHypercubeLinear_UnitTest.cpp
 *
 *  Created on: Apr 2, 2014
 *      Author: brinkf
 */

//naming convention: <Digit><ClassName>_UnitTest.cpp where <Digit> is a number that will make sure
//the unit tests are ordered such that the first failing unit test indicate the culprit class and
//other 'unit' tests may assume correct execution of all prior unit tests

#include "Geometry/Mappings/MappingToPhysHypercubeLinear.hpp"
#include "cassert"

#include "Geometry/ReferenceLine.hpp"
#include "Geometry/PhysicalLine.hpp"
#include "Geometry/ReferenceSquare.hpp"
#include "Geometry/PhysicalQuadrilateral.hpp"
#include "Geometry/ReferenceCube.hpp"
#include "Geometry/PhysicalHexahedron.hpp"
//transformations should map internal points to internal points, external points to external points
//and nodes to nodes so construct the physical geometries such that this can be checked :(

bool isInternal1D(const Geometry::PointPhysical& p) {
	return p[0]>1.4&&p[0]<1.7;
}

bool isInternal2D(const Geometry::PointPhysical& p) {
	return (p[1]-p[0])>1.&&(p[1]+p[0])<7.&&(2.4*p[0]-0.8*p[1])>1.44&&(-1.5*p[0]+1.1*p[1])>.42;
}

bool isInternal3D(const Geometry::PointPhysical& p) {
	return (p[1]-p[0])>1.&&(p[1]+p[0])<7.&&(2.4*p[0]-0.8*p[1])>1.44&&(-1.5*p[0]+1.1*p[1])>.42&&p[2]>0.&&p[2]<4.;
}

int main() {

	//dim 1

	std::vector<unsigned int> pointIndexes;
	std::vector<Geometry::PointPhysical> nodes1D;

	Geometry::PointPhysical point1D(1),compare1D(1);
	Geometry::PointReference refPoint1D(1);

	pointIndexes.push_back(4);
	pointIndexes.push_back(6);

	for(double i=0.;i<10;++i){
		point1D[0]=1.+i/10.;
		nodes1D.push_back(point1D);
	}

	Geometry::ReferenceLine& rGeom = Geometry::ReferenceLine::Instance();

	Geometry::PhysicalLine oops(pointIndexes,nodes1D,&rGeom);
	pointIndexes[1]=7;
	Geometry::PhysicalLine pGeom(pointIndexes,nodes1D,&rGeom);

	Geometry::MappingToPhysHypercubeLinear<1> mapping1D(&pGeom),reinit1D(&oops);
	reinit1D.reinit(&pGeom);

	Geometry::Jacobian jac1D(1,1);

	for(refPoint1D[0]=-2.8189;refPoint1D[0]<3.141;refPoint1D[0]+=0.1) {
		mapping1D.transform(refPoint1D,point1D);
		assert(("transform",rGeom.isInternalPoint(refPoint1D)==isInternal1D(point1D)));
		reinit1D.transform(refPoint1D,point1D);
		assert(("reinit",rGeom.isInternalPoint(refPoint1D)==isInternal1D(point1D)));

		refPoint1D[0]+=-1.e-8;
		mapping1D.transform(refPoint1D,compare1D);
		refPoint1D[0]+=2.e-8;
		mapping1D.transform(refPoint1D,point1D);

		refPoint1D[0]+=-1e-8;
		mapping1D.calcJacobian(refPoint1D,jac1D);
		assert(("jacobian",fabs(jac1D[0]-5.e7*(point1D[0]-compare1D[0]))<1e-5));//estimate is a bit rough, but should work for most mappings
	}

	for(int i=0;i<rGeom.getNumberOfNodes();++i){
		rGeom.getNode(i,refPoint1D);
		pGeom.getNodeCoordinates(i,compare1D);
		mapping1D.transform(refPoint1D,point1D);
		assert(("transform",fabs(point1D[0]-compare1D[0])<1e-12));
	}

	assert(("getTargetDimension",mapping1D.getTargetDimension()==1));

	for(int i=0;i<10;++i){
		mapping1D.getNodeCoordinates(i,compare1D);
		pGeom.getGlobalNodeCoordinates(i,point1D);
		assert(("getNodeCoordinates",compare1D==point1D));
	}

	//dim2

	std::vector<Geometry::PointPhysical> nodes2D;

	Geometry::PointPhysical point2D(2),compare2D(2);
	Geometry::PointReference refPoint2D(2);

	pointIndexes.push_back(12);
	pointIndexes.push_back(13);

	for(double i=0.;i<10;++i){
		point2D[0]=1.+i/10.;
		point2D[1]=2.+i/10.;
		nodes2D.push_back(point2D);
	}
	for(double i=0.;i<10;++i){
		point2D[0]=2.+i/10.;
		point2D[1]=5.-i/10.;
		nodes2D.push_back(point2D);
	}

	Geometry::ReferenceSquare& rGeom2D = Geometry::ReferenceSquare::Instance();

	Geometry::PhysicalQuadrilateral oops2D(pointIndexes,nodes2D,&rGeom2D);
	pointIndexes[3]=18;
	Geometry::PhysicalQuadrilateral pGeom2D(pointIndexes,nodes2D,&rGeom2D);

	Geometry::MappingToPhysHypercubeLinear<2> mapping2D(&pGeom2D),reinit2D(&oops2D);
	reinit2D.reinit(&pGeom2D);

	Geometry::Jacobian jac2D(2,2);

	for(refPoint2D[0]=-1.5189;refPoint2D[0]<1.541;refPoint2D[0]+=0.1) {
		for(refPoint2D[1]=-1.5188;refPoint2D[1]<1.541;refPoint2D[1]+=0.1){
			mapping2D.transform(refPoint2D,point2D);
			assert(("transform"&&(rGeom2D.isInternalPoint(refPoint2D)&&isInternal2D(point2D))||(!rGeom2D.isInternalPoint(refPoint2D)&&!isInternal2D(point2D))));
			reinit2D.transform(refPoint2D,point2D);
			assert(("reinit",rGeom2D.isInternalPoint(refPoint2D)==isInternal2D(point2D)));

			refPoint2D[0]+=-1.e-8;
			mapping2D.transform(refPoint2D,compare2D);
			refPoint2D[0]+=2.e-8;
			mapping2D.transform(refPoint2D,point2D);

			refPoint2D[0]+=-1e-8;
			mapping2D.calcJacobian(refPoint2D,jac2D);
			assert(("jacobian",fabs(jac2D[0]-5.e7*(point2D[0]-compare2D[0]))<1e-5));//estimate is a bit rough, but should work for most mappings
			assert(("jacobian",fabs(jac2D[1]-5.e7*(point2D[1]-compare2D[1]))<1e-5));//implementations are strongly recommended to be more accurate

			refPoint2D[1]+=-1.e-8;
			mapping2D.transform(refPoint2D,compare2D);
			refPoint2D[1]+=2.e-8;
			mapping2D.transform(refPoint2D,point2D);

			refPoint2D[1]+=-1e-8;
			mapping2D.calcJacobian(refPoint2D,jac2D);
			assert(("jacobian",fabs(jac2D[2]-5.e7*(point2D[0]-compare2D[0]))<1e-5));
			assert(("jacobian",fabs(jac2D[3]-5.e7*(point2D[1]-compare2D[1]))<1e-5));
		}
	}

	for(int i=0;i<rGeom2D.getNumberOfNodes();++i){
		rGeom2D.getNode(i,refPoint2D);
		pGeom2D.getNodeCoordinates(i,compare2D);
		mapping2D.transform(refPoint2D,point2D);
		assert(("transform",fabs(point2D[0]-compare2D[0])<1e-12)&&fabs(point2D[1]-compare2D[1])<1e-12);
	}

	assert(("getTargetDimension",mapping2D.getTargetDimension()==2));

	for(int i=0;i<20;++i){
		mapping2D.getNodeCoordinates(i,compare2D);
		pGeom2D.getGlobalNodeCoordinates(i,point2D);
		assert(("getNodeCoordinates",compare2D==point2D));
	}

	//dim3

	std::vector<Geometry::PointPhysical> nodes3D;

	Geometry::PointPhysical point3D(3),compare3D(3);
	Geometry::PointReference refPoint3D(3);

	pointIndexes.push_back(24);
	pointIndexes.push_back(27);
	pointIndexes.push_back(32);
	pointIndexes.push_back(13);

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

	Geometry::ReferenceCube& rGeom3D = Geometry::ReferenceCube::Instance();

	Geometry::PhysicalHexahedron oops3D(pointIndexes,nodes3D,&rGeom3D);
	pointIndexes[7]=38;
	Geometry::PhysicalHexahedron pGeom3D(pointIndexes,nodes3D,&rGeom3D);

	Geometry::MappingToPhysHypercubeLinear<3> mapping3D(&pGeom3D),reinit3D(&oops3D);
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
