/*
 This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
 http://hpgem.org/about-the-code/team
 
 This code is distributed using BSD 3-Clause License. A copy of which can found below.
 
 
 Copyright (c) 2014, University of Twente
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

//naming convention: <Digit><ClassName>_UnitTest.cpp where <Digit> is a number that will make sure
//the unit tests are ordered such that the first failing unit test indicate the culprit class and
//other 'unit' tests may assume correct execution of all prior unit tests

#include "Geometry/ElementGeometry.hpp"
#include <Logger.h>

#include "Geometry/ReferenceCube.hpp"
#include "Geometry/ReferenceHypercube.hpp"
#include "Geometry/ReferenceLine.hpp"
#include "Geometry/ReferencePyramid.hpp"
#include "Geometry/ReferenceSquare.hpp"
#include "Geometry/ReferenceTetrahedron.hpp"
#include "Geometry/ReferenceTriangle.hpp"
#include "Geometry/ReferenceTriangularPrism.hpp"
#include "Geometry/PhysicalHexahedron.hpp"
#include "Geometry/PhysicalLine.hpp"
#include "Geometry/PhysicalOctachoron.hpp"
#include "Geometry/PhysicalPyramid.hpp"
#include "Geometry/PhysicalQuadrilateral.hpp"
#include "Geometry/PhysicalTetrahedron.hpp"
#include "Geometry/PhysicalTriangle.hpp"
#include "Geometry/PhysicalTriangularPrism.hpp"
#include "Geometry/Mappings/MappingToPhysHypercubeLinear.hpp"
#include "Geometry/Mappings/MappingToPhysPyramid.hpp"
#include "Geometry/Mappings/MappingToPhysSimplexLinear.hpp"
#include "Geometry/Mappings/MappingToPhysTriangularPrism.hpp"

#include <cmath>
#include <typeinfo>

int main(){

	//dim 1

	std::vector<std::size_t> pointIndexes;
	std::vector<Geometry::PointPhysical> nodes1D;

	Geometry::PointPhysical point1D(1),compare1D(1);
	Geometry::PointReference orig1D(1);

	Geometry::Jacobian jac(1,1),jaccompare(1,1);

	pointIndexes.push_back(4);
	pointIndexes.push_back(7);

	for(double i=0.;i<1;i+= 0.1){
		point1D[0]=1.+i/10.;
		nodes1D.push_back(point1D);
	}

	Geometry::ElementGeometry* test=new Geometry::ElementGeometry(pointIndexes,nodes1D);
	Geometry::ElementGeometry copy(*test);

	logger.assert_always((typeid(Geometry::MappingToPhysHypercubeLinear<1>)==typeid(*test->getReferenceToPhysicalMap())),"getReferenceToPhysicalMap");
	logger.assert_always((typeid(Geometry::PhysicalLine)==typeid(*test->getPhysicalGeometry())),"getPhysicalGeometry");
	logger.assert_always((typeid(Geometry::ReferenceLine)==typeid(*test->getReferenceGeometry())),"getReferenceGeometry");
	logger.assert_always((test->getNrOfNodes()==2),"getNrOfNodes");

	for(orig1D[0]=-2.8189;orig1D[0]<3.141;orig1D[0]+=0.1){
		compare1D = test->getReferenceToPhysicalMap()->transform(orig1D);
		point1D = test->referenceToPhysical(orig1D);
		logger.assert_always((std::abs(point1D[0]-compare1D[0])<1e-12),"referenceToPhysical");

		jaccompare = test->getReferenceToPhysicalMap()->calcJacobian(orig1D);
		jac = test->calcJacobian(orig1D);
		logger.assert_always((std::abs(jac[0]-jaccompare[0])<1e-12),"calcJacobian");
	}
	std::cout<<*test<<std::endl;

	//dim2
	std::vector<Geometry::PointPhysical> nodes2D;

	Geometry::PointPhysical point2D(2),compare2D(2);
	Geometry::PointReference orig2D(2);
	jac.resize(2,2);
	jaccompare.resize(2,2);

	pointIndexes.push_back(10);

	for(double i=0.;i<1;i+= 0.1){
		point2D[0]=1.+i;
		point2D[1]=2.+i;
		nodes2D.push_back(point2D);
	}

	point2D[0]=3.5;
	point2D[1]=4.6;
	nodes2D.push_back(point2D);
	point2D[0]=6.7;
	point2D[1]=2.8;
	nodes2D.push_back(point2D);

	delete test;
	test = new Geometry::ElementGeometry(pointIndexes,nodes2D);

	logger.assert_always((typeid(Geometry::MappingToPhysSimplexLinear<2>)==typeid(*test->getReferenceToPhysicalMap())),"getReferenceToPhysicalMap");
	logger.assert_always((typeid(Geometry::PhysicalTriangle)==typeid(*test->getPhysicalGeometry())),"getPhysicalGeometry");
	logger.assert_always((typeid(Geometry::ReferenceTriangle)==typeid(*test->getReferenceGeometry())),"getReferenceGeometry");
	logger.assert_always((test->getNrOfNodes()==3),"getNrOfNodes");

	for(orig2D[0]=-2.8189;orig2D[0]<3.141;orig2D[0]+=0.1){
		for(orig2D[1]=-2.8189;orig2D[1]<3.141;orig2D[1]+=0.1){
			compare2D = test->getReferenceToPhysicalMap()->transform(orig2D);
			point2D = test->referenceToPhysical(orig2D);
			logger.assert_always((std::abs(point2D[0]-compare2D[0])<1e-12),"referenceToPhysical");
			logger.assert_always((std::abs(point2D[1]-compare2D[1])<1e-12),"referenceToPhysical");

			jaccompare = test->getReferenceToPhysicalMap()->calcJacobian(orig2D);
			jac = test->calcJacobian(orig2D);
			logger.assert_always((std::abs(jac[0]-jaccompare[0])<1e-12),"calcJacobian");
			logger.assert_always((std::abs(jac[1]-jaccompare[1])<1e-12),"calcJacobian");
			logger.assert_always((std::abs(jac[2]-jaccompare[2])<1e-12),"calcJacobian");
			logger.assert_always((std::abs(jac[3]-jaccompare[3])<1e-12),"calcJacobian");
		}
	}
	std::cout<<*test<<std::endl;

	pointIndexes.push_back(11);

	delete test;
	test = new Geometry::ElementGeometry(pointIndexes,nodes2D);

	logger.assert_always((typeid(Geometry::MappingToPhysHypercubeLinear<2>)==typeid(*test->getReferenceToPhysicalMap())),"getReferenceToPhysicalMap");
	logger.assert_always((typeid(Geometry::PhysicalQuadrilateral)==typeid(*test->getPhysicalGeometry())),"getPhysicalGeometry");
	logger.assert_always((typeid(Geometry::ReferenceSquare)==typeid(*test->getReferenceGeometry())),"getReferenceGeometry");
	logger.assert_always((test->getNrOfNodes()==4),"getNrOfNodes");

	for(orig2D[0]=-2.8189;orig2D[0]<3.141;orig2D[0]+=0.1){
		for(orig2D[1]=-2.8189;orig2D[1]<3.141;orig2D[1]+=0.1){
			compare2D = test->getReferenceToPhysicalMap()->transform(orig2D);
			point2D = test->referenceToPhysical(orig2D);
			logger.assert_always((std::abs(point2D[0]-compare2D[0])<1e-12),"referenceToPhysical");
			logger.assert_always((std::abs(point2D[1]-compare2D[1])<1e-12),"referenceToPhysical");

			jaccompare = test->getReferenceToPhysicalMap()->calcJacobian(orig2D);
			jac = test->calcJacobian(orig2D);
			logger.assert_always((std::abs(jac[0]-jaccompare[0])<1e-12),"calcJacobian");
			logger.assert_always((std::abs(jac[1]-jaccompare[1])<1e-12),"calcJacobian");
			logger.assert_always((std::abs(jac[2]-jaccompare[2])<1e-12),"calcJacobian");
			logger.assert_always((std::abs(jac[3]-jaccompare[3])<1e-12),"calcJacobian");
		}
	}
	std::cout<<*test<<std::endl;

	//dim 3

	std::vector<Geometry::PointPhysical> nodes3D;

	Geometry::PointPhysical point3D(3),compare3D(3);
	Geometry::PointReference orig3D(3);
	jac.resize(3,3);
	jaccompare.resize(3,3);

	for(double i=0.;i<1;i+=0.1){
		point3D[0]=1.+i/10.;
		point3D[1]=2.+i/10.;
		point3D[2]=3.+i/10.;
		nodes3D.push_back(point3D);
	}

	point3D[0]=3.5;
	point3D[1]=4.6;
	point3D[2]=5.4;
	nodes3D.push_back(point3D);
	point3D[0]=6.7;
	point3D[1]=2.8;
	point3D[2]=5.7;
	nodes3D.push_back(point3D);
	point3D[0]=1.4;
	point3D[1]=2.4;
	point3D[2]=5.4;
	nodes3D.push_back(point3D);
	point3D[0]=1.7;
	point3D[1]=2.7;
	point3D[2]=5.7;
	nodes3D.push_back(point3D);
	point3D[0]=3.5;
	point3D[1]=4.6;
	point3D[2]=7.4;
	nodes3D.push_back(point3D);
	point3D[0]=6.7;
	point3D[1]=2.8;
	point3D[2]=7.7;
	nodes3D.push_back(point3D);

	delete test;
	test = new Geometry::ElementGeometry(pointIndexes,nodes3D);

	logger.assert_always((typeid(Geometry::MappingToPhysSimplexLinear<3>)==typeid(*test->getReferenceToPhysicalMap())),"getReferenceToPhysicalMap");
	logger.assert_always((typeid(Geometry::PhysicalTetrahedron)==typeid(*test->getPhysicalGeometry())),"getPhysicalGeometry");
	logger.assert_always((typeid(Geometry::ReferenceTetrahedron)==typeid(*test->getReferenceGeometry())),"getReferenceGeometry");
	logger.assert_always((test->getNrOfNodes()==4),"getNrOfNodes");

	for(orig3D[0]=-2.8189;orig3D[0]<3.141;orig3D[0]+=0.2){
		for(orig3D[1]=-2.8189;orig3D[1]<3.141;orig3D[1]+=0.25){
			for(orig3D[2]=-2.8189;orig3D[2]<3.141;orig3D[2]+=0.3){
				compare3D = test->getReferenceToPhysicalMap()->transform(orig3D);
				point3D = test->referenceToPhysical(orig3D);
				logger.assert_always((std::abs(point3D[0]-compare3D[0])<1e-12),"referenceToPhysical");
				logger.assert_always((std::abs(point3D[1]-compare3D[1])<1e-12),"referenceToPhysical");
				logger.assert_always((std::abs(point3D[2]-compare3D[2])<1e-12),"referenceToPhysical");

				jaccompare = test->getReferenceToPhysicalMap()->calcJacobian(orig3D);
				jac = test->calcJacobian(orig3D);
				logger.assert_always((std::abs(jac[0]-jaccompare[0])<1e-12),"calcJacobian");
				logger.assert_always((std::abs(jac[1]-jaccompare[1])<1e-12),"calcJacobian");
				logger.assert_always((std::abs(jac[2]-jaccompare[2])<1e-12),"calcJacobian");
				logger.assert_always((std::abs(jac[3]-jaccompare[3])<1e-12),"calcJacobian");
				logger.assert_always((std::abs(jac[4]-jaccompare[4])<1e-12),"calcJacobian");
				logger.assert_always((std::abs(jac[5]-jaccompare[5])<1e-12),"calcJacobian");
				logger.assert_always((std::abs(jac[6]-jaccompare[6])<1e-12),"calcJacobian");
				logger.assert_always((std::abs(jac[7]-jaccompare[7])<1e-12),"calcJacobian");
				logger.assert_always((std::abs(jac[8]-jaccompare[8])<1e-12),"calcJacobian");
			}
		}
	}
	std::cout<<*test<<std::endl;

	pointIndexes.push_back(12);

	delete test;
	test = new Geometry::ElementGeometry(pointIndexes,nodes3D);

	logger.assert_always((typeid(Geometry::MappingToPhysPyramid)==typeid(*test->getReferenceToPhysicalMap())),"getReferenceToPhysicalMap");
	logger.assert_always((typeid(Geometry::PhysicalPyramid)==typeid(*test->getPhysicalGeometry())),"getPhysicalGeometry");
	logger.assert_always((typeid(Geometry::ReferencePyramid)==typeid(*test->getReferenceGeometry())),"getReferenceGeometry");
	logger.assert_always((test->getNrOfNodes()==5),"getNrOfNodes");

	for(orig3D[0]=-2.8189;orig3D[0]<3.141;orig3D[0]+=0.2){
		for(orig3D[1]=-2.8189;orig3D[1]<3.141;orig3D[1]+=0.25){
			for(orig3D[2]=-2.8189;orig3D[2]<3.141;orig3D[2]+=0.3){
				compare3D = test->getReferenceToPhysicalMap()->transform(orig3D);
				point3D = test->referenceToPhysical(orig3D);
				logger.assert_always((std::abs(point3D[0]-compare3D[0])<1e-12),"referenceToPhysical");
				logger.assert_always((std::abs(point3D[1]-compare3D[1])<1e-12),"referenceToPhysical");
				logger.assert_always((std::abs(point3D[2]-compare3D[2])<1e-12),"referenceToPhysical");

				jaccompare = test->getReferenceToPhysicalMap()->calcJacobian(orig3D);
				jac = test->calcJacobian(orig3D);
				logger.assert_always((std::abs(jac[0]-jaccompare[0])<1e-12),"calcJacobian");
				logger.assert_always((std::abs(jac[1]-jaccompare[1])<1e-12),"calcJacobian");
				logger.assert_always((std::abs(jac[2]-jaccompare[2])<1e-12),"calcJacobian");
				logger.assert_always((std::abs(jac[3]-jaccompare[3])<1e-12),"calcJacobian");
				logger.assert_always((std::abs(jac[4]-jaccompare[4])<1e-12),"calcJacobian");
				logger.assert_always((std::abs(jac[5]-jaccompare[5])<1e-12),"calcJacobian");
				logger.assert_always((std::abs(jac[6]-jaccompare[6])<1e-12),"calcJacobian");
				logger.assert_always((std::abs(jac[7]-jaccompare[7])<1e-12),"calcJacobian");
				logger.assert_always((std::abs(jac[8]-jaccompare[8])<1e-12),"calcJacobian");
			}
		}
	}
	std::cout<<*test<<std::endl;

	pointIndexes.push_back(13);

	delete test;
	test = new Geometry::ElementGeometry(pointIndexes,nodes3D);

	logger.assert_always((typeid(Geometry::MappingToPhysTriangularPrism)==typeid(*test->getReferenceToPhysicalMap())),"getReferenceToPhysicalMap");
	logger.assert_always((typeid(Geometry::PhysicalTriangularPrism)==typeid(*test->getPhysicalGeometry())),"getPhysicalGeometry");
	logger.assert_always((typeid(Geometry::ReferenceTriangularPrism)==typeid(*test->getReferenceGeometry())),"getReferenceGeometry");
	logger.assert_always((test->getNrOfNodes()==6),"getNrOfNodes");

	for(orig3D[0]=-2.8189;orig3D[0]<3.141;orig3D[0]+=0.2){
		for(orig3D[1]=-2.8189;orig3D[1]<3.141;orig3D[1]+=0.25){
			for(orig3D[2]=-2.8189;orig3D[2]<3.141;orig3D[2]+=0.3){
				compare3D = test->getReferenceToPhysicalMap()->transform(orig3D);
				point3D = test->referenceToPhysical(orig3D);
				logger.assert_always((std::abs(point3D[0]-compare3D[0])<1e-12),"referenceToPhysical");
				logger.assert_always((std::abs(point3D[1]-compare3D[1])<1e-12),"referenceToPhysical");
				logger.assert_always((std::abs(point3D[2]-compare3D[2])<1e-12),"referenceToPhysical");

				jaccompare = test->getReferenceToPhysicalMap()->calcJacobian(orig3D);
				jac = test->calcJacobian(orig3D);
				logger.assert_always((std::abs(jac[0]-jaccompare[0])<1e-12),"calcJacobian");
				logger.assert_always((std::abs(jac[1]-jaccompare[1])<1e-12),"calcJacobian");
				logger.assert_always((std::abs(jac[2]-jaccompare[2])<1e-12),"calcJacobian");
				logger.assert_always((std::abs(jac[3]-jaccompare[3])<1e-12),"calcJacobian");
				logger.assert_always((std::abs(jac[4]-jaccompare[4])<1e-12),"calcJacobian");
				logger.assert_always((std::abs(jac[5]-jaccompare[5])<1e-12),"calcJacobian");
				logger.assert_always((std::abs(jac[6]-jaccompare[6])<1e-12),"calcJacobian");
				logger.assert_always((std::abs(jac[7]-jaccompare[7])<1e-12),"calcJacobian");
				logger.assert_always((std::abs(jac[8]-jaccompare[8])<1e-12),"calcJacobian");
			}
		}
	}
	std::cout<<*test<<std::endl;

	pointIndexes.push_back(14);
	pointIndexes.push_back(15);

	delete test;
	test = new Geometry::ElementGeometry(pointIndexes,nodes3D);

	logger.assert_always((typeid(Geometry::MappingToPhysHypercubeLinear<3>)==typeid(*test->getReferenceToPhysicalMap())),"getReferenceToPhysicalMap");
	logger.assert_always((typeid(Geometry::PhysicalHexahedron)==typeid(*test->getPhysicalGeometry())),"getPhysicalGeometry");
	logger.assert_always((typeid(Geometry::ReferenceCube)==typeid(*test->getReferenceGeometry())),"getReferenceGeometry");
	logger.assert_always((test->getNrOfNodes()==8),"getNrOfNodes");

	for(orig3D[0]=-2.8189;orig3D[0]<3.141;orig3D[0]+=0.2){
		for(orig3D[1]=-2.8189;orig3D[1]<3.141;orig3D[1]+=0.25){
			for(orig3D[2]=-2.8189;orig3D[2]<3.141;orig3D[2]+=0.3){
				compare3D = test->getReferenceToPhysicalMap()->transform(orig3D);
				point3D = test->referenceToPhysical(orig3D);
				logger.assert_always((std::abs(point3D[0]-compare3D[0])<1e-12),"referenceToPhysical");
				logger.assert_always((std::abs(point3D[1]-compare3D[1])<1e-12),"referenceToPhysical");
				logger.assert_always((std::abs(point3D[2]-compare3D[2])<1e-12),"referenceToPhysical");

				jaccompare = test->getReferenceToPhysicalMap()->calcJacobian(orig3D);
				jac = test->calcJacobian(orig3D);
				logger.assert_always((std::abs(jac[0]-jaccompare[0])<1e-12),"calcJacobian");
				logger.assert_always((std::abs(jac[1]-jaccompare[1])<1e-12),"calcJacobian");
				logger.assert_always((std::abs(jac[2]-jaccompare[2])<1e-12),"calcJacobian");
				logger.assert_always((std::abs(jac[3]-jaccompare[3])<1e-12),"calcJacobian");
				logger.assert_always((std::abs(jac[4]-jaccompare[4])<1e-12),"calcJacobian");
				logger.assert_always((std::abs(jac[5]-jaccompare[5])<1e-12),"calcJacobian");
				logger.assert_always((std::abs(jac[6]-jaccompare[6])<1e-12),"calcJacobian");
				logger.assert_always((std::abs(jac[7]-jaccompare[7])<1e-12),"calcJacobian");
				logger.assert_always((std::abs(jac[8]-jaccompare[8])<1e-12),"calcJacobian");
			}
		}
	}
	std::cout<<*test<<std::endl;

	std::vector<Geometry::PointPhysical> nodes4D;

	Geometry::PointPhysical point4D(4),compare4D(4);
	Geometry::PointReference orig4D(4);
	jac.resize(4,4),jaccompare.resize(4,4);

	pointIndexes.push_back(16);
	pointIndexes.push_back(17);
	pointIndexes.push_back(18);
	pointIndexes.push_back(19);
	pointIndexes.push_back(20);
	pointIndexes.push_back(21);
	pointIndexes.push_back(22);
	pointIndexes.push_back(23);

	for(double i=0.;i<1;i+=0.1){
		point4D[0]=1.+i;
		point4D[1]=2.+i;
		point4D[2]=3.+i;
		point4D[3]=1.;
		nodes4D.push_back(point4D);
	}

	point4D[0]=3.5;
	point4D[1]=4.6;
	point4D[2]=5.4;
	point4D[3]=1.;
	nodes4D.push_back(point4D);
	point4D[0]=6.7;
	point4D[1]=2.8;
	point4D[2]=5.7;
	point4D[3]=1.;
	nodes4D.push_back(point4D);
	point4D[0]=1.4;
	point4D[1]=2.4;
	point4D[2]=5.4;
	point4D[3]=2.;
	nodes4D.push_back(point4D);
	point4D[0]=1.7;
	point4D[1]=2.7;
	point4D[2]=5.7;
	point4D[3]=2.;
	nodes4D.push_back(point4D);
	point4D[0]=3.5;
	point4D[1]=4.6;
	point4D[2]=7.4;
	point4D[3]=2.;
	nodes4D.push_back(point4D);
	point4D[0]=6.7;
	point4D[1]=2.8;
	point4D[2]=7.7;
	point4D[3]=2.;
	nodes4D.push_back(point4D);
	point4D[0]=1.4;
	point4D[1]=2.4;
	point4D[2]=5.4;
	point4D[3]=3.;
	nodes4D.push_back(point4D);
	point4D[0]=1.7;
	point4D[1]=2.7;
	point4D[2]=5.7;
	point4D[3]=3.;
	nodes4D.push_back(point4D);
	point4D[0]=3.5;
	point4D[1]=4.6;
	point4D[2]=5.4;
	point4D[3]=3.;
	nodes4D.push_back(point4D);
	point4D[0]=6.7;
	point4D[1]=2.8;
	point4D[2]=5.7;
	point4D[3]=3.;
	nodes4D.push_back(point4D);
	point4D[0]=1.4;
	point4D[1]=2.4;
	point4D[2]=5.4;
	point4D[3]=4.;
	nodes4D.push_back(point4D);
	point4D[0]=1.7;
	point4D[1]=2.7;
	point4D[2]=5.7;
	point4D[3]=4.;
	nodes4D.push_back(point4D);
	point4D[0]=3.5;
	point4D[1]=4.6;
	point4D[2]=7.4;
	point4D[3]=4.;
	nodes4D.push_back(point4D);
	point4D[0]=6.7;
	point4D[1]=2.8;
	point4D[2]=7.7;
	point4D[3]=4.;
	nodes4D.push_back(point4D);

	delete test;
	test = new Geometry::ElementGeometry(pointIndexes,nodes4D);

	logger.assert_always((typeid(Geometry::MappingToPhysHypercubeLinear<4>)==typeid(*test->getReferenceToPhysicalMap())),"getReferenceToPhysicalMap");
	logger.assert_always((typeid(Geometry::PhysicalOctachoron)==typeid(*test->getPhysicalGeometry())),"getPhysicalGeometry");
	logger.assert_always((typeid(Geometry::ReferenceHypercube)==typeid(*test->getReferenceGeometry())),"getReferenceGeometry");
	logger.assert_always((test->getNrOfNodes()==16),"getNrOfNodes");

	for(orig4D[0]=-1.5189;orig4D[0]<1.541;orig4D[0]+=0.25){
		for(orig4D[1]=-1.5189;orig4D[1]<1.541;orig4D[1]+=0.25){
			for(orig4D[2]=-1.5189;orig4D[2]<1.541;orig4D[2]+=0.25){
				for(orig4D[3]=-1.5189;orig4D[3]<1.541;orig4D[3]+=0.3){
					compare4D = test->getReferenceToPhysicalMap()->transform(orig4D);
					point4D = test->referenceToPhysical(orig4D);
					logger.assert_always((std::abs(point4D[0]-compare4D[0])<1e-12),"referenceToPhysical");
					logger.assert_always((std::abs(point4D[1]-compare4D[1])<1e-12),"referenceToPhysical");
					logger.assert_always((std::abs(point4D[2]-compare4D[2])<1e-12),"referenceToPhysical");
					logger.assert_always((std::abs(point4D[3]-compare4D[3])<1e-12),"referenceToPhysical");

					jaccompare = test->getReferenceToPhysicalMap()->calcJacobian(orig4D);
					jac = test->calcJacobian(orig4D);
					logger.assert_always((std::abs(jac[0]-jaccompare[0])<1e-12),"calcJacobian");
					logger.assert_always((std::abs(jac[1]-jaccompare[1])<1e-12),"calcJacobian");
					logger.assert_always((std::abs(jac[2]-jaccompare[2])<1e-12),"calcJacobian");
					logger.assert_always((std::abs(jac[3]-jaccompare[3])<1e-12),"calcJacobian");
					logger.assert_always((std::abs(jac[4]-jaccompare[4])<1e-12),"calcJacobian");
					logger.assert_always((std::abs(jac[5]-jaccompare[5])<1e-12),"calcJacobian");
					logger.assert_always((std::abs(jac[6]-jaccompare[6])<1e-12),"calcJacobian");
					logger.assert_always((std::abs(jac[7]-jaccompare[7])<1e-12),"calcJacobian");
					logger.assert_always((std::abs(jac[8]-jaccompare[8])<1e-12),"calcJacobian");
					logger.assert_always((std::abs(jac[9]-jaccompare[9])<1e-12),"calcJacobian");
					logger.assert_always((std::abs(jac[10]-jaccompare[10])<1e-12),"calcJacobian");
					logger.assert_always((std::abs(jac[11]-jaccompare[11])<1e-12),"calcJacobian");
					logger.assert_always((std::abs(jac[12]-jaccompare[12])<1e-12),"calcJacobian");
					logger.assert_always((std::abs(jac[13]-jaccompare[13])<1e-12),"calcJacobian");
					logger.assert_always((std::abs(jac[14]-jaccompare[14])<1e-12),"calcJacobian");
					logger.assert_always((std::abs(jac[15]-jaccompare[15])<1e-12),"calcJacobian");
				}
			}
		}
	}
	std::cout<<*test<<std::endl;
	std::cout<<copy<<std::endl;

	return 0;
}

