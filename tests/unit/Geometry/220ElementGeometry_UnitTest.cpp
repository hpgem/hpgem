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
#include <cassert>

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

using std::fabs;
int main(){

	//dim 1

	std::vector<unsigned int> pointIndexes;
	std::vector<Geometry::PointPhysical> nodes1D;

	Geometry::PointPhysical point1D(1),compare1D(1);
	Geometry::PointReference orig1D(1);

	Geometry::Jacobian jac(1,1),jaccompare(1,1);

	pointIndexes.push_back(4);
	pointIndexes.push_back(7);

	for(double i=0.;i<10;++i){
		point1D[0]=1.+i/10.;
		nodes1D.push_back(point1D);
	}

	Geometry::ElementGeometry* test=new Geometry::ElementGeometry(pointIndexes,nodes1D);
	Geometry::ElementGeometry copy(*test);

	assert(("getReferenceToPhysicalMap",typeid(Geometry::MappingToPhysHypercubeLinear<1>)==typeid(*test->getReferenceToPhysicalMap())));
	assert(("getPhysicalGeometry",typeid(Geometry::PhysicalLine)==typeid(*test->getPhysicalGeometry())));
	assert(("getReferenceGeometry",typeid(Geometry::ReferenceLine)==typeid(*test->getReferenceGeometry())));
	assert(("getNrOfNodes",test->getNrOfNodes()==2));

	for(orig1D[0]=-2.8189;orig1D[0]<3.141;orig1D[0]+=0.1){
		test->getReferenceToPhysicalMap()->transform(orig1D,compare1D);
		test->referenceToPhysical(orig1D,point1D);
		assert(("referenceToPhysical",fabs(point1D[0]-compare1D[0])<1e-12));

		test->getReferenceToPhysicalMap()->calcJacobian(orig1D,jaccompare);
		test->calcJacobian(orig1D,jac);
		assert(("calcJacobian",fabs(jac[0]-jaccompare[0])<1e-12));
	}
	std::cout<<*test<<std::endl;

	//dim2
	std::vector<Geometry::PointPhysical> nodes2D;

	Geometry::PointPhysical point2D(2),compare2D(2);
	Geometry::PointReference orig2D(2);
	jac.resize(2,2);
	jaccompare.resize(2,2);

	pointIndexes.push_back(10);

	for(double i=0.;i<10;++i){
		point2D[0]=1.+i/10.;
		point2D[1]=2.+i/10.;
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

	assert(("getReferenceToPhysicalMap",typeid(Geometry::MappingToPhysSimplexLinear<2>)==typeid(*test->getReferenceToPhysicalMap())));
	assert(("getPhysicalGeometry",typeid(Geometry::PhysicalTriangle)==typeid(*test->getPhysicalGeometry())));
	assert(("getReferenceGeometry",typeid(Geometry::ReferenceTriangle)==typeid(*test->getReferenceGeometry())));
	assert(("getNrOfNodes",test->getNrOfNodes()==3));

	for(orig2D[0]=-2.8189;orig2D[0]<3.141;orig2D[0]+=0.1){
		for(orig2D[1]=-2.8189;orig2D[1]<3.141;orig2D[1]+=0.1){
			test->getReferenceToPhysicalMap()->transform(orig2D,compare2D);
			test->referenceToPhysical(orig2D,point2D);
			assert(("referenceToPhysical",fabs(point2D[0]-compare2D[0])<1e-12));
			assert(("referenceToPhysical",fabs(point2D[1]-compare2D[1])<1e-12));

			test->getReferenceToPhysicalMap()->calcJacobian(orig2D,jaccompare);
			test->calcJacobian(orig2D,jac);
			assert(("calcJacobian",fabs(jac[0]-jaccompare[0])<1e-12));
			assert(("calcJacobian",fabs(jac[1]-jaccompare[1])<1e-12));
			assert(("calcJacobian",fabs(jac[2]-jaccompare[2])<1e-12));
			assert(("calcJacobian",fabs(jac[3]-jaccompare[3])<1e-12));
		}
	}
	std::cout<<*test<<std::endl;

	pointIndexes.push_back(11);

	delete test;
	test = new Geometry::ElementGeometry(pointIndexes,nodes2D);

	assert(("getReferenceToPhysicalMap",typeid(Geometry::MappingToPhysHypercubeLinear<2>)==typeid(*test->getReferenceToPhysicalMap())));
	assert(("getPhysicalGeometry",typeid(Geometry::PhysicalQuadrilateral)==typeid(*test->getPhysicalGeometry())));
	assert(("getReferenceGeometry",typeid(Geometry::ReferenceSquare)==typeid(*test->getReferenceGeometry())));
	assert(("getNrOfNodes",test->getNrOfNodes()==4));

	for(orig2D[0]=-2.8189;orig2D[0]<3.141;orig2D[0]+=0.1){
		for(orig2D[1]=-2.8189;orig2D[1]<3.141;orig2D[1]+=0.1){
			test->getReferenceToPhysicalMap()->transform(orig2D,compare2D);
			test->referenceToPhysical(orig2D,point2D);
			assert(("referenceToPhysical",fabs(point2D[0]-compare2D[0])<1e-12));
			assert(("referenceToPhysical",fabs(point2D[1]-compare2D[1])<1e-12));

			test->getReferenceToPhysicalMap()->calcJacobian(orig2D,jaccompare);
			test->calcJacobian(orig2D,jac);
			assert(("calcJacobian",fabs(jac[0]-jaccompare[0])<1e-12));
			assert(("calcJacobian",fabs(jac[1]-jaccompare[1])<1e-12));
			assert(("calcJacobian",fabs(jac[2]-jaccompare[2])<1e-12));
			assert(("calcJacobian",fabs(jac[3]-jaccompare[3])<1e-12));
		}
	}
	std::cout<<*test<<std::endl;

	//dim 3

	std::vector<Geometry::PointPhysical> nodes3D;

	Geometry::PointPhysical point3D(3),compare3D(3);
	Geometry::PointReference orig3D(3);
	jac.resize(3,3);
	jaccompare.resize(3,3);

	for(double i=0.;i<10;++i){
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

	assert(("getReferenceToPhysicalMap",typeid(Geometry::MappingToPhysSimplexLinear<3>)==typeid(*test->getReferenceToPhysicalMap())));
	assert(("getPhysicalGeometry",typeid(Geometry::PhysicalTetrahedron)==typeid(*test->getPhysicalGeometry())));
	assert(("getReferenceGeometry",typeid(Geometry::ReferenceTetrahedron)==typeid(*test->getReferenceGeometry())));
	assert(("getNrOfNodes",test->getNrOfNodes()==4));

	for(orig3D[0]=-2.8189;orig3D[0]<3.141;orig3D[0]+=0.2){
		for(orig3D[1]=-2.8189;orig3D[1]<3.141;orig3D[1]+=0.25){
			for(orig3D[2]=-2.8189;orig3D[2]<3.141;orig3D[2]+=0.3){
				test->getReferenceToPhysicalMap()->transform(orig3D,compare3D);
				test->referenceToPhysical(orig3D,point3D);
				assert(("referenceToPhysical",fabs(point3D[0]-compare3D[0])<1e-12));
				assert(("referenceToPhysical",fabs(point3D[1]-compare3D[1])<1e-12));
				assert(("referenceToPhysical",fabs(point3D[2]-compare3D[2])<1e-12));

				test->getReferenceToPhysicalMap()->calcJacobian(orig3D,jaccompare);
				test->calcJacobian(orig3D,jac);
				assert(("calcJacobian",fabs(jac[0]-jaccompare[0])<1e-12));
				assert(("calcJacobian",fabs(jac[1]-jaccompare[1])<1e-12));
				assert(("calcJacobian",fabs(jac[2]-jaccompare[2])<1e-12));
				assert(("calcJacobian",fabs(jac[3]-jaccompare[3])<1e-12));
				assert(("calcJacobian",fabs(jac[4]-jaccompare[4])<1e-12));
				assert(("calcJacobian",fabs(jac[5]-jaccompare[5])<1e-12));
				assert(("calcJacobian",fabs(jac[6]-jaccompare[6])<1e-12));
				assert(("calcJacobian",fabs(jac[7]-jaccompare[7])<1e-12));
				assert(("calcJacobian",fabs(jac[8]-jaccompare[8])<1e-12));
			}
		}
	}
	std::cout<<*test<<std::endl;

	pointIndexes.push_back(12);

	delete test;
	test = new Geometry::ElementGeometry(pointIndexes,nodes3D);

	assert(("getReferenceToPhysicalMap",typeid(Geometry::MappingToPhysPyramid)==typeid(*test->getReferenceToPhysicalMap())));
	assert(("getPhysicalGeometry",typeid(Geometry::PhysicalPyramid)==typeid(*test->getPhysicalGeometry())));
	assert(("getReferenceGeometry",typeid(Geometry::ReferencePyramid)==typeid(*test->getReferenceGeometry())));
	assert(("getNrOfNodes",test->getNrOfNodes()==5));

	for(orig3D[0]=-2.8189;orig3D[0]<3.141;orig3D[0]+=0.2){
		for(orig3D[1]=-2.8189;orig3D[1]<3.141;orig3D[1]+=0.25){
			for(orig3D[2]=-2.8189;orig3D[2]<3.141;orig3D[2]+=0.3){
				test->getReferenceToPhysicalMap()->transform(orig3D,compare3D);
				test->referenceToPhysical(orig3D,point3D);
				assert(("referenceToPhysical",fabs(point3D[0]-compare3D[0])<1e-12));
				assert(("referenceToPhysical",fabs(point3D[1]-compare3D[1])<1e-12));
				assert(("referenceToPhysical",fabs(point3D[2]-compare3D[2])<1e-12));

				test->getReferenceToPhysicalMap()->calcJacobian(orig3D,jaccompare);
				test->calcJacobian(orig3D,jac);
				assert(("calcJacobian",fabs(jac[0]-jaccompare[0])<1e-12));
				assert(("calcJacobian",fabs(jac[1]-jaccompare[1])<1e-12));
				assert(("calcJacobian",fabs(jac[2]-jaccompare[2])<1e-12));
				assert(("calcJacobian",fabs(jac[3]-jaccompare[3])<1e-12));
				assert(("calcJacobian",fabs(jac[4]-jaccompare[4])<1e-12));
				assert(("calcJacobian",fabs(jac[5]-jaccompare[5])<1e-12));
				assert(("calcJacobian",fabs(jac[6]-jaccompare[6])<1e-12));
				assert(("calcJacobian",fabs(jac[7]-jaccompare[7])<1e-12));
				assert(("calcJacobian",fabs(jac[8]-jaccompare[8])<1e-12));
			}
		}
	}
	std::cout<<*test<<std::endl;

	pointIndexes.push_back(13);

	delete test;
	test = new Geometry::ElementGeometry(pointIndexes,nodes3D);

	assert(("getReferenceToPhysicalMap",typeid(Geometry::MappingToPhysTriangularPrism)==typeid(*test->getReferenceToPhysicalMap())));
	assert(("getPhysicalGeometry",typeid(Geometry::PhysicalTriangularPrism)==typeid(*test->getPhysicalGeometry())));
	assert(("getReferenceGeometry",typeid(Geometry::ReferenceTriangularPrism)==typeid(*test->getReferenceGeometry())));
	assert(("getNrOfNodes",test->getNrOfNodes()==6));

	for(orig3D[0]=-2.8189;orig3D[0]<3.141;orig3D[0]+=0.2){
		for(orig3D[1]=-2.8189;orig3D[1]<3.141;orig3D[1]+=0.25){
			for(orig3D[2]=-2.8189;orig3D[2]<3.141;orig3D[2]+=0.3){
				test->getReferenceToPhysicalMap()->transform(orig3D,compare3D);
				test->referenceToPhysical(orig3D,point3D);
				assert(("referenceToPhysical",fabs(point3D[0]-compare3D[0])<1e-12));
				assert(("referenceToPhysical",fabs(point3D[1]-compare3D[1])<1e-12));
				assert(("referenceToPhysical",fabs(point3D[2]-compare3D[2])<1e-12));

				test->getReferenceToPhysicalMap()->calcJacobian(orig3D,jaccompare);
				test->calcJacobian(orig3D,jac);
				assert(("calcJacobian",fabs(jac[0]-jaccompare[0])<1e-12));
				assert(("calcJacobian",fabs(jac[1]-jaccompare[1])<1e-12));
				assert(("calcJacobian",fabs(jac[2]-jaccompare[2])<1e-12));
				assert(("calcJacobian",fabs(jac[3]-jaccompare[3])<1e-12));
				assert(("calcJacobian",fabs(jac[4]-jaccompare[4])<1e-12));
				assert(("calcJacobian",fabs(jac[5]-jaccompare[5])<1e-12));
				assert(("calcJacobian",fabs(jac[6]-jaccompare[6])<1e-12));
				assert(("calcJacobian",fabs(jac[7]-jaccompare[7])<1e-12));
				assert(("calcJacobian",fabs(jac[8]-jaccompare[8])<1e-12));
			}
		}
	}
	std::cout<<*test<<std::endl;

	pointIndexes.push_back(14);
	pointIndexes.push_back(15);

	delete test;
	test = new Geometry::ElementGeometry(pointIndexes,nodes3D);

	assert(("getReferenceToPhysicalMap",typeid(Geometry::MappingToPhysHypercubeLinear<3>)==typeid(*test->getReferenceToPhysicalMap())));
	assert(("getPhysicalGeometry",typeid(Geometry::PhysicalHexahedron)==typeid(*test->getPhysicalGeometry())));
	assert(("getReferenceGeometry",typeid(Geometry::ReferenceCube)==typeid(*test->getReferenceGeometry())));
	assert(("getNrOfNodes",test->getNrOfNodes()==8));

	for(orig3D[0]=-2.8189;orig3D[0]<3.141;orig3D[0]+=0.2){
		for(orig3D[1]=-2.8189;orig3D[1]<3.141;orig3D[1]+=0.25){
			for(orig3D[2]=-2.8189;orig3D[2]<3.141;orig3D[2]+=0.3){
				test->getReferenceToPhysicalMap()->transform(orig3D,compare3D);
				test->referenceToPhysical(orig3D,point3D);
				assert(("referenceToPhysical",fabs(point3D[0]-compare3D[0])<1e-12));
				assert(("referenceToPhysical",fabs(point3D[1]-compare3D[1])<1e-12));
				assert(("referenceToPhysical",fabs(point3D[2]-compare3D[2])<1e-12));

				test->getReferenceToPhysicalMap()->calcJacobian(orig3D,jaccompare);
				test->calcJacobian(orig3D,jac);
				assert(("calcJacobian",fabs(jac[0]-jaccompare[0])<1e-12));
				assert(("calcJacobian",fabs(jac[1]-jaccompare[1])<1e-12));
				assert(("calcJacobian",fabs(jac[2]-jaccompare[2])<1e-12));
				assert(("calcJacobian",fabs(jac[3]-jaccompare[3])<1e-12));
				assert(("calcJacobian",fabs(jac[4]-jaccompare[4])<1e-12));
				assert(("calcJacobian",fabs(jac[5]-jaccompare[5])<1e-12));
				assert(("calcJacobian",fabs(jac[6]-jaccompare[6])<1e-12));
				assert(("calcJacobian",fabs(jac[7]-jaccompare[7])<1e-12));
				assert(("calcJacobian",fabs(jac[8]-jaccompare[8])<1e-12));
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

	for(double i=0.;i<10;++i){
		point4D[0]=1.+i/10.;
		point4D[1]=2.+i/10.;
		point4D[2]=3.+i/10.;
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

	assert(("getReferenceToPhysicalMap",typeid(Geometry::MappingToPhysHypercubeLinear<4>)==typeid(*test->getReferenceToPhysicalMap())));
	assert(("getPhysicalGeometry",typeid(Geometry::PhysicalOctachoron)==typeid(*test->getPhysicalGeometry())));
	assert(("getReferenceGeometry",typeid(Geometry::ReferenceHypercube)==typeid(*test->getReferenceGeometry())));
	assert(("getNrOfNodes",test->getNrOfNodes()==16));

	for(orig4D[0]=-1.5189;orig4D[0]<1.541;orig4D[0]+=0.25){
		for(orig4D[1]=-1.5189;orig4D[1]<1.541;orig4D[1]+=0.25){
			for(orig4D[2]=-1.5189;orig4D[2]<1.541;orig4D[2]+=0.25){
				for(orig4D[3]=-1.5189;orig4D[3]<1.541;orig4D[3]+=0.3){
					test->getReferenceToPhysicalMap()->transform(orig4D,compare4D);
					test->referenceToPhysical(orig4D,point4D);
					assert(("referenceToPhysical",fabs(point4D[0]-compare4D[0])<1e-12));
					assert(("referenceToPhysical",fabs(point4D[1]-compare4D[1])<1e-12));
					assert(("referenceToPhysical",fabs(point4D[2]-compare4D[2])<1e-12));
					assert(("referenceToPhysical",fabs(point4D[3]-compare4D[3])<1e-12));

					test->getReferenceToPhysicalMap()->calcJacobian(orig4D,jaccompare);
					test->calcJacobian(orig4D,jac);
					assert(("calcJacobian",fabs(jac[0]-jaccompare[0])<1e-12));
					assert(("calcJacobian",fabs(jac[1]-jaccompare[1])<1e-12));
					assert(("calcJacobian",fabs(jac[2]-jaccompare[2])<1e-12));
					assert(("calcJacobian",fabs(jac[3]-jaccompare[3])<1e-12));
					assert(("calcJacobian",fabs(jac[4]-jaccompare[4])<1e-12));
					assert(("calcJacobian",fabs(jac[5]-jaccompare[5])<1e-12));
					assert(("calcJacobian",fabs(jac[6]-jaccompare[6])<1e-12));
					assert(("calcJacobian",fabs(jac[7]-jaccompare[7])<1e-12));
					assert(("calcJacobian",fabs(jac[8]-jaccompare[8])<1e-12));
					assert(("calcJacobian",fabs(jac[9]-jaccompare[9])<1e-12));
					assert(("calcJacobian",fabs(jac[10]-jaccompare[10])<1e-12));
					assert(("calcJacobian",fabs(jac[11]-jaccompare[11])<1e-12));
					assert(("calcJacobian",fabs(jac[12]-jaccompare[12])<1e-12));
					assert(("calcJacobian",fabs(jac[13]-jaccompare[13])<1e-12));
					assert(("calcJacobian",fabs(jac[14]-jaccompare[14])<1e-12));
					assert(("calcJacobian",fabs(jac[15]-jaccompare[15])<1e-12));
				}
			}
		}
	}
	std::cout<<*test<<std::endl;
	std::cout<<copy<<std::endl;

	return 0;
}

