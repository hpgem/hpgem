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


//Validates that the mesh contains no gaps or overlaps and that integration works properly by integrating a
//series of (non-)linear functions over the entire domain.

#include "Base/MeshManipulator.hpp"
#include "Base/RectangularMeshDescriptor.hpp"
#include "Integration/ElementIntegral.hpp"
#include "Base/ShortTermStorageElementH1.hpp"
#include "Base/ConfigurationData.hpp"
#include "Base/ElementCacheData.hpp"

#include "Base/CommandLineOptions.hpp"
#include "unordered_set"
#include "cassert"

void testMesh(Base::MeshManipulator* test) {
	class:public Integration::ElementIntegrandBase<LinearAlgebra::NumericalVector>{
		void elementIntegrand(const Base::Element* el, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret){
			ret.resize(1);
			ret[0]=1;
		}
	}one;
	class:public Integration::ElementIntegrandBase<LinearAlgebra::NumericalVector>{
		void elementIntegrand(const Base::Element* el, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret){
			ret.resize(1);
			Geometry::PointPhysical pPhys(p.size());
			ret[0]=0;
			el->referenceToPhysical(p,pPhys);
			for(int i=0;i<p.size();++i){
				ret[0]+=pPhys[i];
			}
		}
	}linear;
	class:public Integration::ElementIntegrandBase<LinearAlgebra::NumericalVector>{
		void elementIntegrand(const Base::Element* el, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret){
			ret.resize(1);
			Geometry::PointPhysical pPhys(p.size());
			ret[0]=1;
			el->referenceToPhysical(p,pPhys);
			for(int i=0;i<p.size();++i){
				ret[0]*=pPhys[i];
			}
		}
	}trilinear;
	Integration::ElementIntegral elIntegral(false);
	elIntegral.setStorageWrapper(new Base::ShortTermStorageElementH1(test->dimension()));
	double total=0;
	LinearAlgebra::NumericalVector result(1);
	for(Base::Element* element:test->getElementsList()){
		elIntegral.integrate(element,&one,result);
		total+=result[0];
	}
	assert(("total mesh volume",fabs(total-1.)<1e-12));
	total=0;
	for(Base::Element* element:test->getElementsList()){
		elIntegral.integrate(element,&linear,result);
		total+=result[0];
	}
	assert(("linear function",fabs(total-.5*test->dimension())<1e-12));
	total=0;
	for(Base::Element* element:test->getElementsList()){
		elIntegral.integrate(element,&trilinear,result);
		total+=result[0];
	}
	assert(("trilinear function",fabs(total-pow(0.5,test->dimension()))<1e-12));
}


int main(int argc, char** argv){
    Base::parse_options(argc,argv);
	// dim 1
	Base::RectangularMeshDescriptor description1D(1),description2D(2),description3D(3);
	description1D.bottomLeft_[0]=0;
	description2D.bottomLeft_[0]=0;
	description2D.bottomLeft_[1]=0;
	description3D.bottomLeft_[0]=0;
	description3D.bottomLeft_[1]=0;
	description3D.bottomLeft_[2]=0;
	description1D.topRight_[0]=1;
	description2D.topRight_[0]=1;
	description2D.topRight_[1]=1;
	description3D.topRight_[0]=1;
	description3D.topRight_[1]=1;
	description3D.topRight_[2]=1;
	description1D.boundaryConditions_[0]=Base::RectangularMeshDescriptor::SOLID_WALL;
	description2D.boundaryConditions_[0]=Base::RectangularMeshDescriptor::SOLID_WALL;
	description2D.boundaryConditions_[1]=Base::RectangularMeshDescriptor::SOLID_WALL;
	description3D.boundaryConditions_[0]=Base::RectangularMeshDescriptor::SOLID_WALL;
	description3D.boundaryConditions_[1]=Base::RectangularMeshDescriptor::SOLID_WALL;
	description3D.boundaryConditions_[2]=Base::RectangularMeshDescriptor::SOLID_WALL;

	description1D.numElementsInDIM_[0]=2;

	Base::MeshManipulator *test = new Base::MeshManipulator(new Base::ConfigurationData(1,1,2,0),false,false,false,2,0);
	test->createTriangularMesh(description1D.bottomLeft_,description1D.topRight_,description1D.numElementsInDIM_);

	testMesh(test);

	delete test;
	test = new Base::MeshManipulator(new Base::ConfigurationData(1,1,2,0),false,false,false,2,0);
	test->createRectangularMesh(description1D.bottomLeft_,description1D.topRight_,description1D.numElementsInDIM_);
	testMesh(test);

	delete test;
	description1D.numElementsInDIM_[0]=3;

	test = new Base::MeshManipulator(new Base::ConfigurationData(1,1,2,0),false,false,false,2,0);
	test->createTriangularMesh(description1D.bottomLeft_,description1D.topRight_,description1D.numElementsInDIM_);

	testMesh(test);

	delete test;
	test = new Base::MeshManipulator(new Base::ConfigurationData(1,1,2,0),false,false,false,2,0);
	test->createRectangularMesh(description1D.bottomLeft_,description1D.topRight_,description1D.numElementsInDIM_);
	testMesh(test);

	// dim 2

	delete test;
	description2D.numElementsInDIM_[0]=2;
	description2D.numElementsInDIM_[1]=3;

	test = new Base::MeshManipulator(new Base::ConfigurationData(2,1,2,0),false,false,false,2,0);
	test->createTriangularMesh(description2D.bottomLeft_,description2D.topRight_,description2D.numElementsInDIM_);

	testMesh(test);

	delete test;
	test = new Base::MeshManipulator(new Base::ConfigurationData(2,1,2,0),false,false,false,2,0);
	test->createRectangularMesh(description2D.bottomLeft_,description2D.topRight_,description2D.numElementsInDIM_);
	testMesh(test);

	delete test;
	description2D.numElementsInDIM_[0]=3;
	description2D.numElementsInDIM_[1]=2;

	test = new Base::MeshManipulator(new Base::ConfigurationData(2,1,2,0),false,false,false,2,0);
	test->createTriangularMesh(description2D.bottomLeft_,description2D.topRight_,description2D.numElementsInDIM_);

	testMesh(test);

	delete test;
	test = new Base::MeshManipulator(new Base::ConfigurationData(2,1,2,0),false,false,false,2,0);
	test->createRectangularMesh(description2D.bottomLeft_,description2D.topRight_,description2D.numElementsInDIM_);
	testMesh(test);

	// dim 3

	delete test;
	description3D.numElementsInDIM_[0]=2;
	description3D.numElementsInDIM_[1]=2;
	description3D.numElementsInDIM_[2]=3;

	test = new Base::MeshManipulator(new Base::ConfigurationData(3,1,2,0),false,false,false,2,0);
	test->createTriangularMesh(description3D.bottomLeft_,description3D.topRight_,description3D.numElementsInDIM_);

	testMesh(test);

	delete test;
	test = new Base::MeshManipulator(new Base::ConfigurationData(3,1,2,0),false,false,false,2,0);
	test->createRectangularMesh(description3D.bottomLeft_,description3D.topRight_,description3D.numElementsInDIM_);
	testMesh(test);

	delete test;
	description3D.numElementsInDIM_[0]=2;
	description3D.numElementsInDIM_[1]=3;
	description3D.numElementsInDIM_[2]=2;

	test = new Base::MeshManipulator(new Base::ConfigurationData(3,1,2,0),false,false,false,2,0);
	test->createTriangularMesh(description3D.bottomLeft_,description3D.topRight_,description3D.numElementsInDIM_);

	testMesh(test);

	delete test;
	test = new Base::MeshManipulator(new Base::ConfigurationData(3,1,2,0),false,false,false,2,0);
	test->createRectangularMesh(description3D.bottomLeft_,description3D.topRight_,description3D.numElementsInDIM_);
	testMesh(test);

	delete test;
	description3D.numElementsInDIM_[0]=3;
	description3D.numElementsInDIM_[1]=2;
	description3D.numElementsInDIM_[2]=2;

	test = new Base::MeshManipulator(new Base::ConfigurationData(3,1,2,0),false,false,false,2,0);
	test->createTriangularMesh(description3D.bottomLeft_,description3D.topRight_,description3D.numElementsInDIM_);

	testMesh(test);

	delete test;
	test = new Base::MeshManipulator(new Base::ConfigurationData(3,1,2,0),false,false,false,2,0);
	test->createRectangularMesh(description3D.bottomLeft_,description3D.topRight_,description3D.numElementsInDIM_);
	testMesh(test);
}



