/*
 * 030Integration_SelfTest.cpp
 *
 *  Created on: Apr 17, 2014
 *      Author: brinkf
 */

//Validates that the mesh contains no gaps or overlaps and that integration works properly by integrating a
//series of (non-)linear functions over the entire domain.

#include "Base/MeshManipulator.hpp"
#include "Base/RectangularMeshDescriptor.hpp"
#include "Integration/ElementIntegral.hpp"

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

int main(){

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

	Base::MeshManipulator *test = new Base::MeshManipulator(new Base::ConfigurationData(1,1,2,1),false,false,false,2,0);
	test->createTriangularMesh(description1D.bottomLeft_,description1D.topRight_,description1D.numElementsInDIM_);

	testMesh(test);

	delete test;
	test = new Base::MeshManipulator(new Base::ConfigurationData(1,1,2,1),false,false,false,2,0);
	test->createRectangularMesh(description1D.bottomLeft_,description1D.topRight_,description1D.numElementsInDIM_);
	testMesh(test);

	delete test;
	description1D.numElementsInDIM_[0]=3;

	test = new Base::MeshManipulator(new Base::ConfigurationData(1,1,2,1),false,false,false,2,0);
	test->createTriangularMesh(description1D.bottomLeft_,description1D.topRight_,description1D.numElementsInDIM_);

	testMesh(test);

	delete test;
	test = new Base::MeshManipulator(new Base::ConfigurationData(1,1,2,1),false,false,false,2,0);
	test->createRectangularMesh(description1D.bottomLeft_,description1D.topRight_,description1D.numElementsInDIM_);
	testMesh(test);

	// dim 2

	delete test;
	description2D.numElementsInDIM_[0]=2;
	description2D.numElementsInDIM_[1]=3;

	test = new Base::MeshManipulator(new Base::ConfigurationData(2,1,2,1),false,false,false,2,0);
	test->createTriangularMesh(description2D.bottomLeft_,description2D.topRight_,description2D.numElementsInDIM_);

	testMesh(test);

	delete test;
	test = new Base::MeshManipulator(new Base::ConfigurationData(2,1,2,1),false,false,false,2,0);
	test->createRectangularMesh(description2D.bottomLeft_,description2D.topRight_,description2D.numElementsInDIM_);
	testMesh(test);

	delete test;
	description2D.numElementsInDIM_[0]=3;
	description2D.numElementsInDIM_[1]=2;

	test = new Base::MeshManipulator(new Base::ConfigurationData(2,1,2,1),false,false,false,2,0);
	test->createTriangularMesh(description2D.bottomLeft_,description2D.topRight_,description2D.numElementsInDIM_);

	testMesh(test);

	delete test;
	test = new Base::MeshManipulator(new Base::ConfigurationData(2,1,2,1),false,false,false,2,0);
	test->createRectangularMesh(description2D.bottomLeft_,description2D.topRight_,description2D.numElementsInDIM_);
	testMesh(test);

	// dim 3

	delete test;
	description3D.numElementsInDIM_[0]=2;
	description3D.numElementsInDIM_[1]=2;
	description3D.numElementsInDIM_[2]=3;

	test = new Base::MeshManipulator(new Base::ConfigurationData(3,1,2,1),false,false,false,2,0);
	test->createTriangularMesh(description3D.bottomLeft_,description3D.topRight_,description3D.numElementsInDIM_);

	testMesh(test);

	delete test;
	test = new Base::MeshManipulator(new Base::ConfigurationData(3,1,2,1),false,false,false,2,0);
	test->createRectangularMesh(description3D.bottomLeft_,description3D.topRight_,description3D.numElementsInDIM_);
	testMesh(test);

	delete test;
	description3D.numElementsInDIM_[0]=2;
	description3D.numElementsInDIM_[1]=3;
	description3D.numElementsInDIM_[2]=2;

	test = new Base::MeshManipulator(new Base::ConfigurationData(3,1,2,1),false,false,false,2,0);
	test->createTriangularMesh(description3D.bottomLeft_,description3D.topRight_,description3D.numElementsInDIM_);

	testMesh(test);

	delete test;
	test = new Base::MeshManipulator(new Base::ConfigurationData(3,1,2,1),false,false,false,2,0);
	test->createRectangularMesh(description3D.bottomLeft_,description3D.topRight_,description3D.numElementsInDIM_);
	testMesh(test);

	delete test;
	description3D.numElementsInDIM_[0]=3;
	description3D.numElementsInDIM_[1]=2;
	description3D.numElementsInDIM_[2]=2;

	test = new Base::MeshManipulator(new Base::ConfigurationData(3,1,2,1),false,false,false,2,0);
	test->createTriangularMesh(description3D.bottomLeft_,description3D.topRight_,description3D.numElementsInDIM_);

	testMesh(test);

	delete test;
	test = new Base::MeshManipulator(new Base::ConfigurationData(3,1,2,1),false,false,false,2,0);
	test->createRectangularMesh(description3D.bottomLeft_,description3D.topRight_,description3D.numElementsInDIM_);
	testMesh(test);
}



