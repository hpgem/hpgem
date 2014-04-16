/*
 * 040Element_UnitTest.cpp
 *
 *  Created on: Apr 15, 2014
 *      Author: brinkf
 */

//naming convention: <Digit><ClassName>_UnitTest.cpp where <Digit> is a number that will make sure
//the unit tests are ordered such that the first failing unit test indicate the culprit class and
//other 'unit' tests may assume correct execution of all prior unit tests

//A lot of the functionality of Element, Face and Edge only makes sense if they are embedded in a larger mesh that is completely set up
//Setting up a testing case to validate this functionality is as hard as generating a mesh. As such this functionality will be tested as
//part of the mesh generation self tests

#include "Base/Element.hpp"
#include "cassert"

#include "Base/AssembleBasisFunctionSet.hpp"
#include "Integration/QuadratureRules/GaussQuadratureRulesForCube.hpp"

int main() {

	std::vector<unsigned int> pointIndexes;
	std::vector<Geometry::PointPhysical> nodes;

	Geometry::PointPhysical point(3);

	pointIndexes.push_back(4);
	pointIndexes.push_back(7);
	pointIndexes.push_back(10);
	pointIndexes.push_back(11);
	pointIndexes.push_back(12);
	pointIndexes.push_back(13);
	pointIndexes.push_back(14);
	pointIndexes.push_back(15);

	for(double i=0.;i<10;++i){
		point[0]=1.+i/10.;
		point[1]=2.+i/10.;
		point[2]=3.+i/10.;
		nodes.push_back(point);
	}

	point[0]=3.5;
	point[1]=4.6;
	point[2]=5.4;
	nodes.push_back(point);
	point[0]=6.7;
	point[1]=2.8;
	point[2]=5.7;
	nodes.push_back(point);
	point[0]=1.4;
	point[1]=2.4;
	point[2]=5.4;
	nodes.push_back(point);
	point[0]=1.7;
	point[1]=2.7;
	point[2]=5.7;
	nodes.push_back(point);
	point[0]=3.5;
	point[1]=4.6;
	point[2]=7.4;
	nodes.push_back(point);
	point[0]=6.7;
	point[1]=2.8;
	point[2]=7.7;
	nodes.push_back(point);

	Base::BasisFunctionSet basisFunctions(3);

	Base::AssembleBasisFunctionSet_3D_Ord3_A1(basisFunctions);

	std::vector<const Base::BasisFunctionSet*> vectorOfFunctions(1,&basisFunctions);

	Base::Element test(pointIndexes,&vectorOfFunctions,nodes,3,14,basisFunctions.size(),18),copy(test);

	assert(("quadrature rule",test.getGaussQuadratureRule()!=NULL));

	test.setQuadratureRulesWithOrder(8);
	copy.setGaussQuadratureRule(&QuadratureRules::Cn3_3_4::Instance());

	assert(("setQuadratureRule",test.getGaussQuadratureRule()->order()>=8));
	assert(("setQuadratureRule",typeid(*copy.getGaussQuadratureRule())==typeid(QuadratureRules::Cn3_3_4)));

	//check set*BasisFunctionSet without breaking preconditions...

	Geometry::PointReference point3D(3);
	for(int i=0;i<basisFunctions.size();++i){
		for(point[0]=-1.5;point[0]<1.51;point[0]+=0.1){
			for(point[1]=-1.5;point[1]<1.51;point[1]+=0.1){
				for(point[2]=-1.5;point[2]<1.51;point[2]+=0.1){
					assert(("basisFunctions",test.basisFunction(i,point3D)==basisFunctions[i]->eval(point3D)));
					assert(("basisFunctions",test.basisFunctionDeriv(i,0,point3D)==basisFunctions[i]->evalDeriv0(point3D)));
					assert(("basisFunctions",test.basisFunctionDeriv(i,1,point3D)==basisFunctions[i]->evalDeriv1(point3D)));
					assert(("basisFunctions",test.basisFunctionDeriv(i,2,point3D)==basisFunctions[i]->evalDeriv2(point3D)));
				}
			}
		}
	}

	//check setFace and setEdge without breaking preconditions...

	cout<<test<<endl;

	return 0;
}
