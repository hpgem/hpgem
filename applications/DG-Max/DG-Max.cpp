/*
 This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
 http://hpgem.org/about-the-code/team
 
 This code is distributed using BSD 3-Clause License. A copy of which can found below.
 
 
 Copyright (c) 2014, Univesity of Twenete
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

//this file should contain all relevant information about how the integrands look like and what problem is solved



#define _USE_MATH_DEFINES
#include <cstdlib>
#include "Output/TecplotDiscontinuousSolutionWriter.hpp"
#include "Output/TecplotPhysicalGeometryIterator.hpp"
#include "LinearAlgebra/NumericalVector.hpp"
#include "BasisFunctionCollection_Curl.hpp"
#include <iostream>
#include "fillMatrices.hpp"
#include "Base/L2Norm.hpp"
#include "BaseExtended.hpp"
#include "math.h"
#include <ctime>
#include "ElementInfos.hpp"
#include "Integration/ElementIntegrandBase.hpp"
#include "Integration/FaceIntegrandBase.hpp"
#include "Output/TecplotSingleElementWriter.hpp"

typedef Base::MeshManipulator BaseMeshManipulatorT;
typedef Base::threeDBasisFunction basisFunctionT;
//typedef std::list<Base::Face<3> >::iterator FaceIteratorT;
	typedef Base::MeshManipulator::FaceIterator               FaceIterator;


/**
 * This class should provide problem specific information about the maxwell equations. 
 */
class DGMax : public hpGemUIExtentions
{
private:

        typedef Geometry::PointPhysical                           PointPhysicalT;
	typedef Base::Element                                     ElementT;
        typedef Base::Face                                         FaceT;
	
	typedef Base::MeshManipulator::FaceIterator               FaceIterator;

public:

    DGMax(int argc,char** argv, MaxwellData* globalConfig, Base::ConfigurationData* elementConfig,matrixFiller* fill):hpGemUIExtentions(argc,argv,globalConfig,elementConfig,fill){}

    /**
     * set up the mesh and complete initialisation of the global data and the configuration data
     */
    bool initialise(){
        int n=getData()->NumberOfIntervals_;
        Geometry::PointPhysical bottomLeft(3), topRight(3);
        std::vector<unsigned int> numElementsOneD(3);
        bottomLeft[0] = 0;
        bottomLeft[1] = 0;
        bottomLeft[2] = 0;
        topRight[0] = 1;
        topRight[1] = 1;
        topRight[2] = 1;
        numElementsOneD[0] = n;
        numElementsOneD[1] = n;
        numElementsOneD[2] = n;
	
	BaseMeshManipulatorT* mesh = new MyMeshManipulator(getConfigData(),getData()->PolynomialOrder_,true,true,true);
	
	//mesh->readCentaurMesh("Cylinder3.hyb");
	//mesh->readCentaurMesh("input_basic2.hyb");
	//mesh->readCentaurMesh("FicheraGlobRefine3.hyb");
	mesh->createTriangularMesh(bottomLeft,topRight,numElementsOneD);
	const_cast<MaxwellData*>(getData())->numberOfUnknowns_=(*mesh->elementColBegin())->getNrOfBasisFunctions();
	setConfigData();
	addMesh(mesh);
        //addMesh("Triangular",bottomLeft, topRight, numElementsOneD);
	
	//actually store the full number of unknowns
	const_cast<MaxwellData*>(getData())->numberOfUnknowns_*=mesh->getElementsList().size();
	for(Base::MeshManipulator::ElementIterator it=mesh->elementColBegin();it!=mesh->elementColEnd();++it){
	    (*it)->setUserData(new ElementInfos(**it));
	}
	
	return true;
    }

    /**
     * Computes element contributions to the mass matrix i.e. phi_i * phi_j
     * returns the contibutions at this gauss point to the entire element matrix in one go
     */
    void elementIntegrand(const ElementT* element, const PointElementReferenceT& p, LinearAlgebra::Matrix& ret){
        ret.resize(element->getNrOfBasisFunctions(),element->getNrOfBasisFunctions());
        //cout<<"\nIn the element integrand for the mass matrix for element id: "<<element->getID();
        ElementInfos* info = static_cast<ElementInfos*> (const_cast<ElementT*>(element)->getUserData());
		NumericalVector phi_i(3),phi_j(3);
		std::vector<NumericalVector> functionValues;
		info->makeFunctionValuesVector(element,p,functionValues);
		for(int i=0;i<element->getNrOfBasisFunctions();++i){
			phi_i=functionValues[i];
				for(int j=0;j<element->getNrOfBasisFunctions();++j){
			phi_j=functionValues[j];
				ret(i,j)=(phi_i[0]*phi_j[0]+phi_i[1]*phi_j[1]+phi_i[2]*phi_j[2])*info->epsilon_;
				ret(j,i)=ret(i,j);
			}
		}
    }

    /**
     * Computes the bits of the face contributions that are common to both the IP method and the brezzi formulation
     * i.e. -.5( (nabla x phi_i) * phi_j + phi_i * (nabla x phi_j) )
     * returns the contibutions at this gauss point to the entire face matrix in one go
     */
    void faceIntegrand(const FaceT* face, const NumericalVector& normal, const PointFaceReferenceT& p, LinearAlgebra::Matrix& ret){
        //cout<<"\nIn the face integrand for the stiffness matrix for element id: "<<face->getPtrElementLeft()->getID();
		ElementT* right;
		ElementT* left=const_cast<ElementT*>(face->getPtrElementLeft());
		ElementInfos* leftInfo = static_cast<ElementInfos*> (left->getUserData());
		ElementInfos* rightInfo;
		PointElementReferenceT pLeft(3),pRight(3);
		face->mapRefFaceToRefElemL(p,pLeft);
		std::vector<NumericalVector> leftValues,leftCurls,rightValues,rightCurls;
		leftInfo->makeFunctionValuesVector(left,pLeft,leftValues);
		leftInfo->makeFunctionCurlsVector(left,pLeft,leftCurls);
		NumericalVector normedNormal(3);
		normedNormal[0] = (normal*(1/Base::L2Norm(normal)))[0];
		normedNormal[1] = (normal*(1/Base::L2Norm(normal)))[1];
		normedNormal[2] = (normal*(1/Base::L2Norm(normal)))[2];
		int leftSize=left->getNrOfBasisFunctions();
		int dimension=left->getNrOfBasisFunctions();
		if(face->isInternal()){
			//cout<<" and element id: "<<face->getPtrElementRight()->getID();
			right=const_cast<ElementT*>(face->getPtrElementRight());
			face->mapRefFaceToRefElemR(p,pRight);
			rightInfo = static_cast<ElementInfos*> (right->getUserData());
			dimension+=right->getNrOfBasisFunctions();
			rightInfo->makeFunctionValuesVector(right,pRight,rightValues);
			rightInfo->makeFunctionCurlsVector(right,pRight,rightCurls);
		}
		ret.resize(dimension,dimension);
		NumericalVector phi_i(3),phi_j(3),phi_i_curl(3),phi_j_curl(3),dummy(3);
		for(int i=0;i<dimension;++i){
			if(i<leftSize){
				dummy=leftValues[i];
			phi_i_curl=leftCurls[i];
			}else{
				dummy=rightValues[i-leftSize];
			phi_i_curl=rightCurls[i-leftSize];
			dummy*=-1;
			}
			OuterProduct(normedNormal,dummy,phi_i);
			for(int j=i;j<dimension;++j){
			if(j<leftSize){
				dummy=leftValues[j];
				phi_j_curl=leftCurls[j];
			}else{
				dummy=rightValues[j-leftSize];
				phi_j_curl=rightCurls[j-leftSize];
				dummy*=-1;
			}
			OuterProduct(normedNormal,dummy,phi_j);
			ret(i,j)=-(face->isInternal()?0.5:1.)*(phi_i[0]*phi_j_curl[0]+phi_i[1]*phi_j_curl[1]+phi_i[2]*phi_j_curl[2]+
											  phi_j[0]*phi_i_curl[0]+phi_j[1]*phi_i_curl[1]+phi_j[2]*phi_i_curl[2]);
				ret(j,i)=ret(i,j);
			}
		}
    }
    
    /**
     * this is where you specify an initial condition
     */
    void initialConditions(const PointPhysicalT& p, NumericalVector& ret){
        exactSolution(p,0,ret);
    }
    
    /**
     * this is where you specify an initial time derivative of the solution
     *
     * support temporarily dropped
     */
    void initialConditionsDeriv(const PointPhysicalT& p, NumericalVector& ret){
        ret[0]=0;
	ret[1]=0;
	ret[2]=0;
    }
    
    /**
     * this is where you specify the spatial part of the source Term
     * assumes that the source term can be split is a spatial part and a time part
     */
    void sourceTerm(const PointPhysicalT& p, NumericalVector& ret){
//         ret[0]=0;
// 	ret[1]=0;
// 	ret[2]=0;
        exactSolution(p,0,ret);
// 	ret*=-1;
	ret*=M_PI*M_PI*8-1;
    }
    
    /**
     * this is where you specify the time part of the source Term
     * assumes that the source term can be split is a spatial part and a time part
     */
    double sourceTermTime(const double t){
        return 1.;
    }
    
    /**
     * this is where you choose the solution of your problem
     * this will only have an effect on the accuracy of your error estimates
     * as a temporary solution remember to also update the exact solution in fillMatrices.cpp
     */
    void exactSolution(const PointPhysicalT& p, const double t, NumericalVector &ret){
         ret[0]=sin(M_PI*2*p[1])*sin(M_PI*2*p[2]);
         ret[1]=sin(M_PI*2*p[2])*sin(M_PI*2*p[0]);
         ret[2]=sin(M_PI*2*p[0])*sin(M_PI*2*p[1]);
 	ret*=cos(sqrt(2)*2*M_PI*t);
//        ret[0]=sin(M_PI*p[1])*sin(M_PI*p[2]);
//        ret[1]=sin(M_PI*p[2])*sin(M_PI*p[0]);
//        ret[2]=sin(M_PI*p[0])*sin(M_PI*p[1]);
//	ret*=cos(sqrt(2)*M_PI*t);
//            ret[0]=p[0]*(1-p[0]);
//            ret[1]=0;
// 	   ret[2]=0;
    }
    
    /**
     * this is where you choose the curl of the solution of your problem
     * this will only have an effect on the accuracy of your error estimates
     */  
    void exactSolutionCurl(const PointPhysicalT& p, const double t, NumericalVector &ret){
         ret[0]=sin(M_PI*2*p[0])*(cos(M_PI*2*p[1])-cos(M_PI*2*p[2]));
         ret[1]=sin(M_PI*2*p[1])*(cos(M_PI*2*p[2])-cos(M_PI*2*p[0]));
         ret[2]=sin(M_PI*2*p[2])*(cos(M_PI*2*p[0])-cos(M_PI*2*p[1]));
 	ret*=cos(sqrt(2)*2*M_PI*t)*2*M_PI;
//        ret[0]=sin(M_PI*p[0])*(cos(M_PI*p[1])-cos(M_PI*p[2]));
//        ret[1]=sin(M_PI*p[1])*(cos(M_PI*p[2])-cos(M_PI*p[0]));
//        ret[2]=sin(M_PI*p[2])*(cos(M_PI*p[0])-cos(M_PI*p[1]));
//	ret*=cos(sqrt(2)*M_PI*t)*M_PI;
//          ret[0]=0;ret[1]=0;ret[2]=0;
    }
};

/**
 * main routine. Suggested usage: make a problem, solve it, tell the user your results.
 */
int main(int argc,char** argv){
  
    //set up timings
    time_t start,end,initialised,solved;
    time(&start);
    
    //read input data
    int elements = 1;
    if(argc>1){
        elements=std::atoi(argv[1]);
	cout<<"using "<<elements*elements*elements*5<<" elements"<<endl;
    }
    int order = 1;
    if(argc>2){
        order=std::atoi(argv[2]);
	cout<<"using polynomial order: "<<order<<endl;
    }else{
        cout<<"usage:./Maxwell.out <elements> <order> [<petsc-args>]";
	exit(1);
    }
    //set up problem and decide flux type
    DGMax problem(argc-2,&argv[2],new MaxwellData(elements,order),new Base::ConfigurationData(3,1,0,2),new matrixFillerBR(125));
    try{
        problem.initialise();
	time(&initialised);
	
	//choose what problem to solve
        problem.solveHarmonic();
	time(&solved);
	char filename[]="output.dat";
	problem.makeOutput(filename);
	time(&end);
	
	//display timing data
	cout<<"Initialisation took "<<difftime(initialised,start)<<" seconds."<<endl;
	cout<<"Solving the problem took "<<difftime(solved,initialised)<<" seconds."<<endl;
	cout<<"The rest took "<<difftime(end,solved)<<" seconds."<<endl;
	
    }catch(const char* message){
        std::cout << message;
    }
    return 0;
}
