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

#define hpGEM_INCLUDE_PETSC_SUPPORT

#include "Base/HpgemUISimplified.hpp"
#include "petscksp.h"
#include "Utilities/GlobalMatrix.hpp"
#include "Utilities/GlobalVector.hpp"
#include "Base/Norm2.hpp"

//If this test ever breaks it is not a bad thing per se.
//If the results are still readable by tecplot, and you are convinced that your changes improved the code,
//you should update the data file to reflect the updated result. Always confer with other developers if you do this.

class Laplace : public Base::HpgemUISimplified,public Output::TecplotSingleElementWriter{

	class:public Integration::ElementIntegrandBase<LinearAlgebra::NumericalVector>{

		//has an almost exact copy near the element integrand for the source term
	    double sourceTerm(const PointPhysicalT& p){
	    	double ret=sin(2*M_PI*p[0]);
	    	if(p.size()>1){
	    		ret*=cos(2*M_PI*p[1])/2.;
	    	}
	    	if(p.size()>2){
	    		ret*=cos(2*M_PI*p[2])*2.;
	    	}
	    	return ret;
	    }

		void elementIntegrand(const Base::Element* el, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret){
			ret.resize(1);
			Geometry::PointPhysical pPhys(p.size());
			el->referenceToPhysical(p,pPhys);
			el->getSolution(0,p,ret);
			ret[0]-=sourceTerm(pPhys);
			ret[0]*=ret[0];
		}
	}errors;

public:
	int n_,p_,DIM_;
	Base::MeshType type_;
	double penaltyParameter_;
	Laplace(int n, int p,int DIM, Base::MeshType type):
		HpgemUISimplified(DIM,p),n_(n),p_(p),DIM_(DIM),type_(type),penaltyParameter_(3*n_*p_*(p_+DIM_-1)+1){ }

	bool initialise(){
    	RectangularMeshDescriptorT description(DIM_);
    	for(int i=0;i<DIM_;++i){
    		description.bottomLeft_[i]=0;
    		description.topRight_[i]=1;
    		description.numElementsInDIM_[i]=n_;
    		description.boundaryConditions_[i]=RectangularMeshDescriptorT::SOLID_WALL;
    	}
    	addMesh(description,type_,1,1,1,1);
		return true;
	}

	//stiffness matrix
	void elementIntegrand(const Base::Element* el, const Geometry::PointReference& p, LinearAlgebra::Matrix& ret){
    	int n=el->getNrOfBasisFunctions();
    	LinearAlgebra::NumericalVector phiDerivI(DIM_),phiDerivJ(DIM_);
    	ret.resize(n,n);
    	for(int i=0;i<n;++i){
			el->basisFunctionDeriv(i,p,phiDerivI);
    		for(int j=0;j<n;++j){
    			el->basisFunctionDeriv(j,p,phiDerivJ);
    			ret(j,i)=phiDerivI*phiDerivJ;
    		}
    	}
	}

	//has an exact copy in the error class
    double sourceTerm(const PointPhysicalT& p){
    	double ret=sin(2*M_PI*p[0])*(4*M_PI*M_PI);
    	if(DIM_>1){
    		ret*=cos(2*M_PI*p[1]);
    	}
    	if(DIM_>2){
    		ret*=cos(2*M_PI*p[2])*3;
    	}
    	return ret;
    }

	//RHS
	void elementIntegrand(const Base::Element* el, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret){
    	PointPhysicalT pPhys(DIM_);
    	el->referenceToPhysical(p,pPhys);
    	ret.resize(el->getNrOfBasisFunctions());
    	for(int i=0;i<el->getNrOfBasisFunctions();++i){
    		ret[i]=el->basisFunction(i,p)*sourceTerm(pPhys);
    	}
	}

	//face discontinuities
	void faceIntegrand(const Base::Face* fa, const LinearAlgebra::NumericalVector& normal, const Geometry::PointReference& p, LinearAlgebra::Matrix& ret){
    	int n=fa->getNrOfBasisFunctions();
    	ret.resize(n,n);
    	LinearAlgebra::NumericalVector phiNormalI(DIM_),phiNormalJ(DIM_),phiDerivI(DIM_),phiDerivJ(DIM_);
		PointPhysicalT pPhys(DIM_);
		fa->referenceToPhysical(p,pPhys);
    	for(int i=0;i<n;++i){
			fa->basisFunctionNormal(i,normal,p,phiNormalI);
			fa->basisFunctionDeriv(i,p,phiDerivI);
    		for(int j=0;j<n;++j){
    			fa->basisFunctionNormal(j,normal,p,phiNormalJ);
    			fa->basisFunctionDeriv(j,p,phiDerivJ);
    			if(fa->isInternal()){
    				ret(j,i)=-(phiNormalI*phiDerivJ+phiNormalJ*phiDerivI)/2+
    							penaltyParameter_*phiNormalI*phiNormalJ;
    			}else if(fabs(pPhys[0])<1e-12||fabs(pPhys[0]-1.)<1e-12){//Dirichlet
    				ret(j,i)=-(phiNormalI*phiDerivJ+phiNormalJ*phiDerivI)+
    							penaltyParameter_*phiNormalI*phiNormalJ*2;
    			}else{//homogeneous Neumann
    				ret(j,i)=0;
    			}
    		}
    	}
	}

	//boundary conditions
	void faceIntegrand(const Base::Face* fa, const LinearAlgebra::NumericalVector& normal, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret){
		int n=fa->getNrOfBasisFunctions();
		ret.resize(n);
		PointPhysicalT pPhys(DIM_);
		fa->referenceToPhysical(p,pPhys);
		if(fabs(pPhys[0])<1e-9||fabs(pPhys[0]-1)<1e-9){//Dirichlet
			LinearAlgebra::NumericalVector phiDeriv(DIM_);
			for(int i=0;i<n;++i){
				fa->basisFunctionDeriv(i,p,phiDeriv);
				ret[i]=(-normal*phiDeriv/Utilities::norm2(normal)+penaltyParameter_*fa->basisFunction(i,p))*0;
			}
		}else{
			for(int i=0;i<n;++i){
				ret[i]=0;
			}
		}
	}

	double initialConditions(const Geometry::PointPhysical& p){
		return 0;
	}

    void writeToTecplotFile(const ElementT* element,const  PointReferenceT& p, ostream& out){
    	LinearAlgebra::NumericalVector value(1);
    	element->getSolution(0,p,value);
    	out<<value[0];
    }

	double solve(){
    	doAllElementIntegration();
    	doAllFaceIntegration();
    	Utilities::GlobalPetscMatrix A(HpgemUI::meshes_[0],0,0);
    	Utilities::GlobalPetscVector b(HpgemUI::meshes_[0],0,0),x(HpgemUI::meshes_[0]);

    	b.assemble();
    	VecSet(x,0);

    	KSP ksp;
    	KSPCreate(MPI_COMM_WORLD,&ksp);
    	KSPSetTolerances(ksp,1e-12,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
    	KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN);
    	KSPSetFromOptions(ksp);
    	KSPSolve(ksp,b,x);

    	x.writeTimeLevelData(0);

    	std::ofstream outFile("030TecplotOutput_SelfTest_output.dat");
		Output::TecplotDiscontinuousSolutionWriter writeFunc(outFile,"test","01","value");
		writeFunc.write(meshes_[0],"monomial solution",false,this);
	}
};

int main(){
	PetscInitializeNoArguments();
	Laplace test8(8,5,2,Base::TRIANGULAR);
	test8.initialise();
	test8.solve();
	//actual test is done by comparing output files
	PetscFinalize();
	return 0;
}


