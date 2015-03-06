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


//#define hpGEM_INCLUDE_PETSC_SUPPORT//temporarily activating this definition makes development easier on some IDEs

#include "Base/HpgemUISimplified.hpp"
#include "petscksp.h"
#include "Utilities/GlobalMatrix.hpp"
#include "Utilities/GlobalVector.hpp"
#include "Base/Norm2.hpp"
#include "Base/ShortTermStorageElementH1.hpp"
#include "Base/FaceCacheData.hpp"
#include "Base/ElementCacheData.hpp"
#include "Geometry/PointPhysical.hpp"
#include "Base/RectangularMeshDescriptor.hpp"
#include "Integration/ElementIntegral.hpp"
#include <cmath>
#include "Logger.h"

//If this test ever breaks it is not a bad thing per se. However, once this breaks a thorough convergence analysis needs to be done.
//If the results still show the theoretically optimal order of convergence, and you are convinced that your changes improved the code,
//you should change the numbers in this test to reflect the updated result. Always confer with other developers if you do this.

class Laplace : public Base::HpgemUISimplified{

	class:public Integration::ElementIntegrandBase<LinearAlgebra::NumericalVector>{

		//has an almost exact copy near the element integrand for the source term
	    double sourceTerm(const PointPhysicalT& p){
	    	double ret=std::sin(2*M_PI*p[0]);
	    	if(p.size()>1){
	    		ret*=std::cos(2*M_PI*p[1])/2.;
	    	}
	    	if(p.size()>2){
	    		ret*=std::cos(2*M_PI*p[2])*2.;
	    	}
	    	return ret;
	    }

		void elementIntegrand(const Base::Element* el, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret){
			ret.resize(1);
			Geometry::PointPhysical pPhys = el->referenceToPhysical(p);
			ret = el->getSolution(0,p);
			ret[0]-=sourceTerm(pPhys);
			ret[0]*=ret[0];
		}
	}errors;

public:
	std::size_t n_,p_,DIM_;
	Base::MeshType type_;
	double penaltyParameter_;
	Laplace(std::size_t numBasisFuns, std::size_t p,std::size_t DIM, Base::MeshType type):
		HpgemUISimplified(DIM,p),n_(numBasisFuns),p_(p),DIM_(DIM),type_(type),penaltyParameter_(3*n_*p_*(p_+DIM_-1)+1){ }

	bool initialise(){
    	RectangularMeshDescriptorT description(DIM_);
    	for(std::size_t i=0;i<DIM_;++i){
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
	    std::size_t numBasisFuns=el->getNrOfBasisFunctions();
    	LinearAlgebra::NumericalVector phiDerivI(DIM_),phiDerivJ(DIM_);
    	ret.resize(numBasisFuns,numBasisFuns);
    	for(std::size_t i=0;i<numBasisFuns;++i){
			phiDerivI = el->basisFunctionDeriv(i,p);
    		for(std::size_t j=0;j<numBasisFuns;++j){
    			phiDerivJ = el->basisFunctionDeriv(j,p);
    			ret(j,i)=phiDerivI*phiDerivJ;
    		}
    	}
	}

	//has an exact copy in the error class
    double sourceTerm(const PointPhysicalT& p){
    	double ret=std::sin(2*M_PI*p[0])*(4*M_PI*M_PI);
    	if(DIM_>1){
    		ret*=std::cos(2*M_PI*p[1]);
    	}
    	if(DIM_>2){
    		ret*=std::cos(2*M_PI*p[2])*3;
    	}
    	return ret;
    }

	//RHS
	void elementIntegrand(const Base::Element* el, const Geometry::PointReference& p, LinearAlgebra::NumericalVector& ret){
    	PointPhysicalT pPhys = el->referenceToPhysical(p);
    	ret.resize(el->getNrOfBasisFunctions());
    	for(std::size_t i=0;i<el->getNrOfBasisFunctions();++i){
    		ret[i]=el->basisFunction(i,p)*sourceTerm(pPhys);
    	}
	}

	//face discontinuities
	void faceIntegrand(const Base::Face* fa, const LinearAlgebra::NumericalVector& normal, const Geometry::PointReference& p, LinearAlgebra::Matrix& ret){
	    std::size_t numBasisFuns=fa->getNrOfBasisFunctions();
    	ret.resize(numBasisFuns,numBasisFuns);
    	LinearAlgebra::NumericalVector phiNormalI(DIM_),phiNormalJ(DIM_),phiDerivI(DIM_),phiDerivJ(DIM_);
		PointPhysicalT pPhys = fa->referenceToPhysical(p);
    	for(std::size_t i=0;i<numBasisFuns;++i){
			phiNormalI = fa->basisFunctionNormal(i,normal,p);
			phiDerivI = fa->basisFunctionDeriv(i,p);
    		for(std::size_t j=0;j<numBasisFuns;++j){
    			phiNormalJ = fa->basisFunctionNormal(j,normal,p);
    			phiDerivJ = fa->basisFunctionDeriv(j,p);
    			if(fa->isInternal()){
    				ret(j,i)=-(phiNormalI*phiDerivJ+phiNormalJ*phiDerivI)/2+
    							penaltyParameter_*phiNormalI*phiNormalJ;
    			}else if(std::abs(pPhys[0])<1e-12||std::abs(pPhys[0]-1.)<1e-12){//Dirichlet
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
	    std::size_t numBasisFuns=fa->getNrOfBasisFunctions();
		ret.resize(numBasisFuns);
		PointPhysicalT pPhys = fa->referenceToPhysical(p);
		if(std::abs(pPhys[0])<1e-9||std::abs(pPhys[0]-1)<1e-9){//Dirichlet
			LinearAlgebra::NumericalVector phiDeriv(DIM_);
			for(std::size_t i=0;i<numBasisFuns;++i){
				phiDeriv = fa->basisFunctionDeriv(i,p);
				ret[i]=(-normal*phiDeriv/Utilities::norm2(normal)+penaltyParameter_*fa->basisFunction(i,p))*0;
			}
		}else{
			for(std::size_t i=0;i<numBasisFuns;++i){
				ret[i]=0;
			}
		}
	}

	double initialConditions(const Geometry::PointPhysical& p){
		return 0;
	}

	double solve(){
    	doAllElementIntegration();
    	doAllFaceIntegration();
    	Utilities::GlobalPetscMatrix A(HpgemUI::meshes_[0],0,0);
    	Utilities::GlobalPetscVector b(HpgemUI::meshes_[0],0,0),x(HpgemUI::meshes_[0]);

    	b.assemble();
    	VecSet(x,0);

    	KSP ksp;
    	KSPCreate(PETSC_COMM_WORLD,&ksp);
    	KSPSetTolerances(ksp,1e-12,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
    	KSPSetOperators(ksp,A,A);
    	KSPSetFromOptions(ksp);
    	KSPSolve(ksp,b,x);

    	x.writeTimeLevelData(0);

    	double sum(0);
    	LinearAlgebra::NumericalVector ret(1);
    	Integration::ElementIntegral elIntegral(false);
    	elIntegral.setStorageWrapper(new Base::ShortTermStorageElementH1(meshes_[0]->dimension()));
    	for(Base::Element* el:meshes_[0]->getElementsList()){
    		ret = elIntegral.integrate(el,&errors);
    		sum+=ret[0];
    	}
    	return std::sqrt(sum);
	}
};

int main(int argc,char** argv){
    Base::parse_options(argc,argv);
	//no 3D testing due to speed related issues
	Laplace test0(1,2,1,Base::RECTANGULAR);
	test0.initialise();
        std::cout.precision(10);
        std::cout<<test0.solve()<<std::endl;
	logger.assert_always((std::abs(test0.solve()-0.35188045)<1e-8),"comparison to old results");
	Laplace test1(2,3,1,Base::RECTANGULAR);
	test1.initialise();
        std::cout<<test1.solve()<<std::endl;
	logger.assert_always((std::abs(test1.solve()-0.01607749)<1e-8),"comparison to old results");
	Laplace test2(4,4,1,Base::RECTANGULAR);
	test2.initialise();
        std::cout<<test2.solve()<<std::endl;
	logger.assert_always((std::abs(test2.solve()-0.00007200)<1e-8),"comparison to old results");
	Laplace test3(8,5,1,Base::RECTANGULAR);
	test3.initialise();
	logger.assert_always((std::abs(test3.solve()-0.00000008)<1e-8),"comparison to old results");
	Laplace test4(16,1,1,Base::RECTANGULAR);
	test4.initialise();
        std::cout<<test4.solve()<<std::endl;
	logger.assert_always((std::abs(test4.solve()-0.00880380)<1e-8),"comparison to old results");
	Laplace test5(1,2,2,Base::TRIANGULAR);
	test5.initialise();
        std::cout<<test5.solve()<<std::endl;
	//logger.assert_always(("comparison to old results",std::abs(test5.solve()-0.18046613)<1e-8));
	logger.assert_always((std::abs(test5.solve()-0.17226144)<1e-8),"comparison to old results");
	Laplace test6(2,3,2,Base::TRIANGULAR);
	test6.initialise();
        std::cout<<test6.solve()<<std::endl;
	logger.assert_always((std::abs(test6.solve()-0.01782337)<1e-8),"comparison to old results");
	Laplace test7(4,4,2,Base::TRIANGULAR);
	test7.initialise();
	logger.assert_always((std::abs(test7.solve()-0.00035302)<1e-8),"comparison to old results");
	Laplace test8(8,5,2,Base::TRIANGULAR);
	test8.initialise();
	logger.assert_always((std::abs(test8.solve()-0.00000061)<1e-8),"comparison to old results");
	Laplace test9(16,1,2,Base::TRIANGULAR);
	test9.initialise();
        std::cout<<test9.solve()<<std::endl;
	logger.assert_always((std::abs(test9.solve()-0.00448270)<1e-8),"comparison to old results");
	return 0;
}
