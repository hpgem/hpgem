/*
 * main.cpp
 *
 * Laplace equations (test problem)
 *
 *  Created on: Jan 8, 2014
 *      Author: brinkf
 */

#include "Base/HpgemUISimplified.hpp"
#include "Base/Norm2.hpp"
#include "GlobalAssembly.hpp"
#include "petscksp.h"
#include "Output/TecplotSingleElementWriter.hpp"
#include <fstream>

const unsigned int DIM=1;
double penaltyParameter=256;

class Laplace : public Base::HpgemUISimplified,Output::TecplotSingleElementWriter{

public:
	Laplace(int n,int p):HpgemUISimplified(DIM),n_(n),p_(p){}

	///set up the mesh
    bool virtual initialise(){
    	RectangularMeshDescriptorT description(DIM);
    	for(int i=0;i<DIM;++i){
    		description.bottomLeft_[i]=0;
    		description.topRight_[i]=1;
    		description.numElementsInDIM_[i]=n_;
    		description.boundaryConditions_[i]=RectangularMeshDescriptorT::SOLID_WALL;
    	}
    	addMesh(description,Base::RECTANGULAR,1,1,1,1);
    	return true;
    }

    ///You pass the reference point to the basisfunctions. Internally the basisfunctions will be mapped to the physical element
    ///so you wont have to do any transformations yourself
    virtual void elementIntegrand(const ElementT* element, const PointReferenceT& p, LinearAlgebra::Matrix& ret){
    	int n=element->getNrOfBasisFunctions();
    	NumericalVector phiDerivI(DIM),phiDerivJ(DIM);
    	for(int i=0;i<n;++i){
    		for(int j=0;j<n;++j){
    			element->basisFunctionDeriv(i,p,phiDerivI);
    			element->basisFunctionDeriv(j,p,phiDerivJ);
    			ret(j,i)=phiDerivI*phiDerivJ;
    		}
    	}
    }

    ///You pass the reference point to the basisfunctions. Internally the basisfunctions will be mapped to the physical element
    ///so you wont have to do any transformations yourself. If you expect 4 matrices here, you can assume that ret is block structured such
    ///that basisfunctions belonging to the left element are indexed first
    virtual void faceIntegrand(const FaceT* face, const PointPhysicalT& normal, const PointReferenceOnTheFaceT& p,  LinearAlgebra::Matrix& ret){
    	NumericalVector unitNormal(&normal[0],normal.size());
    	unitNormal/=Utilities::norm2(normal);
    	const ElementT *elLeft(face->getPtrElementLeft()),*elRight(face->getPtrElementRight());
    	PointReferenceT pLeft(DIM),pRight(DIM);
    	int n(elLeft->getNrOfBasisFunctions()),nLeft(n);
    	face->mapRefFaceToRefElemL(p,pLeft);
    	PointPhysicalT pPhys(DIM);
    	elLeft->referenceToPhysical(pLeft,pPhys);

    	if(face->isInternal()){
    		n+=elRight->getNrOfBasisFunctions();
    		face->mapRefFaceToRefElemR(p,pRight);
    	}
    	ret.resize(n,n);

    	double phiI,phiJ;
    	NumericalVector phiDerivI(DIM),phiDerivJ(DIM);
    	for(int i=0;i<n;++i){
    		for(int j=0;j<n;++j){
    			if(i<nLeft){
    				phiI=elLeft->basisFunction(i,pLeft);
    				elLeft->basisFunctionDeriv(i,pLeft,phiDerivI);
    			}else{
    				phiI=-elRight->basisFunction(i-nLeft,pRight);//more accurately the '-' should be with the normal vector
    				elRight->basisFunctionDeriv(i-nLeft,pRight,phiDerivI);
    			}
    			if(j<nLeft){
    				phiJ=elLeft->basisFunction(j,pLeft);
    				elLeft->basisFunctionDeriv(j,pLeft,phiDerivJ);
    			}else{
    				phiJ=-elRight->basisFunction(j-nLeft,pRight);//more accurately the '-' should be with the normal vector
    				elRight->basisFunctionDeriv(j-nLeft,pRight,phiDerivJ);
    			}
    			if(pPhys[0]==1||pPhys[0]==0){//Dirichlet
					ret(j,i)=-unitNormal*((phiI)*phiDerivJ+(phiJ)*phiDerivI)
							  +penaltyParameter*(phiI)*(phiJ);
    			}else if(face->isInternal()){
					ret(j,i)=-unitNormal*(phiI*phiDerivJ+phiJ*phiDerivI)/2
							  +penaltyParameter*phiI*phiJ;
    		//	}else if(pPhys[0]==0){//Robin
    		//		ret(j,i)=-phiI*phiJ*0;
    			}else{//homogeneous Neumann
    				ret(j,i)=0;
    			}
    		}
    	}
    }

    ///The vector edition of the face integrand is meant to support implementation of boundary conditions
    virtual void faceIntegrand(const FaceT* face, const PointPhysicalT& normal, const PointReferenceOnTheFaceT& p,  LinearAlgebra::NumericalVector& ret){
    	if(face->isInternal()){
    		int n=face->getPtrElementLeft()->getNrOfBasisFunctions()+face->getPtrElementRight()->getNrOfBasisFunctions();
    		ret.resize(n);
    		for(int i=0;i<n;++i){
    			ret[i]=0;
    		}
    	}else{
    		int n=face->getPtrElementLeft()->getNrOfBasisFunctions();
    		ret.resize(n);
    		const ElementT* element=face->getPtrElementLeft();
    		PointReferenceT pElement(DIM);
    		PointPhysicalT pPhys(DIM);
    		face->mapRefFaceToRefElemL(p,pElement);
    		element->referenceToPhysical(pElement,pPhys);
    		//if(face->faceType_==(...))
    		if(pPhys[0]==2000){//Neumann and robin(toggled off)
    			for(int i=0;i<n;++i){//be carefull in 1D; this boundary condition always concerns df\dn
    				ret[i]=element->basisFunction(i,pElement)*-1;
    			}
    		}else if(pPhys[0]==1||pPhys[0]==0){//Dirichlet
    	    	NumericalVector unitNormal(&normal[0],normal.size());
    	    	unitNormal/=Utilities::norm2(normal);
				LinearAlgebra::NumericalVector phiDeriv(DIM);
    			for(int i=0;i<n;++i){
    				element->basisFunctionDeriv(i,pElement,phiDeriv);
    				ret[i]=(-unitNormal*phiDeriv+penaltyParameter*element->basisFunction(i,pElement))*0;
    			}
    		}else{
    			for(int i=0;i<n;++i){
    				ret[i]=0;
    			}
    		}
		}
    }

    virtual void initialConditions(const PointPhysicalT& p){
    	// initial conditions are not needed for a steady-state problem
    }

    double sourceTerm(const PointPhysicalT& p){
    	return sin(2*M_PI*p[0])*(4*M_PI*M_PI)*cos(2*M_PI*p[1])*cos(2*M_PI*p[2])*3;
    }

    void elementIntegrand(const ElementT* element, const PointReferenceT& p, LinearAlgebra::NumericalVector& ret){
    	PointPhysicalT pPhys(DIM);
    	element->referenceToPhysical(p,pPhys);
    	for(int i=0;i<element->getNrOfBasisFunctions();++i){
    		ret[i]=element->basisFunction(i,p)*sourceTerm(pPhys);
    	}
    }

    void writeToTecplotFile(const ElementT* element,const  PointReferenceT& p, ostream& out){
    	double value(0);
    	Geometry::PointPhysical pPhys(DIM);
    	element->referenceToPhysical(p,pPhys);
    	for(int i=0;i<element->getNrOfBasisFunctions();++i){
    		value+=element->basisFunction(i,p)*element->getData(0,0,i);
    	}
    	out<<value-sin(2*M_PI*pPhys[0]);
    }

    bool solve(){
    	doAllElementIntegration(0);
    	doAllFaceIntegration(0);
    	Utilities::GlobalPetscMatrix A(HpgemUI::meshes_[0],0,0);
    	Utilities::GlobalPetscVector b(HpgemUI::meshes_[0],0,0),x(HpgemUI::meshes_[0]);
    	b.assemble();

    	KSP ksp;
    	KSPCreate(MPI_COMM_WORLD,&ksp);
    	KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN);
    	KSPSetFromOptions(ksp);

    	KSPSolve(ksp,b,x);

    	x.writeTimeLevelData(0);
    	std::ofstream outFile("output.dat");
		Output::TecplotDiscontinuousSolutionWriter writeFunc(outFile,"test","0","value");
		writeFunc.write(meshes_[0],"Steady state solution",false,this);
		return true;
    }

private:

    //number of elements per cardinal direction
    int n_;

    //polynomial order of the approximation
    int p_;
};

int main(int argc, char **argv){
	PetscInitialize(&argc,&argv,NULL,NULL);
	Laplace test(8,2);
	test.initialise();
	test.solve();
	PetscFinalize();
	return 0;
}


