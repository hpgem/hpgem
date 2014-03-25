/*
 * main.cpp
 *
 * Laplace equations (test problem)
 *
 *  Created on: Jan 8, 2014
 *      Author: brinkf
 */

#define hpGEM_INCLUDE_PETSC_SUPPORT

#include "Base/HpgemUISimplified.hpp"
#include "Base/Norm2.hpp"
#include "Utilities/GlobalMatrix.hpp"
#include "Utilities/GlobalVector.hpp"
#include "petscksp.h"
#include "Output/TecplotSingleElementWriter.hpp"
#include <fstream>
#include "Geometry/ReferenceTetrahedron.hpp"
#include "Utilities/BasisFunctions1DH1ConformingLine.hpp"
#include "Utilities/BasisFunctions2DH1ConformingSquare.hpp"
#include "Utilities/BasisFunctions2DH1ConformingTriangle.hpp"
#include "Utilities/BasisFunctions3DH1ConformingCube.hpp"
#include "Utilities/BasisFunctions3DH1ConformingTetrahedron.hpp"

const unsigned int DIM=2;
//penaltyParameter>3np(p+dim-1) - set in the code
double penaltyParameter=-3.141;

class Laplace : public Base::HpgemUISimplified,Output::TecplotSingleElementWriter{

public:
	Laplace(int n,int p):HpgemUISimplified(DIM,p),n_(n),p_(p){
		penaltyParameter=3*n_*p_*(p_+DIM-1)+1;
	}

	///set up the mesh
    bool virtual initialise(){
    	RectangularMeshDescriptorT description(DIM);
    	for(int i=0;i<DIM;++i){
    		description.bottomLeft_[i]=0;
    		description.topRight_[i]=1;
    		description.numElementsInDIM_[i]=n_;
    		description.boundaryConditions_[i]=RectangularMeshDescriptorT::SOLID_WALL;
    	}
    	addMesh(description,Base::TRIANGULAR,1,1,1,1);
    	meshes_[0]->setDefaultBasisFunctionSet(Utilities::createInteriorBasisFunctionSet2DH1Triangle(p_));
    	std::vector<const Base::BasisFunctionSet*> bFsets;
    	Utilities::createVertexBasisFunctionSet2DH1Triangle(p_,bFsets);
    	meshes_[0]->addVertexBasisFunctionSet(bFsets);
    	std::vector<const Base::OrientedBasisFunctionSet*> oBFsets;
    	Utilities::createFaceBasisFunctionSet2DH1Triangle(p_,oBFsets);
    	meshes_[0]->addFaceBasisFunctionSet(oBFsets);
    	//oBFsets.clear();
    	//Utilities::createEdgeBasisFunctionSet3DH1Tetrahedron(p_,oBFsets);
    	//meshes_[0]->addEdgeBasisFunctionSet(oBFsets);
    	return true;
    }

    ///You pass the reference point to the basisfunctions. Internally the basisfunctions will be mapped to the physical element
    ///so you wont have to do any transformations yourself
    virtual void elementIntegrand(const ElementT* element, const PointReferenceT& p, LinearAlgebra::Matrix& ret){
    	int n=element->getNrOfBasisFunctions();
    	NumericalVector phiDerivI(DIM),phiDerivJ(DIM);
    	ret.resize(n,n);
    	for(int i=0;i<n;++i){
			element->basisFunctionDeriv(i,p,phiDerivI);
    		for(int j=0;j<n;++j){
    			element->basisFunctionDeriv(j,p,phiDerivJ);
    			ret(j,i)=phiDerivI*phiDerivJ;
    		}
    	}
    }

    ///You pass the reference point to the basisfunctions. Internally the basisfunctions will be mapped to the physical element
    ///so you wont have to do any transformations yourself. If you expect 4 matrices here, you can assume that ret is block structured such
    ///that basisfunctions belonging to the left element are indexed first
    ///note that using a consistent flux has no effect if you also use conforming basisfunctions
    virtual void faceIntegrand(const FaceT* face, const NumericalVector& normal, const PointReferenceOnTheFaceT& p,  LinearAlgebra::Matrix& ret){
    	int n=face->getNrOfBasisFunctions();
    	ret.resize(n,n);
    	NumericalVector phiNormalI(DIM),phiNormalJ(DIM),phiDerivI(DIM),phiDerivJ(DIM);
		PointPhysicalT pPhys(DIM);
		face->referenceToPhysical(p,pPhys);
    	for(int i=0;i<n;++i){
			face->basisFunctionNormal(i,normal,p,phiNormalI);
			face->basisFunctionDeriv(i,p,phiDerivI);
    		for(int j=0;j<n;++j){
    			face->basisFunctionNormal(j,normal,p,phiNormalJ);
    			face->basisFunctionDeriv(j,p,phiDerivJ);
        		//if(face->faceType_==(...))
    			if(face->isInternal()){
    				ret(j,i)=-(phiNormalI*phiDerivJ+phiNormalJ*phiDerivI)/2+
    							penaltyParameter*phiNormalI*phiNormalJ;
    			}else if(fabs(pPhys[0])<1e-9||fabs(pPhys[0]-1.)<1e-9){//Dirichlet
    				ret(j,i)=-(phiNormalI*phiDerivJ+phiNormalJ*phiDerivI)+
    							penaltyParameter*phiNormalI*phiNormalJ*2;
    			//}else if(fabs(pPhys[0]-1)<1e-9){//Robin
    			//	ret(j,i)=-phiNormalI*phiNormalJ*0;
    			}else{//homogeneous Neumann
    				ret(j,i)=0;
    			}
    		}
    	}
    }

    ///The vector edition of the face integrand is meant for implementation of boundary conditions
    virtual void faceIntegrand(const FaceT* face, const NumericalVector& normal, const PointReferenceOnTheFaceT& p,  LinearAlgebra::NumericalVector& ret){
		int n=face->getNrOfBasisFunctions();
		ret.resize(n);
		PointPhysicalT pPhys(DIM);
		face->referenceToPhysical(p,pPhys);
		//if(face->faceType_==(...))
		/*if(fabs(pPhys[0]-1)<1e-9){//Neumann and robin
			for(int i=0;i<n;++i){//be carefull in 1D; this boundary condition always concerns df\dn
				ret[i]=face->basisFunction(i,p)*-1;
			}
		}else */if(fabs(pPhys[0])<1e-9||fabs(pPhys[0]-1)<1e-9){//Dirichlet
			LinearAlgebra::NumericalVector phiDeriv(DIM);
			for(int i=0;i<n;++i){
				face->basisFunctionDeriv(i,p,phiDeriv);
				ret[i]=(-normal*phiDeriv/Utilities::norm2(normal)+penaltyParameter*face->basisFunction(i,p))*0;
			}
		}else{
			for(int i=0;i<n;++i){
				ret[i]=0;
			}
		}
    }

    virtual double initialConditions(const PointPhysicalT& p){
    	// initial conditions are not needed for a steady-state problem
    	return 0;
    }

    double sourceTerm(const PointPhysicalT& p){
    	return sin(2*M_PI*p[0])*(4*M_PI*M_PI)*cos(2*M_PI*p[1]);//*cos(2*M_PI*p[2])*3;
    }

    ///interpolates the source term
    void elementIntegrand(const ElementT* element, const PointReferenceT& p, LinearAlgebra::NumericalVector& ret){
    	PointPhysicalT pPhys(DIM);
    	element->referenceToPhysical(p,pPhys);
    	ret.resize(element->getNrOfBasisFunctions());
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
    	out<<value;
    }

    ///guarantees the linear system keeps the dirichlet boundary conditions in place. Assumes x already contains expansion coefficients for
    ///the boundaries in question and noise in all other entries
    void insertDirichletBoundary(Utilities::GlobalPetscMatrix& A, Utilities::GlobalPetscVector& b, Utilities::GlobalPetscVector& x){
    	int numberOfRows(0);
    	std::vector<int> rows(0);
    	Geometry::PointPhysical pPhys(DIM);
    	Geometry::PointReference centre(DIM-1);
    	for(Base::Face* face:meshes_[0]->getFacesList()){
    		face->getReferenceGeometry()->getCenter(centre);
    		face->referenceToPhysical(centre,pPhys);
    		//if(face->faceType_=(...))
    		if(fabs(pPhys[0])<1e-9||fabs(pPhys[0]-1)<1e-9){
    			A.getMatrixBCEntries(face,numberOfRows,rows);
    		}
    	}
		int ierr=MatZeroRows(A,numberOfRows,&rows[0],1.0,x,b);
		CHKERRV(ierr);
    }

    bool solve(){
    	doAllElementIntegration();
    	//doAllFaceIntegration();
    	Utilities::GlobalPetscMatrix A(HpgemUI::meshes_[0],0,-1);
    	Utilities::GlobalPetscVector b(HpgemUI::meshes_[0],0,-1),x(HpgemUI::meshes_[0]);

    	b.assemble();
    	VecSet(x,0);
    	insertDirichletBoundary(A,b,x);

    	KSP ksp;
    	KSPCreate(MPI_COMM_WORLD,&ksp);
    	KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN);
    	KSPSetFromOptions(ksp);
    	KSPSolve(ksp,b,x);
    	KSPConvergedReason conferge;
    	KSPGetConvergedReason(ksp,&conferge);
    	int iterations;
    	KSPGetIterationNumber(ksp,&iterations);
    	cout<<"KSP solver ended because of "<<KSPConvergedReasons[conferge]<<" in "<<iterations<<" iterations."<<endl;

    	x.writeTimeLevelData(0);

    	std::ofstream outFile("output.dat");
		Output::TecplotDiscontinuousSolutionWriter writeFunc(outFile,"test","01","value");
		writeFunc.write(meshes_[0],"continuous solution",false,this);
		return true;
    }

private:

    //number of elements per cardinal direction
    int n_;

    //polynomial order of the approximation
    int p_;
};

int main(int argc, char **argv){
	try{
		int n,p;
		if(argc>2){
			n=std::atoi(argv[1]);
			p=std::atoi(argv[2]);
			argv[2]=argv[0];
			argc-=2;
			argv+=2;
		}else{
			throw "usage: Laplace.out n p [petsc-args]";
		}

		PetscInitialize(&(argc),&(argv),NULL,NULL);
		Laplace test(n,p);
		test.initialise();
		test.solve();
		PetscFinalize();
		return 0;
	}
	catch(const char* e){
		cout<<e;
	}
}


