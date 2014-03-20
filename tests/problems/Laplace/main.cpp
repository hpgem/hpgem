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
#include <chrono>

const unsigned int DIM=2;
//penaltyParameter>2p^2/h=2np^2
double penaltyParameter=-3.141;

void viewBasisFunctionInTecplot(const Base::BaseBasisFunction* function, const Geometry::ReferenceGeometry& geometry, ostream& file)
{
//	file<<"TITLE= \"BasisFunction\"\nVARIABLES= x,y,z,value"<<endl;
	file<<"ZONE I="<<41 <<" J= "<<41<<" K= "<<41<<endl;
	Geometry::PointReference p(3);
	for(p[0]=-1.;p[0]<1.001;p[0]+=0.05){
		for(p[1]=-1.;p[1]<1.001;p[1]+=0.05){
			for(p[2]=-1.;p[2]<1.001;p[2]+=0.05){
				file<<p[0]<<" "<<p[1]<<" "<<p[2]<<" ";
				if(geometry.isInternalPoint(p)){
					file<<function->eval(p)<<" "<<function->evalDeriv0(p)<<" "<<function->evalDeriv1(p)<<" "<<function->evalDeriv2(p)<<endl;
				}else{
					file<<"0.0 0.0 0.0 0.0"<<endl;
				}
			}
		}
	}
}

class dummyDouble{
public:
	dummyDouble(double d):d_(d){}
	void axpy(double a,dummyDouble x){
		d_+=a*x.d_;
	}
	void operator*=(dummyDouble d){
		d_*=d.d_;
	}
	void operator+=(dummyDouble d){
		d_+=d.d_;
	}
	operator double(){
		return d_;
	}
private:
	double d_;
};

class Laplace : public Base::HpgemUISimplified,Output::TecplotSingleElementWriter,public Integration::ElementIntegrandBase<dummyDouble>{

public:	//I want to override the basisfunctions so for now just claim the polynomial order should be 2 (so the basis function set always exists)
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
    	addMesh(description,Base::TRIANGULAR,1,1,1,1);//not exactly 'typesafe' but this works for sane users
    	/*meshes_[0]->setDefaultBasisFunctionSet(Utilities::createInteriorBasisFunctionSet3DH1Tetrahedron(p_));
    	std::vector<const Base::BasisFunctionSet*> bFsets;
    	Utilities::createVertexBasisFunctionSet3DH1Tetrahedron(p_,bFsets);
    	meshes_[0]->addVertexBasisFunctionSet(bFsets);
    	std::vector<const Base::OrientedBasisFunctionSet*> oBFsets;
    	Utilities::createFaceBasisFunctionSet3DH1Tetrahedron(p_,oBFsets);
    	meshes_[0]->addFaceBasisFunctionSet(oBFsets);
    	oBFsets.clear();
    	Utilities::createEdgeBasisFunctionSet3DH1Tetrahedron(p_,oBFsets);
    	meshes_[0]->addEdgeBasisFunctionSet(oBFsets);*/
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
    			//}else if(pPhys[0]==1){//Robin
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
		/*if(pPhys[0]==2000){//Neumann and robin(toggled off)
			for(int i=0;i<n;++i){//be carefull in 1D; this boundary condition always concerns df\dn
				ret[i]=face->basisFunction(i,p)*-1;
			}
		}else */if(fabs(pPhys[0])<1e-9||fabs(pPhys[0]-1)<1e-9){//Dirichlet
			//NumericalVector unitNormal(&normal[0],normal.size());
			//unitNormal/=Utilities::norm2(normal);
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
    	//return 2*p[1]*p[1]*(1-p[1])*(1-p[1])+(8*p[1]*(1-p[1])-2*p[1]*p[1]-2*(1-p[1])*(1-p[1]))*p[0]*(1-p[0]);
    }

    void elementIntegrand(const ElementT* element, const PointReferenceT& p, dummyDouble& ret){
    	PointPhysicalT pPhys(DIM);
    	element->referenceToPhysical(p,pPhys);
    	ret=0;
    	for(int i=0;i<element->getNrOfBasisFunctions();++i){
    		ret+=element->basisFunction(i,p)*element->getData(0,0,i);
    	}
    	ret+=-sourceTerm(pPhys)/(4*M_PI*M_PI)/(DIM);
    	//ret +=-pPhys[0]*(1-pPhys[0])*pPhys[1]*pPhys[1]*(1-pPhys[1])*(1-pPhys[1]);
    	ret*=ret;
    }

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
    	std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    	doAllElementIntegration();
    	doAllFaceIntegration();
    	Utilities::GlobalPetscMatrix A(HpgemUI::meshes_[0],0,0);
    	Utilities::GlobalPetscVector b(HpgemUI::meshes_[0],0,0),x(HpgemUI::meshes_[0]);

    	b.assemble();
    	VecSet(x,0);
    	//insertDirichletBoundary(A,b,x);

    	//b.writeTimeLevelData(0);
    	//std::ofstream testFile("testOutput.dat");
		//Output::TecplotDiscontinuousSolutionWriter testWriteFunc(testFile,"test","012","value");
		//testWriteFunc.write(meshes_[0],"source term",false,this);

    	KSP ksp;
    	KSPCreate(MPI_COMM_WORLD,&ksp);
    	KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN);
    	KSPSetFromOptions(ksp);
    	KSPSolve(ksp,b,x);
    	KSPConvergedReason conferge;
    	KSPGetConvergedReason(ksp,&conferge);
    	int iterations;
    	KSPGetIterationNumber(ksp,&iterations);
    	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    	cout<<"KSP solver ended because of "<<KSPConvergedReasons[conferge]<<" in "<<iterations<<" iterations."<<endl;

    	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double> >(end-start);
    	cout<<"This took "<<time_span.count()<<" seconds."<<endl;
    	x.writeTimeLevelData(0);

    	std::ofstream outFile("output.dat");
		Output::TecplotDiscontinuousSolutionWriter writeFunc(outFile,"test","01","value");
		writeFunc.write(meshes_[0],"Monomial solution",false,this);
    	//viewBasisFunctionInTecplot((*elementColBegin())->getBasisFunction(0),Geometry::ReferenceLine::Instance(),outFile);
    	//viewBasisFunctionInTecplot((*elementColBegin())->getBasisFunction(1),Geometry::ReferenceLine::Instance(),outFile);

    	Integration::ElementIntegral elIntegral(false);
    	dummyDouble error(0),elError(0);
    	for(Base::Element* el:meshes_[0]->getElementsList()){
    		elIntegral.integrate<dummyDouble>(el,this,elError);
    		error+=elError;
    		elError=0;
    	}
    	cout<<"The L2-error in the monomial case is: "<<sqrt(error)<<endl;

    	/*Base::BasisFunctionSet* DG= new Base::BasisFunctionSet(p_);
    	switch(p_){
    	case 1:
    		Base::AssembleBasisFunctionSet_2D_Ord1_A1(*DG);
    		break;
    	case 2:
    		Base::AssembleBasisFunctionSet_2D_Ord2_A1(*DG);
    		break;
    	case 3:
    		Base::AssembleBasisFunctionSet_2D_Ord3_A1(*DG);
    		break;
    	case 4:
    		Base::AssembleBasisFunctionSet_2D_Ord4_A1(*DG);
    		break;
    	case 5:
    		Base::AssembleBasisFunctionSet_2D_Ord5_A1(*DG);
    		break;
    	}

    	meshes_[0]->setDefaultBasisFunctionSet(Utilities::createDGBasisFunctionSet3DH1Tetrahedron(p_));

    	VecSet(x,0);//no cheating with pre-converged solutions

    	start = std::chrono::steady_clock::now();

    	doAllElementIntegration();
    	doAllFaceIntegration();
    	Utilities::GlobalPetscMatrix DGA(HpgemUI::meshes_[0],0,0);
    	Utilities::GlobalPetscVector DGb(HpgemUI::meshes_[0],0,0);
    	DGb.assemble();
    	x.reset();

    	KSPDestroy(&ksp);//stupid stubborn ksp...
    	KSPCreate(MPI_COMM_WORLD,&ksp);
    	KSPSetOperators(ksp,DGA,DGA,DIFFERENT_NONZERO_PATTERN);
    	KSPSetFromOptions(ksp);
    	KSPSolve(ksp,DGb,x);
    	KSPGetConvergedReason(ksp,&conferge);
    	KSPGetIterationNumber(ksp,&iterations);
    	end=std::chrono::steady_clock::now();
    	cout<<"KSP solver ended because of "<<KSPConvergedReasons[conferge]<<" in "<<iterations<<" iterations."<<endl;
    	time_span=std::chrono::duration_cast<std::chrono::duration<double> >(end-start);
    	cout<<"This took "<<time_span.count()<<" seconds."<<endl;
    	x.writeTimeLevelData(0);

		writeFunc.write(meshes_[0],"Discontinuous solution",false,this);
		for(int i=0;i<(*elementColBegin())->getNrOfBasisFunctions();++i){
			//viewBasisFunctionInTecplot((*elementColBegin())->getBasisFunction(i),Geometry::ReferenceSquare::Instance(),outFile);
		}
    	//viewBasisFunctionInTecplot((*elementColBegin())->getBasisFunction(13),Geometry::ReferenceTetrahedron::Instance(),outFile);

		error=0;
		elError=0;
    	for(Base::Element* el:meshes_[0]->getElementsList()){
    		elIntegral.integrate<dummyDouble>(el,this,elError);
    		error+=elError;
    		elError=0;
    	}
    	cout<<"The L2-error in the DG case is: "<<sqrt(error)<<endl;*/
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


