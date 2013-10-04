//
//  HEuler.cpp
//  
//
//  Created by Shavarsh Nurijanyan on 9/3/13.
//
//

#include "HEuler.hpp"


HEuler::HEuler(HEulerGlobalVariables* global, const HEulerConfigurationData* config):
    Base::HpgemUI<DIM>(global, config),
    P_(),
    Q_()
{
    switch(config->solutionType_)
	{
		case HEulerConfigurationData::INCOMPRESSIBLE_WALLS:
			exactSolution_ = new InitialVelocityConstructorTaylor (config->lx_,config->ly_,config->lz_);
            break;
		case HEulerConfigurationData::INCOMPRESSIBLE_PERIODIC:
			exactSolution_ = new Incompressible3DPeriodic();
            break;
		case HEulerConfigurationData::COMPRESSIBLE_PERIODIC:
			exactSolution_ = new Compressible3DPeriodic();
            break;
		case HEulerConfigurationData::INCOMPRESSIBLE_ONETHIRDPERIODIC:
			exactSolution_ = new Compressible3DOneThirdPeriodic(config->lx_,config->ly_,config->lz_);
			break;
		case HEulerConfigurationData::COMPRESSIBLE_WALLS:
			exactSolution_ = new Compressible3DWalls();
            break;
		default:
			cout<<"Can't happen. Major error!"<<endl;
            break;
	}

        //exactSolution_= new Compressible3DPeriodic();
}

HEuler::~HEuler()
{
    
        MatDestroy(&P_);
        MatDestroy(&Q_);
        PetscFinalize();
}

void
HEuler::outputMatrix(Mat& matrix, const string& name)
{
    cout << "Mat Assembly for "<< name<<endl;
    MatAssemblyBegin(matrix,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(matrix,MAT_FINAL_ASSEMBLY);
    
    PetscViewer viewer;
    PetscViewerASCIIOpen(PETSC_COMM_SELF,name.c_str(),&viewer);
    MatView(matrix,viewer);
}
void
HEuler::outputMatrix(const Mat& matrix, const string& name)const 
{
    cout << "Mat Assembly for "<< name<<endl;
    MatAssemblyBegin(matrix,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(matrix,MAT_FINAL_ASSEMBLY);
    
    PetscViewer viewer;
    PetscViewerASCIIOpen(PETSC_COMM_SELF,name.c_str(),&viewer);
    MatView(matrix,viewer);
}
void
HEuler::outputVectorMatlab(const Vec& vec, const string& name)const
{
    cout << "Vec Assembly for "<< name<<endl;
    VecAssemblyBegin(vec);
    VecAssemblyEnd(vec);
    
    PetscViewer viewer;
        //PetscViewerBinaryOpen(PETSC_COMM_WORLD, name.c_str(), FILE_MODE_WRITE, &viewer);
     PetscViewerASCIIOpen(PETSC_COMM_SELF,name.c_str(),&viewer);
    VecView(vec, viewer);
}

void
HEuler::outputVectorMatlab(Vec& vec, const string& name)
{
    cout << "Vec Assembly for "<< name<<endl;
    VecAssemblyBegin(vec);
    VecAssemblyEnd(vec);
    
    PetscViewer viewer;
        //    PetscViewerBinaryOpen(PETSC_COMM_WORLD, name.c_str(), FILE_MODE_WRITE, &viewer);
    PetscViewerASCIIOpen(PETSC_COMM_SELF,name.c_str(),&viewer);

    VecView(vec, viewer);
}

void
HEuler::printFullMatrixInfo(Mat& matrix, const string& name)
{
    PetscInt m=0,n=0;
    MatGetSize(matrix,&m,&n);
    
    MatInfo info;
    MatGetInfo(matrix,MAT_LOCAL, &info);
    cout<<name<<endl;
    printf("N = %d, N = %d, block_size = %d, memory = %d, assemblies = %d, mallocs = %d, matrix nonzeros (SeqBAIJ format) = %d, allocated nonzeros= %d\n", m, n, (int)info.block_size, (int)info.memory, (int)info.assemblies, (int)info.mallocs,(int)info.nz_used,(int)info.nz_allocated);
}

bool
HEuler::initialiseMesh()
{
    RectangularMeshDescriptor<DIM> rectangularMesh;
    
    const HEulerConfigurationData* config = static_cast<const HEulerConfigurationData*>(configData_);
    
   
    if(config->solutionType_==HEulerConfigurationData::COMPRESSIBLE_WALLS || config->solutionType_==HEulerConfigurationData::INCOMPRESSIBLE_WALLS)
    {
        rectangularMesh.boundaryConditions_[0] = RectangularMeshDescriptor<DIM>::SOLID_WALL;
        rectangularMesh.boundaryConditions_[1] = RectangularMeshDescriptor<DIM>::SOLID_WALL;
        rectangularMesh.boundaryConditions_[2] = RectangularMeshDescriptor<DIM>::SOLID_WALL;
    }
    else if(config->solutionType_==HEulerConfigurationData::COMPRESSIBLE_PERIODIC || config->solutionType_==HEulerConfigurationData::INCOMPRESSIBLE_PERIODIC)
    {
        rectangularMesh.boundaryConditions_[0] = RectangularMeshDescriptor<DIM>::PERIODIC;
        rectangularMesh.boundaryConditions_[1] = RectangularMeshDescriptor<DIM>::PERIODIC;
        rectangularMesh.boundaryConditions_[2] = RectangularMeshDescriptor<DIM>::PERIODIC;
    }
    else if(config->solutionType_==HEulerConfigurationData::INCOMPRESSIBLE_ONETHIRDPERIODIC)
    {
        rectangularMesh.boundaryConditions_[0] = RectangularMeshDescriptor<DIM>::SOLID_WALL;
        rectangularMesh.boundaryConditions_[1] = RectangularMeshDescriptor<DIM>::PERIODIC;
        rectangularMesh.boundaryConditions_[2] = RectangularMeshDescriptor<DIM>::SOLID_WALL;
    }
    rectangularMesh.bottomLeft_[0]       = 0;
    rectangularMesh.bottomLeft_[1]       = 0;
    rectangularMesh.bottomLeft_[2]       = 0;
    rectangularMesh.topLeft_[0]          = config->lx_;
    rectangularMesh.topLeft_[1]          = config->ly_;
    rectangularMesh.topLeft_[2]          = config->lz_;
    rectangularMesh.numElementsInDIM_[0] = config->nx_;
    rectangularMesh.numElementsInDIM_[1] = config->ny_;
    rectangularMesh.numElementsInDIM_[2] = config->nz_;
    
    MeshId id = addMesh(rectangularMesh);
    
    HEulerGlobalVariables* globalData = static_cast<HEulerGlobalVariables*>(globalData_);
    
    globalData->nElements_ = getNumberOfElements(id);
    globalData->dt_         = config->onePeriod_/config->numOfTimeStepInOnePeriod_;
 
    
    return true;
}

void
HEuler::calculateMassMatrix(const ElementT* element, const PointReferenceT& p, LinearAlgebra::Matrix& massMatrix)
{
    unsigned int numOfDegreesOfFreedom = element->getNrOfBasisFunctions();
    
    massMatrix.resize(numOfDegreesOfFreedom, numOfDegreesOfFreedom);
    
        // cout << "POINT="<<p<<endl;
    
    for (unsigned int i=0; i < numOfDegreesOfFreedom; ++i)
    {
        for (unsigned int j=0; j < numOfDegreesOfFreedom; ++j)
        {
            massMatrix(i,j) = element->basisFunction(j,p) * element->basisFunction(i,p);
            
                //cout << "i="<<i<<", j="<<j<<"phiI"<<element->basisFunction(i,p)<<", phiJ"<<element->basisFunction(j,p)<<endl;
        }
    }
}

void
HEuler::calculateLocalEnergy(const ElementT& element, const PointReferenceT& p, double& returnValue)
{
    double extra=0;
    ElementT::SolutionVector  solution;
    
    element.getSolution(0, p, solution);
    
    HEulerConfigurationData::SolutionType sType = static_cast<const HEulerConfigurationData*>(configData_)->solutionType_;
    
    if (sType  == HEulerConfigurationData::COMPRESSIBLE_WALLS || sType  == HEulerConfigurationData::COMPRESSIBLE_PERIODIC)
    {
        extra =0.5*(solution[0]*solution[0]);
    }
    returnValue = 0.5*(solution[1]*solution[1] + solution[2]*solution[2] + solution[3]*solution[3])+extra;
}

void
HEuler::elementIntegrand(const ElementT* element, const PointReferenceT& p, ElementIntegralData& returnObject)
{
    
    unsigned int numberOfDegreesOfFreedom = element->getNrOfBasisFunctions();
    
    LinearAlgebra::Matrix& xDerReturnData = returnObject.xGrad_;
    LinearAlgebra::Matrix& yDerReturnData = returnObject.yGrad_;
    LinearAlgebra::Matrix& zDerReturnData = returnObject.zGrad_;
    
    xDerReturnData.resize(numberOfDegreesOfFreedom, numberOfDegreesOfFreedom);
    yDerReturnData.resize(numberOfDegreesOfFreedom, numberOfDegreesOfFreedom);
    zDerReturnData.resize(numberOfDegreesOfFreedom, numberOfDegreesOfFreedom);
    
    
    NumericalVector grads(3);
    
    for (unsigned int i=0; i < numberOfDegreesOfFreedom; ++i)
    {
        Utilities::PhysGradientOfBasisFunction<DIM, ElementT> obj(element, i);
        obj(p, grads);
        
        for (unsigned int j=0; j < numberOfDegreesOfFreedom; ++j)
        {
            xDerReturnData(i,j) = element->basisFunction(j,p) * grads[0];
            yDerReturnData(i,j) = element->basisFunction(j,p) * grads[1];
            zDerReturnData(i,j) = element->basisFunction(j,p) * grads[2];
        }
    }
}

void
HEuler::faceIntegrand(const FaceT* face,          const PointPhysicalT& normal,
                   const PointReferenceOnTheFaceT& p,  FluxData& ret)
{
    if (face->isInternal())
    {
        const double magn                     = Utilities::norm2<DIM>(normal);
        unsigned int numberOfDegreesOfFreedom = face->getPtrElementLeft()->getNrOfBasisFunctions();
        
        PointReferenceT 	pL, pR;
        double              bFL, bFR, BFevalL, BFevalR;
        
        double theta = 0.5;//static_cast<const HEulerConfigurationData*>(configData_)->theta_;
        
        double nx    = normal[0]/magn;
        double ny    = normal[1]/magn;
        double nz    = normal[2]/magn;
        
        const ElementT* const left   = face->getPtrElementLeft();
        const ElementT* const right  = face->getPtrElementRight();
        
        face->mapRefFaceToRefElemL(p, pL);
        face->mapRefFaceToRefElemR(p, pR);
        
        
        for (int j = 0; j< numberOfDegreesOfFreedom; ++j)
        {
            bFL 	= 	left->basisFunction(j,pL);
            bFR 	=   right->basisFunction(j, pR);
            
            LinearAlgebra::Matrix& leftReturnData   = ret.left_[j];
            
            LinearAlgebra::Matrix& rightReturnData  = ret.right_[j];
            

            for (unsigned int i=0; i < numberOfDegreesOfFreedom; ++i)
            {
                
                BFevalL = left->basisFunction(i,pL);
                
                BFevalR = right->basisFunction(i,pR);
                
                leftReturnData(0,i) = 	BFevalL*nx*(1-theta)*bFL;
                leftReturnData(1,i) = 	BFevalL*ny*(1-theta)*bFL;
                leftReturnData(2,i) = 	BFevalL*nz*(1-theta)*bFL;
                leftReturnData(3,i) = 	BFevalR*nx*(theta)*bFL;
                leftReturnData(4,i) = 	BFevalR*ny*(theta)*bFL;
                leftReturnData(5,i) = 	BFevalR*nz*(theta)*bFL;
                
                leftReturnData(6,i) = 	bFL*nx*(1-theta)*BFevalL;
                leftReturnData(7,i) = 	bFL*ny*(1-theta)*BFevalL;
                leftReturnData(8,i) = 	bFL*nz*(1-theta)*BFevalL;
                leftReturnData(9,i) = 	bFL*nx*(theta)*BFevalR;
                leftReturnData(10,i) = 	bFL*ny*(theta)*BFevalR;
                leftReturnData(11,i) = 	bFL*nz*(theta)*BFevalR;
                
                
                rightReturnData(0,i) = 	BFevalL*nx*(1-theta)*bFR;
                rightReturnData(1,i) = 	BFevalL*ny*(1-theta)*bFR;
                rightReturnData(2,i) = 	BFevalL*nz*(1-theta)*bFR;
                rightReturnData(3,i) = 	BFevalR*nx*(theta)*bFR;
                rightReturnData(4,i) = 	BFevalR*ny*(theta)*bFR;
                rightReturnData(5,i) = 	BFevalR*nz*(theta)*bFR;
                
                rightReturnData(6,i) = 	bFR*nx*(1-theta)*BFevalL;
                rightReturnData(7,i) = 	bFR*ny*(1-theta)*BFevalL;
                rightReturnData(8,i) = 	bFR*nz*(1-theta)*BFevalL;
                rightReturnData(9,i) = 	bFR*nx*(theta)*BFevalR;
                rightReturnData(10,i) = bFR*ny*(theta)*BFevalR;
                rightReturnData(11,i) = bFR*nz*(theta)*BFevalR;
                
            }
        }
            //ret.print();
    }
    
}

void
HEuler::correctInitialProjectionOfVelocity(const Vec& UInit, Vec& UCorrected)const
{
	unsigned int nb 	= configData_->numberOfBasisFunctions_;
    unsigned int Nl     = static_cast<HEulerGlobalVariables*>(globalData_)->nElements_*nb;

	Mat& globalDIV      = static_cast<HEulerGlobalVariables*>(globalData_)->DivergenceFreeMatrix_;
	
	PetscInt a(10);
	PetscReal max;
    
	double reltol = 1.e-13;
	double abstol = 1.e-13;
	
    Vec RHS, MU;
	
	VecCreateSeq(PETSC_COMM_SELF, Nl, &RHS);
	VecCreateSeq(PETSC_COMM_SELF, 3*Nl, &UCorrected);
	
	MatMult(globalDIV, UInit, RHS);
        // 	outputVectorMatlab(RHS,"rhs.dat");
	
    
	VecScale(RHS,-1);
	VecAssemblyBegin(RHS);
	VecAssemblyEnd(RHS);
	VecMax(RHS, &a, &max);
	cout << "DIV U Error ="<<max<<endl;
	
	if (max > abstol)
	{
			//[A=DIV*DIV'/6;]
            //    [mu=gmres(A,b,10, 1.e-15, 1.e+8);]
		KSP ksp;
		
            // Preconditioner
		PC pc;
		KSPCreate(PETSC_COMM_SELF, &ksp);
		KSPSetOperators(ksp, globalDIV, globalDIV, SAME_NONZERO_PATTERN);
		cout << "2"<<endl;
        
            //	        KSPSetType(ksp, KSPPREONLY);
		KSPSetType(ksp, KSPLSQR);
		cout << "3"<<endl;
		KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
		cout << "4"<<endl;
		
		
		KSPGetPC(ksp,&pc);
		cout << "5"<<endl;
        
		
		cout << "6"<<endl;
		
		KSPSetFromOptions(ksp);
            //PCFactorSetLevels(pc, 1);
            //PCFactorSetFill(pc, 2.945);
		
		cout << "7"<<endl;
		
		PCSetType(pc, PCNONE);
		
        
		KSPSetTolerances(ksp, reltol, abstol, PETSC_DEFAULT, PETSC_DEFAULT);
		KSPSetUp(ksp);
		
		/*! Solving linear system.*/
		cout << "8"<<endl;
		KSPSolve(ksp, RHS, UCorrected);
            //KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
		
		KSPConvergedReason reason;
		KSPGetConvergedReason(ksp, &reason);
			//CHKERRQ(ierr);
		if(reason < 0)
		{
			PetscInt its;
			KSPGetIterationNumber(ksp, &its);
			PetscReal rnorm;
			KSPGetResidualNorm(ksp, &rnorm);
			cout << "\t\tPetSc: Solving the linear system has failed. Reason code: "
            << reason << endl << "Check KSPConvergedReason for the reason" << endl
            << "\t\tPetSc: Residual after " << int(its) << " iterations : ";
			cout.setf(ios::scientific, ios::floatfield);
			cout.precision(4);
			cout << rnorm << endl;
			cout.setf(ios::fixed, ios::floatfield);
			cout.precision(5);
		}
        
        else
		{
			if(reason > 0)
			{
				PetscInt its;
				KSPGetIterationNumber(ksp, &its);
				PetscReal rnorm;
				KSPGetResidualNorm(ksp, &rnorm);
				cout << "\t\tPetSc: Solving the linear system has succeeded. Reason code: "
                << reason << endl << "Check KSPConvergedReason for the reason" << endl
                << "\t\tPetsc: Convergence in " << int(its) << " iterations : ";
				cout.setf(ios::scientific, ios::floatfield);
				cout.precision(4);
				cout << rnorm << endl;
				cout.setf(ios::fixed, ios::floatfield);
				cout.precision(5);
			}
			else
			{
				cout << "\t\tPetSc: Solving the linear system is still under way" << endl;
			}
		}
        
        
        cout << "9"<<endl;
			//[u=velocity+DIV'*mu;]
        
            //  			outputVectorMatlab(UCorrected,"mu.dat");
        
        cout << "UC is ready!"<<endl;
        
        a	=10;
        max	=0;
        VecMax(UCorrected, &a, &max);
        cout << "Difference between UInit and UCorrected "<<max<<endl;
			//[maxE=max(u+velocity);]
        VecAXPY(UCorrected, 1, UInit);
            //   [error=DIV*u;]
        KSPDestroy(&ksp);
	}
	else
	{
		VecCopy(UInit,UCorrected); 
	}
	
    cout << "*************************************************************************"<<endl;
    
    VecDestroy(&RHS);
}

void
HEuler::calculatePressure(const Mat& A, const Mat& Ah,const Vec& UCorrected)
{
	unsigned int nb 	= configData_->numberOfBasisFunctions_;
    unsigned int Nl     = static_cast<HEulerGlobalVariables*>(globalData_)->nElements_*nb;
    
    unsigned int Nu = Nl;
    unsigned int Nv = Nl;
    unsigned int Nw = Nl;
    
	Mat& globalDIV      = static_cast<HEulerGlobalVariables*>(globalData_)->DivergenceFreeMatrix_;
    
    unsigned int N      = Nu+Nv+Nw;
	
	Vec Lambda, RHS;
	
	VecCreateSeq(PETSC_COMM_SELF, Nl, &RHS);
	VecCreateSeq(PETSC_COMM_SELF, Nl, &Lambda);
	
	
	MatMult(Ah, UCorrected, RHS);
	PC pcL;
	KSP ksp;
	KSPCreate(PETSC_COMM_SELF, &ksp);
    
//    outputMatrix(A, "A.dat");
//    outputVectorMatlab(RHS, "RHS.dat");
    
	KSPSetOperators(ksp, A, A, SAME_NONZERO_PATTERN);
	cout << "2L"<<endl;
	
        //KSPSetType(ksp, KSPPREONLY);
	KSPSetType(ksp, KSPGMRES);
	KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
	cout << "3L"<<endl;
	
	
	cout << "4L"<<endl;
    
    
	KSPGetPC(ksp,&pcL);
	cout << "5L"<<endl;
	
	
	cout << "6L"<<endl;
	
	
	KSPSetFromOptions(ksp);
	
	
	
        // 	PCFactorSetFill(pcL, 6.45759);
        // 	PCFactorSetLevels(pcL,1.18);
    
	
        //PCSetType(pcL, PCNONE);
	
    cout << "7L"<<endl;
	double reltolL = 1.e-13;
	double abstolL = 1.e-13;
	KSPSetTolerances(ksp, reltolL, abstolL, PETSC_DEFAULT, PETSC_DEFAULT);
	KSPSetUp(ksp);
	
	/*! Solving linear system.*/
	cout << "8L"<<endl;
        //KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
    KSPSolve(ksp, RHS, Lambda);
    
	KSPConvergedReason reason;
	KSPGetConvergedReason(ksp, &reason);
	if(reason < 0)
	{
		PetscInt its;
		KSPGetIterationNumber(ksp, &its);
		PetscReal rnorm;
		KSPGetResidualNorm(ksp, &rnorm);
		cout << "\t\tPetSc: Solving the linear system has failed. Reason code: "
        << reason << endl << "Check KSPConvergedReason for the reason" << endl
        << "\t\tPetSc: Residual after " << int(its) << " iterations : ";
		cout.setf(ios::scientific, ios::floatfield);
		cout.precision(4);
		cout << rnorm << endl;
		cout.setf(ios::fixed, ios::floatfield);
		cout.precision(5);
	}
    
	else
	{
		if(reason > 0)
		{
			PetscInt its;
			KSPGetIterationNumber(ksp, &its);
			PetscReal rnorm;
			KSPGetResidualNorm(ksp, &rnorm);
			cout << "\t\tPetSc: Solving the linear system has succeeded. Reason code: "
            << reason << endl << "Check KSPConvergedReason for the reason" << endl
            << "\t\tPetsc: Convergence in " << int(its) << " iterations : ";
			cout.setf(ios::scientific, ios::floatfield);
			cout.precision(4);
			cout << rnorm << endl;
			cout.setf(ios::fixed, ios::floatfield);
			cout.precision(5);
		}
		else
		{
			cout << "\t\tPetSc: Solving the linear system is still under way" << endl;
		}
	}
        //*********************************************************************************************************************************
        //*************************************Calculate lambda directly END******************************************************************
	
	PetscScalar* XTEMP = new PetscScalar [N];
	PetscScalar* XTEMP1 = new PetscScalar [Nl];
	VecGetArray(UCorrected, &XTEMP);
    
//    outputVectorMatlab(Lambda, "Lambda.txt");
 	VecGetArray(Lambda, &XTEMP1);
	int pos=0;
	int pos1=0;
    
    const HEulerConfigurationData* config = static_cast<const HEulerConfigurationData*>(configData_);
   
    
    for (ConstElementIterator cit = elementColBegin(); cit != elementColEnd(); ++cit)
    {
        Base::Element<DIM>* element= *cit;
        int k = element->getID();
        
        unsigned int pos = k*nb;
        
        for (int i =0; i<nb;++i)
        {
                //cout << "p="<<XTEMP1[pos+i]<<endl;
            element->setData(0,3,i, -XTEMP1[pos+i]);//set lambda
            
            
            element->setData(0,0,i, XTEMP[pos+i]);//set U
            
            element->setData(0,1,i, XTEMP[Nu+pos+i]);//set V
            
            element->setData(0,2,i, XTEMP[Nu+Nv+pos+i]);//set W
            
            
        }
    }
    
            //update the energy.
            //const Element<dim>& el = *itM;
		
    
            //ExactEnergyWrapper  integrandEx(gv, el, data, 0, init);
            //integrateOverElement(el, integrandEx, energyOverElement);
            //element.energyExact_= energyOverElement[0];
        
//		EnergyWrapper  integrand(gv, el, data);
//		integrateOverElement(el, integrand, energyOverElement);
//		element.setEnergy(energyOverElement[0]);
//        
//		L2ErrorNorm  error(gv, el, data, 0, exact);
//		integrateOverElement(el, error, errorVector);
//		element.setUError(errorVector[0]);
//		element.setVError(errorVector[1]);
//		element.setWError(errorVector[2]);
//		element.setPError(errorVector[3]);
	
	VecRestoreArray(UCorrected, &XTEMP);
 	VecRestoreArray(Lambda, &XTEMP1);
	
	KSPDestroy(&ksp);
	VecDestroy(&Lambda);
	VecDestroy(&RHS);
	
}


void
HEuler::createIncompressibleSystem()
{
 	PetscInfoAllow(PETSC_TRUE, "history.txt");
        // 	DumpElementsFaces(mesh,data);return;
	
	Mat C;
	Mat BF;
    
	Mat C1;
	Mat BF1;
	Mat MInv;
	Mat sMInv;
  	Mat Ml;
	Mat A;
	Mat Ah;
    Mat& globalDIV      = static_cast<HEulerGlobalVariables*>(globalData_)->DivergenceFreeMatrix_;
    
    unsigned int nb 	= configData_->numberOfBasisFunctions_;
    unsigned int Nu     = static_cast<HEulerGlobalVariables*>(globalData_)->nElements_*nb;
    unsigned int Nv     = Nu;
    unsigned int Nw     = Nu;
    unsigned int Nl     = Nu;
    unsigned int N      = Nu+Nv+Nw;
    
	
    double dt           = static_cast<HEulerGlobalVariables*>(globalData_)->dt_;
    
    cout << "dt=" <<dt<<endl;
    
	Vec 												UInit, RHS, UCorrected;
	PetscViewer 										viewer;
	
	MatInfo info;
    
	double omega1		= 0.0;
    double omega2		= 0.0;
    double omega3		= 1.0;
    double theta 		= 0.5; //later on read from config file or something else
	
    double u		    = 0;
	double v            = 0;
	double w            = 0;
	double lambda		= 0;
	double massij		= 0;
	
	
    
    

//	cout<< "nElements_="<<gv.nElements_<<endl;
//	cout<< "N="<<N<<", nb="<<nb<<endl;
	
	cout << "**************Starting create Matrices*************"<<endl;
	MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nu+Nv+Nw, 	Nu+Nv+Nw,    1*nb, 				PETSC_NULL, &MInv);
	MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nl, 		Nl,          1*nb, 				PETSC_NULL, &sMInv);
 	MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nl, 		Nl,          1*nb, 				PETSC_NULL, &Ml);
    
	    //MatCreateSeqAIJ(PETSC_COMM_SELF, 	N, 		N,    	1*nb, 			PETSC_NULL, &M);
        //MatCreateSeqAIJ(PETSC_COMM_SELF, 	N1, 	N1,    	(3*7)*nb, 			PETSC_NULL, &A);
        //MatCreateSeqAIJ(PETSC_COMM_SELF, 	N1, 	3*N,    (3*7+7)*nb, 			PETSC_NULL, &Ah);
	
	
	MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nu+Nv+Nw, 	Nu+Nv+Nw,    12*nb, 			PETSC_NULL, &C);
	
    
	MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nu+Nv+Nw, 	Nl,         (7)*nb, 			PETSC_NULL, &BF);//number of possible nonzero blocks are 7: element and his 6 neighbours
	
	MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nl, 	Nu+Nv+Nw,    (3*7)*nb, 			PETSC_NULL, &globalDIV);//number of possible nonzero blocks are 7: element and his 6 neighbours
    
	VecCreateSeq(PETSC_COMM_SELF, Nu+Nv+Nw, &UInit);
	VecCreateSeq(PETSC_COMM_SELF, Nu+Nv+Nw, &UCorrected);
	VecCreateSeq(PETSC_COMM_SELF, Nl, &RHS);
//	VecCreateSeq(PETSC_COMM_SELF, Nl, &gv.LambdaConstraint_);
	
	
	
	
	
	cout << "**************Mass Matrix Calculated*************"<<endl;
	double intMassMatrix;
//	timer t0;
//	std::cout << t0.elapsed() << " - : Element Integration started\n";
    ElementIntegralData gradMass;
    bool useCache = false;
    ElementIntegralT   elIntegral(useCache);
    typedef void  (HEuler::*Integrand)(const ElementT* , const PointReferenceT&, ElementIntegralData&);
    Integrand gradMassInteg = &HEuler::elementIntegrand;
    
    
    std::cout << " - : Element Integration started\n";
    
    for (ConstElementIterator cit = elementColBegin(); cit != elementColEnd(); ++cit)
    {

        int m  = (*cit)->getID();
        unsigned int pos = (*cit)->getID()*nb;
        
        HEulerElementData* elementData = static_cast<HEulerElementData*>((*cit)->getUserData());
        
        ElementT* element =(*cit);
        
        elIntegral.integrate(element, gradMassInteg, gradMass, this);
		
		for (unsigned int i = 0; i<nb;++i)
 		{
			u = element->getData(0,0,i);
				
			VecSetValue(UInit, pos+i,     u, INSERT_VALUES);
			
             v = element->getData(0,1,i);
			VecSetValue(UInit, Nu+pos+i,   v, INSERT_VALUES);
			
             w = element->getData(0,2,i);;
			VecSetValue(UInit, Nu+Nv+pos+i, w, INSERT_VALUES);
			
            lambda	= element->getData(0,3,i);;
            
//            VecSetValue(Lambda, pos+i, lambda, INSERT_VALUES);
//
//        
//			Lambda phi_i(i);
//			integrateOverElement(*itM, phi_i, tempInteg);
//			
//			VecSetValue(gv.LambdaConstraint_, pos1, tempInteg[0], INSERT_VALUES);///int_/Omega{Lambda}=0 constraint
		}
        
		for (unsigned int j = 0; j<nb;++j)
		{
			for (unsigned int i=0; i<nb;++i)//can be optimized later!!!
 			{
                    //integrateOverElement(*itM, scalarProduct<dim>(ExpansionServer::basisFunction(i),ExpansionServer::basisFunction(j)), 					massVij=calculatePhiIDotPhiJForUWithV(i,j, elementM, itM);
                massij=elementData->massMatrix_(j,i);
                MatSetValue(C, pos+j,      Nu+pos+i,            massij*(omega3),  ADD_VALUES);//V coefficient
//                MatSetValue(C, pos+j,       Nu+Nv+pos+i,   		massij*(-omega2), ADD_VALUES);   //W coefficient  omega2=0 done for memory shortage
                
                MatSetValue(C, Nu+pos+j,   pos+i,               massij*(-omega3), ADD_VALUES);//U coefficient
//                MatSetValue(C, Nu+pos+j,    Nu+Nv+pos+i,        massij*(omega1),  ADD_VALUES);//W coefficient  omega1=0 done for memory shortage
//                MatSetValue(C, Nu+Nv+pos+j, pos+i,              massij*(omega2),  ADD_VALUES);//U coefficient    omega2=0 done for memory shortage
//                MatSetValue(C, Nu+Nv+pos+j, Nu+pos+i,           massij*(-omega1), ADD_VALUES);   //V coefficient omega1=0 done for memory shortage
                
                MatSetValue(Ml,    pos+j, 		pos+i,     massij , ADD_VALUES);
            
                MatSetValue(MInv, pos+j,       pos+i,    	elementData->invMassMatrix_(j,i) , ADD_VALUES);
                MatSetValue(MInv, Nu+pos+j,    Nu+pos+i,    elementData->invMassMatrix_(j,i) , ADD_VALUES);
                MatSetValue(MInv, Nu+Nv+pos+j, Nu+Nv+pos+i, elementData->invMassMatrix_(j,i) , ADD_VALUES);
                
                MatSetValue(sMInv, pos+j, pos+i,           elementData->invMassMatrix_(j,i) , ADD_VALUES);
                
                MatSetValue(BF, pos+j, 		 pos+i,   	-gradMass.xGrad_(i,j),  ADD_VALUES);//U coefficient
                MatSetValue(BF, Nu+pos+j, 	 pos+i,   	-gradMass.yGrad_(i,j),  ADD_VALUES);//V coefficient
                MatSetValue(BF, Nu+Nv+pos+j, pos+i,   	-gradMass.zGrad_(i,j),  ADD_VALUES);//W coefficient
                
                
                MatSetValue(globalDIV, pos+j, pos+i,   	    gradMass.xGrad_(j,i),  ADD_VALUES);//U coefficient
                MatSetValue(globalDIV, pos+j, Nu+pos+i,   	gradMass.yGrad_(j,i),  ADD_VALUES);//V coefficient
                MatSetValue(globalDIV, pos+j, Nu+Nv+pos+i,	gradMass.zGrad_(j,i),  ADD_VALUES);//W coefficient
                                                                                          
 			}
		}
    }
	

        //        outputMatrix(globalDIV, "globalDIV.dat");
        //        outputMatrix(BF, "BF.dat");
    
    printFullMatrixInfo(globalDIV, "globalDIV");
    printFullMatrixInfo(BF, "BF");
    
    cout << "Mat Assembly for "<< "MINV"<<endl;
    
    MatAssemblyBegin(MInv,MAT_FINAL_ASSEMBLY);
        //CHKERRQ(ierr);
    
    MatAssemblyEnd(MInv,MAT_FINAL_ASSEMBLY);
    cout << "Mat Assembly for "<< "MINVSm"<<endl;
        ////
    MatAssemblyBegin(sMInv,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(sMInv,MAT_FINAL_ASSEMBLY);
    
    unsigned int posR=0;
    unsigned int posL=0;
    unsigned int eR=0;
    unsigned int eL=0;
    
    std::cout << " - : Face Integration started\n";
        //
        //
    FluxData fData(nb);
    typedef void  (HEuler::*FaceIntegrand)(const FaceT*, const PointPhysicalT& normal , const PointReferenceOnTheFaceT&, FluxData&);
    FaceIntegrand faceInteg = &HEuler::faceIntegrand;
    FaceIntegralT   faceIntegral(useCache);
    

	for (ConstFaceIterator citFe = faceColBegin(); citFe != faceColEnd(); ++citFe)
    {
        
        if ((*citFe)->getPtrElementRight()== NULL) // boundary face
        {
        }
        else
        {
            eR = (*citFe)->getPtrElementRight()->getID();
            eL = (*citFe)->getPtrElementLeft()->getID();
            
            posR  = eR*nb;
            
            posL  = eL*nb;
            
            faceIntegral.integrate((*citFe), faceInteg, fData, this);
            
            for (unsigned int j=0;j<nb;++j)
            {
                for (unsigned int i=0;i<nb;++i)
                {
                    MatSetValue(BF, posR+j, 		posL+i,     fData.right_[j](0,i),  	ADD_VALUES);//U coefficient
                    MatSetValue(BF, posR+j, 		posR+i,    -fData.right_[j](3,i),  	ADD_VALUES);//U coefficient
                    MatSetValue(BF, Nu+posR+j,      posL+i,     fData.right_[j](1,i),  	ADD_VALUES);//V coefficient
                    MatSetValue(BF, Nu+posR+j,      posR+i,    -fData.right_[j](4,i),  	ADD_VALUES);//V coefficient
                    MatSetValue(BF, posL+j, 		posL+i,   	fData.left_[j](0,i),  	ADD_VALUES);//U coefficient
                    MatSetValue(BF, posL+j, 		posR+i,    -fData.left_[j](3,i),  	ADD_VALUES);//U coefficient
                    MatSetValue(BF, Nu+posL+j,      posL+i,   	fData.left_[j](1,i),  	ADD_VALUES);//V coefficient
                    MatSetValue(BF, Nu+posL+j,      posR+i,    -fData.left_[j](4,i),  	ADD_VALUES);//V coefficient
                    MatSetValue(BF, Nu+Nv+posL+j, 	posR+i,    -fData.left_[j](5,i),  	ADD_VALUES);//W coefficient
                    MatSetValue(BF, Nu+Nv+posL+j, 	posL+i,   	fData.left_[j](2,i),  	ADD_VALUES);//W coefficient
                    MatSetValue(BF, Nu+Nv+posR+j, 	posR+i,    -fData.right_[j](5,i),  	ADD_VALUES);//W coefficient
                    MatSetValue(BF, Nu+Nv+posR+j, 	posL+i,     fData.right_[j](2,i),  	ADD_VALUES);//W coefficient
                    
                    MatSetValue(globalDIV, posR+j, posL+i,    	  fData.right_[j](6,i),  	ADD_VALUES);//U coefficient
                    MatSetValue(globalDIV, posL+j, posL+i,   	 -fData.left_[j](6,i),  	ADD_VALUES);//U coefficient
                    MatSetValue(globalDIV, posL+j, Nu+posL+i,    -fData.left_[j](7,i),  	ADD_VALUES);//V coefficient
                    MatSetValue(globalDIV, posR+j, Nu+posL+i,     fData.right_[j](7,i),  	ADD_VALUES);//V coefficient
                    MatSetValue(globalDIV, posR+j, Nu+posR+i,     fData.right_[j](10,i),  	ADD_VALUES);//V coefficient
                    MatSetValue(globalDIV, posL+j, Nu+posR+i,    -fData.left_[j](10,i),  	ADD_VALUES);//V coefficient
                    MatSetValue(globalDIV, posL+j, posR+i,   	 -fData.left_[j](9,i),  	ADD_VALUES);//U coefficient
                    MatSetValue(globalDIV, posR+j, posR+i,    	  fData.right_[j](9,i),  	ADD_VALUES);//U coefficient
                    MatSetValue(globalDIV, posL+j, Nu+Nv+posL+i, -fData.left_[j](8,i),  	ADD_VALUES);//W coefficient
                    MatSetValue(globalDIV, posR+j, Nu+Nv+posL+i,  fData.right_[j](8,i),  	ADD_VALUES);//W coefficient
                    MatSetValue(globalDIV, posL+j, Nu+Nv+posR+i, -fData.left_[j](11,i),  	ADD_VALUES);//W coefficient
                    MatSetValue(globalDIV, posR+j, Nu+Nv+posR+i,  fData.right_[j](11,i),  	ADD_VALUES);//W coefficient
                }
            }
                //                cout <<"***********************************"<<endl;
        }
    }
    
    MatAssemblyBegin(globalDIV,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(globalDIV,MAT_FINAL_ASSEMBLY);
    
    MatAssemblyBegin(MInv,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(MInv,MAT_FINAL_ASSEMBLY);
    
    MatAssemblyBegin(BF,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(BF,MAT_FINAL_ASSEMBLY);
    
    
//    outputMatrix(globalDIV, "globalDIV.dat");
//    outputMatrix(BF, "BF.dat");
    
    printFullMatrixInfo(globalDIV, "globalDIV");
    printFullMatrixInfo(BF, "BF");
    
    printFullMatrixInfo(MInv, "MInv");

	
// 	outputMatrix(gv.DivergenceFreeMatrix_, "DIV.dat");
        //  	outputMatrixMatlab(BF, "BF.dat");
        // 	outputMatrixMatlab(C, "C.dat");
        // 	outputMatrixMatlab(MInv, "MInv.dat");
    
        // 	ierr = MatAssemblyBegin(MInv,MAT_FINAL_ASSEMBLY);
        //     ierr = MatAssemblyEnd(MInv,MAT_FINAL_ASSEMBLY);
	
	cout << "Mat Assembly for "<< "C"<<endl;
 	MatAssemblyBegin(Ml,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Ml,MAT_FINAL_ASSEMBLY);
	
    MatAssemblyBegin(C,MAT_FINAL_ASSEMBLY);
   	MatAssemblyEnd(C,MAT_FINAL_ASSEMBLY);
	
	cout << "Mat Assembly for "<< "BF"<<endl;
	MatAssemblyBegin(BF,MAT_FINAL_ASSEMBLY);
   	MatAssemblyEnd(BF,MAT_FINAL_ASSEMBLY);
	
	
	printFullMatrixInfo(MInv, "MInv");
        //********************DIV=[Ex+F4u, Ey+F4v, Ez+F4w]*******
    cout << "Mat Assembly for "<< "DivergenceFreeMatrix_"<<endl;
	
    MatAssemblyBegin(globalDIV, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(globalDIV, MAT_FINAL_ASSEMBLY);
	cout << "Mat Assembly for "<< "blubla"<<endl;
        //
	
        //*******************************************************
	cout << "Mat Assembly for "<< "before"<<endl;
	
	std::cout << " - : Create rotational matrix!\n";
	cout << "Mat Assembly for "<< "after"<<endl;
	
	double fillC=1;//(double)(3*N*2*nb)/(3*N*2*nb +3*N*nb);
	cout << "Mat Assembly for "<< "touch C"<<endl;
	
  	MatMatMult(MInv, C, MAT_INITIAL_MATRIX, 1., &C1);
	
    cout << "Mat Assembly for "<< "after touch C"<<endl;
	std::cout << " - : Create rotational matrix1!\n";
	MatDestroy(&C);
	printFullMatrixInfo(C1, "C1");
	double fillBF=1.0;//(double)(3*N*7*nb)/(3*N*7*nb +3*N*nb);
	MatMatMult(MInv, BF, MAT_INITIAL_MATRIX, fillBF, &BF1);
    
        //	MatMatMult(sMInv, gv.DivergenceFreeMatrix_, MAT_INITIAL_MATRIX, fillBF, &gv.DivergenceFreeMatrix_);
        //MatMatMult(BF1, Ml, MAT_INITIAL_MATRIX, fillBF, &BF1);
	
	std::cout<< " - : Create rotational matrix2!\n";
	MatDestroy(&BF);
	printFullMatrixInfo(BF1, "BF1");
 	
 	
    
        //MatMatMult(BF, M, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &BF);
    
	MatDestroy(&MInv);
	MatDestroy(&Ml);
    
	MatMatMult(globalDIV, BF1, MAT_INITIAL_MATRIX, 1, &A);
	std::cout << " - : Create rotational matrix3!\n";
	MatMatMult(globalDIV, C1, MAT_INITIAL_MATRIX, 1, &Ah);
	std::cout << " - : Finished rotational matrix4!\n";
	
	printFullMatrixInfo(Ah, "Ah");
	printFullMatrixInfo(A, "A");
        //outputMatrix(BF, "BF.txt");
	cout << "Mat Assembly for "<< "C1"<<endl;
	MatAssemblyBegin(C1,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(C1,MAT_FINAL_ASSEMBLY);
	
	cout << "Mat Assembly for "<< "BF1"<<endl;
	MatAssemblyBegin(BF1,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(BF1,MAT_FINAL_ASSEMBLY);
	
	cout << "Mat Assembly for "<< "A"<<endl;
	MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
	
    
	cout << "Mat Assembly for "<< "Ah"<<endl;
	MatAssemblyBegin(Ah,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(Ah,MAT_FINAL_ASSEMBLY);
	
	
	printFullMatrixInfo(Ah, "Ah");
	printFullMatrixInfo(A, "A");
	
        //outputMatrix(Ah, "Ah,dat");
        //outputMatrix(A, "A.dat");
	
	printFullMatrixInfo(C1, "C1");
	printFullMatrixInfo(BF1, "BF1");
	
    
        //outputMatrix(Ah, "Ah.txt");
	
    correctInitialProjectionOfVelocity(UInit, UCorrected);
        //UCorrected=UInit;
	cout << "Correction is done!"<<endl;
	
	MatMult(globalDIV, UCorrected, RHS);
	
    PetscInt a(10);
    PetscReal max;
    VecMax(RHS, &a, &max);
	
	cout << "Discrete divergemceNEW="<<max<<endl;

 	double dummy=0;
	const PetscInt* cols;
	const PetscScalar* values;
	PetscInt 	numberOfNonZeros;
	std::cout<< " - : Started Creating P matrix!\n";
	
	PetscInt* nnzP= new PetscInt[Nu+Nv+Nw+Nl];
	PetscInt* nnzQ= new PetscInt[Nu+Nv+Nw+Nl];
   
    for (int i =0;i<Nu;++i)
    {
        nnzP[i]=(3*nb+1);
        nnzQ[i]=(3*nb+1);
        nnzP[Nu+i]=(3*nb+1);
        nnzQ[Nu+i]=(3*nb+1);
        nnzP[Nu+Nv+i]=(3*nb+1);
        nnzQ[Nu+Nv+i]=(3*nb+1);
        nnzP[Nu+Nv+Nw+i]=(21*nb+1);
        nnzQ[Nu+Nv+Nw+i]=(21*nb+1);
    }
//        // 	nnzP[Nu+Nv+Nw+Nl]=Nl;
//        // 	nnzQ[Nu+Nv+Nw+Nl]=Nl;
//	
	std::cout<< " stage1"<<endl;
        // 	MatCreateSeqAIJ(PETSC_COMM_SELF, Nu+Nv+Nw+Nl+1, Nu+Nv+Nw+Nl, PETSC_NULL, nnzP,&gv.P_);
    
	MatCreateSeqAIJ(PETSC_COMM_SELF, Nu+Nv+Nw+Nl, Nu+Nv+Nw+Nl, PETSC_NULL, nnzP,&P_);
    
    std::cout<< " stage2"<<endl;
	cout <<endl<<Nu+Nv+Nw<<endl;
	MatCreateSeqAIJ(PETSC_COMM_SELF, Nu+Nv+Nw+Nl, Nu+Nv+Nw, PETSC_NULL, nnzQ,&Q_);
        // 	MatCreateSeqAIJ(PETSC_COMM_SELF, Nu+Nv+Nw+Nl+1, Nu+Nv+Nw, PETSC_NULL, nnzQ,&gv.Q_);
	delete [] nnzP;
	delete [] nnzQ;
    std::cout << " stage3"<<endl;
    
    std::cout << " stage4"<<endl;
    cout <<endl<<Nu+Nv+Nw<<endl;
    
    printFullMatrixInfo(P_, "P_");
    printFullMatrixInfo(Q_, "Q_");
    
    for (unsigned int i =0;i<N;++i)
    {
        MatSetValue(P_, i,     i,     1, ADD_VALUES);//creating 2*I matrix in the upper part
        
        MatSetValue(Q_, i,     i,     1, ADD_VALUES);//creating 2*I matrix in the upper part
    }
    
    printFullMatrixInfo(P_, "P_");
    printFullMatrixInfo(Q_, "Q_");
    
    cout <<endl<<N<<endl;
//	
//        // 	PetscScalar* lambdaCons;
//        //
//        // 	VecGetArray(gv.LambdaConstraint_, &lambdaCons);
//        //
//        // 	for (int i =0; i<Nl;++i)
//        // 	{
//        // 		MatSetValue(gv.P_, 	N+Nl,  N+i,    lambdaCons[i], 	ADD_VALUES);
//        // 		MatSetValue(gv.Q_, 	N+Nl,  N+i,    0, 	ADD_VALUES);
//        // 	}
//        // 	VecRestoreArray(gv.LambdaConstraint_, &lambdaCons);
//	cout <<"asdasdasd"<<endl;
//
    numberOfNonZeros=0;
	for (unsigned int i=0;i<N;++i)
	{
        
            //add extra lambda equation
		
		
        MatGetRow(C1, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0)
            {
                MatSetValue(P_, 	i,     cols[j],     -0.5*dt*dummy, 	ADD_VALUES);
                MatSetValue(Q_, 	i,     cols[j],     0.5*dt*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(C1, i, &numberOfNonZeros, &cols, &values);
        
        MatGetRow(BF1, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0)
            {
                MatSetValue(P_, 	i,     N+cols[j],     -dt*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(BF1, i, &numberOfNonZeros, &cols, &values);
        if (i < Nl)
        {
            MatGetRow(Ah, i, &numberOfNonZeros, &cols, &values);
            for (int j=0;j<numberOfNonZeros;++j)
            {
                dummy = (values[j]);
                if (dummy!=0)
                {
                    MatSetValue(P_, 	N+i,     cols[j],     0.5*dummy, 		ADD_VALUES);
                    MatSetValue(Q_, 	N+i,     cols[j],     -0.5*dummy, 	ADD_VALUES);
                }
            }
            MatRestoreRow(Ah, i, &numberOfNonZeros, &cols, &values);
        }
        if (i < Nl)
        {
            
            MatGetRow(A, i, &numberOfNonZeros, &cols, &values);
            for (int j=0;j<numberOfNonZeros;++j)
            {
                dummy = (values[j]);
                
                
                
                if (dummy!=0)
                {
//                    cout << "N+i"<<N+i<<",N+cols[j]="<<N+cols[j]<<endl;

                    MatSetValue(P_, 	N+i,     	N+cols[j],     dummy, 		ADD_VALUES);
                }
            }
            MatRestoreRow(A, i, &numberOfNonZeros, &cols, &values);
        }
        
    }
    cout << "Mat Assembly for "<< "P"<<endl;
    MatAssemblyBegin(P_,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(P_,MAT_FINAL_ASSEMBLY);
    
    cout << "Mat Assembly for "<< "Q"<<endl;
    MatAssemblyBegin(Q_,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Q_,MAT_FINAL_ASSEMBLY);
    
    printFullMatrixInfo(P_, "P_");
    printFullMatrixInfo(Q_, "Q_");
    
    
        //outputMatrix(P_, "P.dat");
        //outputMatrix(Q_, "Q.dat");//
    std::cout<< " - : Finished Creating matrix P !\n";
        //*********************************************************************************************************************************
        //*************************************Calculate lambda directly BEGIN******************************************************************
    
        //outputMatrix(gv.DivergenceFreeMatrix_, "DIV.txt");
        //outputMatrix(gv.P_, "P.txt");
        //outputMatrix(gv.Q_, "Q.txt");
    
    calculatePressure(A, Ah, UCorrected);
        //*********************************************************************************************************************** 
        //    	 outputMatrixMatlab(gv.P_, "P.dat");
	
	MatDestroy(&BF1);
	MatDestroy(&A);
	MatDestroy(&Ah);
	MatDestroy(&C1);
	
	cout<<"12"<<endl;
	VecDestroy(&UInit);
	VecDestroy(&UCorrected);
	VecDestroy(&RHS);
    cout << "**************END of contructing matrices*************"<<endl;
}

void
HEuler::createCompressibleSystem()
{
    PetscInfoAllow(PETSC_TRUE, "history.txt");
        // 	DumpElementsFaces(mesh,data);return;
    
    Mat C;
    Mat BF;
    
    double dt = static_cast<HEulerGlobalVariables*>(globalData_)->dt_;
    
    cout << "dt=" <<dt<<endl;
    
    Mat BF1;
    Mat MInv;
    Mat MInvSm;
    Mat DIV;
    Mat A;
    Mat Ah;
    Mat& globalDIV =  static_cast<HEulerGlobalVariables*>(globalData_)->DivergenceFreeMatrix_;
    
    unsigned int nb 	= configData_->numberOfBasisFunctions_;
    unsigned int Nu     = static_cast<HEulerGlobalVariables*>(globalData_)->nElements_*nb;
    unsigned int Nv     = Nu;
    unsigned int Nw     = Nu;
    unsigned int Nl     = Nu;
    unsigned int N      = Nu+Nv+Nw;
    
    
    
    MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nu+Nv+Nw, 	Nu+Nv+Nw,    	1*nb, 				PETSC_NULL, &MInv);
    
    MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nu,         Nu,             1*nb, 				PETSC_NULL, &MInvSm);
    
    cout << "**************Starting create Matrices*************"<<endl;
    
    
        //        MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nu+Nv+Nw, 	Nu+Nv+Nw,    nb, 				PETSC_NULL, &C);
    
    
    MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nu+Nv+Nw, 	Nl,     (7)*nb, 			PETSC_NULL, &BF);//number of possible nonzero blocks are 7: element and his 6 neighbours
    
    MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nl, 	Nu+Nv+Nw,    (3*7)*nb, 			PETSC_NULL, &globalDIV);//number of possible nonzero blocks are 7: element and his 6 neighbours
                                                                                                     
    
    
    printFullMatrixInfo(BF, "BF");
    
    
    
    cout << "**************Mass Matrix Calculated*************"<<endl;
    
    cout.flush();
    
    
    
    ElementIntegralData gradMass;
    bool useCache = false;
    ElementIntegralT   elIntegral(useCache);
    typedef void  (HEuler::*Integrand)(const ElementT* , const PointReferenceT&, ElementIntegralData&);
    Integrand gradMassInteg = &HEuler::elementIntegrand;
    
    
    std::cout << " - : Element Integration started\n";
    
    for (ConstElementIterator cit = elementColBegin(); cit != elementColEnd(); ++cit)
    {
        
        int m  = (*cit)->getID();
        unsigned int pos = (*cit)->getID()*nb;
        
        HEulerElementData* elementData = static_cast<HEulerElementData*>((*cit)->getUserData());

        ElementT* elem =(*cit);
        elIntegral.integrate(elem, gradMassInteg, gradMass, this);
        
        for (unsigned int j = 0; j<nb;++j)
        {
            for (unsigned int i=0; i<nb;++i)//can be optimized later!!!
            {
                MatSetValue(MInv, pos+j,       pos+i,    		elementData->invMassMatrix_(j,i) , ADD_VALUES);
                MatSetValue(MInv, Nu+pos+j,    Nu+pos+i,    	elementData->invMassMatrix_(j,i) , ADD_VALUES);
                MatSetValue(MInv, Nu+Nv+pos+j, Nu+Nv+pos+i,    elementData->invMassMatrix_(j,i) , ADD_VALUES);
                
                MatSetValue(MInvSm, pos+j, pos+i,         elementData->invMassMatrix_(j,i) , ADD_VALUES);
               
                MatSetValue(BF, pos+j, 		 pos+i,   	-gradMass.xGrad_(i,j),  ADD_VALUES);//U coefficient
                MatSetValue(BF, Nu+pos+j, 	 pos+i,   	-gradMass.yGrad_(i,j),  ADD_VALUES);//V coefficient
                MatSetValue(BF, Nu+Nv+pos+j, pos+i,   	-gradMass.zGrad_(i,j),  ADD_VALUES);//W coefficient
                
                
                MatSetValue(globalDIV, pos+j, pos+i,   	    gradMass.xGrad_(j,i),  ADD_VALUES);//U coefficient
                MatSetValue(globalDIV, pos+j, Nu+pos+i,   	gradMass.yGrad_(j,i),  ADD_VALUES);//V coefficient
                MatSetValue(globalDIV, pos+j, Nu+Nv+pos+i,	gradMass.zGrad_(j,i),  ADD_VALUES);//W coefficient
            }
        }
    }
    
    
//        outputMatrix(globalDIV, "globalDIV.dat");
//        outputMatrix(BF, "BF.dat");

    printFullMatrixInfo(globalDIV, "globalDIV");
    printFullMatrixInfo(BF, "BF");
    
    cout << "Mat Assembly for "<< "MINV"<<endl;
    
    MatAssemblyBegin(MInv,MAT_FINAL_ASSEMBLY);
        //CHKERRQ(ierr);
    
    MatAssemblyEnd(MInv,MAT_FINAL_ASSEMBLY);
    cout << "Mat Assembly for "<< "MINVSm"<<endl;
        ////
    MatAssemblyBegin(MInvSm,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(MInvSm,MAT_FINAL_ASSEMBLY);

    unsigned int posR=0;
    unsigned int posL=0;
    unsigned int eR=0;
    unsigned int eL=0;

    std::cout << " - : Face Integration started\n";
        //
        //
    FluxData fData(nb);
    typedef void  (HEuler::*FaceIntegrand)(const FaceT*, const PointPhysicalT& normal , const PointReferenceOnTheFaceT&, FluxData&);
    FaceIntegrand faceInteg = &HEuler::faceIntegrand;
    FaceIntegralT   faceIntegral(useCache);
    
    for (ConstFaceIterator citFe = faceColBegin(); citFe != faceColEnd(); ++citFe)
    {
        
        if ((*citFe)->getPtrElementRight()== NULL) // boundary face
        {
        }
        else
        {
            eR = (*citFe)->getPtrElementRight()->getID();
            eL = (*citFe)->getPtrElementLeft()->getID();
            
            posR  = eR*nb;
            
            posL  = eL*nb;
            
            faceIntegral.integrate((*citFe), faceInteg, fData, this);

            for (unsigned int j=0;j<nb;++j)
            {
                for (unsigned int i=0;i<nb;++i)
                {
                    MatSetValue(BF, posR+j, 		posL+i,     fData.right_[j](0,i),  	ADD_VALUES);//U coefficient
                    MatSetValue(BF, posR+j, 		posR+i,    -fData.right_[j](3,i),  	ADD_VALUES);//U coefficient
                    MatSetValue(BF, Nu+posR+j,      posL+i,     fData.right_[j](1,i),  	ADD_VALUES);//V coefficient
                    MatSetValue(BF, Nu+posR+j,      posR+i,    -fData.right_[j](4,i),  	ADD_VALUES);//V coefficient
                    MatSetValue(BF, posL+j, 		posL+i,   	fData.left_[j](0,i),  	ADD_VALUES);//U coefficient
                    MatSetValue(BF, posL+j, 		posR+i,    -fData.left_[j](3,i),  	ADD_VALUES);//U coefficient
                    MatSetValue(BF, Nu+posL+j,      posL+i,   	fData.left_[j](1,i),  	ADD_VALUES);//V coefficient
                    MatSetValue(BF, Nu+posL+j,      posR+i,    -fData.left_[j](4,i),  	ADD_VALUES);//V coefficient
                    MatSetValue(BF, Nu+Nv+posL+j, 	posR+i,    -fData.left_[j](5,i),  	ADD_VALUES);//W coefficient
                    MatSetValue(BF, Nu+Nv+posL+j, 	posL+i,   	fData.left_[j](2,i),  	ADD_VALUES);//W coefficient
                    MatSetValue(BF, Nu+Nv+posR+j, 	posR+i,    -fData.right_[j](5,i),  	ADD_VALUES);//W coefficient
                    MatSetValue(BF, Nu+Nv+posR+j, 	posL+i,     fData.right_[j](2,i),  	ADD_VALUES);//W coefficient
                    
                    MatSetValue(globalDIV, posR+j, posL+i,    	  fData.right_[j](6,i),  	ADD_VALUES);//U coefficient
                    MatSetValue(globalDIV, posL+j, posL+i,   	 -fData.left_[j](6,i),  	ADD_VALUES);//U coefficient
                    MatSetValue(globalDIV, posL+j, Nu+posL+i,    -fData.left_[j](7,i),  	ADD_VALUES);//V coefficient
                    MatSetValue(globalDIV, posR+j, Nu+posL+i,     fData.right_[j](7,i),  	ADD_VALUES);//V coefficient
                    MatSetValue(globalDIV, posR+j, Nu+posR+i,     fData.right_[j](10,i),  	ADD_VALUES);//V coefficient
                    MatSetValue(globalDIV, posL+j, Nu+posR+i,    -fData.left_[j](10,i),  	ADD_VALUES);//V coefficient
                    MatSetValue(globalDIV, posL+j, posR+i,   	 -fData.left_[j](9,i),  	ADD_VALUES);//U coefficient
                    MatSetValue(globalDIV, posR+j, posR+i,    	  fData.right_[j](9,i),  	ADD_VALUES);//U coefficient
                    MatSetValue(globalDIV, posL+j, Nu+Nv+posL+i, -fData.left_[j](8,i),  	ADD_VALUES);//W coefficient
                    MatSetValue(globalDIV, posR+j, Nu+Nv+posL+i,  fData.right_[j](8,i),  	ADD_VALUES);//W coefficient
                    MatSetValue(globalDIV, posL+j, Nu+Nv+posR+i, -fData.left_[j](11,i),  	ADD_VALUES);//W coefficient
                    MatSetValue(globalDIV, posR+j, Nu+Nv+posR+i,  fData.right_[j](11,i),  	ADD_VALUES);//W coefficient
                }
            }
                //                cout <<"***********************************"<<endl;
        }
    }
    
    
    MatAssemblyBegin(globalDIV,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(globalDIV,MAT_FINAL_ASSEMBLY);
    
    MatAssemblyBegin(MInv,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(MInv,MAT_FINAL_ASSEMBLY);
   
    MatAssemblyBegin(BF,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(BF,MAT_FINAL_ASSEMBLY);
    
    
    
    
        //outputMatrix(globalDIV, "globalDIV.dat");
        //outputMatrix(BF, "BF.dat");
    
    printFullMatrixInfo(globalDIV, "globalDIV");
    printFullMatrixInfo(BF, "BF");
    
    printFullMatrixInfo(MInv, "MInv");
    
        //********************DIV=[Ex+F4u, Ey+F4v, Ez+F4w]*******
    cout << "Mat Assembly for "<< "DivergenceFreeMatrix_"<<endl;
    MatAssemblyBegin(globalDIV,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(globalDIV,MAT_FINAL_ASSEMBLY);
    
        //*******************************************************
    std::cout << " - : Create rotational matrix!\n";
 
    double fillBF=1;//(double)(3*N*7*nb)/(3*N*7*nb +3*N*nb);
    MatMatMult(MInv, BF, MAT_INITIAL_MATRIX, fillBF, &BF1);
    
    std::cout << " - : Create rotational matrix2!\n";
    MatDestroy(&BF);
    printFullMatrixInfo(BF1, "BF1");
    PetscErrorCode ierr;
    
    
    double fillDIV=1;//(double)(3*N*7*nb)/(3*N*7*nb +3*N*nb);
    MatMatMult(MInvSm, globalDIV, MAT_INITIAL_MATRIX, fillBF, &DIV);
    
    std::cout << " - : Create rotational matrix3!\n";
   
    MatDestroy(&globalDIV);
    printFullMatrixInfo(DIV, "DIV");
    
    MatDestroy(&MInv);
    MatDestroy(&MInvSm);
    
    
    cout << "Mat Assembly for "<< "C1"<<endl;
 
    cout << "Mat Assembly for "<< "BF1"<<endl;
    ierr = MatAssemblyBegin(BF1,MAT_FINAL_ASSEMBLY);
    ierr = MatAssemblyEnd(BF1,MAT_FINAL_ASSEMBLY);
    
    cout << "Mat Assembly for "<< "DIV"<<endl;
    ierr = MatAssemblyBegin(DIV,MAT_FINAL_ASSEMBLY);
    ierr = MatAssemblyEnd(DIV,MAT_FINAL_ASSEMBLY);
    
    
    
    printFullMatrixInfo(DIV, "DIV");
        //        printFullMatrixInfo(C1, "C1");
    printFullMatrixInfo(BF1, "BF1");
    

    
    double dummy=0;
    const PetscInt* cols;
    const PetscScalar* values;
    PetscInt 	numberOfNonZeros;
    std::cout << " - : Started Creating P matrix!\n";
    PetscInt* nnzP= new PetscInt[Nu+Nv+Nw+Nl];
    PetscInt* nnzQ= new PetscInt[Nu+Nv+Nw+Nl];

    
    int noFV=5;
    for (int i =0;i<Nu;++i)
    {
        nnzP[i]=(3*nb+1);
        nnzQ[i]=(3*nb+1);
        nnzP[Nu+i]=(3*nb+1);
        nnzQ[Nu+i]=(3*nb+1);
        nnzP[Nu+Nv+i]=(3*nb+1);
        nnzQ[Nu+Nv+i]=(3*nb+1);
        nnzP[Nu+Nv+Nw+i]=(9*nb+1);
        nnzQ[Nu+Nv+Nw+i]=(9*nb+1);
    }

    std::cout << " stage1"<<endl;
    
    MatCreateSeqAIJ(PETSC_COMM_SELF, Nu+Nv+Nw+Nl, Nu+Nv+Nw+Nl, PETSC_NULL, nnzP,&P_);
    std::cout << " stage2"<<endl;
    cout <<endl<<Nu+Nv+Nw<<endl;
    
    MatCreateSeqAIJ(PETSC_COMM_SELF, Nu+Nv+Nw+Nl, Nu+Nv+Nw+Nl, PETSC_NULL, nnzQ,&Q_);
    delete [] nnzP;
    delete [] nnzQ;
    
    std::cout << " stage3"<<endl;
    
    std::cout << " stage4"<<endl;
    cout <<endl<<Nu+Nv+Nw<<endl;
    printFullMatrixInfo(P_, "P_");
    printFullMatrixInfo(Q_, "Q_");
    for (unsigned int i =0;i<N+Nl;++i)
    {
        MatSetValue(P_, i,     i,     1, ADD_VALUES);//creating 2*I matrix in the upper part
        
        MatSetValue(Q_, i,     i,     1, ADD_VALUES);//creating 2*I matrix in the upper part
    }
    printFullMatrixInfo(P_, "P_");
    printFullMatrixInfo(Q_, "Q_");
    
    cout <<endl<<N<<endl;
    
    for (unsigned int i=0;i<N;++i)
    {
        MatGetRow(BF1, i, &numberOfNonZeros, &cols, &values);
        for (int j=0;j<numberOfNonZeros;++j)
        {
            dummy = (values[j]);
            if (dummy!=0.0)
            {
                 MatSetValue(P_, 	i,     N+cols[j],    -0.5*dt*dummy, 	ADD_VALUES);
                MatSetValue(Q_, 	i,     N+cols[j],     0.5*dt*dummy, 	ADD_VALUES);
            }
        }
        MatRestoreRow(BF1, i, &numberOfNonZeros, &cols, &values);
        
        if (i<Nl)
        {
            MatGetRow(DIV, i, &numberOfNonZeros, &cols, &values);
            for (int j=0;j<numberOfNonZeros;++j)
            {
                dummy = (values[j]);
                if (dummy!=0.0)
                {
                    MatSetValue(P_, 	N+i,     cols[j],     -0.5*dt*dummy, 		ADD_VALUES);
                    MatSetValue(Q_, 	N+i,     cols[j],     0.5*dt*dummy, 	ADD_VALUES);
                }
            }
            MatRestoreRow(DIV, i, &numberOfNonZeros, &cols, &values);
        }
    }
    cout << "Mat Assembly for "<< "DIV"<<endl;
    ierr = MatAssemblyBegin(P_,MAT_FINAL_ASSEMBLY);
    ierr = MatAssemblyEnd(P_,MAT_FINAL_ASSEMBLY);
    
    cout << "Mat Assembly for "<< "DIV"<<endl;
    ierr = MatAssemblyBegin(Q_,MAT_FINAL_ASSEMBLY);
    ierr = MatAssemblyEnd(Q_,MAT_FINAL_ASSEMBLY);
    
    printFullMatrixInfo(P_, "P_");
    printFullMatrixInfo(Q_, "Q_");
    
    
        //outputMatrix(P_, "P.dat");
        //    outputMatrix(Q_, "Q.dat");
    
        //
}

    /// create Mass Matrices, store them as a User Defined Element Data
    /// calculate projection of the every unknown on the FEM spaces.
void
HEuler::initialConditions()
{
    unsigned int ldof = static_cast<const HEulerConfigurationData*>(configData_)->numberOfBasisFunctions_;
    
    HEulerElementData* elemData;
    
    bool useCache = false;
    
    ElementIntegralT   elIntegral(useCache);
    
    typedef void  (HEuler::*Integrand)(const ElementT* , const PointReferenceT&, LinearAlgebra::Matrix&);
    Integrand massMatrixIntegrand = &HEuler::calculateMassMatrix;
    
    LinearAlgebra::NumericalVector rightHand(ldof);
    
    
    ElementT::SolutionVector       numerical(ldof);
    
    LinearAlgebra::Matrix          invMassM(ldof, ldof);
    
    cout << "ldof="<<ldof<<endl;
    
    InitCondU       uEx(exactSolution_);
    InitCondV       vEx(exactSolution_);
    InitCondW       wEx(exactSolution_);
    InitCondLambda  pOrRhoEx(exactSolution_);
    
    cout <<"start calculations of inital condition!"<<endl;
    unsigned int count =0;
    for (ConstElementIterator cit = elementColBegin(); cit != elementColEnd(); ++cit)
    {
        ElementT* elem =(*cit);
        
        elemData = new HEulerElementData(ldof);
        
        LinearAlgebra::Matrix& massMatrix = elemData->massMatrix_;
        LinearAlgebra::Matrix& invMassM   = elemData->invMassMatrix_;
        
        elIntegral.integrate(elem, massMatrixIntegrand, massMatrix, this);
     
        massMatrix.inverse(invMassM);
        
        elem->setUserData(elemData);
        
        elIntegral.integrate(elem, uEx, rightHand);
        
        numerical = invMassM*rightHand;// projection of U
        
        
        elem->setTimeLevelData(0,0,numerical);
        
        elIntegral.integrate(elem, vEx, rightHand);
        
        numerical = invMassM*rightHand;// projection of V
        
        elem->setTimeLevelData(0,1,numerical);
        
        
        elIntegral.integrate(elem, wEx, rightHand);
        
        numerical = invMassM*rightHand;// projection of W
        elem->setTimeLevelData(0,2,numerical);
        
        elIntegral.integrate(elem, pOrRhoEx, rightHand);
        
        numerical = invMassM*rightHand;// projection of P
        
            //cout << "P="<<numerical<<endl;
        elem->setTimeLevelData(0,3,numerical);
    }
    
    cout <<"finish calculations of inital condition!"<<endl;
    
}

void
HEuler::solve()
{
    const HEulerConfigurationData* config = static_cast<const HEulerConfigurationData*>(configData_);
    HEulerGlobalVariables* globalData = static_cast<HEulerGlobalVariables*>(globalData_);
    
    
    unsigned int nb 	= configData_->numberOfBasisFunctions_;
    unsigned int Nu     =  static_cast<HEulerGlobalVariables*>(globalData_)->nElements_*nb;
    unsigned int Nv     = Nu;
    unsigned int Nw     = Nu;
    unsigned int Nl     = Nu;
    unsigned int N      = Nu+Nv+Nw;
    
    double endTime      = config->numOfPeriods_*config->onePeriod_ ;
    int nplotPerPeriod  = config->numOfPeriodsInOnePlotStep_;
    
    
    cout << "endTime="<<endTime<<endl;
    Vec UExact, Lambda, RHS, RH;
    
    VecCreateSeq(PETSC_COMM_SELF, N+Nl, &Lambda);
    VecCreateSeq(PETSC_COMM_SELF, N+Nl, &RHS);
    
    
    if(config->solutionType_==HEulerConfigurationData::COMPRESSIBLE_WALLS || config->solutionType_==HEulerConfigurationData::COMPRESSIBLE_PERIODIC)
    {
        VecCreateSeq(PETSC_COMM_SELF, N+Nl, &UExact);
    }
    else
    {
        VecCreateSeq(PETSC_COMM_SELF, N, &UExact);
    }
    
    VecCreateSeq(PETSC_COMM_SELF, Nl, &RH);
    
    double u, v, w, uEx, vEx, wEx, pres, presLap, rh,l;
    
    cout<<"Initializing the solver with velocity..."<<endl;
    
    for (ConstElementIterator cit = elementColBegin(); cit != elementColEnd(); ++cit)
    {
        Base::Element<DIM>* element= *cit;
        int m = element->getID();
        
        unsigned int pos = m*nb;
        
        
        for (int i =0; i<nb;++i)
        {
            u = element->getData(0,0,i);
            VecSetValue(Lambda, pos+i,     u, INSERT_VALUES);
            VecSetValue(UExact, pos+i,     u, INSERT_VALUES);
            
            v = element->getData(0,1,i);;
            VecSetValue(Lambda, Nu+pos+i,   v, INSERT_VALUES);
            VecSetValue(UExact, Nu+pos+i,   v, INSERT_VALUES);
            
            w = element->getData(0,2,i);;
            VecSetValue(UExact, Nu+Nv+pos+i, w, INSERT_VALUES);
            VecSetValue(Lambda, Nu+Nv+pos+i, w, INSERT_VALUES);
            
            l = element->getData(0,3,i);;
            if(config->solutionType_==HEulerConfigurationData::COMPRESSIBLE_WALLS ||config->solutionType_==HEulerConfigurationData::COMPRESSIBLE_PERIODIC)
                VecSetValue(UExact, Nu+Nv+Nw+pos+i, l, INSERT_VALUES);
            VecSetValue(Lambda, Nu+Nv+Nw+pos+i, l, INSERT_VALUES);
        }
    }
    VecAssemblyBegin(UExact);
    VecAssemblyEnd(UExact);
    
    double currTime=0;
    
    int iterations=1;
    int iterations1=0;
    PetscScalar* XTEMP = new PetscScalar [Nu+Nv+Nw+Nl];
    
    
    KSP ksp;
        // Preconditioner
    PC pc;
    cout << "Solve"<<endl;
        // Create a solver
    KSPCreate(PETSC_COMM_SELF, &ksp);
    
    cout << "ksp create"<<endl;
    printFullMatrixInfo(P_, "P_");
    printFullMatrixInfo(Q_, "Q_");
    
        //! Direct solver: only a preconditioner is needed. To change
        //this to an iterative method, you can change the solver to
        //KSPMINRES, KSPGMRES, etc. Similarly, you can choose PCNONE
        //(instead of PCLU) for the preconditioner below if you want to
        //use a nonpreconditioned iterative method.
    
    KSPSetOperators(ksp, P_, P_, SAME_NONZERO_PATTERN);
    
    cout << "Setup operators"<<endl;
    
    KSPSetType(ksp, KSPGMRES);
    cout << "Setup solver"<<endl;
    
        // 	KSPSetType(ksp, KSPLSQR);
    KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
    
    cout << "Setup Initial guess"<<endl;
    
    KSPGetPC(ksp,&pc);
    
    cout << "Setup pc"<<endl;
    
    
    
    PCSetType(pc, PCNONE);
    
        PCSetType(pc, PCILU);
        //PCSetType(pc, PCCHOLESKY);
        //PCILUSetFill(pc, 4);
    
    KSPSetFromOptions(ksp);
        //	PCFactorSetFill(pc, 3.83016);
        //	PCFactorSetLevels(pc,2);
    
    cout << "Setup options"<<endl;
    
    
    double reltol = 1.0e-14;
    double abstol = 1.0e-14;
    KSPSetTolerances(ksp, reltol, abstol, PETSC_DEFAULT, PETSC_DEFAULT);
    cout << "Setup tolerance"<<endl;
    
    
    
    KSPSetUp(ksp);
    cout << "Setup ksp"<<endl;
    
    
    
    
        //! Solving linear system.
    
    while(currTime<=endTime)
    {
            //Lambda                         integrandLambda(&Lambda);
        
        cout<<"construct proper rightHandside"<<endl;
        
            ///##########################construct proper rightHandside#####################
            ///***********************************************************************************
        MatMult(Q_, UExact, RHS);
        
            ///***********************************************************************************
            ///***************************construct proper rightHandside**************************
        
        
        cout << "Finalizing vector creation"<<endl;
        
        
            // KSP Solver initializer
            // 	    outputVectorMatlab(RHS, "RHS.dat");
        KSPSolve(ksp, RHS, Lambda);
            //KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
        
        KSPConvergedReason reason;
        KSPGetConvergedReason(ksp, &reason);
        if(reason < 0)
        {
            PetscInt its;
            KSPGetIterationNumber(ksp, &its);
            PetscReal rnorm;
            KSPGetResidualNorm(ksp, &rnorm);
            cout << "\t\tPetSc: Solving the linear system has failed. Reason code: "
            << reason << endl << "Check KSPConvergedReason for the reason" << endl
            << "\t\tPetSc: Residual after " << int(its) << " iterations : ";
            cout.setf(ios::scientific, ios::floatfield);
            cout.precision(4);
            cout << rnorm << endl;
            cout.setf(ios::fixed, ios::floatfield);
            cout.precision(5);
        }
        
        else
        {
            if(reason > 0)
            {
                PetscInt its;
                KSPGetIterationNumber(ksp, &its);
                PetscReal rnorm;
                KSPGetResidualNorm(ksp, &rnorm);
                cout << "\t\tPetSc: Solving the linear system has succeeded. Reason code: "
                << reason << endl << "Check KSPConvergedReason for the reason" << endl
                << "\t\tPetsc: Convergence in " << int(its) << " iterations : ";
                cout.setf(ios::scientific, ios::floatfield);
                cout.precision(4);
                cout << rnorm << endl;
                cout.setf(ios::fixed, ios::floatfield);
                cout.precision(5);
            }
            else
            {
                cout << "\t\tPetSc: Solving the linear system is still under way" << endl;
            }
        }
        
        cout << "Solved a timestep" << endl;
        
        
        VecGetArray(Lambda, &XTEMP);
        int pos=0;
        int k=0;
        currTime+=globalData->dt_;
        for (ConstElementIterator cit = elementColBegin(); cit != elementColEnd(); ++cit)
        {
            Base::Element<DIM>* element= *cit;
            int k = element->getID();
            
            unsigned int pos = k*nb;
            
            for (int i =0; i<nb;++i)
            {
                element->setData(0,3,i, XTEMP[Nu+Nv+Nw+pos+i]);//set rho
                
                if(config->solutionType_==HEulerConfigurationData::COMPRESSIBLE_WALLS ||config->solutionType_==HEulerConfigurationData::COMPRESSIBLE_PERIODIC)
                {
                    VecSetValue(UExact, Nu+Nv+Nw+pos+i, XTEMP[Nu+Nv+Nw+pos+i], INSERT_VALUES);
                }
                element->setData(0,0,i, XTEMP[pos+i]);//set U
                
                VecSetValue(UExact, pos+i,     		XTEMP[pos+i], INSERT_VALUES);
                
                element->setData(0,1,i, XTEMP[Nu+pos+i]);//set V
                
                VecSetValue(UExact, Nu+pos+i,   	XTEMP[Nu+pos+i], INSERT_VALUES);
                
                element->setData(0,2,i, XTEMP[Nu+Nv+pos+i]);//set W
                
                VecSetValue(UExact, Nu+Nv+pos+i, 	XTEMP[Nu+Nv+pos+i], INSERT_VALUES);
            }
        }
        
        VecRestoreArray(Lambda, &XTEMP);
        VecAssemblyBegin(UExact);
        VecAssemblyEnd(UExact);
        if (iterations%nplotPerPeriod==0)
            output(currTime);
        
        cout<<"currTime="<<currTime<<endl;
        iterations++;
    }
    
    KSPDestroy(&ksp);
    VecDestroy(&UExact);
    VecDestroy(&RHS);
    VecDestroy(&Lambda);
    VecDestroy(&RH);
    
    cout << "The end of the time loop"<<endl;
    
}

void
HEuler::output(double time)
{
    
    string dxName  = "x.dat";
    string dName   = "u.dat";
    string dExName = "uEx.dat";
    char outF[20];
    ofstream    	l2ErrorFile("l2error.dat");
    ofstream    	energyFile("energy.dat");
    ofstream    	divFile("div.dat");
    ofstream    	energyExfile("energyEx.dat");
    
    string          outputFormat("u,uExact, uError, v, vExact, vError, w, wExact, wError, Lambda,LambdaExact, LambdaError, Energy");
    
    sprintf(outF,"outputTec%f.dat",time);
    
    string          outFile(outF);
    
    ofstream    file1Dx0(dxName.c_str());
    ofstream    file1D0(dName.c_str());
    ofstream    file1DEx0(dExName.c_str());
    std::ofstream file3D;
    file3D.open (outFile.c_str());
    
    Output::TecplotDiscontinuousSolutionWriter<DIM> out(file3D,"RectangularMesh","012",outputFormat);
    
    
    TecplotWriteFunction userWriteFunction(&energyFile, &divFile, &energyExfile, exactSolution_, &file1Dx0, &file1D0,&file1DEx0, &l2ErrorFile);
    userWriteFunction.time_=time;
    
    cout <<"OUTPUTING...."<<endl;
    ostringstream ostr;
    ostr << "t = " << time;
    out.write(meshes_[0],std::string(ostr.str()),false, userWriteFunction);
    cout <<"FINISH OUTPUTING...."<<endl;
 
    file1Dx0.close();
    file1D0.close();
    file1DEx0.close();
    divFile.close();
    energyFile.close();
    energyExfile.close();
    l2ErrorFile.close();
    
    file3D.close();

    cout <<"FINISH...."<<endl;
}