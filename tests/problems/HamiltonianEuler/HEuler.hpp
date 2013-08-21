#include "cstring"
using std::string;

#include <petscmat.h>
#include <petscksp.h>


#include "Output/TecplotDiscontinuousSolutionWriter.hpp"
#include "Output/TecplotPhysicalGeometryIterator.hpp"
#include "LinearAlgebra/NumericalVector.hpp"
#include "Base/GlobalData.hpp"
#include "Base/ConfigurationData.hpp"
#include "Base/HpgemUI.hpp"
#include "Base/Norm2.hpp"
#include "Base/PhysGradientOfBasisFunction.hpp"
#include "InitialConditions.hpp"
#include "TecplotOutput.hpp"
using Base::RectangularMeshDescriptor;
using Base::HpgemUI;
using Base::GlobalData;
using Base::ConfigurationData;


const unsigned int  DIM = 3;
typedef std::vector<LinearAlgebra::Matrix>   VectorOfMatrices;
struct ElementIntegralData
{
        //optimize later!
    ElementIntegralData operator*= (const double& scalar){xGrad_*=scalar; yGrad_*=scalar; zGrad_*=scalar; return *this;}
    void axpy(double a, const ElementIntegralData& x){ xGrad_.axpy(a, x.xGrad_); yGrad_.axpy(a, x.yGrad_); zGrad_.axpy(a, x.zGrad_);}
    
    LinearAlgebra::Matrix xGrad_;
    LinearAlgebra::Matrix yGrad_;
    LinearAlgebra::Matrix zGrad_;
};

struct FluxData
{
    FluxData(unsigned int nb)
    {
        left_.resize(nb);
        right_.resize(nb);
        
        for (unsigned int i=0; i<nb; ++i)
        {
            LinearAlgebra::Matrix& left= left_[i];
            LinearAlgebra::Matrix& right= right_[i];
            left.resize(12,nb);
            right.resize(12,nb);
        }
    }
    void print()
    {
        cout <<"left="<<endl;
        for (unsigned int n=0; n < left_.size();++n)
        {
            cout <<"n="<<n<<", "<<left_[n]<<endl;
        }
        
        cout <<"right="<<endl;
        for (unsigned int n=0; n < right_.size();++n)
        {
            cout <<"n="<<n<<", "<<right_[n]<<endl;
        }
        
    }
    FluxData operator*= (const double& scalar)
    {
        for (unsigned int n=0; n < left_.size();++n)
        {
            left_[n]*=scalar;
            right_[n]*=scalar;
        }
        return *this;
    }
    void axpy(double a, const FluxData& x)
    {
        for (unsigned int n=0; n < left_.size();++n)
        {
            left_[n].axpy(a, x.left_[n]);
            right_[n].axpy(a, x.right_[n]);
        }
    }
    
    VectorOfMatrices    left_;
    VectorOfMatrices    right_;
};



struct HEulerElementData: public UserElementData
{
    HEulerElementData(unsigned int ndof):
        massMatrix_(ndof, ndof),
        invMassMatrix_(ndof, ndof)
    {
    }
    
    LinearAlgebra::Matrix   massMatrix_;
    LinearAlgebra::Matrix   invMassMatrix_;
};

struct HEulerGlobalVariables: public GlobalData
{
    unsigned int    nElements_;
    Mat             DivergenceFreeMatrix_;
    
    double          dt_;
};



struct HEulerConfigurationData: public ConfigurationData
{
    enum  SolutionType 		{INCOMPRESSIBLE_WALLS, INCOMPRESSIBLE_ONETHIRDPERIODIC, INCOMPRESSIBLE_PERIODIC, COMPRESSIBLE_WALLS, COMPRESSIBLE_PERIODIC};

    HEulerConfigurationData(unsigned int numberOfUnknowns, unsigned int numberOfBasisFunctions, unsigned int  numberOfTimeLevels=1, const string& fileName="in.txt", SolutionType type=COMPRESSIBLE_PERIODIC):
        ConfigurationData(numberOfUnknowns, numberOfBasisFunctions, numberOfTimeLevels=1),
        solutionType_(type)
    {
            /// reading from a file
        theta_=0.5;
        numOfPeriods_=1;
        numOfTimeStepInOnePeriod_=10;
        numOfPeriodsInOnePlotStep_=1;
        onePeriod_=1;
        
    }
    
public:
    SolutionType    solutionType_;
    
    unsigned int    nx_;
    unsigned int    ny_;
    unsigned int    nz_;
    
    double          lx_;
    double          ly_;
    double          lz_;
    
    double          theta_;
    
    double          numOfPeriods_;
    double          numOfTimeStepInOnePeriod_;
    double          numOfPeriodsInOnePlotStep_;
    double          onePeriod_;
};

class HEuler: public HpgemUI<DIM>
{
public:
    typedef HpgemUI<DIM>                        HpgemUIT;
    typedef Integration::ElementIntegral<DIM>   ElementIntegralT;
    typedef Integration::FaceIntegral<DIM>      FaceIntegralT;
    typedef ExactSolutionBase<DIM>              ExactSolutionT;
    typedef PointReference<DIM-1>               PointReferenceOnTheFaceT;
    
        //using HpgemUIT::ElementT;
        //    using HpgemUIT::PointReferenceT;
    
public:
    HEuler(HEulerGlobalVariables* global, const HEulerConfigurationData* config):
        HpgemUI(global, config),
        P_(),
        Q_()
    {
        exactSolution_= new Compressible3DPeriodic();
    }
    
    ~HEuler()
    {
        
        MatDestroy(&P_);
        MatDestroy(&Q_);
    }
    
public:
    
    void outputMatrix(Mat& matrix, const string& name)
    {
        cout << "Mat Assembly for "<< name<<endl;
        MatAssemblyBegin(matrix,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(matrix,MAT_FINAL_ASSEMBLY);
        
        PetscViewer viewer;
        PetscViewerASCIIOpen(PETSC_COMM_SELF,name.c_str(),&viewer);
        MatView(matrix,viewer);
    }
    
    void outputVectorMatlab(Vec& vec, const string& name)
    {
        cout << "Vec Assembly for "<< name<<endl;
        VecAssemblyBegin(vec);
        VecAssemblyEnd(vec);
        
        PetscViewer viewer;
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, name.c_str(), FILE_MODE_WRITE, &viewer);
        VecView(vec, viewer);
    }
    
    void printFullMatrixInfo(Mat& matrix, const string& name)
    {
        PetscInt m=0,n=0;
        MatGetSize(matrix,&m,&n);
        
        MatInfo info;
        MatGetInfo(matrix,MAT_LOCAL, &info);
        cout<<name<<endl;
        printf("N = %d, N = %d, block_size = %d, memory = %d, assemblies = %d, mallocs = %d, matrix nonzeros (SeqBAIJ format) = %d, allocated nonzeros= %d\n", m, n, (int)info.block_size, (int)info.memory, (int)info.assemblies, (int)info.mallocs,(int)info.nz_used,(int)info.nz_allocated);
    }
    
    bool initialiseMesh()
    {
        RectangularMeshDescriptor<DIM> rectangularMesh;
       
        rectangularMesh.bottomLeft_[0]       = 0;
        rectangularMesh.bottomLeft_[1]       = 0;
        rectangularMesh.bottomLeft_[2]       = 0;
        rectangularMesh.topLeft_[0]          = 1;
        rectangularMesh.topLeft_[1]          = 1;
        rectangularMesh.topLeft_[2]          = 1;
        rectangularMesh.numElementsInDIM_[0] = 10;
        rectangularMesh.numElementsInDIM_[1] = 10;
        rectangularMesh.numElementsInDIM_[2] = 10;
        
        MeshId id = addMesh(rectangularMesh);
        
        HEulerGlobalVariables* globalData = static_cast<HEulerGlobalVariables*>(globalData_);
        
        globalData->nElements_ = getNumberOfElements(id);
        
        const HEulerConfigurationData* config = static_cast<const HEulerConfigurationData*>(configData_);
        
        globalData->dt_         = config->onePeriod_/config->numOfTimeStepInOnePeriod_;


        return true;
    }

    void calculateMassMatrix(const ElementT* element, const PointReferenceT& p, LinearAlgebra::Matrix& massMatrix)
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
        
//        cout <<"*******************"<<endl;
    }
    
    void calculateLocalEnergy(const ElementT& element, const PointReferenceT& p, double& returnValue)
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
    
    void elementIntegrand(const ElementT* element, const PointReferenceT& p, ElementIntegralData& returnObject)
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
    
    void faceIntegrand(const FaceT* face,          const PointPhysicalT& normal,
                       const PointReferenceOnTheFaceT& p,  FluxData& ret)
    {
        if (face->isInternal())
		{
            const double magn                     = Utilities::norm2<DIM>(normal);
            unsigned int numberOfDegreesOfFreedom = face->getPtrElementLeft()->getNrOfBasisFunctions();
            
            PointReferenceT 	pL, pR;
            double              bFL, bFR, BFevalL, BFevalR;
            
            double theta = 0.5;//static_cast<const HEulerConfigurationData*>(configData_)->theta_;
            
//            cout << "theta="<<theta<<endl;
//cout <<"normal="<<normal<<endl;
            double nx    = normal[0]/magn;
            double ny    = normal[1]/magn;
            double nz    = normal[2]/magn;
            
            const ElementT* const left   = face->getPtrElementLeft();
            const ElementT* const right  = face->getPtrElementRight();

            face->mapRefFaceToRefElemL(p, pL);
            face->mapRefFaceToRefElemR(p, pR);
            
//            cout << "p="<<p<<endl;
//            cout << "pL="<<pL<<endl;
//            
     		
			for (int j = 0; j< numberOfDegreesOfFreedom; ++j)
			{
                bFL 	= 	left->basisFunction(j,pL);
                bFR 	=   right->basisFunction(j, pR);
                
                LinearAlgebra::Matrix& leftReturnData   = ret.left_[j];
                
                LinearAlgebra::Matrix& rightReturnData  = ret.right_[j];

//                cout <<"leftReturnData(11,i)"<<leftReturnData<<endl;
//                
//                cout <<"rightReturnData(11,i)"<<rightReturnData<<endl;
//                
        //                /// provide numerical fluxes for
                for (unsigned int i=0; i < numberOfDegreesOfFreedom; ++i)
                {
                    
                    BFevalL = left->basisFunction(i,pL);
//                    cout << "BFevalL="<<BFevalL<<endl;
//                    
//                    cout << "bFL="<<bFL<<endl;
//                    
//                    cout << "nx="<<nx<<endl;
//                    
                    
                    
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
    void createTheCompressibleSystem()
    {
        PetscInfoAllow(PETSC_TRUE, "history.txt");
            // 	DumpElementsFaces(mesh,data);return;
        
        Mat C;
        Mat BF;
        
        double dt = static_cast<HEulerGlobalVariables*>(globalData_)->dt_;
        
        cout << "dt=" <<dt<<endl;
        
//        Mat C1;
        Mat BF1;
        Mat MInv;
        Mat MInvSm;
        Mat DIV;
        Mat A;
        Mat Ah;
        Mat& globalDIV =  static_cast<HEulerGlobalVariables*>(globalData_)->DivergenceFreeMatrix_;
        unsigned int nb 	= configData_->numberOfBasisFunctions_;
        unsigned int Nu     =  static_cast<HEulerGlobalVariables*>(globalData_)->nElements_*nb;
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
//        MatCreateSeqAIJ(PETSC_COMM_SELF, 	Nl, 	Nu+Nv+Nw,    (3*7)*nb, 			PETSC_NULL, &DIV);//number of possible nonzero blocks are 7: element and his 6 neighbours
        
        
//        VecCreateSeq(PETSC_COMM_SELF, Nu+Nv+Nw, &UInit);
//        VecCreateSeq(PETSC_COMM_SELF, Nl, &RHS);
        
        
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
////
////            for (unsigned int i = 0; i<nb;++i)
////            {
////                u 		= elementM.getU(i);
////                VecSetValue(UInit, pos+i,     u, INSERT_VALUES);
////                v 		= elementM.getV(i);
////                VecSetValue(UInit, Nu+pos+i,   v, INSERT_VALUES);
////                w 		= elementM.getW(i);
////                VecSetValue(UInit, Nu+Nv+pos+i, w, INSERT_VALUES);
////                lambda	= elementM.getL(i);
////                //VecSetValue(Lambda, pos1+i, lambda, INSERT_VALUES);
////            }
////
//           
            ElementT* elem =(*cit);
            elIntegral.integrate(elem, gradMassInteg, gradMass, this);

            for (unsigned int j = 0; j<nb;++j)
            {
                for (unsigned int i=0; i<nb;++i)//can be optimized later!!!
                {
                        MatSetValue(MInv, pos+j,       pos+i,    		elementData->invMassMatrix_(j,i) , ADD_VALUES);
                        MatSetValue(MInv, Nu+pos+j,    Nu+pos+i,    	elementData->invMassMatrix_(j,i) , ADD_VALUES);
                        MatSetValue(MInv, Nu+Nv+pos+j, Nu+Nv+pos+i,    elementData->invMassMatrix_(j,i) , ADD_VALUES);

                        //massVij=calculatePhiIDotPhiJForUWithV(i,j, elementM, itM);
                        //MatSetValue(C, posu+j, 	Nu+posv+i,   		massVij*(omega3),  ADD_VALUES);//V coefficient

                        //massUij=calculatePhiIDotPhiJForUWithV(j, i, elementM, itM);
                        //MatSetValue(C, Nu+posv+j, posu+i,     		massUij*(-omega3), ADD_VALUES);//U coefficient

                  
                    
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
//        
        printFullMatrixInfo(globalDIV, "globalDIV");
        printFullMatrixInfo(BF, "BF");

        cout << "Mat Assembly for "<< "MINV"<<endl;

        MatAssemblyBegin(MInv,MAT_FINAL_ASSEMBLY);
            //CHKERRQ(ierr);
        MatAssemblyEnd(MInv,MAT_FINAL_ASSEMBLY);
            //outputMatrix(MInv, "MInv.dat");

        cout << "Mat Assembly for "<< "MINVSm"<<endl;
////
        MatAssemblyBegin(MInvSm,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(MInvSm,MAT_FINAL_ASSEMBLY);
////
////        std::cout << t0.elapsed() << " - : ElementIntegrals done\n";
////        
//        PetscViewer viewer1;
////
        unsigned int posR=0;
        unsigned int posL=0;
        unsigned int eR=0;
        unsigned int eL=0;
//
//        printFullMatrixInfo(BF, "BF");
//
//        printFullMatrixInfo(gv.DivergenceFreeMatrix_, "gv.DivergenceFreeMatrix_");
//
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
                
//                cout << "posL="<<posL<<endl;
//                cout << "posR="<<posR<<endl;

                faceIntegral.integrate((*citFe), faceInteg, fData, this);
//                cout <<"***********************************"<<endl;
                
//                cout <<"Face is="<<*citFe<<endl;
                
//                fData.print();
                
                
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
                        
//                            cout <<"fData.right_[j](8,i)="<<fData.right_[j](8,i)<<endl;
                            MatSetValue(globalDIV, posL+j, Nu+Nv+posR+i, -fData.left_[j](11,i),  	ADD_VALUES);//W coefficient
//                            cout <<"-fData.left_[j](11,i)="<<-fData.left_[j](11,i)<<endl;
                            MatSetValue(globalDIV, posR+j, Nu+Nv+posR+i,  fData.right_[j](11,i),  	ADD_VALUES);//W coefficient
                    }
                }
//                cout <<"***********************************"<<endl;
            }
        }
        
        
        MatAssemblyBegin(globalDIV,MAT_FINAL_ASSEMBLY);
            //CHKERRQ(ierr);
        MatAssemblyEnd(globalDIV,MAT_FINAL_ASSEMBLY);
        
        MatAssemblyBegin(MInv,MAT_FINAL_ASSEMBLY);
            //CHKERRQ(ierr);
        MatAssemblyEnd(MInv,MAT_FINAL_ASSEMBLY);
        MatAssemblyBegin(BF,MAT_FINAL_ASSEMBLY);
            //CHKERRQ(ierr);
        MatAssemblyEnd(BF,MAT_FINAL_ASSEMBLY);
        
        

        
        outputMatrix(globalDIV, "globalDIV.dat");
        outputMatrix(BF, "BF.dat");
        
        printFullMatrixInfo(globalDIV, "globalDIV");
        printFullMatrixInfo(BF, "BF");
        
        printFullMatrixInfo(MInv, "MInv");
        
            //********************DIV=[Ex+F4u, Ey+F4v, Ez+F4w]*******
        cout << "Mat Assembly for "<< "DivergenceFreeMatrix_"<<endl;
        MatAssemblyBegin(globalDIV,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(globalDIV,MAT_FINAL_ASSEMBLY);
            //
        
            //*******************************************************
        std::cout << " - : Create rotational matrix!\n";
            //double fillC=1;//(double)(3*N*2*nb)/(3*N*2*nb +3*N*nb);
            //MatMatMult(MInv, C, MAT_INITIAL_MATRIX, fillC, &C1);
            //std::cout << t0.elapsed() << " - : Create rotational matrix1!\n";
            //MatDestroy(&C);
            //printFullMatrixInfo(C1, "C1");
        double fillBF=1;//(double)(3*N*7*nb)/(3*N*7*nb +3*N*nb);
        MatMatMult(MInv, BF, MAT_INITIAL_MATRIX, fillBF, &BF1);
        
        std::cout << " - : Create rotational matrix2!\n";
        MatDestroy(&BF);
        printFullMatrixInfo(BF1, "BF1");
        PetscErrorCode ierr; 
        
            //MatMatMult(BF, M, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &BF);
        
        double fillDIV=1;//(double)(3*N*7*nb)/(3*N*7*nb +3*N*nb);
        MatMatMult(MInvSm, globalDIV, MAT_INITIAL_MATRIX, fillBF, &DIV);
        
        std::cout << " - : Create rotational matrix3!\n";
            //     outputMatrixMatlab(gv.DivergenceFreeMatrix_, "DIV.dat");
        
        MatDestroy(&globalDIV);
        printFullMatrixInfo(DIV, "DIV");
        
        MatDestroy(&MInv);
        MatDestroy(&MInvSm);
        
        
            //outputMatrix(BF, "BF.txt");
        cout << "Mat Assembly for "<< "C1"<<endl;
            //ierr = MatAssemblyBegin(C1,MAT_FINAL_ASSEMBLY);
//        ierr = MatAssemblyEnd(C1,MAT_FINAL_ASSEMBLY);
        
        cout << "Mat Assembly for "<< "BF1"<<endl;
        ierr = MatAssemblyBegin(BF1,MAT_FINAL_ASSEMBLY);
        ierr = MatAssemblyEnd(BF1,MAT_FINAL_ASSEMBLY);
        
        cout << "Mat Assembly for "<< "DIV"<<endl;
        ierr = MatAssemblyBegin(DIV,MAT_FINAL_ASSEMBLY);
        ierr = MatAssemblyEnd(DIV,MAT_FINAL_ASSEMBLY);
        
        
        
        printFullMatrixInfo(DIV, "DIV");
//        printFullMatrixInfo(C1, "C1");
        printFullMatrixInfo(BF1, "BF1");
        
        
            //outputMatrix(Ah, "Ah.txt");
        
//        UCorrected=UInit;
        

        
        double dummy=0;
        const PetscInt* cols;
        const PetscScalar* values;
        PetscInt 	numberOfNonZeros;
        std::cout << " - : Started Creating P matrix!\n";
        PetscInt* nnzP= new PetscInt[Nu+Nv+Nw+Nl];
        PetscInt* nnzQ= new PetscInt[Nu+Nv+Nw+Nl];
            //int noFV=(nb==1? 10:0);
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
        
        /*
         for (int i =0;i<Nw;++i)
         {
         if (i<Nu)
         {
         nnzP[i]=(1+2+noFV)*nb+1;
         nnzQ[i]=(1+2+noFV)*nb+1;
         }
         if (i<Nv)
         {
         nnzP[Nu+i]=(1+2+noFV)*nb+1;;
         nnzQ[Nu+i]=(1+2+noFV)*nb+1;
         }
         
         if (i<Nw)
         {
         nnzP[Nu+Nv+i]=(2+noFV)*nb+1;
         nnzQ[Nu+Nv+i]=(2+noFV)*nb+1;
         }
         if (i<Nl)
         {
         nnzP[Nu+Nv+Nw+i]=(18+noFV)*nb;
         nnzQ[Nu+Nv+Nw+i]=(18+noFV)*nb;
         }
         }*/
        
            // 	std::cout << t0.elapsed() << " stage1"<<endl;
            //
            // 	MatCreateSeqAIJ(PETSC_COMM_SELF, Nu+Nv+Nw+Nl, Nu+Nv+Nw+Nl, PETSC_NULL, nnzP,&P_);
            // 	std::cout << t0.elapsed() << " stage2"<<endl;
            // 	cout <<endl<<Nu+Nv+Nw<<endl;
            //
            // 	MatCreateSeqAIJ(PETSC_COMM_SELF, Nu+Nv+Nw+Nl, Nu+Nv+Nw+Nl, 	PETSC_NULL, nnzQ,&Q_);
            // 	delete [] nnzP;
            // 	delete [] nnzQ;
            //
            //
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
        
            //        outputMatrix(globalDIV, "globalDIV.dat");
//        outputMatrix(BF1, "BF1.dat");
        for (unsigned int i=0;i<N;++i)
        {
//            MatGetRow(C1, i, &numberOfNonZeros, &cols, &values);
//            for (int j=0;j<numberOfNonZeros;++j)
//            {
//                dummy = (values[j]);
//                if (dummy!=0)
//                {
//                    MatSetValue(P_, 	i,     cols[j],     -0.5*dt*dummy, 	ADD_VALUES);
//                    MatSetValue(Q_, 	i,     cols[j],     0.5*dt*dummy, 	ADD_VALUES);
//                }
//            }
//            MatRestoreRow(C1, i, &numberOfNonZeros, &cols, &values);
            
            MatGetRow(BF1, i, &numberOfNonZeros, &cols, &values);
            for (int j=0;j<numberOfNonZeros;++j)
            {
                dummy = (values[j]);
                if (dummy!=0.0)
                {
                        //cout << "FUCKL"<<endl;
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

        
        outputMatrix(P_, "P.dat");
        outputMatrix(Q_, "Q.dat");

//
    }
        /// create Mass Matrices, store them as a User Defined Element Data
        /// calculate projection of the every unknown on the FEM spaces.
    void initialConditions()
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
//            cout <<"count="<<++count<<endl;
//            cout << "ldof= "<<ldof<<endl;
            
//            cout << "Element"<<**cit<<endl;;
            ElementT* elem =(*cit);
            
            elemData = new HEulerElementData(ldof);
            
            LinearAlgebra::Matrix& massMatrix = elemData->massMatrix_;
            LinearAlgebra::Matrix& invMassM   = elemData->invMassMatrix_;
            
            elIntegral.integrate(elem, massMatrixIntegrand, massMatrix, this);
//
//            cout <<"MASS=" <<massMatrix<<endl;
//
            massMatrix.inverse(invMassM);
            
            elem->setUserData(elemData);
//            
//            cout <<"INVERSEMASS=" <<invMassM<<endl;
//
            elIntegral.integrate(elem, uEx, rightHand);
//
//            cout <<"right=" <<rightHand<<endl;
            
            numerical = invMassM*rightHand;// projection of U
         

//            cout <<"numericalU=" <<numerical<<endl;

            elem->setTimeLevelData(0,0,numerical);

            elIntegral.integrate(elem, vEx, rightHand);
//            cout <<"right=" <<rightHand<<endl;
            numerical = invMassM*rightHand;// projection of V
            
            elem->setTimeLevelData(0,1,numerical);
//            cout <<"numericalV=" <<numerical<<endl;

            
            elIntegral.integrate(elem, wEx, rightHand);
//            cout <<"right=" <<rightHand<<endl;
            numerical = invMassM*rightHand;// projection of W
            elem->setTimeLevelData(0,2,numerical);
//            cout <<"numericalW=" <<numerical<<endl;
            
            elIntegral.integrate(elem, pOrRhoEx, rightHand);
//            cout <<"right=" <<rightHand<<endl;
            numerical = invMassM*rightHand;// projection of P
            elem->setTimeLevelData(0,3,numerical);
//            cout <<"numericalP=" <<numerical<<endl;

//            cout << "MESH";
//            cout.flush();
//            MeshManipulatorT* mesh=(meshes_[0]);
//            (mesh)->outputMesh(cout)
        
        }
            
               cout <<"finish calculations of inital condition!"<<endl;
        
    }
    
    void solve()
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
            // 	VecCreateSeq(PETSC_COMM_SELF, N+Nl+1, &RHS);
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
        
//        cout<<"Checking the divergence free condition..."<<endl;
//        if(gdata->solutionType_==HEulerConfigurationData::INCOMPRESSIBLE_WALLS || data->solutionType_==HEulerConfigurationData::INCOMPRESSIBLE_PERIODIC || data->solutionType_==HEulerConfigurationData::INCOMPRESSIBLE_ONETHIRDPERIODIC)
//        {
//            MatMult(gv.DivergenceFreeMatrix_, UExact, RH);
//            
//            PetscInt a1(10);
//            PetscReal maxE;
//            
//            VecMax(RH, &a1, &maxE);
//            userWriteFunction.time_ = currTime;
//            userWriteFunction.writeDiv(maxE);
//        }
        double currTime=0;
        
        int iterations=1;
        int iterations1=0;
        PetscScalar* XTEMP = new PetscScalar [Nu+Nv+Nw+Nl];
//        FixedVector<NumType, 1> energyOverElement;
//        FixedVector<NumType, 4> errorVector;
        
        
        
        
        KSP ksp;
            // Preconditioner
        PC pc;
        cout << "Solve"<<endl;
            // Create a solver
        KSPCreate(PETSC_COMM_SELF, &ksp);
        
        cout << "ksp create"<<endl;
        printFullMatrixInfo(P_, "P_");
        printFullMatrixInfo(Q_, "Q_");
        
        /*! Direct solver: only a preconditioner is needed. To change
         this to an iterative method, you can change the solver to
         KSPMINRES, KSPGMRES, etc. Similarly, you can choose PCNONE
         (instead of PCLU) for the preconditioner below if you want to
         use a nonpreconditioned iterative method. !*/
        
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
		
            //PCSetType(pc, PCICC);
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

        
		
        
        /*! Solving linear system.*/
        
        while(currTime<=endTime)
        {
                //Lambda                         integrandLambda(&Lambda);
            
            cout<<"construct proper rightHandside"<<endl;
			
                ///***************************construct proper rightHandside**************************
                ///***********************************************************************************
            MatMult(Q_, UExact, RHS);
            
                ///***********************************************************************************
                ///***************************construct proper rightHandside**************************
			
            
            /*		VecAssemblyBegin(gv.Q);
             VecAssemblyEnd(gv.Q);*/
            
            
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
//                EnergyWrapper  integrand(gv, el, allData);
//                integrateOverElement(el, integrand, energyOverElement);
//                element.setEnergy(energyOverElement[0]);
//                
//                PressureIntegral  presInt(allData, el);
//                integrateOverElement(el, presInt, energyOverElement);
//                element.setPressureIntegral(energyOverElement[0]);
//                
//                L2ErrorNorm  error(gv, el, allData, currTime,exactSolution);
//                integrateOverElement(el, error, errorVector);
//                element.setUError(errorVector[0]);
//                element.setVError(errorVector[1]);
//                element.setWError(errorVector[2]);
//                element.setPError(errorVector[3]);
            }
            
            VecRestoreArray(Lambda, &XTEMP);
            VecAssemblyBegin(UExact);
            VecAssemblyEnd(UExact);
            
            
            
//            userWriteFunction.time_ = currTime;
//            userWriteFunction.writeEnergy(mesh, allData); 
//			
//            if(gv.solution_==GlobalVariables::INCOMPRESSIBLE_WALLS ||gv.solution_==HEulerConfigurationData::INCOMPRESSIBLE_PERIODIC ||gv.solution_==HEulerConfigurationData::INCOMPRESSIBLE_ONETHIRDPERIODIC)
//            {
//                MatMult(gv.DivergenceFreeMatrix_, UExact, RH);
//                PetscInt a(10);	
//                PetscReal max;
//                VecMax(RH, &a, &max);
//                userWriteFunction.writeDiv(max);
//            }
			
            
            
            
            
           
			output(currTime);
            
////            if (iterations ==1 || iterations%nplotPerPeriod==0) 
////            {
////                sprintf(str2,"x.dat%d",++iterations1);
////                sprintf(str,"u.dat%d",iterations1);
////                sprintf(str1,"uEx.dat%d",iterations1);
////                
////                    // 		   ofstream    file1Dx(str2);
////                    // 		   ofstream    file1D(str);
////                    // 		   ofstream    file1DEx(str1);
////                    // 		   
////                userWriteFunction.time_ = currTime;
////                    // 		   userWriteFunction.setOutputFiles(&file1Dx, &file1D, &file1DEx);
////                tempFileName= formFileName(outFile, currTime, endTime);
////                
////                ofstream tempFile(tempFileName.c_str());
////                
////                TecplotDiscontinuousSolutionWriter3D tecWriter(tempFile, "wbdata", "012", outputFormat);
////                tecWriter.write(mesh, std::string(ostr.str()), false, userWriteFunction);
////                tempFile.close();
////                    //            file1Dx.close();
////                    // 		   file1D.close();
////                    // 		   file1DEx.close();
////                cout << "t = " << currTime<<endl;
////                userWriteFunction.outputErrors();
//                    //userWriteFunction.writeEnergy(mesh, allData);
////        }
////output(currTime);
            cout<<"currTime="<<currTime<<endl;
            iterations++;
        }
            //*************************************************************************************        
        
        KSPDestroy(&ksp);
        VecDestroy(&UExact);
        VecDestroy(&RHS);
        VecDestroy(&Lambda); 
        VecDestroy(&RH); 
//        if(gv.solution_==HEulerConfigurationData::INCOMPRESSIBLE_WALLS ||gv.solution_==HEulerConfigurationData::INCOMPRESSIBLE_PERIODIC ||gv.solution_==HEulerConfigurationData::INCOMPRESSIBLE_ONETHIRDPERIODIC)
//            MatDestroy(&gv.DivergenceFreeMatrix_);
            //PCDestroy(pc);
            //delete[] XTEMP;
        
        PetscFinalize();
        cout << "The end of the time loop"<<endl;

    }
    
    void output(double time=0.0)
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
            //meshes_[0]->outputMesh(cout);
        
        file1Dx0.close();
        file1D0.close();
        file1DEx0.close();
        divFile.close();
        energyFile.close();
        energyExfile.close();
        l2ErrorFile.close();

        
        file3D.close();
            //userWriteFunction.outputErrors();
        
            //userWriteFunction.writeEnergy();
        
        cout <<"FINISH...."<<endl;
        
    }
    
private:
    ExactSolutionT*          exactSolution_;
    Mat                      P_;
    Mat                      Q_;
};