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

typedef Base::MeshManipulator<3> BaseMeshManipulatorT;
typedef Base::threeDBasisFunction basisFunctionT;
//typedef std::list<Base::Face<3> >::iterator FaceIteratorT;
	typedef Base::MeshManipulator<3u>::FaceIterator               FaceIterator;


/**
 * This class should provide problem specific information about the maxwell equations. 
 */
class DomokosProblem : public hpGemUIExtentions<3>
{
private:

        typedef Geometry::PointPhysical<3>                            PointPhysicalT;
	typedef Base::Element<3>                                      ElementT;
        typedef Base::Face<3>                                         FaceT;
	
	typedef Base::MeshManipulator<3u>::FaceIterator               FaceIterator;

public:

    DomokosProblem(int argc,char** argv, MaxwellData* globalConfig, Base::ConfigurationData* elementConfig,matrixFiller* fill):hpGemUIExtentions< 3 >(argc,argv,globalConfig,elementConfig,fill){}

    /**
     * set up the mesh and complete initialisation of the global data and the configuration data
     */
    bool initialise(){
        int n=getData()->NumberOfIntervals_;
        Geometry::PointPhysical<3> bottomLeft, topRight;
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
	
	//FIXME workaround voor het maken van eigen basisfuncties
	BaseMeshManipulatorT* mesh = new MyMeshManipulator(getConfigData(),getData()->PolynomialOrder_,true,true,true);
	
	//The number of unknowns per element in the maxwell equations is equal to the number of basis functions.
	const_cast<MaxwellData*>(getData())->numberOfUnknowns_=mesh->collBasisFSet_[0]->size();
	setConfigData();
	mesh->readCentaurMesh("Cylinder2.hyb");
	//mesh->readCentaurMesh("input_basic2.hyb");
	//mesh->createTriangularMesh(bottomLeft,topRight,numElementsOneD);
	addMesh(mesh);
        //addMesh("Triangular",bottomLeft, topRight, numElementsOneD);
	
	//actually store the full number of unknowns
	const_cast<MaxwellData*>(getData())->numberOfUnknowns_*=mesh->getElementsList().size();
	for(Base::MeshManipulator<3>::ElementIterator it=mesh->elementColBegin();it!=mesh->elementColEnd();++it){
	    (*it)->setUserData(new ElementInfos(**it));
	}
	
	return true;
    }
    
    /**
     * Computes element contributions to the stiffness matrix i.e. (nabla x phi_i) * (nabla x phi_j)
     * returns the contibutions at this gauss point to the entire element matrix in one go
     */
    void elementStiffnessIntegrand(const ElementT* element, const PointElementReferenceT& p, LinearAlgebra::Matrix& ret){
        ElementInfos* info = static_cast<ElementInfos*> (const_cast<ElementT*>(element)->getUserData());
	NumericalVector phi_i(3),phi_j(3);
	std::vector<NumericalVector> functionCurls;
	info->makeFunctionCurlsVector(element,p,functionCurls);
	for(int i=0;i<element->getNrOfBasisFunctions();++i){
	    phi_i=functionCurls[i];
	    for(int j=i;j<element->getNrOfBasisFunctions();++j){
		phi_j=functionCurls[j];
		ret(i,j)=phi_i[0]*phi_j[0]+phi_i[1]*phi_j[1]+phi_i[2]*phi_j[2];
	        ret(j,i)=ret(i,j);
	    }
	}
    }

    /**
     * Computes element contributions to the mass matrix i.e. phi_i * phi_j
     * returns the contibutions at this gauss point to the entire element matrix in one go
     */
    void elementMassIntegrand(const ElementT* element, const PointElementReferenceT& p, LinearAlgebra::Matrix& ret){
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
    void faceIntegrand(const FaceT* face, const PointPhysicalT& normal, const PointFaceReferenceT& p, LinearAlgebra::Matrix& ret){	
	ElementT* right;
	ElementT* left=const_cast<ElementT*>(face->getPtrElementLeft());
	ElementInfos* leftInfo = static_cast<ElementInfos*> (left->getUserData());
	ElementInfos* rightInfo;
	PointElementReferenceT pLeft,pRight;
	face->mapRefFaceToRefElemL(p,pLeft);
	std::vector<NumericalVector> leftValues,leftCurls,rightValues,rightCurls;
	leftInfo->makeFunctionValuesVector(left,pLeft,leftValues);
	leftInfo->makeFunctionCurlsVector(left,pLeft,leftCurls);
	NumericalVector normedNormal(3);
	normedNormal[0] = (normal*(1/Base::L2Norm<3>(normal)))[0];
	normedNormal[1] = (normal*(1/Base::L2Norm<3>(normal)))[1];
	normedNormal[2] = (normal*(1/Base::L2Norm<3>(normal)))[2];
	int leftSize=left->getNrOfBasisFunctions();
	int dimension=left->getNrOfBasisFunctions();
	if(face->isInternal()){
	    right=const_cast<ElementT*>(face->getPtrElementRight());
	    face->mapRefFaceToRefElemR(p,pRight);
	    rightInfo = static_cast<ElementInfos*> (right->getUserData());
	    dimension+=right->getNrOfBasisFunctions();
	    rightInfo->makeFunctionValuesVector(right,pRight,rightValues);
	    rightInfo->makeFunctionCurlsVector(right,pRight,rightCurls);
	}
//	ret.resize(dimension,dimension);
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
	
// 	if(face->isInternal()){
// 	    PointPhysicalT pPhysLeft,pPhysRight;
// 	    left->referenceToPhysical(pLeft,pPhysLeft);
// 	    right->referenceToPhysical(pRight,pPhysRight);
// 	    if(Base::L2Norm<3>(pPhysLeft-pPhysRight)>1e-9){
// 	       cout<<"WARNING: the left element thinks the current integration point is located at :"<<endl<<pPhysLeft<<endl<<"but the right element thinks it is at "<<endl<<pPhysRight<<endl;
// 	    }
// 	}
	
    }
    
    /**
     * Computes the bits of the face contributions that are only used in the IP method
     * i.e. eta_F( (n x phi_i) * (n x phi_j) )
     * returns the contibutions at this gauss point to the entire face matrix in one go
     */
    void faceIntegrandIPPart(const FaceT *face, const PointPhysicalT &normal, const PointFaceReferenceT &p, LinearAlgebra::Matrix &ret){
        ElementT* right;
	ElementT* left=const_cast<ElementT*>(face->getPtrElementLeft());
	ElementInfos* leftInfo = static_cast<ElementInfos*> (left->getUserData());
	ElementInfos* rightInfo;
	PointElementReferenceT pLeft,pRight;
	face->mapRefFaceToRefElemL(p,pLeft);
	std::vector<NumericalVector> leftValues,rightValues;
	leftInfo->makeFunctionValuesVector(left,pLeft,leftValues);
	NumericalVector normedNormal(3);
	normedNormal[0] = (normal*(1/Base::L2Norm<3>(normal)))[0];
	normedNormal[1] = (normal*(1/Base::L2Norm<3>(normal)))[1];
	normedNormal[2] = (normal*(1/Base::L2Norm<3>(normal)))[2];
	int leftSize=left->getNrOfBasisFunctions();
	int dimension=left->getNrOfBasisFunctions();
	NumericalVector phi_i(3),phi_j(3),dummy(3);
	if(face->isInternal()){
	    right=const_cast<ElementT*>(face->getPtrElementRight());
	    face->mapRefFaceToRefElemR(p,pRight);
	    rightInfo = static_cast<ElementInfos*> (right->getUserData());
	    dimension+=right->getNrOfBasisFunctions();
	    rightInfo->makeFunctionValuesVector(right,pRight,rightValues);
	}
//	ret.resize(dimension,dimension);
	for(int i=0;i<dimension;++i){
	    if(i<leftSize){
	        dummy=leftValues[i];
	    }else{
	        dummy=rightValues[i-leftSize];
		dummy*=-1;
	    }
	    OuterProduct(normedNormal,dummy,phi_i);
	    for(int j=i;j<dimension;++j){
		if(j<leftSize){
		    dummy=leftValues[j];
		}else{
		    dummy=rightValues[j-leftSize];
		    dummy*=-1;
		}
		OuterProduct(normedNormal,dummy,phi_j);
		ret(i,j)=getData()->StabCoeff_*(phi_i[0]*phi_j[0]+phi_i[1]*phi_j[1]+phi_i[2]*phi_j[2]);
	        ret(j,i)=ret(i,j);
	    }
	}
    }    
    
    /**
     * Computes the bits of the face contributions that are only used in the BR formulation
     * more accurately only returns phi_i * (n x phi_j) the matrix product should be done elsewhere
     * returns the contibutions at this gauss point to the entire face matrix in one go
     */
    void faceIntegrandBRPart(const FaceT *face, const PointPhysicalT &normal, const PointFaceReferenceT &p, LinearAlgebra::Matrix &ret){
        ElementT* right;
        ElementT* left=const_cast<ElementT*>(face->getPtrElementLeft());
	ElementInfos* leftInfo = static_cast<ElementInfos*> (left->getUserData());
	ElementInfos* rightInfo;
	PointElementReferenceT pLeft,pRight;
	double localepsilon;
	face->mapRefFaceToRefElemL(p,pLeft);
	std::vector<NumericalVector> leftValues,rightValues;
	leftInfo->makeFunctionValuesVector(left,pLeft,leftValues);
	NumericalVector normedNormal(3);
	normedNormal[0] = (normal*(1/Base::L2Norm<3>(normal)))[0];
	normedNormal[1] = (normal*(1/Base::L2Norm<3>(normal)))[1];
	normedNormal[2] = (normal*(1/Base::L2Norm<3>(normal)))[2];
	int leftSize=left->getNrOfBasisFunctions();
	int dimension=left->getNrOfBasisFunctions();
	NumericalVector phi_i(3),phi_j(3),dummy(3);
	if(face->isInternal()){
	    right=const_cast<ElementT*>(face->getPtrElementRight());
	    face->mapRefFaceToRefElemR(p,pRight);
	    rightInfo = static_cast<ElementInfos*> (right->getUserData());
	    dimension+=right->getNrOfBasisFunctions();
	    rightInfo->makeFunctionValuesVector(right,pRight,rightValues);
	}
//	ret.resize(dimension,dimension);
	for(int i=0;i<dimension;++i){
	    if(i<leftSize){
	        phi_i=leftValues[i];
		localepsilon=leftInfo->epsilon_;
	    }else{
	        phi_i=rightValues[i-leftSize];
		localepsilon=rightInfo->epsilon_;
	    }
	    for(int j=0;j<dimension;++j){
		if(j<leftSize){
		    dummy=leftValues[j];
		}else{
		    dummy=rightValues[j-leftSize];
		    dummy*=-1;
		}
		OuterProduct(normedNormal,dummy,phi_j);
		ret(j,i)=(face->isInternal()?1:2)*(phi_i[0]*phi_j[0]+phi_i[1]*phi_j[1]+phi_i[2]*phi_j[2])*sqrt(localepsilon);
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
     */
    void exactSolution(const PointPhysicalT& p, const double t, NumericalVector &ret){
        ret[0]=sin(M_PI*2*p[1])*sin(M_PI*2*p[2]);
        ret[1]=sin(M_PI*2*p[2])*sin(M_PI*2*p[0]);
        ret[2]=sin(M_PI*2*p[0])*sin(M_PI*2*p[1]);
	ret*=cos(sqrt(2)*2*M_PI*t);
//         ret[0]=sin(M_PI*p[1])*sin(M_PI*p[2]);
//         ret[1]=sin(M_PI*p[2])*sin(M_PI*p[0]);
//         ret[2]=sin(M_PI*p[0])*sin(M_PI*p[1]);
// 	ret*=cos(sqrt(2)*M_PI*t);
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
//         ret[0]=sin(M_PI*p[0])*(cos(M_PI*p[1])-cos(M_PI*p[2]));
//         ret[1]=sin(M_PI*p[1])*(cos(M_PI*p[2])-cos(M_PI*p[0]));
//         ret[2]=sin(M_PI*p[2])*(cos(M_PI*p[0])-cos(M_PI*p[1]));
// 	ret*=cos(sqrt(2)*M_PI*t)*M_PI;
//          ret[0]=0;ret[1]=0;ret[2]=0;
    }
};

/**
 * main routine. Suggested usage: make a problem, solve it, tell the user your results.
 */
int main(int argc,char** argv){
  
  
    time_t start,end,initialised,solved;
    time(&start);
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
    DomokosProblem problem(argc-2,&argv[2],new MaxwellData(elements,order),new Base::ConfigurationData(0,0,1),new matrixFillerBR);
    try{
        problem.initialise();
	time(&initialised);
        problem.solveHarmonic();
	time(&solved);
	char filename[]="output.dat";
	problem.makeOutput(filename);
	time(&end);
	cout<<"Initialisation took "<<difftime(initialised,start)<<" seconds."<<endl;
	cout<<"Solving the problem took "<<difftime(solved,initialised)<<" seconds."<<endl;
	cout<<"The rest took "<<difftime(end,solved)<<" seconds."<<endl;
	
    }catch(const char* message){
        std::cout << message;
    }
    return 0;
}
