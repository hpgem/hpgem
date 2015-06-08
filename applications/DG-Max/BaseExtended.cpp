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

#include "BaseExtended.h"
#include "kspaceData.h"
#include "Utilities/GlobalMatrix.h"
#include "Utilities/GlobalVector.h"
#include "petscksp.h"
#include "Base/ConfigurationData.h"
#include "Geometry/PointPhysical.h"
#include "Integration/ElementIntegral.h"
#include "Integration/FaceIntegral.h"
#include "Base/ShortTermStorageElementHcurl.h"
#include "Base/ShortTermStorageFaceHcurl.h"
#include "Output/TecplotDiscontinuousSolutionWriter.h"
#include <valarray>
#include <cmath>
#include "Base/ElementCacheData.h"
#include "Base/FaceCacheData.h"

void hpGemUIExtentions::setConfigData()
{
    // const_cast<Base::ConfigurationData*>(configData_)->numberOfBasisFunctions_=globalData_->numberOfUnknowns_;
    const_cast<Base::GlobalData*>(globalData_)->numberOfTimeLevels_ = configData_->numberOfTimeLevels_;
    const_cast<Base::ConfigurationData*>(configData_)->numberOfUnknowns_ = 1;
}

const Base::ConfigurationData* hpGemUIExtentions::getConfigData()
{
    return configData_;
}

const MaxwellData* hpGemUIExtentions::getData() const
{
    return dynamic_cast<const MaxwellData*>(globalData_);
}

/*
 void hpGemUIExtentions::elementIntegrand(const Base::HpgemUI::ElementT* element, const PointElementReferenceT& p, LinearAlgebra::NumericalVector& ret) {
 ret.resize(element->getNrOfBasisFunctions());
 //ElementInfos* info = static_cast<ElementInfos*>(const_cast<ElementT*>(element)->getUserData());
 PointPhysicalT pPhys(3);
 element->referenceToPhysical(p,pPhys);
 LinearAlgebra::NumericalVector val(3),phi(3);
 //std::vector<LinearAlgebra::NumericalVector> functionValues;
 //info->makeFunctionValuesVector(element,p,functionValues);
 initialConditions(pPhys,val);
 for(int i=0; i<element->getNrOfBasisFunctions(); ++i) {
 //phi=functionValues[i];
 element->basisFunction(i, p, phi);
 ret(i)=phi[0]*val[0]+phi[1]*val[1]+phi[2]*val[2];
 }
 }
 */

/*
 void hpGemUIExtentions::initialConditionsDeriv(const Base::HpgemUI::ElementT* element, const PointElementReferenceT& p, LinearAlgebra::Matrix& ret){
 ret.resize(element->getNrOfBasisFunctions(),1);
 //ElementInfos* info = static_cast<ElementInfos*>(const_cast<ElementT*>(element)->getUserData());
 PointPhysicalT pPhys(3);
 element->referenceToPhysical(p,pPhys);
 //std::vector<LinearAlgebra::NumericalVector> functionValues;
 //info->makeFunctionValuesVector(element,p,functionValues);
 LinearAlgebra::NumericalVector val(3),phi(3);
 initialConditionsDeriv(pPhys,val);
 for(int i=0; i<element->getNrOfBasisFunctions(); ++i) {
 //phi=functionValues[i];
 element->basisFunction(i, p, phi);
 ret(i,0)=phi[0]*val[0]+phi[1]*val[1]+phi[2]*val[2];
 }
 }
 
 */
double hpGemUIExtentions::sourceTermTime(const double t)
{
    return 1.;
}

void hpGemUIExtentions::exactSolution(const PointPhysicalT& p, const double t, LinearAlgebra::MiddleSizeVector &ret)
{
    //ret[0]=sin(M_PI*2*p[1])*sin(M_PI*2*p[2]);
    //ret[1]=sin(M_PI*2*p[2])*sin(M_PI*2*p[0]);
    //ret[2]=sin(M_PI*2*p[0])*sin(M_PI*2*p[1]);
    //ret*=cos(sqrt(2)*2*M_PI*t);
    
    ret[0] = sin(M_PI * p[1]) * sin(M_PI * p[2]);
    ret[1] = sin(M_PI * p[2]) * sin(M_PI * p[0]);
    ret[2] = sin(M_PI * p[0]) * sin(M_PI * p[1]);
    ret *= cos(sqrt(2) * M_PI * t);
    
    //     ret[0]=p[0]*(1-p[0]);
    //     ret[1]=0;
    // 	   ret[2]=0;
}

void hpGemUIExtentions::exactSolutionCurl(const PointPhysicalT& p, const double t, LinearAlgebra::MiddleSizeVector &ret)
{
    //ret[0]=sin(M_PI*2*p[0])*(cos(M_PI*2*p[1])-cos(M_PI*2*p[2]));
    //ret[1]=sin(M_PI*2*p[1])*(cos(M_PI*2*p[2])-cos(M_PI*2*p[0]));
    //ret[2]=sin(M_PI*2*p[2])*(cos(M_PI*2*p[0])-cos(M_PI*2*p[1]));
    //ret*=cos(sqrt(2)*2*M_PI*t)*2*M_PI;
    
    ret[0] = sin(M_PI * p[0]) * (cos(M_PI * p[1]) - cos(M_PI * p[2]));
    ret[1] = sin(M_PI * p[1]) * (cos(M_PI * p[2]) - cos(M_PI * p[0]));
    ret[2] = sin(M_PI * p[2]) * (cos(M_PI * p[0]) - cos(M_PI * p[1]));
    ret *= cos(sqrt(2) * M_PI * t) * M_PI;
    
    //          ret[0]=0;ret[1]=0;ret[2]=0;
}

void hpGemUIExtentions::writeToTecplotFile(const Base::HpgemUI::ElementT* element, const PointElementReferenceT& p, std::ostream& output)
{
    //ElementInfos* info = static_cast<ElementInfos*>(const_cast<ElementT*>(element)->getUserData());
    LinearAlgebra::MiddleSizeVector data = const_cast<ElementT*>(element)->getTimeLevelData(timelevel_);
    LinearAlgebra::MiddleSizeVector results1(3), results(3), curls1(3), curls(3), waveDirection(3), Eorth(3), Horth(3);
    //std::vector<LinearAlgebra::NumericalVector> values,curlsVec;
    //info->makeFunctionValuesVector(element,p,values);
    //info->makeFunctionCurlsVector(element,p,curlsVec);
    for (int i = 0; i < element->getNrOfBasisFunctions(); ++i)
    {
        
        //results+=values[i]*data[i];
        //curls+=curlsVec[i]*data[i];
        
        element->basisFunction(i, p, results1);
        element->basisFunctionCurl(i, p, curls1);
        results += results1 * data[i];
        curls += curls1 * data[i];
    }
    output << results[0] << " " << results[1] << " " << results[2] << " " << curls[0] << " " << curls[1] << " " << curls[2] << std::endl;
}
/*
 template<>
 void hpGemUIExtentions<3>::writeTecplotFile(const Base::HpgemUI< 3 >::MeshManipulatorT& mesh, const char* zonetitle, const int timelevel, std::ofstream& file, const bool existingFile){
 int numberOfNodes;
 int numberOfElements;//tecplot somewhy messes up the ordering of the elements for linked zones so dont reuse this
 long int numberOfNodesPosition;
 long int numberOfElementsPosition;
 if(!existingFile) {
 file<<"TITLE = \"The electric field\"\nVARIABLES = \"x0\", \"x1\", \"x2\", \"E0\", \"E1\", \"E2\", \"H0\", \"H1\", \"H2\"\n";
 }
 file<<"ZONE T = \""<<zonetitle<<"\", ZONETYPE = FETETRAHEDRON, DATAPACKING = POINT, N = ";
 numberOfNodesPosition = file.tellp();
 file<<"         , E = ";
 numberOfElementsPosition = file.tellp();
 file<<"        \n";
 numberOfElements=0;
 numberOfNodes=0;
 PointElementReferenceT p;
 PointPhysicalT pPhys;
 for(ConstElementIterator it=mesh.elementColBegin(); it!=mesh.elementColEnd(); ++it) {
 ++numberOfElements;
 numberOfNodes+=4;
 for(int i=0; i<4; ++i) {
 (*it)->getReferenceGeometry()->getNode(i,p);
 (*it)->referenceToPhysical(p,pPhys);
 file<<pPhys[0]<<" "<<pPhys[1]<<" "<<pPhys[2]<<" ";
 timelevel_=timelevel;
 writeFieldValues(*(*it),p,file);
 }
 file<<"\n";
 }
 long int currentPosition=file.tellp();
 file.seekp(numberOfNodesPosition);
 file<<numberOfNodes;
 file.seekp(numberOfElementsPosition);
 file<<numberOfElements;
 file.seekp(currentPosition);
 for(int i=1; i<=numberOfNodes; ++i) {
 file<<i<<" ";
 if(i%1024==0){//some fileviewers cant handle very long lines
 file<<std::endl;
 }
 }
 file<<std::endl;
 }
 */

void hpGemUIExtentions::elementIntegrand(const Base::HpgemUI::ElementT* element, const PointElementReferenceT& p, errorData& ret)
{
    //ret.resize(2,1);
    ElementInfos* info = static_cast<ElementInfos*>(const_cast<ElementT*>(element)->getUserData());
    PointPhysicalT pPhys(3);
    element->referenceToPhysical(p, pPhys);
    //std::vector<LinearAlgebra::NumericalVector> functionValues,functionCurls;
    //info->makeFunctionValuesVector(element,p,functionValues);
    //info->makeFunctionCurlsVector(element,p,functionCurls);
    
    LinearAlgebra::MiddleSizeVector phi(3), phiCurl(3), error(3), errorCurl(3);
    
    exactSolution(pPhys, measureTimes_[timelevel_], error); //maybe move exactSolution from DGMax to here
    exactSolutionCurl(pPhys, measureTimes_[timelevel_], errorCurl); //mayebe move exactSolutionCurl from DGMax to here
    LinearAlgebra::MiddleSizeVector data = element->getTimeLevelData(timelevel_);
    
    for (int i = 0; i < element->getNrOfBasisFunctions(); ++i)
    {
        //phi=functionValues[i];
        //phiCurl=functionCurls[i];
        element->basisFunction(i, p, phi);
        element->basisFunctionCurl(i, p, phiCurl);
        error -= data[i] * phi;
        errorCurl -= data[i] * phiCurl;
    }
    ret[0] = Base::L2Norm(error) * Base::L2Norm(error);
    ret[1] = Base::L2Norm(errorCurl) * Base::L2Norm(errorCurl);
}

/*
 void hpGemUIExtentions::faceIntegrand(const Base::HpgemUI::FaceT* face, const LinearAlgebra::NumericalVector& normal, const PointFaceReferenceT& p, errorData& ret) {
 ElementT* element=const_cast<ElementT*>(face->getPtrElementLeft());
 ElementInfos* info = static_cast<ElementInfos*>(element->getUserData());
 PointElementReferenceT pElement(3);
 PointPhysicalT PPhys(3);
 face->mapRefFaceToRefElemL(p,pElement);
 element->referenceToPhysical(pElement,PPhys);
 LinearAlgebra::NumericalVector error(3),phi(3),dummy(3),normedNormal(3);
 normedNormal[0] = (normal*(1/Base::L2Norm(normal)))[0];
 normedNormal[1] = (normal*(1/Base::L2Norm(normal)))[1];
 normedNormal[2] = (normal*(1/Base::L2Norm(normal)))[2];
 exactSolution(PPhys,measureTimes_[timelevel_],dummy);
 OuterProduct(normedNormal,dummy,error);
 //std::vector<LinearAlgebra::NumericalVector> functionValues;
 //info->makeFunctionValuesVector(element,pElement,functionValues);
 LinearAlgebra::Matrix data = element->getTimeLevelData(timelevel_);
 for(int i=0; i<element->getNrOfBasisFunctions(); ++i) {
 //dummy=functionValues[i];
 element->basisFunction(i, pElement, dummy);
 OuterProduct(normedNormal,dummy,phi);
 error-=data(i,0)*phi;
 }
 if(face->isInternal()) {
 element=const_cast<ElementT*>(face->getPtrElementRight());
 ElementInfos* infoAlso = static_cast<ElementInfos*>(element->getUserData());
 face->mapRefFaceToRefElemR(p,pElement);
 element->referenceToPhysical(pElement,PPhys);
 exactSolution(PPhys,measureTimes_[timelevel_],dummy);
 OuterProduct(normedNormal,dummy,phi);
 error-=phi;
 //std::vector<LinearAlgebra::NumericalVector> functionValuesToo;
 //infoAlso->makeFunctionValuesVector(element,pElement,functionValuesToo);
 data=element->getTimeLevelData(timelevel_);
 for(int i=0; i<element->getNrOfBasisFunctions();++i) {
 //dummy=functionValuesToo[i];
 element->basisFunction(i, pElement, dummy);
 OuterProduct(normedNormal,dummy,phi);
 error+=data(i,0)*phi;//iets van een procentje te weinig
 }
 }
 ret[0]=Base::L2Norm(error)*Base::L2Norm(error);
 }
 */

void hpGemUIExtentions::faceIntegrand(const Base::HpgemUI::FaceT* face, const LinearAlgebra::MiddleSizeVector& normal, const PointFaceReferenceT& p, errorData& ret)
{
    ElementT* element = const_cast<ElementT*>(face->getPtrElementLeft());
    PointElementReferenceT pElement(3);
    PointPhysicalT PPhys(3);
    
    face->mapRefFaceToRefElemL(p, pElement);
    element->referenceToPhysical(pElement, PPhys);
    
    LinearAlgebra::MiddleSizeVector error(3), phi(3), dummy(3), dummy2(3), dummy3(3), normedNormal(3);
    normedNormal[0] = (normal * (1 / Base::L2Norm(normal)))[0];
    normedNormal[1] = (normal * (1 / Base::L2Norm(normal)))[1];
    normedNormal[2] = (normal * (1 / Base::L2Norm(normal)))[2];
    
    exactSolution(PPhys, measureTimes_[timelevel_], dummy); //maybe move exactSolution from DGMax to here
    OuterProduct(normedNormal, dummy, error);
    
    //ElementInfos* info = static_cast<ElementInfos*>(element->getUserData());
    //std::vector<LinearAlgebra::NumericalVector> functionValues;
    //info->makeFunctionValuesVector(element,pElement,functionValues);
    
    int n = face->getPtrElementLeft()->getNrOfBasisFunctions();
    int M = face->getNrOfBasisFunctions();
    //int M = face->getPtrElementLeft()->getNrOfBasisFunctions() + face->getPtrElementRight()->getNrOfBasisFunctions();
    
    LinearAlgebra::MiddleSizeVector data = element->getTimeLevelData(timelevel_);
    for (int i = 0; i < n; ++i)
    {
        
        //dummy=functionValues[i];
        face->basisFunctionNormal(i, normal, p, phi);
        error -= data[i] * phi;
    }
    if (face->isInternal())
    {
        element = const_cast<ElementT*>(face->getPtrElementRight());
        face->mapRefFaceToRefElemR(p, pElement);
        element->referenceToPhysical(pElement, PPhys);
        exactSolution(PPhys, measureTimes_[timelevel_], dummy2); //maybe move exactSolution from DGMax to here
        OuterProduct(normedNormal, dummy2, dummy3);
        
        error -= dummy3;
        
        //ElementInfos* infoAlso = static_cast<ElementInfos*>(element->getUserData());
        //std::vector<LinearAlgebra::NumericalVector> functionValuesToo;
        //infoAlso->makeFunctionValuesVector(element,pElement,functionValuesToo);
        
        data = element->getTimeLevelData(timelevel_);
        n = face->getPtrElementLeft()->getNrOfBasisFunctions();
        for (int i = n; i < M; ++i)
        {
            //dummy=functionValuesToo[i];
            
            face->basisFunctionNormal(i, normal, p, phi);
            error += data[i - n] * phi; //iets van een procentje te weinig
        }
    }
    ret[0] = Base::L2Norm(error) * Base::L2Norm(error);
}

void hpGemUIExtentions::computeErrors(Base::HpgemUI::MeshManipulatorT& Mesh, int timelevel, double& L2Norm, double& InfNorm, double& HCurlNorm, double& DGNorm)
{
    //this should probalby be toggled off when the exact solution is not known
    //alternatively something based on richardson extrapolation can be done
    int TotalAmountOfProcessors, localProcessorNumber, numberOfElements;
    MPI_Comm_size(PETSC_COMM_WORLD, &TotalAmountOfProcessors);
    MPI_Comm_rank(PETSC_COMM_WORLD, &localProcessorNumber);
    numberOfElements = meshes_[0]->getElementsList().size();
    
    //ElementFunction elF = &hpGemUIExtentions::elementErrorIntegrand;
    this->timelevel_ = timelevel;
    //LinearAlgebra::Matrix integrationResults(2,1);
    errorData integrationResults;
    
    Base::ShortTermStorageElementBase* localElement_;
    localElement_ = new Base::ShortTermStorageElementHcurl(3);
    
    Integration::ElementIntegral elIntegral(false);
    elIntegral.setStorageWrapper(localElement_);
    
    L2Norm = 0;
    HCurlNorm = 0;
    for (ElementIterator it = Mesh.elementColBegin(); it != Mesh.elementColEnd(); ++it)
    {
        if (((*it)->getID()) / (numberOfElements / TotalAmountOfProcessors + 1) == localProcessorNumber)
        {
            elIntegral.integrate<errorData>((*it), this, integrationResults);
            L2Norm += std::abs(integrationResults[0]);
            HCurlNorm += std::abs(integrationResults[0]) + std::abs(integrationResults[1]);
            integrationResults.reset();
        }
    }
    DGNorm = HCurlNorm;
    //FaceFunction faF = & hpGemUIExtentions::faceErrorIntegrand;
    //integrationResults.resize(1,1);
    
    Base::ShortTermStorageFaceBase* localFace_;
    localFace_ = new Base::ShortTermStorageFaceHcurl(3); // creating object for H curl transformation
            
    Integration::FaceIntegral faIntegral(false);
    faIntegral.setStorageWrapper(localFace_);
    
    for (FaceIterator it = Mesh.faceColBegin(); it != Mesh.faceColEnd(); ++it)
    {
        if (((*it)->getPtrElementLeft()->getID()) / (numberOfElements / TotalAmountOfProcessors + 1) == localProcessorNumber)
        {
            faIntegral.integrate<errorData>(*it, this, integrationResults);
            DGNorm += std::abs(integrationResults[0]);
            integrationResults.reset();
        }
    }
}

void hpGemUIExtentions::makeOutput(char* filename)
{
    //set up containers for the errors
    Vec L2Norm, InfNorm, HCurlNorm, DGNorm;
    double L2NormEntry, InfNormEntry, HCurlNormEntry, DGNormEntry;
    ierr_ = VecCreate(PETSC_COMM_WORLD, &L2Norm);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecSetSizes(L2Norm, PETSC_DECIDE, getData()->numberOfTimeLevels_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecSetUp(L2Norm);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecDuplicate(L2Norm, &HCurlNorm);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecDuplicate(L2Norm, &DGNorm);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    
    //do tecplot related stuff
    std::ofstream fileWriter;
    fileWriter.open(filename);
    int dimensions[3] = {0, 1, 2};
    Output::TecplotDiscontinuousSolutionWriter tecplotWriter(fileWriter, "The electric field", "012", "E0,E1,E2,H0,H1,H2");
    std::stringstream string;
    char parsed[20] = "";
    string.readsome(parsed, 20); //clean storage
    string << "t=" << measureTimes_[0];
    int read = string.readsome(parsed, 20);
    string.readsome(&parsed[read], 20 - read); //readsome doen't want to read in one go on the first try :(
    //writeFunction writeDataFunction=&hpGemUIExtentions::writeToTecplotFile;
    timelevel_ = 0;
    tecplotWriter.write(meshes_[0], parsed, false, this);
    //writeTecplotFile(*meshes_[0],parsed,0,fileWriter,false);
    
    //compute errors
    computeErrors(*meshes_[0], 0, L2NormEntry, InfNormEntry, HCurlNormEntry, DGNormEntry);
    VecSetValue(L2Norm, 0, L2NormEntry, ADD_VALUES);
    VecSetValue(HCurlNorm, 0, HCurlNormEntry, ADD_VALUES);
    VecSetValue(DGNorm, 0, DGNormEntry, ADD_VALUES);
    for (int i = 0; i < getData()->numberOfTimeLevels_ - 1; ++i)
    {
        timelevel_ = i;
        computeErrors(*meshes_[0], i + 1, L2NormEntry, InfNormEntry, HCurlNormEntry, DGNormEntry); //also updates the timelevel_
        VecSetValue(L2Norm, i + 1, L2NormEntry, ADD_VALUES);
        VecSetValue(HCurlNorm, i + 1, HCurlNormEntry, ADD_VALUES);
        VecSetValue(DGNorm, i + 1, DGNormEntry, ADD_VALUES);
        string << "t=" << measureTimes_[i + 1];
        string.readsome(parsed, 20);
        tecplotWriter.write(meshes_[0], parsed, true, this);
        //writeTecplotFile(*meshes_[0],parsed,i+1,fileWriter,true);
    }
    ierr_ = VecAssemblyBegin(L2Norm);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecAssemblyBegin(HCurlNorm);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecAssemblyBegin(DGNorm);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecAssemblyEnd(L2Norm);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecAssemblyEnd(HCurlNorm);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecAssemblyEnd(DGNorm);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecSqrtAbs(L2Norm);
    ierr_ = VecSqrtAbs(HCurlNorm);
    ierr_ = VecSqrtAbs(DGNorm);
    VecView(L2Norm, 0);
    //	VecView(L2Norm,PETSC_VIEWER_DRAW_WORLD);
    //VecView(InfNorm,0);//FIXME not computed yet
    //	VecView(InfNorm,PETSC_VIEWER_DRAW_WORLD);
    VecView(HCurlNorm, 0);
    //	VecView(HCurlNorm,PETSC_VIEWER_DRAW_WORLD);
    VecView(DGNorm, 0);
    //	VecView(DGNorm,PETSC_VIEWER_DRAW_WORLD);
    ierr_ = VecDestroy(&HCurlNorm);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecDestroy(&DGNorm);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecDestroy(&L2Norm);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
}

hpGemUIExtentions::hpGemUIExtentions(int argc, char** argv, MaxwellData* globalConfig, Base::ConfigurationData* elementConfig, MatrixAssembly* fluxType)
        : HpgemAPIBase(globalConfig, elementConfig), assembler(fluxType)
{
    ierr_ = SlepcInitialize(&argc, &argv, (char*) 0, "PETSc help");
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = PetscLogBegin();
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecCreate(PETSC_COMM_WORLD, &x_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecCreate(PETSC_COMM_WORLD, &RHS_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecCreate(PETSC_COMM_WORLD, &derivative_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = MatCreate(PETSC_COMM_WORLD, &M_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = MatCreate(PETSC_COMM_WORLD, &S_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = MatSetType(M_, "mpibaij");
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = MatSetType(S_, "mpibaij");
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = KSPCreate(PETSC_COMM_WORLD, &solver_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = EPSCreate(PETSC_COMM_WORLD, &eigenSolver_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    storage_ = new const PetscScalar*;
    measureTimes_ = new double[getConfigData()->numberOfTimeLevels_];
}

hpGemUIExtentions::~hpGemUIExtentions()
{
    delete storage_;
    delete[] measureTimes_;
    ierr_ = VecDestroy(&x_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecDestroy(&RHS_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecDestroy(&derivative_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = MatDestroy(&M_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = MatDestroy(&S_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = KSPDestroy(&solver_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = EPSDestroy(&eigenSolver_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    PetscViewer log;
    ierr_ = PetscViewerASCIIOpen(PETSC_COMM_WORLD, "maxwell.log", &log);
    ierr_ = PetscLogView(log);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = PetscViewerDestroy(&log);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = SlepcFinalize();
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
}

Base::HpgemUI::MeshId hpGemUIExtentions::addMesh(Base::HpgemUI::MeshManipulatorT* mesh)
{
    meshes_.push_back(mesh);
    return meshes_.size() - 1;
}

void hpGemUIExtentions::makeShiftMatrix(LinearAlgebra::MiddleSizeVector& direction, Vec& waveVecMatrix)
{
    for (ElementIterator it = meshes_[0]->elementColBegin(); it != meshes_[0]->elementColEnd(); ++it)
    {
        for (int j = 0; j < (*it)->getNrOfBasisFunctions(); ++j)
        {
            PointPhysicalT centerPhys(3);
            PointReferenceT center(3);
            (*it)->getReferenceGeometry()->getCenter(center);
            ;
            (*it)->referenceToPhysical(center, centerPhys);
            //this extra accuracy is probably irrelevant and a lot of extra ugly to get it working
            //static_cast<Base::threeDBasisFunction*>((*it)->basisFunctionSet_->vecOfBasisFcn_[j])->getReasonableNode(*(*it),centerPhys);
            
            PetscScalar value = exp(std::complex<double>(0, direction[0] * centerPhys[0] + direction[1] * centerPhys[1] + direction[2] * centerPhys[2]));
            VecSetValue(waveVecMatrix, ((*it)->getID()) * (*it)->getNrOfBasisFunctions() + j, value, INSERT_VALUES);
        }
    }
}

void hpGemUIExtentions::findBoundaryBlocks(std::vector<IS>& xRow, std::vector<IS>& xCol, std::vector<IS>& yRow, std::vector<IS>& yCol, std::vector<IS>& zRow, std::vector<IS>& zCol)
{
    int nx(0), ny(0), nz(0);
    int places[] = {0};
    int blocksize = 0;
    for (FaceIterator it = meshes_[0]->faceColBegin(); it != meshes_[0]->faceColEnd(); ++it)
    {
        assert((*it)->isInternal());
        PointFaceReferenceT p(3);
        (*it)->getReferenceGeometry()->getCenter(p);
        PointElementReferenceT pLeft(3), pRight(3);
        PointPhysicalT pLeftPhys(3), pRightPhys(3);
        (*it)->mapRefFaceToRefElemL(p, pLeft);
        (*it)->mapRefFaceToRefElemR(p, pRight);
        (*it)->getPtrElementLeft()->referenceToPhysical(pLeft, pLeftPhys);
        (*it)->getPtrElementRight()->referenceToPhysical(pRight, pRightPhys);
        //if the left coordinate is not close to the right coordinate it is a boundary face
        if (Base::L2Norm(pLeftPhys - pRightPhys) > 1e-3)
        { //pretty lousy tolerance, but this norm should be either 1 or 0
            //cout<<pLeftPhys<<std::endl<<pRightPhys<<std::endl<<pLeftPhys-pRightPhys<<std::endl<<(pLeftPhys[0]-pRightPhys[0])*(pLeftPhys[0]-pRightPhys[0])<<std::endl;;
            if ((pLeftPhys[0] - pRightPhys[0]) * (pLeftPhys[0] - pRightPhys[0]) > 1e-3)
            {
                //cout<<"X!"<<std::endl;
                xRow.resize(nx + 2);
                xCol.resize(nx + 2);
                places[0] = (*it)->getPtrElementLeft()->getID();
                ISCreateBlock(PETSC_COMM_WORLD, (*it)->getPtrElementLeft()->getNrOfBasisFunctions(), 1, places, PETSC_COPY_VALUES, &xRow[nx]);
                ISCreateBlock(PETSC_COMM_WORLD, (*it)->getPtrElementLeft()->getNrOfBasisFunctions(), 1, places, PETSC_COPY_VALUES, &xCol[nx + 1]);
                places[0] = (*it)->getPtrElementRight()->getID();
                ISCreateBlock(PETSC_COMM_WORLD, (*it)->getPtrElementRight()->getNrOfBasisFunctions(), 1, places, PETSC_COPY_VALUES, &xRow[nx + 1]);
                ISCreateBlock(PETSC_COMM_WORLD, (*it)->getPtrElementRight()->getNrOfBasisFunctions(), 1, places, PETSC_COPY_VALUES, &xCol[nx]);
                blocksize = (*it)->getPtrElementLeft()->getNrOfBasisFunctions();
                nx += 2;
            }
            else if ((pLeftPhys[1] - pRightPhys[1]) * (pLeftPhys[1] - pRightPhys[1]) > 1e-3)
            {
                //cout<<"Y!"<<std::endl;
                yRow.resize(ny + 2);
                yCol.resize(ny + 2);
                places[0] = (*it)->getPtrElementLeft()->getID();
                ISCreateBlock(PETSC_COMM_WORLD, (*it)->getPtrElementLeft()->getNrOfBasisFunctions(), 1, places, PETSC_COPY_VALUES, &yRow[ny]);
                ISCreateBlock(PETSC_COMM_WORLD, (*it)->getPtrElementLeft()->getNrOfBasisFunctions(), 1, places, PETSC_COPY_VALUES, &yCol[ny + 1]);
                places[0] = (*it)->getPtrElementRight()->getID();
                ISCreateBlock(PETSC_COMM_WORLD, (*it)->getPtrElementRight()->getNrOfBasisFunctions(), 1, places, PETSC_COPY_VALUES, &yRow[ny + 1]);
                ISCreateBlock(PETSC_COMM_WORLD, (*it)->getPtrElementRight()->getNrOfBasisFunctions(), 1, places, PETSC_COPY_VALUES, &yCol[ny]);
                ny += 2;
            }
            else
            {
                //cout<<"Z!"<<std::endl;
                zRow.resize(nz + 2);
                zCol.resize(nz + 2);
                assert((pLeftPhys[2] - pRightPhys[2]) * (pLeftPhys[2] - pRightPhys[2]) > (1.e-3));
                places[0] = (*it)->getPtrElementLeft()->getID();
                ISCreateBlock(PETSC_COMM_WORLD, (*it)->getPtrElementLeft()->getNrOfBasisFunctions(), 1, places, PETSC_COPY_VALUES, &zRow[nz]);
                ISCreateBlock(PETSC_COMM_WORLD, (*it)->getPtrElementLeft()->getNrOfBasisFunctions(), 1, places, PETSC_COPY_VALUES, &zCol[nz + 1]);
                places[0] = (*it)->getPtrElementRight()->getID();
                ISCreateBlock(PETSC_COMM_WORLD, (*it)->getPtrElementRight()->getNrOfBasisFunctions(), 1, places, PETSC_COPY_VALUES, &zRow[nz + 1]);
                ISCreateBlock(PETSC_COMM_WORLD, (*it)->getPtrElementRight()->getNrOfBasisFunctions(), 1, places, PETSC_COPY_VALUES, &zCol[nz]);
                nz += 2;
            }
        }
    }
}

void hpGemUIExtentions::LDOSIntegrand(const Base::HpgemUI::ElementT* element, const PointElementReferenceT& p, double& ret)
{ //currently LDOS is computed by evaluation at a point so integrand is a bit of a misnomer
    ElementInfos* info = static_cast<ElementInfos*>(const_cast<ElementT*>(element)->getUserData());
    LinearAlgebra::MiddleSizeVector Phi(3), PhiRealI(3), PhiRealJ(3), PhiImagI(3), PhiImagJ(3);
    //std::vector<LinearAlgebra::NumericalVector> functionValues(element->getNrOfBasisFunctions());
    //info->makeFunctionValuesVector(element,p,functionValues);
    for (int i = 0; i < element->getNrOfBasisFunctions(); ++i)
    {
        
        //PhiRealI+=functionValues[i]*element->getData(0,0,i);
        //PhiImagI+=functionValues[i]*element->getData(1,0,i);
        //PhiRealJ+=functionValues[i]*element->getData(2,0,i);
        //PhiImagJ+=functionValues[i]*element->getData(3,0,i);
        
        element->basisFunction(i, p, Phi);
        
        PhiRealI += Phi * element->getData(0, 0, i);
        PhiImagI += Phi * element->getData(1, 0, i);
        PhiRealJ += Phi * element->getData(2, 0, i);
        PhiImagJ += Phi * element->getData(3, 0, i);
    }
    //current implementation: orientation averaged LDOS; modify the next expression for other orientations
    ret = PhiRealI[0] * PhiRealJ[0] + PhiRealI[1] * PhiRealJ[1] + PhiRealI[2] * PhiRealJ[2] + PhiImagI[0] * PhiImagJ[0] + PhiImagI[1] * PhiImagJ[1] + PhiImagI[2] * PhiImagJ[2];
    ret *= 48; //assume computation is done only in the irreducible brillouin zone
    //cout<<ret<<std::endl;
}

void hpGemUIExtentions::makeFunctionValue(Vec eigenVector, LinearAlgebra::MiddleSizeVector& result)
{
    Vec scaledVec;
    double partialResult;
    int index(0);
    ierr_ = VecDuplicate(eigenVector, &scaledVec);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = MatMult(M_, eigenVector, scaledVec);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    result.resize(4 * meshes_[0]->getNumberOfElements());
    for (ElementIterator it = meshes_[0]->elementColBegin(); it != meshes_[0]->elementColEnd(); ++it)
    {
        ierr_ = VecGetArrayRead(eigenVector, storage_);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        for (int i = 0; i < (*it)->getNrOfBasisFunctions(); ++i)
        {
            (*it)->setData(0, 0, i, (*storage_)[(*it)->getNrOfBasisFunctions() * ((*it)->getID()) + i].real());
            (*it)->setData(1, 0, i, (*storage_)[(*it)->getNrOfBasisFunctions() * ((*it)->getID()) + i].imag());
        }
        ierr_ = VecRestoreArrayRead(eigenVector, storage_);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        ierr_ = VecGetArrayRead(scaledVec, storage_);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        for (int i = 0; i < (*it)->getNrOfBasisFunctions(); ++i)
        {
            (*it)->setData(2, 0, i, (*storage_)[(*it)->getNrOfBasisFunctions() * ((*it)->getID()) + i].real());
            (*it)->setData(3, 0, i, (*storage_)[(*it)->getNrOfBasisFunctions() * ((*it)->getID()) + i].imag());
        }
        ierr_ = VecRestoreArrayRead(scaledVec, storage_);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        for (int i = 0; i < 4; ++i)
        {
            PointElementReferenceT p(3);
            (*it)->getReferenceGeometry()->getCenter(p);
            LDOSIntegrand(*it, p, partialResult);
            result[index] = partialResult;
            ++index;
        }
    }
    ierr_ = VecDestroy(&scaledVec);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
}

void hpGemUIExtentions::solveTimeDependant()
{
    std::cout << "doing a time dependent simulation" << std::endl;
    const MaxwellData* actualdata = getData();
    MHasToBeInverted_ = true;
    assembler->fillMatrices(this);
    
    Utilities::GlobalPetscMatrix M_(meshes_[0], 0, -1), S_(meshes_[0], 1, 0);
    std::cout << "GlobalPetscMatrix initialised" << std::endl;
    Utilities::GlobalPetscVector x_(meshes_[0], 0, -1), derivative_(meshes_[0], 1, -1), RHS_(meshes_[0], 2, 0);
    std::cout << "GlobalPetscVector initialised" << std::endl;
    x_.assemble();
    std::cout << "x_ assembled" << std::endl;
    RHS_.assemble();
    std::cout << "RHS_ assembled" << std::endl;
    derivative_.assemble();
    std::cout << "derivative_ assembled" << std::endl;
    
    Vec dummy, dummy2;
    ierr_ = VecDuplicate(derivative_, &dummy);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecDuplicate(derivative_, &dummy2);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = MatMult(M_, derivative_, dummy);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecCopy(dummy, derivative_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = MatMult(M_, x_, dummy);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecCopy(dummy, x_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    
    double tau = 0.05 / ((2 * actualdata->PolynomialOrder_ + 1) * actualdata->NumberOfIntervals_);
    double Nsteps = floor((actualdata->EndTime_ - actualdata->StartTime_) / tau) + 1;
    tau = (actualdata->EndTime_ - actualdata->StartTime_) / Nsteps;
    Nsteps = (actualdata->EndTime_ - actualdata->StartTime_) / tau;
    tau = (actualdata->EndTime_ - actualdata->StartTime_) / Nsteps;
    Nsteps = (actualdata->EndTime_ - actualdata->StartTime_) / tau;
    
    int measureStep = double(Nsteps) / double(getData()->numberOfTimeLevels_ - 1) + 1; //round up
            
    double t = actualdata->StartTime_;
    measureTimes_[0] = t;
    
    double eta = sourceTermTime(t);
    ierr_ = VecScale(RHS_, 1 / eta);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    
    double scale0 = (1 - actualdata->Sigma_ * tau / 2);
    double scale1 = 1 / (1 + actualdata->Sigma_ * tau / 2);
    
    int measureAmount = 0;
    
    std::cout << tau << " " << Nsteps << std::endl;
    
    for (int i = 0; i < Nsteps; ++i)
    {
        if (i % measureStep == 0)
        {
            //getArrayRead only gets local values, so make the entire std::vector local
            /* Vec dummy4;
             IS ISallNumbers;
             int zero[]= {0};
             ierr_=ISCreateBlock(PETSC_COMM_SELF,getData()->numberOfUnknowns_,1,zero,PETSC_COPY_VALUES,&ISallNumbers);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
             ierr_=VecCreateSeq(PETSC_COMM_SELF,getData()->numberOfUnknowns_,&dummy4);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
             ierr_=VecGetSubVector(x_,ISallNumbers,&dummy4);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
             ierr_=VecGetArrayRead(dummy4,storage_);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
             for(ElementIterator it=meshes_[0]->elementColBegin(); it!=meshes_[0]->elementColEnd(); ++it) {
             for(int j=0; j<(*it)->getNrOfBasisFunctions(); ++j) {
             //remnant from before timelevelData was nice and allowed you to set things, should be made less ugly
             const_cast<LinearAlgebra::NumericalVector*>(&(*it)->getTimeLevelData(measureAmount))->operator[](j)=(*storage_)[(*it)->getNrOfBasisFunctions()*((*it)->getID())+j].real();
             }
             */
            x_.writeTimeLevelData(measureAmount);
            //}
            measureTimes_[measureAmount] = t;
            measureAmount++;
            //ierr_=VecRestoreArrayRead(dummy4,storage_);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
        }
        
        //leap-frog sceme (Yee)
        ierr_ = VecAXPY(x_, tau / 2, derivative_);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        ierr_ = MatMult(S_, x_, dummy);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_); //starting here:dummy contaings partial update of x
                
        eta = sourceTermTime(t);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        ierr_ = VecCopy(RHS_, dummy2);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_); //starting here: dummy2 contains time-scaled version of RHS
        ierr_ = VecScale(dummy2, 0.5 * tau * eta);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        ierr_ = VecAYPX(dummy, -tau, dummy2);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        eta = sourceTermTime(t + tau);
        ierr_ = VecCopy(RHS_, dummy2);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        ierr_ = VecAXPY(dummy, 0.5 * tau * eta, dummy2);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        ierr_ = VecScale(dummy, scale1);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        ierr_ = MatMult(M_, dummy, dummy2);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_); //starting here: dummy2 contiains partial update of x
        ierr_ = VecAYPX(derivative_, scale0 * scale1, dummy2);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        ierr_ = VecAXPY(x_, tau / 2, derivative_);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        t = t + tau;
    }
    //getArrayRead only gets local values, so make the entire std::vector local
    
    /*
     Vec dummy4;
     IS ISallNumbers;
     int zero[]= {0};
     ierr_=ISCreateBlock(PETSC_COMM_SELF,getData()->numberOfUnknowns_,1,zero,PETSC_COPY_VALUES,&ISallNumbers);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
     ierr_=VecCreateSeq(PETSC_COMM_SELF,getData()->numberOfUnknowns_,&dummy4);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
     ierr_=VecGetSubVector(x_,ISallNumbers,&dummy4);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
     ierr_=VecGetArrayRead(dummy4,storage_);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
     for(ElementIterator it=meshes_[0]->elementColBegin(); it!=meshes_[0]->elementColEnd(); ++it) {
     for(int j=0; j<(*it)->getNrOfBasisFunctions(); ++j) {
     const_cast<LinearAlgebra::NumericalVector*>(&(*it)->getTimeLevelData(measureAmount))->operator[](j)=(*storage_)[(*it)->getNrOfBasisFunctions()*((*it)->getID())+j].real();
     }
     }
     */
    x_.writeTimeLevelData(measureAmount);
    measureTimes_[measureAmount] = t;
    //VecRestoreArrayRead(dummy4,storage_);
    //ierr_=VecDestroy(&dummy);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
    //ierr_=VecDestroy(&dummy2);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
}

void hpGemUIExtentions::solveHarmonic()
{
    std::cout << "finding a time-harmonic solution" << std::endl;
    const MaxwellData* actualdata = getData();
    
    MHasToBeInverted_ = false;
    assembler->fillMatrices(this);
    std::cout << "Fill Matrices call completed" << std::endl;
    
    Utilities::GlobalPetscMatrix M_(meshes_[0], 0, -1), S_(meshes_[0], 1, 0);
    std::cout << "GlobalPetscMatrix initialised" << std::endl;
    Utilities::GlobalPetscVector x_(meshes_[0], 0, -1), derivative_(meshes_[0], 1, -1), RHS_(meshes_[0], 2, 0);
    std::cout << "GlobalPetscVector initialised" << std::endl;
    x_.assemble();
    std::cout << "x_ assembled" << std::endl;
    RHS_.assemble();
    std::cout << "RHS_ assembled" << std::endl;
    derivative_.assemble();
    std::cout << "derivative_ assembled" << std::endl;
    
    PC preconditioner;
    ierr_ = KSPGetPC(solver_, &preconditioner);
    ierr_ = PCSetType(preconditioner, "jacobi");
    ierr_ = KSPSetPC(solver_, preconditioner);
    
    ierr_ = MatAXPY(S_, -1, M_, SUBSET_NONZERO_PATTERN);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    
    ierr_ = KSPSetTolerances(solver_, 1e-10, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    
    ierr_ = KSPSetType(solver_, "minres");
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    //	ierr=MatConvert(S,"aij",MAT_REUSE_MATRIX,&S);CHKERRABORT(PETSC_COMM_WORLD,ierr);
    
    //everything that is set in the code, but before this line is overridden by command-line options
    ierr_ = KSPSetFromOptions(solver_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    //everything that is set in the code, but after this line overrides the comand-line options
    //Changed due to PETSC-3.5.3
    //ierr_=KSPSetOperators(solver_,S_,S_,DIFFERENT_NONZERO_PATTERN);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
    ierr_ = KSPSetOperators(solver_, S_, S_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = KSPSetUp(solver_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = KSPSolve(solver_, RHS_, x_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    
    x_.writeTimeLevelData(0); //DOUBTFUL
            
    //getArrayRead only gets local values, so make the entire std::vector local
    /*
     Vec dummy;
     IS ISallNumbers;
     int zero[]= {0};
     ierr_=ISCreateBlock(PETSC_COMM_SELF,getData()->numberOfUnknowns_,1,zero,PETSC_COPY_VALUES,&ISallNumbers);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
     ierr_=VecCreateSeq(PETSC_COMM_SELF,getData()->numberOfUnknowns_,&dummy);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
     ierr_=VecGetSubVector(x_,ISallNumbers,&dummy);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
     ierr_=VecGetArrayRead(dummy,storage_);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
     for(ElementIterator it=elementColBegin(); it!=elementColEnd(); ++it) {
     for(int j=0; j<(*it)->getNrOfBasisFunctions(); ++j) {
     const_cast<LinearAlgebra::NumericalVector*>(&(*it)->getTimeLevelData(0))->operator[](j)=(*storage_)[(*it)->getNrOfBasisFunctions()*((*it)->getID())+j].real();
     }
     }
     */
    //measureTimes_[0]=0;
    //ierr_=VecRestoreArrayRead(x_,storage_);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
}

void hpGemUIExtentions::exportMatrixes()
{
    std::cout << "genereting Matlab scripts to load the matrixes" << std::endl;
    MHasToBeInverted_ = false;
    assembler->fillMatrices(this);
    PetscViewer viewM, viewS;
    PetscViewerASCIIOpen(MPI_COMM_WORLD, "M.m", &viewM);
    PetscViewerASCIIOpen(MPI_COMM_WORLD, "S.m", &viewS);
    PetscViewerSetFormat(viewM, PETSC_VIEWER_ASCII_MATLAB);
    PetscViewerSetFormat(viewS, PETSC_VIEWER_ASCII_MATLAB);
    MatView(M_, viewM);
    MatView(S_, viewS);
    PetscViewerDestroy(&viewM);
    PetscViewerDestroy(&viewS);
}

void hpGemUIExtentions::solveEigenvalues()
{
    std::cout << "finding a bunch of eigenvalues" << std::endl;
    const MaxwellData* actualdata = getData();
    int degreesOfFreedomPerElement = configData_->numberOfBasisFunctions_ * configData_->numberOfUnknowns_;
    std::valarray<PetscScalar> blockvalues(degreesOfFreedomPerElement * degreesOfFreedomPerElement);
    int measureAmount = 0;
    
    std::vector<IS> xboundaryRow, xboundaryCol, yboundaryRow, yboundaryCol, zboundaryRow, zboundaryCol;
    findBoundaryBlocks(xboundaryRow, xboundaryCol, yboundaryRow, yboundaryCol, zboundaryRow, zboundaryCol);
    
    //For IP-DG solving a general eigenproblem is slightly faster,
    //but the Brezzi formulation needs an inverse of each block in the mass matrix anyway
    //so there the general eigenproblem is just more work
    MHasToBeInverted_ = true;
    assembler->fillMatrices(this); //This part has to be rewritten by using the GlobalPetscMatrix and GlobalPetscVector from Utilities.
            
    Utilities::GlobalPetscMatrix M_(meshes_[0], 0, -1), S_(meshes_[0], 1, 0);
    std::cout << "GlobalPetscMatrix initialised" << std::endl;
    Utilities::GlobalPetscVector x_(meshes_[0], 0, -1), derivative_(meshes_[0], 1, -1), RHS_(meshes_[0], 2, 0);
    std::cout << "GlobalPetscVector initialised" << std::endl;
    x_.assemble();
    std::cout << "x_ assembled" << std::endl;
    RHS_.assemble();
    std::cout << "RHS_ assembled" << std::endl;
    derivative_.assemble();
    std::cout << "derivative_ assembled" << std::endl;
    
    // optionally divide the eigenvalues by pi^2 (for debugging / visual error tracking)
    //ierr_=MatScale(M_,1/M_PI/M_PI);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
    
    //The eigensolver doesn't like block structures
    Mat product, dummy;
    ierr_ = MatMatMult(M_, S_, MAT_INITIAL_MATRIX, 1.0, &product);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    
    //ierr_=EPSSetProblemType(eigenSolver_,EPS_HEP);
    ierr_ = EPSSetOperators(eigenSolver_, product, NULL);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = EPSSetUp(eigenSolver_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = EPSSetDimensions(eigenSolver_, 24, PETSC_DECIDE, PETSC_DECIDE);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    // 	ierr_=EPSSetType(eigenSolver_,"jd");CHKERRABORT(PETSC_COMM_WORLD,ierr_);
    // 	ST transformation;
    // 	KSP spectralProblem;
    // 	PC preconditioner;
    // 	ierr_=EPSGetST(eigenSolver_,&transformation);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
    // 	ierr_=STSetType(transformation,"precond");CHKERRABORT(PETSC_COMM_WORLD,ierr_);
    // 	ierr_=STGetKSP(transformation,&spectralProblem);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
    // 	ierr_=KSPGetPC(spectralProblem,&preconditioner);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
    // 	ierr_=PCSetType(preconditioner,"ilu");CHKERRABORT(PETSC_COMM_WORLD,ierr_);
    // 	ierr_=EPSSetST(eigenSolver_,transformation);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
    
    //everything that is set in the code, but before this line is overridden by command-line options
    ierr_ = EPSSetFromOptions(eigenSolver_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    //everything that is set in the code, but after this line overrides the comand-line options
    
    ierr_ = EPSSolve(eigenSolver_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    
    //if you are only interested in the default set of eigenvalues use this bit of code for output
    //         ierr_=EPSPrintSolution(eigenSolver_,0);
    //         CHKERRABORT(PETSC_COMM_WORLD,ierr_);
    //
    // 	ierr_=EPSGetEigenvector(eigenSolver_,1,x_,derivative_);//actually eps,nr,real(x),imag(x)
    // 	ierr_=VecGetArrayRead(x_,storage_);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
    // 	for(ElementIterator it=meshes_[0]->elementColBegin();it!=meshes_[0]->elementColEnd();++it){
    // 	    for(int j=0;j<(*it)->getNrOfBasisFunctions();++j){
    // 		const_cast<LinearAlgebra::Matrix*>(&(*it)->getTimeLevelData(0))->operator[](j)=(*storage_)[(*it)->getNrOfBasisFunctions()*((*it)->getID())+j].real();
    // 	    }
    // 	}
    // 	measureTimes_[0]=0;
    // 	ierr_=VecRestoreArrayRead(x_,storage_);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
    
    //if you want all possible eigenvalues (including using the wave-vector) use this 'bit' of code
    //create some storage std::vectors
    PetscScalar neededOnlyForRealPetsc, eigenvalue;
    Vec *eigenvalues, example, *eigenVectors;
    eigenvalues = new Vec[24];
    eigenVectors = new Vec[40]; //a few extra in case SLEPc finds more than the requested amount of eigenvalues
    ierr_ = VecCreate(PETSC_COMM_WORLD, &example);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecSetSizes(example, PETSC_DECIDE, 61);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecSetUp(example);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecDuplicateVecs(example, 24, &eigenvalues);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecDuplicateVecs(x_, 40, &eigenVectors);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecDestroy(&example);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    int converged;
    ierr_ = EPSGetConverged(eigenSolver_, &converged);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    for (int i = 0; i < converged && i < 24; ++i)
    {
        ierr_ = EPSGetEigenvalue(eigenSolver_, i, &eigenvalue, &neededOnlyForRealPetsc);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        ierr_ = VecSetValue(eigenvalues[i], 0, eigenvalue, ADD_VALUES);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    }
    
    Vec waveVec, waveVecConjugate;
    const int *rows, *columns;
    ierr_ = VecCreate(PETSC_COMM_WORLD, &waveVec);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecSetType(waveVec, "mpi");
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecSetBlockSize(waveVec, configData_->numberOfBasisFunctions_);
    ierr_ = VecSetSizes(waveVec, PETSC_DECIDE, globalData_->numberOfUnknowns_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecSetUp(waveVec);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    LinearAlgebra::MiddleSizeVector k(3);
    k[0] = M_PI / 20.;
    k[1] = 0;
    k[2] = 0;
    makeShiftMatrix(k, waveVec);
    ierr_ = VecAssemblyBegin(waveVec);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecAssemblyEnd(waveVec);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecDuplicate(waveVec, &waveVecConjugate);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecCopy(waveVec, waveVecConjugate);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecConjugate(waveVecConjugate);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    
    PetscLayout distribution;
    int localProcessorNumber, blockProcessorNumber, localsize;
    MPI_Comm_rank(PETSC_COMM_WORLD, &localProcessorNumber);
    ierr_ = PetscLayoutCreate(PETSC_COMM_WORLD, &distribution);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = MatGetLocalSize(product, &localsize, NULL);
    ierr_ = PetscLayoutSetLocalSize(distribution, localsize);
    ierr_ = PetscLayoutSetUp(distribution);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    
    for (int i = 1; i < 61; ++i)
    {
        std::cout << i << std::endl;
        if (i == 21)
        {
            //these are only increments, actually this make the wavevector move from pi,0,0 to pi,pi,0
            k[0] = 0;
            k[1] = M_PI / 20.;
            
            //recompute the shifts
            makeShiftMatrix(k, waveVec);
            ierr_ = VecAssemblyBegin(waveVec);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            ierr_ = VecAssemblyEnd(waveVec);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            ierr_ = VecCopy(waveVec, waveVecConjugate);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            ierr_ = VecConjugate(waveVecConjugate);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        }
        else if (i == 41)
        {
            k[1] = 0;
            k[2] = M_PI / 20.;
            makeShiftMatrix(k, waveVec);
            ierr_ = VecAssemblyBegin(waveVec);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            ierr_ = VecAssemblyEnd(waveVec);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            ierr_ = VecCopy(waveVec, waveVecConjugate);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            ierr_ = VecConjugate(waveVecConjugate);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        }
        ierr_ = EPSGetInvariantSubspace(eigenSolver_, eigenVectors);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        ierr_ = MatDiagonalScale(product, waveVec, waveVecConjugate);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        
        //this is probably not the best way to do this
        //a cleaner set-up can be made using two iterators
        for (int j = 0; j < xboundaryCol.size(); ++j)
        {
            ierr_ = ISGetIndices(xboundaryRow[j], &rows);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            ierr_ = ISGetIndices(xboundaryCol[j], &columns);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            PetscLayoutFindOwner(distribution, rows[0], &blockProcessorNumber);
            //if(blockProcessorNumber==localProcessorNumber){
            ierr_ = MatGetValues(product, degreesOfFreedomPerElement, rows, degreesOfFreedomPerElement, columns, &blockvalues[0]);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            blockvalues *= exp(std::complex<double>(0, k[0] * (j % 2 == 0 ? -1 : 1)));
            ierr_ = MatSetValues(product, degreesOfFreedomPerElement, rows, degreesOfFreedomPerElement, columns, &blockvalues[0], INSERT_VALUES);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            //}
            ierr_ = MatAssemblyBegin(product, MAT_FINAL_ASSEMBLY);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            ierr_ = MatAssemblyEnd(product, MAT_FINAL_ASSEMBLY);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            ierr_ = ISRestoreIndices(xboundaryRow[j], &rows);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            ierr_ = ISRestoreIndices(xboundaryCol[j], &columns);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        }
        for (int j = 0; j < yboundaryCol.size(); ++j)
        {
            ierr_ = ISGetIndices(yboundaryRow[j], &rows);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            ierr_ = ISGetIndices(yboundaryCol[j], &columns);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            PetscLayoutFindOwner(distribution, rows[0], &blockProcessorNumber);
            //if(blockProcessorNumber==localProcessorNumber){
            ierr_ = MatGetValues(product, degreesOfFreedomPerElement, rows, degreesOfFreedomPerElement, columns, &blockvalues[0]);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            blockvalues *= exp(std::complex<double>(0, k[1] * (j % 2 == 0 ? -1 : 1)));
            ierr_ = MatSetValues(product, degreesOfFreedomPerElement, rows, degreesOfFreedomPerElement, columns, &blockvalues[0], INSERT_VALUES);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            //}
            ierr_ = MatAssemblyBegin(product, MAT_FINAL_ASSEMBLY);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            ierr_ = MatAssemblyEnd(product, MAT_FINAL_ASSEMBLY);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            ierr_ = ISRestoreIndices(yboundaryRow[j], &rows);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            ierr_ = ISRestoreIndices(yboundaryCol[j], &columns);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        }
        for (int j = 0; j < zboundaryCol.size(); ++j)
        {
            ierr_ = ISGetIndices(zboundaryRow[j], &rows);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            ierr_ = ISGetIndices(zboundaryCol[j], &columns);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            PetscLayoutFindOwner(distribution, rows[0], &blockProcessorNumber);
            //if(blockProcessorNumber==localProcessorNumber){
            ierr_ = MatGetValues(product, degreesOfFreedomPerElement, rows, degreesOfFreedomPerElement, columns, &blockvalues[0]);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            blockvalues *= exp(std::complex<double>(0, k[2] * (j % 2 == 0 ? -1 : 1)));
            ierr_ = MatSetValues(product, degreesOfFreedomPerElement, rows, degreesOfFreedomPerElement, columns, &blockvalues[0], INSERT_VALUES);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            //}
            ierr_ = MatAssemblyBegin(product, MAT_FINAL_ASSEMBLY);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            ierr_ = MatAssemblyEnd(product, MAT_FINAL_ASSEMBLY);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            ierr_ = ISRestoreIndices(zboundaryRow[j], &rows);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            ierr_ = ISRestoreIndices(zboundaryCol[j], &columns);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        }
        
        //outputs 'old' data
        if (i % 20 == 1)
        {
            //getArrayRead only gets local values, so make the entire std::vector local
            // Vec dummy4;
            // IS ISallNumbers;
            // int zero[]= {0};
            // ierr_=ISCreateBlock(PETSC_COMM_SELF,getData()->numberOfUnknowns_,1,zero,PETSC_COPY_VALUES,&ISallNumbers);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
            // ierr_=VecCreateSeq(PETSC_COMM_SELF,getData()->numberOfUnknowns_,&dummy4);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
            //for(int j=0;j<10&&j<converged;++j){
            //     ierr_=VecGetSubVector(eigenVectors[j],ISallNumbers,&dummy4);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
            //     ierr_=VecGetArrayRead(dummy4,storage_);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
            //     for(ElementIterator it=meshes_[0]->elementColBegin(); it!=meshes_[0]->elementColEnd(); ++it) {
            //         for(int k=0; k<(*it)->getNrOfBasisFunctions(); ++k) {
            //             const_cast<LinearAlgebra::NumericalVector*>(&(*it)->getTimeLevelData(measureAmount))->operator[](k)=(*storage_)[(*it)->getNrOfBasisFunctions()*((*it)->getID())+k].real();
            //         }
            //     }
            //     ierr_=VecRestoreArrayRead(dummy4,storage_);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
            
            x_.writeTimeLevelData(measureAmount);
            measureTimes_[measureAmount] = i / 20;
            measureAmount++;
            //}
        }
        
        ierr_ = EPSSetOperators(eigenSolver_, product, NULL);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        ierr_ = EPSSetInitialSpace(eigenSolver_, converged, eigenVectors);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        ierr_ = EPSSetUp(eigenSolver_);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        ierr_ = EPSSolve(eigenSolver_);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        ierr_ = EPSGetConverged(eigenSolver_, &converged);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        
        for (int j = 0; j < converged && j < 24; ++j)
        {
            ierr_ = EPSGetEigenvalue(eigenSolver_, j, &eigenvalue, &neededOnlyForRealPetsc);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            ierr_ = VecSetValue(eigenvalues[j], i, eigenvalue, INSERT_VALUES);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        }
    }
    for (int i = 0; i < xboundaryCol.size(); ++i)
    {
        ISDestroy(&xboundaryCol[i]);
        ISDestroy(&xboundaryRow[i]);
    }
    for (int i = 0; i < yboundaryCol.size(); ++i)
    {
        ISDestroy(&yboundaryCol[i]);
        ISDestroy(&yboundaryRow[i]);
    }
    for (int i = 0; i < zboundaryCol.size(); ++i)
    {
        ISDestroy(&zboundaryCol[i]);
        ISDestroy(&zboundaryRow[i]);
    }
    for (int i = 0; i < 20; ++i)
    {
        ierr_ = VecAssemblyBegin(eigenvalues[i]);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        ierr_ = VecAssemblyEnd(eigenvalues[i]);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        ierr_ = VecView(eigenvalues[i], 0);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    }
    
    ierr_ = EPSGetInvariantSubspace(eigenSolver_, eigenVectors);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    
    //Vec dummy4;
    //IS ISallNumbers;
    //int zero[]= {0};
    //ierr_=ISCreateBlock(PETSC_COMM_SELF,getData()->numberOfUnknowns_,1,zero,PETSC_COPY_VALUES,&ISallNumbers);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
    //ierr_=VecCreateSeq(PETSC_COMM_SELF,getData()->numberOfUnknowns_,&dummy4);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
    //  for(int j=0;j<10&&j<converged;++j){
    //    ierr_=VecGetSubVector(eigenVectors[j],ISallNumbers,&dummy4);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
    //    ierr_=VecGetArrayRead(dummy4,storage_);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
    //    for(ElementIterator it=meshes_[0]->elementColBegin(); it!=meshes_[0]->elementColEnd(); ++it) {
    //        for(int k=0; k<(*it)->getNrOfBasisFunctions(); ++k) {
    //            const_cast<LinearAlgebra::NumericalVector*>(&(*it)->getTimeLevelData(measureAmount))->operator[](k)=(*storage_)[(*it)->getNrOfBasisFunctions()*((*it)->getID())+k].real();
    //        }
    //    }
    //    ierr_=VecRestoreArrayRead(dummy4,storage_);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
    x_.writeTimeLevelData(measureAmount);
    measureTimes_[measureAmount] = 3;
    measureAmount++;
    // }
    
    ierr_ = VecDestroy(&waveVec);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecDestroy(&waveVecConjugate);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    
    //always clean up after you are done
    ierr_ = MatDestroy(&product);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
}

void hpGemUIExtentions::solveDOS()
{
    std::cout << "finding the local density of states" << std::endl;
    const MaxwellData* actualdata = getData();
    int degreesOfFreedomPerElement = configData_->numberOfBasisFunctions_;
    std::valarray<PetscScalar> blockvalues(degreesOfFreedomPerElement * degreesOfFreedomPerElement);
    int measureAmount = 0;
    std::vector<double> eigenvalues;
    std::vector<LinearAlgebra::MiddleSizeVector> functionValues;
    std::vector<IS> xboundaryRow, xboundaryCol, yboundaryRow, yboundaryCol, zboundaryRow, zboundaryCol;
    findBoundaryBlocks(xboundaryRow, xboundaryCol, yboundaryRow, yboundaryCol, zboundaryRow, zboundaryCol);
    KspaceData brillouinZone(5);
    
    //For IP-DG solving a general eigenproblem is slightly faster,
    //but the Brezzi formulation needs an inverse of each block in the mass matrix anyway
    //so there the general eigenproblem is just more work
    MHasToBeInverted_ = true;
    assembler->fillMatrices(this);
    
    Utilities::GlobalPetscMatrix M_(meshes_[0], 0, -1), S_(meshes_[0], 1, 0);
    std::cout << "GlobalPetscMatrix initialised" << std::endl;
    Utilities::GlobalPetscVector x_(meshes_[0], 0, -1), derivative_(meshes_[0], 1, -1), RHS_(meshes_[0], 2, 0);
    std::cout << "GlobalPetscVector initialised" << std::endl;
    x_.assemble();
    std::cout << "x_ assembled" << std::endl;
    RHS_.assemble();
    std::cout << "RHS_ assembled" << std::endl;
    derivative_.assemble();
    std::cout << "derivative_ assembled" << std::endl;
    
    //Divide the eigenvalues by pi^2 (because some people cant do that by heart and the useful info shouldn't be hidden)
    //ierr_=MatScale(M_,1/M_PI/M_PI);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
    
    //The eigensolver doesn't like block structures
    //ierr_=MatConvert(M_,"aij",MAT_REUSE_MATRIX,&M_);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
    //ierr_=MatConvert(S_,"aij",MAT_REUSE_MATRIX,&S_);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
    Mat product, dummy;
    ierr_ = MatMatMult(M_, S_, MAT_INITIAL_MATRIX, 1.0, &product);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    
    //    ierr_=EPSSetProblemType(eigenSolver_,EPS_HEP);
    ierr_ = EPSSetOperators(eigenSolver_, product, NULL);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = EPSSetUp(eigenSolver_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = EPSSetDimensions(eigenSolver_, 24, PETSC_DECIDE, PETSC_DECIDE);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    // 	ierr_=EPSSetType(eigenSolver_,"jd");CHKERRABORT(PETSC_COMM_WORLD,ierr_);
    // 	ST transformation;
    // 	KSP spectralProblem;
    // 	PC preconditioner;
    // 	ierr_=EPSGetST(eigenSolver_,&transformation);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
    // 	ierr_=STSetType(transformation,"precond");CHKERRABORT(PETSC_COMM_WORLD,ierr_);
    // 	ierr_=STGetKSP(transformation,&spectralProblem);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
    // 	ierr_=KSPGetPC(spectralProblem,&preconditioner);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
    // 	ierr_=PCSetType(preconditioner,"ilu");CHKERRABORT(PETSC_COMM_WORLD,ierr_);
    // 	ierr_=EPSSetST(eigenSolver_,transformation);CHKERRABORT(PETSC_COMM_WORLD,ierr_);
    
    //everything that is set in the code, but before this line is overridden by command-line options
    ierr_ = EPSSetFromOptions(eigenSolver_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    //everything that is set in the code, but after this line overrides the comand-line options
    
    PetscScalar target;
    ierr_ = EPSGetTarget(eigenSolver_, &target);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    
    ierr_ = EPSSolve(eigenSolver_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    
    PetscScalar neededOnlyForRealPetsc, eigenvalue;
    Vec example, *eigenVectors;
    eigenVectors = new Vec[40]; //a few extra in case SLEPc finds more than the requested amount of eigenvalues
    ierr_ = VecCreate(PETSC_COMM_WORLD, &example);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecSetSizes(example, PETSC_DECIDE, 61);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecSetUp(example);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecDuplicateVecs(x_, 40, &eigenVectors);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecDestroy(&example);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    int converged;
    ierr_ = EPSGetConverged(eigenSolver_, &converged);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    
    ierr_ = EPSGetInvariantSubspace(eigenSolver_, eigenVectors);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    for (int i = 0; i < converged && i < 24; ++i)
    {
        ierr_ = EPSGetEigenvalue(eigenSolver_, i, &eigenvalue, &neededOnlyForRealPetsc);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        //SLEPSc sorts the eigenvalues by closest to the target first, this reorders them for the better smallest first
        //if SLEPSc manages to find some of the zero eigenvalues they should be skipped
        if (eigenvalue.real() > 1e-6)
        {
            LinearAlgebra::MiddleSizeVector functionvalue;
            //makeFunctionValue(eigenVectors[i],functionvalue);
            if (eigenvalue.real() < target.real())
            {
                eigenvalues.insert(eigenvalues.begin(), sqrt(eigenvalue.real()));
                functionValues.insert(functionValues.begin(), functionvalue);
            }
            else
            {
                eigenvalues.push_back(sqrt(eigenvalue.real()));
                functionValues.push_back(functionvalue);
            }
        }
    }
    //eigenvalues.insert(eigenvalues.begin(),2,0);//at the infinite wavelength limit the two constant functions also contribute to the DOS
    brillouinZone.setOmega(eigenvalues);
    //brillouinZone.setFunctionValues(functionValues);
    
    Vec waveVec, waveVecConjugate;
    const int *rows, *columns;
    ierr_ = VecCreate(PETSC_COMM_WORLD, &waveVec);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecSetType(waveVec, "mpi");
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecSetBlockSize(waveVec, configData_->numberOfBasisFunctions_);
    ierr_ = VecSetSizes(waveVec, PETSC_DECIDE, globalData_->numberOfUnknowns_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecSetUp(waveVec);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    LinearAlgebra::MiddleSizeVector k(3);
    k[0] = M_PI / 60;
    k[1] = 0;
    k[2] = 0;
    makeShiftMatrix(k, waveVec);
    ierr_ = VecAssemblyBegin(waveVec);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecAssemblyEnd(waveVec);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecDuplicate(waveVec, &waveVecConjugate);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecCopy(waveVec, waveVecConjugate);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecConjugate(waveVecConjugate);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    
    PetscLayout distribution;
    int localProcessorNumber, blockProcessorNumber, localsize;
    MPI_Comm_rank(PETSC_COMM_WORLD, &localProcessorNumber);
    ierr_ = PetscLayoutCreate(PETSC_COMM_WORLD, &distribution);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = MatGetLocalSize(product, &localsize, NULL);
    ierr_ = PetscLayoutSetLocalSize(distribution, localsize);
    ierr_ = PetscLayoutSetUp(distribution);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    
    while (brillouinZone.hasNextPoint())
    {
        k = brillouinZone.nextPoint();
        ierr_ = MatDiagonalScale(product, waveVec, waveVecConjugate);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        
        //find a better strategy -- at least there is no communication in this set-up so MatAssemblyBegin should be reasonably fast
        //a cleaner set-up can be made using two iterators
        for (int j = 0; j < xboundaryCol.size(); ++j)
        {
            ierr_ = ISGetIndices(xboundaryRow[j], &rows);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            ierr_ = ISGetIndices(xboundaryCol[j], &columns);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            PetscLayoutFindOwner(distribution, rows[0], &blockProcessorNumber);
            if (blockProcessorNumber == localProcessorNumber)
            {
                ierr_ = MatGetValues(product, degreesOfFreedomPerElement, rows, degreesOfFreedomPerElement, columns, &blockvalues[0]);
                CHKERRABORT(PETSC_COMM_WORLD, ierr_);
                blockvalues *= exp(std::complex<double>(0, k[0] * (j % 2 == 0 ? -1 : 1)));
                ierr_ = MatSetValues(product, degreesOfFreedomPerElement, rows, degreesOfFreedomPerElement, columns, &blockvalues[0], INSERT_VALUES);
                CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            }
            ierr_ = MatAssemblyBegin(product, MAT_FINAL_ASSEMBLY);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            ierr_ = MatAssemblyEnd(product, MAT_FINAL_ASSEMBLY);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            ierr_ = ISRestoreIndices(xboundaryRow[j], &rows);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            ierr_ = ISRestoreIndices(xboundaryCol[j], &columns);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        }
        for (int j = 0; j < yboundaryCol.size(); ++j)
        {
            ierr_ = ISGetIndices(yboundaryRow[j], &rows);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            ierr_ = ISGetIndices(yboundaryCol[j], &columns);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            PetscLayoutFindOwner(distribution, rows[0], &blockProcessorNumber);
            if (blockProcessorNumber == localProcessorNumber)
            {
                ierr_ = MatGetValues(product, degreesOfFreedomPerElement, rows, degreesOfFreedomPerElement, columns, &blockvalues[0]);
                CHKERRABORT(PETSC_COMM_WORLD, ierr_);
                blockvalues *= exp(std::complex<double>(0, k[1] * (j % 2 == 0 ? -1 : 1)));
                ierr_ = MatSetValues(product, degreesOfFreedomPerElement, rows, degreesOfFreedomPerElement, columns, &blockvalues[0], INSERT_VALUES);
                CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            }
            ierr_ = MatAssemblyBegin(product, MAT_FINAL_ASSEMBLY);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            ierr_ = MatAssemblyEnd(product, MAT_FINAL_ASSEMBLY);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            ierr_ = ISRestoreIndices(yboundaryRow[j], &rows);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            ierr_ = ISRestoreIndices(yboundaryCol[j], &columns);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        }
        for (int j = 0; j < zboundaryCol.size(); ++j)
        {
            ierr_ = ISGetIndices(zboundaryRow[j], &rows);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            ierr_ = ISGetIndices(zboundaryCol[j], &columns);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            PetscLayoutFindOwner(distribution, rows[0], &blockProcessorNumber);
            if (blockProcessorNumber == localProcessorNumber)
            {
                ierr_ = MatGetValues(product, degreesOfFreedomPerElement, rows, degreesOfFreedomPerElement, columns, &blockvalues[0]);
                CHKERRABORT(PETSC_COMM_WORLD, ierr_);
                blockvalues *= exp(std::complex<double>(0, k[2] * (j % 2 == 0 ? -1 : 1)));
                ierr_ = MatSetValues(product, degreesOfFreedomPerElement, rows, degreesOfFreedomPerElement, columns, &blockvalues[0], INSERT_VALUES);
                CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            }
            ierr_ = MatAssemblyBegin(product, MAT_FINAL_ASSEMBLY);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            ierr_ = MatAssemblyEnd(product, MAT_FINAL_ASSEMBLY);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            ierr_ = ISRestoreIndices(zboundaryRow[j], &rows);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            ierr_ = ISRestoreIndices(zboundaryCol[j], &columns);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        }
        
        ierr_ = EPSSetOperators(eigenSolver_, product, NULL);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        ierr_ = EPSSetInitialSpace(eigenSolver_, converged, eigenVectors);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        ierr_ = EPSSetUp(eigenSolver_);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        ierr_ = EPSSolve(eigenSolver_);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        ierr_ = EPSGetInvariantSubspace(eigenSolver_, eigenVectors);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        ierr_ = EPSGetConverged(eigenSolver_, &converged);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
        
        eigenvalues.clear();
        functionValues.clear();
        for (int j = 0; j < converged && j < 24; ++j)
        {
            ierr_ = EPSGetEigenvalue(eigenSolver_, j, &eigenvalue, &neededOnlyForRealPetsc);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            if (eigenvalue.real() > 1e-6)
            {
                LinearAlgebra::MiddleSizeVector functionvalue;
                //makeFunctionValue(eigenVectors[j],functionvalue);
                if (eigenvalue.real() < target.real())
                {
                    eigenvalues.insert(eigenvalues.begin(), sqrt(eigenvalue.real()));
                    functionValues.insert(functionValues.begin(), functionvalue);
                }
                else
                {
                    eigenvalues.push_back(sqrt(eigenvalue.real()));
                    functionValues.push_back(functionvalue);
                }
            }
        }
        brillouinZone.setOmega(eigenvalues);
        //brillouinZone.setFunctionValues(functionValues);
    }
    for (int i = 0; i < xboundaryCol.size(); ++i)
    {
        ISDestroy(&xboundaryCol[i]);
        ISDestroy(&xboundaryRow[i]);
    }
    for (int i = 0; i < yboundaryCol.size(); ++i)
    {
        ISDestroy(&yboundaryCol[i]);
        ISDestroy(&yboundaryRow[i]);
    }
    for (int i = 0; i < zboundaryCol.size(); ++i)
    {
        ISDestroy(&zboundaryCol[i]);
        ISDestroy(&zboundaryRow[i]);
    }
    
    ierr_ = VecDestroy(&waveVec);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecDestroy(&waveVecConjugate);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    
    //always clean up after you are done
    ierr_ = MatDestroy(&product);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    
    LinearAlgebra::MiddleSizeVector result(meshes_[0]->getNumberOfElements());
    for (int i = 0; i < 1001; ++i)
    {
        brillouinZone.getIntegral(double(i) / 100., result);
        std::cout << result[0] << " " << result[6] << " " << result[10] << std::endl;
        result[0] = 0;
        //result[6]=0;
        //result[10]=0;
    }
}

//this is where you specify an initial condition
void hpGemUIExtentions::initialConditions(const PointPhysicalT& p, LinearAlgebra::MiddleSizeVector& ret)
{
    hpGemUIExtentions::exactSolution(p, 0, ret);
}

void hpGemUIExtentions::initialConditionsDeriv(const PointPhysicalT& p, LinearAlgebra::MiddleSizeVector& ret)
{
    ret[0] = 0;
    ret[1] = 0;
    ret[2] = 0;
}

/**
 * this is where you specify the spatial part of the source Term
 * assumes that the source term can be split is a spatial part and a time part
 */

void hpGemUIExtentions::sourceTerm(const PointPhysicalT& p, LinearAlgebra::MiddleSizeVector& ret)
{
    //  ret[0]=0;
    // ret[1]=0;
    // ret[2]=0;
    hpGemUIExtentions::exactSolution(p, 0, ret);
    // 	ret*=-1;
    //ret*=M_PI*M_PI*8-1;
    ret *= M_PI * M_PI * 2 - 1;
}

void hpGemUIExtentions::initialExactSolution(const PointPhysicalT& p, LinearAlgebra::MiddleSizeVector &ret)
{
    //ret[0]=sin(M_PI*2*p[1])*sin(M_PI*2*p[2]);
    //ret[1]=sin(M_PI*2*p[2])*sin(M_PI*2*p[0]);
    //ret[2]=sin(M_PI*2*p[0])*sin(M_PI*2*p[1]);
    ret[0] = sin(M_PI * p[1]) * sin(M_PI * p[2]);
    ret[1] = sin(M_PI * p[2]) * sin(M_PI * p[0]);
    ret[2] = sin(M_PI * p[0]) * sin(M_PI * p[1]);
    //            ret[0]=p[0]*(1-p[0]);
    //            ret[1]=0;
    // 	  		  ret[2]=0;
}

void hpGemUIExtentions::boundaryConditions(const Geometry::PointPhysical &p, LinearAlgebra::MiddleSizeVector &ret)
{
    initialExactSolution(p, ret);
}

void hpGemUIExtentions::anonymous1::elementIntegrand(const Base::HpgemUI::ElementT* element, const PointElementReferenceT& p, LinearAlgebra::MiddleSizeMatrix& ret)
{
    //std::cout<<"Anonymous 1 called for element integration"<<std::endl;
    ret.resize(element->getNrOfBasisFunctions(), element->getNrOfBasisFunctions());
    //cout<<"\nIn the element integrand for the mass matrix for element id: "<<element->getID();
    ElementInfos* info = static_cast<ElementInfos*>(const_cast<ElementT*>(element)->getUserData());
    LinearAlgebra::MiddleSizeVector phi_i(3), phi_j(3);
    //std::vector<LinearAlgebra::NumericalVector> functionValues;
    //info->makeFunctionValuesVector(element,p,functionValues);
    for (int i = 0; i < element->getNrOfBasisFunctions(); ++i)
    {
        //phi_i=functionValues[i];
        element->basisFunction(i, p, phi_i);
        for (int j = 0; j < element->getNrOfBasisFunctions(); ++j)
        {
            //phi_j=functionValues[j];
            element->basisFunction(j, p, phi_j);
            ret(i, j) = (phi_i[0] * phi_j[0] + phi_i[1] * phi_j[1] + phi_i[2] * phi_j[2]) * info->epsilon_;
            ret(j, i) = ret(i, j);
            //std::cout<<ret(i, j)<<std::endl;
        }
    }
}

void hpGemUIExtentions::anonymous2::elementIntegrand(const ElementT* element, const PointElementReferenceT& p, LinearAlgebra::MiddleSizeMatrix& ret)
{
    //cout<<"\nIn the element integrand for the stiffness matrix for element id: "<<element->getID();
    ret.resize(element->getNrOfBasisFunctions(), element->getNrOfBasisFunctions());
    //ElementInfos* info = static_cast<ElementInfos*> (const_cast<ElementT*>(element)->getUserData());
    LinearAlgebra::MiddleSizeVector phi_i(3), phi_j(3);
    //std::vector<LinearAlgebra::NumericalVector> functionCurls;
    //info->makeFunctionCurlsVector(element,p,functionCurls);
    for (int i = 0; i < element->getNrOfBasisFunctions(); ++i)
    {
        //phi_i=functionCurls[i];
        element->basisFunctionCurl(i, p, phi_i);
        for (int j = i; j < element->getNrOfBasisFunctions(); ++j)
        {
            //phi_j=functionCurls[j];
            element->basisFunctionCurl(j, p, phi_j);
            ret(i, j) = phi_i[0] * phi_j[0] + phi_i[1] * phi_j[1] + phi_i[2] * phi_j[2];
            ret(j, i) = ret(i, j);
            //std::cout<<ret(i, j)<<std::endl;
        }
    }
}

void hpGemUIExtentions::anonymous3::elementIntegrand(const ElementT* element, const PointElementReferenceT& p, LinearAlgebra::MiddleSizeVector& ret)
{
    ret.resize(element->getNrOfBasisFunctions());
    //ElementInfos* info = static_cast<ElementInfos*>(const_cast<ElementT*>(element)->getUserData());
    PointPhysicalT pPhys(3);
    element->referenceToPhysical(p, pPhys);
    LinearAlgebra::MiddleSizeVector val(3), phi(3);
    //std::vector<LinearAlgebra::NumericalVector> functionValues;
    //info->makeFunctionValuesVector(element,p,functionValues);
    sourceTerm(pPhys, val);
    for (int i = 0; i < element->getNrOfBasisFunctions(); ++i)
    {
        //phi=functionValues[i];
        element->basisFunction(i, p, phi);
        ret(i) = phi[0] * val[0] + phi[1] * val[1] + phi[2] * val[2];
        //std::cout<<ret(i)<<std::endl;
    }
}

void hpGemUIExtentions::anonymous4::elementIntegrand(const Base::HpgemUI::ElementT* element, const PointElementReferenceT& p, LinearAlgebra::MiddleSizeVector& ret)
{
    
    ret.resize(element->getNrOfBasisFunctions());
    PointPhysicalT pPhys(3);
    element->referenceToPhysical(p, pPhys);
    LinearAlgebra::MiddleSizeVector val(3), phi(3);
    initialConditions(pPhys, val);
    
    for (int i = 0; i < element->getNrOfBasisFunctions(); ++i)
    {
        element->basisFunction(i, p, phi);
        
        ret(i) = phi[0] * val[0] + phi[1] * val[1] + phi[2] * val[2];
        
    }
    //Called twice for fill matrices of both Initial Conditions as well as Initial Conditions Deriv.
}

void hpGemUIExtentions::anonymous5::elementIntegrand(const Base::HpgemUI::ElementT* element, const PointElementReferenceT& p, LinearAlgebra::MiddleSizeVector& ret)
{
    ret.resize(element->getNrOfBasisFunctions());
    //ElementInfos* info = static_cast<ElementInfos*>(const_cast<ElementT*>(element)->getUserData());
    PointPhysicalT pPhys(3);
    element->referenceToPhysical(p, pPhys);
    //std::vector<LinearAlgebra::NumericalVector> functionValues;
    //info->makeFunctionValuesVector(element,p,functionValues);
    LinearAlgebra::MiddleSizeVector val(3), phi(3);
    initialConditionsDeriv(pPhys, val);
    for (int i = 0; i < element->getNrOfBasisFunctions(); ++i)
    {
        //phi=functionValues[i];
        element->basisFunction(i, p, phi);
        ret(i) = phi[0] * val[0] + phi[1] * val[1] + phi[2] * val[2];
        //std::cout<<ret(i)<<std::endl;
    }
}

void hpGemUIExtentions::anonymous6::faceIntegrand(const FaceT* face, const LinearAlgebra::MiddleSizeVector& normal, const PointFaceReferenceT& p, LinearAlgebra::MiddleSizeMatrix& ret)
{
    
    //int n = face->getPtrElementLeft()->getNrOfBasisFunctions();
    //int M = face->getPtrElementLeft()->getNrOfBasisFunctions() + face->getPtrElementRight()->getNrOfBasisFunctions();
    
    int M = face->getNrOfBasisFunctions();
    
    ret.resize(M, M);
    LinearAlgebra::MiddleSizeVector phi_i_normal(3), phi_j_normal(3), phi_i_curl(3), phi_j_curl(3);
    for (int i = 0; i < M; ++i)
    {
        face->basisFunctionCurl(i, p, phi_i_curl);
        face->basisFunctionNormal(i, normal, p, phi_i_normal);
        
        for (int j = i; j < M; ++j)
        {
            face->basisFunctionCurl(j, p, phi_j_curl);
            face->basisFunctionNormal(j, normal, p, phi_j_normal);
            //std::cout<<phi_j_curl<<std::endl;
            //std::cout<<normal<<std::endl;
            
            ret(i, j) = -(face->isInternal() ? 0.5 : 1.) * (phi_i_normal[0] * phi_j_curl[0] + phi_i_normal[1] * phi_j_curl[1] + phi_i_normal[2] * phi_j_curl[2] + phi_j_normal[0] * phi_i_curl[0] + phi_j_normal[1] * phi_i_curl[1] + phi_j_normal[2] * phi_i_curl[2]);
            ret(j, i) = ret(i, j);
            
            //std::cout<<ret(i, j)<<std::endl;
        }
    }
}

void hpGemUIExtentions::anonymous7::faceIntegrand(const FaceT* face, const LinearAlgebra::MiddleSizeVector& normal, const PointElementReferenceT& p, LinearAlgebra::MiddleSizeMatrix& ret)
{
    //cout<<"\nIn the face integrand for the stiffness matrix (IP-only part) for element id: "<<face->getPtrElementLeft()->getID();
    
    //int M = face->getPtrElementLeft()->getNrOfBasisFunctions() + face->getPtrElementRight()->getNrOfBasisFunctions();
    
    int M = face->getNrOfBasisFunctions();
    ret.resize(M, M);
    
    LinearAlgebra::MiddleSizeVector phi_i(3), phi_j(3);
    for (int i = 0; i < M; ++i)
    {
        face->basisFunctionNormal(i, normal, p, phi_i);
        
        for (int j = i; j < M; ++j)
        {
            face->basisFunctionNormal(j, normal, p, phi_j);
            ret(i, j) = 37 * (phi_i[0] * phi_j[0] + phi_i[1] * phi_j[1] + phi_i[2] * phi_j[2]);
            ret(j, i) = ret(i, j);
            //std::cout<<ret(i, j)<<std::endl;
        }
    }
    
}

void hpGemUIExtentions::anonymous8::faceIntegrand(const FaceT* face, const LinearAlgebra::MiddleSizeVector& normal, const PointElementReferenceT& p, LinearAlgebra::MiddleSizeVector& ret)
{
    
    ElementT* left = const_cast<ElementT*>(face->getPtrElementLeft());
    Geometry::PointReference PLeft(3);
    face->mapRefFaceToRefElemL(p, PLeft);
    
    PointPhysicalT PPhys(3);
    left->referenceToPhysical(PLeft, PPhys);
    LinearAlgebra::MiddleSizeVector normedNormal(3);
    
    normedNormal[0] = (normal * (1 / Base::L2Norm(normal)))[0];
    normedNormal[1] = (normal * (1 / Base::L2Norm(normal)))[1];
    normedNormal[2] = (normal * (1 / Base::L2Norm(normal)))[2];
    
    LinearAlgebra::MiddleSizeVector val(3), phi(3), phi_curl(3), dummy(3);
    
    boundaryConditions(PPhys, dummy); //assumes the initial conditions and the boundary conditions match
            
    OuterProduct(normedNormal, dummy, val);
    int n = face->getPtrElementLeft()->getNrOfBasisFunctions();
    ret.resize(n);
    
    for (int i = 0; i < n; ++i)
    {
        face->basisFunctionNormal(i, normal, p, phi);
        face->basisFunctionCurl(i, p, phi_curl);
        
        ret(i) = -(phi_curl[0] * val[0] + phi_curl[1] * val[1] + phi_curl[2] * val[2]) + 37 * (phi[0] * val[0] + phi[1] * val[1] + phi[2] * val[2]);
        //std::cout<<ret(i)<<std::endl;
        
    }
    
}

void hpGemUIExtentions::anonymous9::faceIntegrand(const FaceT* face, const LinearAlgebra::MiddleSizeVector& normal, const PointElementReferenceT& p, LinearAlgebra::MiddleSizeMatrix& ret)
{
    //cout<<"\nIn the face integrand for the stiffness matrix (BR-only part) for element id: "<<face->getPtrElementLeft()->getID();
    ElementT* right;
    ElementT* left = const_cast<ElementT*>(face->getPtrElementLeft());
    ElementInfos* leftInfo = static_cast<ElementInfos*>(left->getUserData());
    ElementInfos* rightInfo;
    
    //PointElementReferenceT pLeft(3),pRight(3);
    //face->mapRefFaceToRefElemL(p,pLeft);
    
    double localepsilon;
    
    if (face->isInternal())
    {
        //cout<<" and element id: "<<face->getPtrElementRight()->getID();
        //face->mapRefFaceToRefElemR(p,pRight);
        right = const_cast<ElementT*>(face->getPtrElementRight());
        rightInfo = static_cast<ElementInfos*>(right->getUserData());
    }
    
    int n = face->getPtrElementLeft()->getNrOfBasisFunctions();
    int M = face->getNrOfBasisFunctions();
    
    //int M = face->getPtrElementLeft()->getNrOfBasisFunctions() + face->getPtrElementRight()->getNrOfBasisFunctions();
    
    ret.resize(M, M);
    
    LinearAlgebra::MiddleSizeVector phi_i(3), phi_j(3);
    
    for (int i = 0; i < M; ++i)
    {
        if (i < n)
        {
            localepsilon = leftInfo->epsilon_;
        }
        else
        {
            localepsilon = rightInfo->epsilon_;
        }
        
        localepsilon = leftInfo->epsilon_;
        face->basisFunction(i, p, phi_i);
        //std::cout<<phi_i<<std::endl;
        for (int j = 0; j < M; ++j)
        {
            face->basisFunctionNormal(j, normal, p, phi_j);
            
            ret(j, i) = (face->isInternal() ? 1 : 2) * (phi_i[0] * phi_j[0] + phi_i[1] * phi_j[1] + phi_i[2] * phi_j[2]) * sqrt(localepsilon);
            //std::cout<<ret<<std::endl;
        }
    }
    
}

void hpGemUIExtentions::anonymous10::faceIntegrand(const FaceT* face, const LinearAlgebra::MiddleSizeVector& normal, const PointElementReferenceT& p, LinearAlgebra::MiddleSizeVector& ret)
{
    
    int n = face->getPtrElementLeft()->getNrOfBasisFunctions();
    
    ElementT* left = const_cast<ElementT*>(face->getPtrElementLeft());
    PointElementReferenceT PLeft(3);
    face->mapRefFaceToRefElemL(p, PLeft);
    
    PointPhysicalT PPhys(3);
    left->referenceToPhysical(PLeft, PPhys);
    LinearAlgebra::MiddleSizeVector normedNormal(3);
    
    normedNormal[0] = (normal * (1 / Base::L2Norm(normal)))[0];
    normedNormal[1] = (normal * (1 / Base::L2Norm(normal)))[1];
    normedNormal[2] = (normal * (1 / Base::L2Norm(normal)))[2];
    
    LinearAlgebra::MiddleSizeVector val(3), dummy(3), phi(3);
    boundaryConditions(PPhys, dummy); //assumes the initial conditions and the boundary conditions match
    OuterProduct(normedNormal, dummy, val);
    //std::cout<<val<<std::endl;
    ret.resize(n);
    for (int i = 0; i < n; ++i)
    {
        face->basisFunction(i, p, phi);
        
        //std::cout<<phi<<std::endl;
        
        ret(i) = 2 * (phi[0] * val[0] + phi[1] * val[1] + phi[2] * val[2]);
    }
}

void hpGemUIExtentions::anonymous11::faceIntegrand(const FaceT* face, const LinearAlgebra::MiddleSizeVector& normal, const PointElementReferenceT& p, LinearAlgebra::MiddleSizeVector& ret)
{
    
    int n = face->getPtrElementLeft()->getNrOfBasisFunctions();
    
    ElementT* left = const_cast<ElementT*>(face->getPtrElementLeft());
    PointElementReferenceT PLeft(3);
    face->mapRefFaceToRefElemL(p, PLeft);
    
    PointPhysicalT PPhys(3);
    left->referenceToPhysical(PLeft, PPhys);
    LinearAlgebra::MiddleSizeVector normedNormal(3);
    
    normedNormal[0] = (normal * (1 / Base::L2Norm(normal)))[0];
    normedNormal[1] = (normal * (1 / Base::L2Norm(normal)))[1];
    normedNormal[2] = (normal * (1 / Base::L2Norm(normal)))[2];
    
    LinearAlgebra::MiddleSizeVector val(3), dummy(3), phi_curl(3);
    boundaryConditions(PPhys, dummy); //assumes the initial conditions and the boundary conditions match
    OuterProduct(normedNormal, dummy, val);
    //std::cout<<val<<std::endl;
    ret.resize(n);
    for (int i = 0; i < n; ++i)
    {
        face->basisFunctionCurl(i, p, phi_curl);
        
        //std::cout<<phi<<std::endl;
        
        ret(i) = -(phi_curl[0] * val[0] + phi_curl[1] * val[1] + phi_curl[2] * val[2]);
    }
    
}

