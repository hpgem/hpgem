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
#include "Base/MpiContainer.h"
#include "Base/ConfigurationData.h"
#include "Geometry/PointPhysical.h"
#include "Integration/ElementIntegral.h"
#include "Integration/FaceIntegral.h"
#include "Base/HCurlConformingTransformation.h"
#include "Output/TecplotDiscontinuousSolutionWriter.h"
#include <valarray>
#include <cmath>
#include <complex>
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
    return static_cast<const MaxwellData*>(globalData_);
}

double hpGemUIExtentions::sourceTermTime(const double t)
{
    return 1.0; // for comparison with Freekjan's report. The 1.0 is not physically meaningful. The source term is actually 0.0, but we want to avoid division by zero in timedependent code.
    
    // for comparison with time-integration paper by Domokos Sarmany
    //return -cos(t)-1.0/4.0*cos(t/2.0)-1.0/9.0*cos(t/3.0)+2*M_PI*M_PI*(cos(t)+cos(.5*t)+cos(1.0/3.0*t));
}

void hpGemUIExtentions::exactSolution(const Geometry::PointPhysical<DIM>& p, const double t, LinearAlgebra::SmallVector<DIM>& ret)
{
      //ret[0]=sin(M_PI*2*p[1])*sin(M_PI*2*p[2]);
      //ret[1]=sin(M_PI*2*p[2])*sin(M_PI*2*p[0]);
      //ret[2]=sin(M_PI*2*p[0])*sin(M_PI*2*p[1]);
      //ret*=cos(sqrt(2)*2*M_PI*t);
    
    // for comparison with Freekjan's report.
    ret[0]=sin(M_PI*p[1])*sin(M_PI*p[2]);
    ret[1]=sin(M_PI*p[2])*sin(M_PI*p[0]);
    ret[2]=sin(M_PI*p[0])*sin(M_PI*p[1]);
    ret*=cos(sqrt(2)*M_PI*t);
    //ret[0] = p[2];
    //ret[1] = p[0];
    //ret[2] = p[1];
    
    //     ret[0]=p[0]*(1-p[0]);
    //     ret[1]=0;
    // 	   ret[2]=0;
    
    // for comparison with the time-integration paper by Domokos Sarmany
    //ret[0]=sin(M_PI*p[1])*sin(M_PI*p[2]);
    //ret[1]=sin(M_PI*p[2])*sin(M_PI*p[0]);
    //ret[2]=sin(M_PI*p[0])*sin(M_PI*p[1]);
    //ret*=(cos(t)+cos(.5*t)+cos(1.0/3.0*t));
}

void hpGemUIExtentions::exactSolutionCurl(const Geometry::PointPhysical<DIM>& p, const double t, LinearAlgebra::SmallVector<DIM>& ret)
{
      //ret[0]=sin(M_PI*2*p[0])*(cos(M_PI*2*p[1])-cos(M_PI*2*p[2]));
      //ret[1]=sin(M_PI*2*p[1])*(cos(M_PI*2*p[2])-cos(M_PI*2*p[0]));
      //ret[2]=sin(M_PI*2*p[2])*(cos(M_PI*2*p[0])-cos(M_PI*2*p[1]));
      //ret*=cos(sqrt(2)*2*M_PI*t)*2*M_PI;
    
    // for comparison with Freekjan's report.
    ret[0] = sin(M_PI * p[0]) * (cos(M_PI * p[1]) - cos(M_PI * p[2]));
    ret[1] = sin(M_PI * p[1]) * (cos(M_PI * p[2]) - cos(M_PI * p[0]));
    ret[2] = sin(M_PI * p[2]) * (cos(M_PI * p[0]) - cos(M_PI * p[1]));
    ret *= cos(sqrt(2) * M_PI * t) * M_PI;
    
    //ret[0] = 1.0;
    //ret[1] = 1.0;
    //ret[2] = 1.0;
    
    //          ret[0]=0;ret[1]=0;ret[2]=0;
    
    // for comparison with the time-integration paper by Domokos Sarmany
    //ret[0] = sin(M_PI * p[0]) * (cos(M_PI * p[1]) - cos(M_PI * p[2]));
    //ret[1] = sin(M_PI * p[1]) * (cos(M_PI * p[2]) - cos(M_PI * p[0]));
    //ret[2] = sin(M_PI * p[2]) * (cos(M_PI * p[0]) - cos(M_PI * p[1]));
    //ret *= (cos(t)+cos(.5*t)+cos(1.0/3.0*t)) * M_PI;
}


 void hpGemUIExtentions::writeToTecplotFile(const Base::Element* element, const PointReferenceT& p, std::ostream& output)
//Needs to be changed in arguments, use PhysicalElement instead of Reference element in the argument list
{
   
   // const Base::Element* element = el.getElement();
    //const Geometry::PointReference<DIM>& p = el.getPointReference();

    
    LinearAlgebra::MiddleSizeVector data;
    data = const_cast<ElementT*>(element)->getTimeIntegrationVector(timelevel_);
    //std::cout<<data.size()<<std::endl;
    LinearAlgebra::SmallVector<DIM> results1, results, curls1, curls;
    
    for (int i = 0; i < element->getNrOfBasisFunctions(); ++i)
    {
        element->basisFunction(i, p, results1); //Needs to be corrected
        curls1 = element->basisFunctionCurl(i, p); //Needs to be corrected
        results1 = results1 * std::real(data[i]);
        results += results1;
        curls1 = curls1 * std::real(data[i]);
        curls += curls1;
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

void hpGemUIExtentions::elementIntegrand(Base::PhysicalElement<DIM>& el, errorData &ret)
{
    
    const Base::Element* element = el.getElement();
    const Geometry::PointReference<DIM>& p = el.getPointReference();
    PointPhysicalT pPhys;
    pPhys = element->referenceToPhysical(p);
    
    
    LinearAlgebra::SmallVector<DIM> phi, phiCurl, error, errorCurl;
    
    exactSolution(pPhys, measureTimes_[timelevel_], error);
    exactSolutionCurl(pPhys, measureTimes_[timelevel_], errorCurl);
    LinearAlgebra::MiddleSizeVector data;
    data = element->getTimeIntegrationVector(timelevel_);
    for (int i = 0; i < element->getNrOfBasisFunctions(); ++i)
    {
        el.basisFunction(i, phi);
        phiCurl = el.basisFunctionCurl(i);
        error -= (std::real(data[i]) * phi);
        errorCurl -= (std::real(data[i]) * phiCurl);
    }
    ret[0] = Base::L2Norm(error) * Base::L2Norm(error);
    ret[1] = Base::L2Norm(errorCurl) * Base::L2Norm(errorCurl);
}

void hpGemUIExtentions::faceIntegrand(Base::PhysicalFace<DIM>& fa, errorData &ret)
{
    const Base::Face* face = fa.getFace();
    LinearAlgebra::SmallVector<DIM> normal = fa.getNormalVector();
    const Geometry::PointReference<2>& p = fa.getPointReference();

    ElementT* element = const_cast<ElementT*>(face->getPtrElementLeft());
    PointPhysicalT PPhys;
    const PointElementReferenceT& pElement = face->mapRefFaceToRefElemL(p);
    
    PPhys = element->referenceToPhysical(pElement);
    LinearAlgebra::SmallVector<DIM> error, phiNormal1, phiNormal2, dummy, dummy2, dummy3, normedNormal;
    normedNormal[0] = (normal * (1 / Base::L2Norm(normal)))[0];
    normedNormal[1] = (normal * (1 / Base::L2Norm(normal)))[1];
    normedNormal[2] = (normal * (1 / Base::L2Norm(normal)))[2];
    
    exactSolution(PPhys, measureTimes_[timelevel_], dummy);
    OuterProduct(normedNormal, dummy, error);
    int n = face->getPtrElementLeft()->getNrOfBasisFunctions();
    LinearAlgebra::MiddleSizeVector data = element->getTimeIntegrationVector(timelevel_); //Issue regarding parallelisation is in this line....it goes out of bound for memory
    
    for (int i = 0; i < n; ++i)
    {
        fa.basisFunctionUnitNormal(i, phiNormal1);
        error -= std::real(data[i]) * phiNormal1;
    }
    if (face->isInternal())
    {
        element = const_cast<ElementT*>(face->getPtrElementRight());
        const PointElementReferenceT& pElement = face->mapRefFaceToRefElemR(p);
        PPhys = element->referenceToPhysical(pElement);
        exactSolution(PPhys, measureTimes_[timelevel_], dummy2);
        OuterProduct(normedNormal, dummy2, dummy3);
        error -= dummy3;
        
        data = element->getTimeIntegrationVector(timelevel_);
        int M = face->getPtrElementLeft()->getNrOfBasisFunctions() + face->getPtrElementRight()->getNrOfBasisFunctions();
        
        for (int i = n; i < M; ++i)
        {
            fa.basisFunctionUnitNormal(i, phiNormal2);
            error -= (std::real(data[i - n]) * phiNormal2);
        }
    }
    ret[0] = Base::L2Norm(error) * Base::L2Norm(error);
    //To remove double contribution of flux computed on the boudary faces by different processors
    if(face->getFaceType() == Geometry::FaceType::SUBDOMAIN_BOUNDARY ||face->getFaceType() == Geometry::FaceType::PERIODIC_SUBDOMAIN_BC)
    {
        ret[0] /= 2.;
    }

}

void hpGemUIExtentions::computeErrors(Base::MeshManipulator<DIM>& Mesh, int timelevel, double& L2Norm, double& InfNorm, double& HCurlNorm, double& DGNorm)
{
    //this should probalby be toggled off when the exact solution is not known
    //alternatively something based on richardson extrapolation can be done
    
    int TotalAmountOfProcessors, localProcessorNumber, numberOfElements;
    MPI_Comm_size(PETSC_COMM_WORLD, &TotalAmountOfProcessors);
    MPI_Comm_rank(PETSC_COMM_WORLD, &localProcessorNumber);
    numberOfElements = meshes_[0]->getElementsList().size();
    
    this->timelevel_ = timelevel;
    errorData integrationResults;
    
    Integration::ElementIntegral<DIM> elIntegral(false);
    elIntegral.setTransformation(std::shared_ptr<Base::CoordinateTransformation<DIM> > (new Base::HCurlConformingTransformation<DIM>()));
    
    L2Norm = 0;
    HCurlNorm = 0;
    
    for (ElementIterator it = Mesh.elementColBegin(); it != Mesh.elementColEnd(); ++it)
    {
        
            integrationResults = elIntegral.integrate<errorData>((*it), this);
            L2Norm += std::abs(integrationResults[0]);
            HCurlNorm += std::abs(integrationResults[0]) + std::abs(integrationResults[1]);
            integrationResults.reset();
    }
    DGNorm = HCurlNorm;
            
    Integration::FaceIntegral<DIM> faIntegral(false);
    faIntegral.setTransformation(std::shared_ptr<Base::CoordinateTransformation<DIM> >(new Base::HCurlConformingTransformation<DIM>()));
    //Error in parallelisation occurs because of the following loop
    
    for (FaceIterator it = Mesh.faceColBegin(); it != Mesh.faceColEnd(); ++it)
    {
            integrationResults = faIntegral.integrate<errorData>(*it, this); //Error in parallelisation occurs because of this integration routine to calculate errors for the faces
            DGNorm += std::abs(integrationResults[0]);
            integrationResults.reset();
    }
    
}

void hpGemUIExtentions::makeOutput(char* filename)
{
    
    //set up containers for the errors
    Vec L2Norm, HCurlNorm, DGNorm;
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
    
    Output::TecplotDiscontinuousSolutionWriter<DIM> tecplotWriter(fileWriter, "The electric field", "012", "E0,E1,E2,H0,H1,H2");
    
    std::stringstream string;
    char parsed[20] = "";
    string.readsome(parsed, 20); //clean storage
    string << "t=" << measureTimes_[0];
    
    int read = string.readsome(parsed, 20);
    string.readsome(&parsed[read], 20 - read);
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
    VecView(HCurlNorm, 0);
    VecView(DGNorm, 0);
    ierr_ = VecDestroy(&HCurlNorm);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecDestroy(&DGNorm);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecDestroy(&L2Norm);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    
}

hpGemUIExtentions::hpGemUIExtentions(MaxwellData* globalConfig, Base::ConfigurationData* elementConfig, MatrixAssembly* fluxType)
        : HpgemAPIBase(globalConfig, elementConfig), assembler(fluxType)
{
    ierr_ = PetscLogBegin();
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = MatCreate(PETSC_COMM_WORLD, &M_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = MatCreate(PETSC_COMM_WORLD, &S_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecCreate(PETSC_COMM_WORLD, &x_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecCreate(PETSC_COMM_WORLD, &derivative_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecCreate(PETSC_COMM_WORLD, &RHS_);
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
    ierr_ = VecDestroy(&RHS_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    
    ierr_ = VecDestroy(&derivative_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    
    ierr_ = VecDestroy(&x_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    
    ierr_ = MatDestroy(&S_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    
    ierr_ = MatDestroy(&M_);
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
    
}

std::size_t hpGemUIExtentions::addMesh(Base::MeshManipulator<DIM>* mesh)
{
    meshes_.push_back(mesh);
    return meshes_.size() - 1;
}

void hpGemUIExtentions::makeShiftMatrix(LinearAlgebra::SmallVector<DIM>& direction, Vec& waveVecMatrix)
{
    for (ElementIterator it = meshes_[0]->elementColBegin(); it != meshes_[0]->elementColEnd(); ++it)
    {
        for (int j = 0; j < (*it)->getNrOfBasisFunctions(); ++j)
        {
            PointPhysicalT centerPhys;
            //PointReferenceT center(3);
            const PointReferenceT& center = (*it)->getReferenceGeometry()->getCenter();
            ;
            centerPhys = (*it)->referenceToPhysical(center);
            //this extra accuracy is probably irrelevant and a lot of extra ugly to get it working
            
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
        logger.assert_always((*it)->isInternal(), "Internal face boundary");
        const PointFaceReferenceT& p = (*it)->getReferenceGeometry()->getCenter();
        PointPhysicalT pLeftPhys, pRightPhys;
        const PointElementReferenceT& pLeft = (*it)->mapRefFaceToRefElemL(p);
        const PointElementReferenceT& pRight = (*it)->mapRefFaceToRefElemR(p);
        pLeftPhys = (*it)->getPtrElementLeft()->referenceToPhysical(pLeft);
        pRightPhys = (*it)->getPtrElementRight()->referenceToPhysical(pRight);
        //if the left coordinate is not close to the right coordinate it is a boundary face
        if (Base::L2Norm(pLeftPhys - pRightPhys) > 1e-3)
        { //pretty lousy tolerance, but this norm should be either 1 or 0
            if ((pLeftPhys[0] - pRightPhys[0]) * (pLeftPhys[0] - pRightPhys[0]) > 1e-3)
            {
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
                zRow.resize(nz + 2);
                zCol.resize(nz + 2);
                logger.assert_always((pLeftPhys[2] - pRightPhys[2]) * (pLeftPhys[2] - pRightPhys[2]) > (1.e-3), "Boundary Block in z direction");
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
/*
void hpGemUIExtentions::LDOSIntegrand(Base::PhysicalElement<DIM>& element, double &ret)
{ //currently LDOS is computed by evaluation at a point so integrand is a bit of a misnomer
    //ElementInfos* info = static_cast<ElementInfos*>(const_cast<ElementT*>(element)->getUserData());
    LinearAlgebra::SmallVector<DIM> Phi, PhiRealI, PhiRealJ, PhiImagI, PhiImagJ;
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
            //PointElementReferenceT p(3);
            const PointElementReferenceT& p = (*it)->getReferenceGeometry()->getCenter();
            LDOSIntegrand(*it, partialResult);
            result[index] = partialResult;
            ++index;
        }
    }
    ierr_ = VecDestroy(&scaledVec);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
}
*/

void hpGemUIExtentions::GetCoeffCO4(LinearAlgebra::SmallVector<6>& alpha,
		 LinearAlgebra::SmallVector<6>& beta,
		 LinearAlgebra::SmallVector<6>& alpha_sum,
		 LinearAlgebra::SmallVector<6>& beta_sum,
		 LinearAlgebra::SmallVector<6>& scale0,
		 LinearAlgebra::SmallVector<6>& scale1,
		 const double& tau)
{
  alpha[0] = 0.0;
  beta[0]  = 0.0;
  scale0[0] = 0.0;
  scale1[0] = 0.0;

  alpha[1] = (146 + 5*std::sqrt(19.0))/540;
  beta[5]  = alpha[1];

  alpha[2] = (-2 + 10*std::sqrt(19.0)) / 135;
  beta[4]  = alpha[2];

  alpha[3] = 1.0/5.0;
  beta[3]  = alpha[3];

  alpha[4] = (-23 - 20*std::sqrt(19.0)) / 270;
  beta[2]  = alpha[4];

  alpha[5] = (14 - std::sqrt(19.0)) / 108;
  beta[1]  = alpha[5];

  for(unsigned int i = 0; i < 6; ++i)
    {
      alpha_sum[i] = 0.0;
      beta_sum[i] = 0.0;
      for(unsigned int j = 0; j <= i; ++j)
	{
	  alpha_sum[i] += alpha[j];
	  beta_sum[i]  += beta[j];
	}
      scale0[i] = 1;
      scale1[i] = 1;
    }

//   for(unsigned int i = 0; i < 6; ++i)
//     cout << scale0[i] << "\t" << scale1[i] << endl;

//   for(unsigned int i = 0; i < 6; ++i)
//     cout << alpha_sum[i] << "\t" << beta_sum[i] << endl;

  return;
}

void hpGemUIExtentions::solveTimeDependent(bool useCO2, bool useCO4)
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
    
    Vec dummy, dummy2, dummy3; 
    ierr_ = VecDuplicate(derivative_, &dummy);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    if (useCO4 == true)
    {
        ierr_ = VecDuplicate(derivative_, &dummy3);
        CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    }
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
    
    double tau;
    if (useCO2 == true)
    {
        tau = 0.05 / ((2 * actualdata->PolynomialOrder_ + 1) * actualdata->NumberOfIntervals_);
    }
        
    if (useCO4 == true)
    {
        tau = 0.3/( (2*actualdata->PolynomialOrder_+1)*actualdata->NumberOfIntervals_ );
    }
    
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
    
    double scale0, scale1;
    if (useCO2 == true)
    {
        scale0 = (1 - actualdata->Sigma_ * tau / 2);
        scale1 = 1 / (1 + actualdata->Sigma_ * tau / 2);
    }
    
    LinearAlgebra::SmallVector<6> alpha, beta, alpha_sum, beta_sum, scale0vector, scale1vector;
    double sourceTime;
    if (useCO4 == true)
    {
        GetCoeffCO4(alpha, beta, alpha_sum, beta_sum, scale0vector, scale1vector, tau);
    }
    int measureAmount = 0;
    
    std::cout << tau << " " << Nsteps << std::endl;

    for (int i = 0; i < Nsteps; ++i)
    {
        if (i % measureStep == 0)
        {
            x_.writeTimeIntegrationVector(measureAmount);
            
            measureTimes_[measureAmount] = t;
            measureAmount++;
        }
        if (useCO2 == true)
        {
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
       
        if (useCO4 == true)
        {
            for(unsigned int k = 1; k < 6; ++k)
            {
                ierr_ = VecAXPY(x_, tau*(alpha[k-1] + beta[k]), derivative_);
                CHKERRABORT(PETSC_COMM_WORLD, ierr_);

                ierr_ = MatMult(S_, x_, dummy);
                CHKERRABORT(PETSC_COMM_WORLD, ierr_);

                sourceTime = t + (alpha_sum[k-1] + beta_sum[k-1])*tau;
                eta = sourceTermTime(sourceTime);
                ierr_ = VecCopy(RHS_, dummy2);
                CHKERRABORT(PETSC_COMM_WORLD, ierr_);
                ierr_ = VecScale(dummy2, beta[k]*tau*eta);
                CHKERRABORT(PETSC_COMM_WORLD, ierr_);
                VecAYPX(dummy, -(beta[k]+alpha[k])*tau, dummy2);
                CHKERRABORT(PETSC_COMM_WORLD, ierr_);

                sourceTime = t + (alpha_sum[k] + beta_sum[k])*tau;
                eta = sourceTermTime(sourceTime);
                ierr_ = VecCopy(RHS_, dummy2);
                CHKERRABORT(PETSC_COMM_WORLD, ierr_);
                ierr_ = VecScale(dummy2, alpha[k]*tau*eta);
                CHKERRABORT(PETSC_COMM_WORLD, ierr_);
                ierr_ = VecAXPY(dummy, 1, dummy2);
                CHKERRABORT(PETSC_COMM_WORLD, ierr_);

                ierr_ = VecScale(dummy, scale1vector[k]);
                CHKERRABORT(PETSC_COMM_WORLD, ierr_);
                ierr_ = MatMult(M_, dummy, dummy3);
                CHKERRABORT(PETSC_COMM_WORLD, ierr_);
                VecAYPX(derivative_, scale0vector[k]*scale1vector[k], dummy3);
                CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            }
            ierr_ = VecAXPY(x_, tau*alpha[5], derivative_);
            CHKERRABORT(PETSC_COMM_WORLD, ierr_);
            t += tau;
        }
    //getArrayRead only gets local values, so make the entire std::vector local
        
    x_.writeTimeIntegrationVector(measureAmount);
    measureTimes_[measureAmount] = t;
    }
}

void hpGemUIExtentions::solveHarmonic()
{
    std::cout << "finding a time-harmonic solution" << std::endl;
    const MaxwellData* actualdata = getData();
    MHasToBeInverted_ = false;
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

    //everything that is set in the code, but before this line is overridden by command-line options
    ierr_ = KSPSetFromOptions(solver_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    
    ierr_ = KSPSetOperators(solver_, S_, S_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = KSPSetUp(solver_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = KSPSolve(solver_, RHS_, x_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    
    x_.writeTimeIntegrationVector(0); //DOUBTFUL
    measureTimes_[0]=0;
    
    synchronize(0);
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
    int degreesOfFreedomPerElement = configData_->numberOfBasisFunctions_ * configData_->numberOfUnknowns_;
    std::valarray<PetscScalar> blockvalues(degreesOfFreedomPerElement * degreesOfFreedomPerElement);
    int measureAmount = 0;
    
    std::vector<IS> xboundaryRow, xboundaryCol, yboundaryRow, yboundaryCol, zboundaryRow, zboundaryCol;
    //findBoundaryBlocks(xboundaryRow, xboundaryCol, yboundaryRow, yboundaryCol, zboundaryRow, zboundaryCol);
    
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
   
    Mat product;
    ierr_ = MatMatMult(M_, S_, MAT_INITIAL_MATRIX, 1.0, &product);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    
    ierr_ = EPSSetOperators(eigenSolver_, product, NULL);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = EPSSetUp(eigenSolver_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = EPSSetDimensions(eigenSolver_, 24, PETSC_DECIDE, PETSC_DECIDE);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    
    //everything that is set in the code, but before this line is overridden by command-line options
    ierr_ = EPSSetFromOptions(eigenSolver_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    //everything that is set in the code, but after this line overrides the comand-line options
    
    ierr_ = EPSSolve(eigenSolver_);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    
    PetscScalar neededOnlyForRealPetsc, eigenvalue;
    Vec *eigenvalues, example, *eigenVectors;
    eigenvalues = new Vec[24];
    eigenVectors = new Vec[40]; //a few extra in case SLEPc finds more than the requested amount of eigenvalues
    ierr_ = VecCreate(PETSC_COMM_WORLD, &example);
    CHKERRABORT(PETSC_COMM_WORLD, ierr_);
    ierr_ = VecSetSizes(example, PETSC_DECIDE, 1);// SH put back to 61 when turning k back on
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
    LinearAlgebra::SmallVector<DIM> k;
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
    
    /*for (int i = 1; i < 61; ++i)
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
            x_.writeTimeIntegrationVector(measureAmount);
            measureTimes_[measureAmount] = i / 20;
            measureAmount++;
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
    }*/
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
    
    x_.writeTimeIntegrationVector(measureAmount);
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
/*
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
    LinearAlgebra::NumericalVector k(3);
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
                LinearAlgebra::NumericalVector functionvalue;
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
    
    LinearAlgebra::NumericalVector result(meshes_[0]->getNumberOfElements());
    for (int i = 0; i < 1001; ++i)
    {
        brillouinZone.getIntegral(double(i) / 100., result);
        std::cout << result[0] << " " << result[6] << " " << result[10] << std::endl;
        result[0] = 0;
        //result[6]=0;
        //result[10]=0;
    }
}
*/
//this is where you specify an initial condition
void hpGemUIExtentions::initialConditions(const PointPhysicalT& p, LinearAlgebra::SmallVector<DIM>& ret)
{
    hpGemUIExtentions::exactSolution(p, 0, ret);
    }

void hpGemUIExtentions::initialConditionsDeriv(const PointPhysicalT& p, LinearAlgebra::SmallVector<DIM>& ret)
{
    ret[0] = 0;
    ret[1] = 0;
    ret[2] = 0;
}

/**
 * this is where you specify the spatial part of the source Term
 * assumes that the source term can be split is a spatial part and a time part
 */

void hpGemUIExtentions::sourceTerm(const PointPhysicalT& p, LinearAlgebra::SmallVector<DIM>& ret)
{
    hpGemUIExtentions::exactSolution(p, 0, ret);
    // 	ret*=-1;
    //ret*=M_PI*M_PI*8-1;
    ret *= M_PI * M_PI * 2 - 1;
    //ret *= -1; 
    // for comparison with the time-dependent code in Freekjan's report.
    //ret[0] = 0.0;
    //ret[1] = 0.0;
    //ret[2] = 0.0;
    
    // for comparison with the time-integration paper by Domokos Sarmany
    //ret[0] = hpGemUIExtentions::sourceTermTime(0)*sin(M_PI * p[1]) * sin(M_PI * p[2]);
    //ret[1] = hpGemUIExtentions::sourceTermTime(0)*sin(M_PI * p[2]) * sin(M_PI * p[0]);
    //ret[2] = hpGemUIExtentions::sourceTermTime(0)*sin(M_PI * p[0]) * sin(M_PI * p[1]);
}

void hpGemUIExtentions::initialExactSolution(const PointPhysicalT& p, LinearAlgebra::SmallVector<DIM>& ret)
{
    //ret[0]=sin(M_PI*2*p[1])*sin(M_PI*2*p[2]);
    //ret[1]=sin(M_PI*2*p[2])*sin(M_PI*2*p[0]);
    //ret[2]=sin(M_PI*2*p[0])*sin(M_PI*2*p[1]);
    
    ret[0] = sin(M_PI * p[1]) * sin(M_PI * p[2]);
    ret[1] = sin(M_PI * p[2]) * sin(M_PI * p[0]);
    ret[2] = sin(M_PI * p[0]) * sin(M_PI * p[1]);
    
    //ret[0] = p[2];
    //ret[1] = p[0];
    //ret[2] = p[1];
    
    //            ret[0]=p[0]*(1-p[0]);
    //            ret[1]=0;
    // 	  		  ret[2]=0;
    
    // for comparison with the time-integration paper by Domokos Sarmany
    //ret[0] = 3.0*sin(M_PI * p[1]) * sin(M_PI * p[2]);
    //ret[1] = 3.0*sin(M_PI * p[2]) * sin(M_PI * p[0]);
    //ret[2] = 3.0*sin(M_PI * p[0]) * sin(M_PI * p[1]);
}

void hpGemUIExtentions::boundaryConditions(const PointPhysicalT &p, LinearAlgebra::SmallVector<DIM>& ret)
{
    initialExactSolution(p, ret);
}

//elementMassIntegrand

void hpGemUIExtentions::anonymous1::elementIntegrand(Base::PhysicalElement<DIM>& el, LinearAlgebra::MiddleSizeMatrix& ret)
{
    const Base::Element* element = el.getElement();
    ret.resize(element->getNrOfBasisFunctions(), element->getNrOfBasisFunctions());
    ElementInfos* info = static_cast<ElementInfos*>(const_cast<ElementT*>(element)->getUserData());
    LinearAlgebra::SmallVector<DIM> phi_i, phi_j;
    for (int i = 0; i < element->getNrOfBasisFunctions(); ++i)
    {
        el.basisFunction(i, phi_i);
        for (int j = 0; j < element->getNrOfBasisFunctions(); ++j)
        {
            el.basisFunction(j, phi_j);
            ret(i, j) = (phi_i[0] * phi_j[0] + phi_i[1] * phi_j[1] + phi_i[2] * phi_j[2]) * info->epsilon_;
            ret(j, i) = ret(i, j);
        }
    }
}

//elementStiffnessIntegrand

void hpGemUIExtentions::anonymous2::elementIntegrand(Base::PhysicalElement<DIM>& el, LinearAlgebra::MiddleSizeMatrix& ret)
{
    const Base::Element* element = el.getElement();
    ret.resize(element->getNrOfBasisFunctions(), element->getNrOfBasisFunctions());
    LinearAlgebra::SmallVector<DIM> phi_i, phi_j;
    for (int i = 0; i < element->getNrOfBasisFunctions(); ++i)
    {
        phi_i = el.basisFunctionCurl(i);
        for (int j = i; j < element->getNrOfBasisFunctions(); ++j)
        {
            phi_j = el.basisFunctionCurl(j);
            ret(i, j) = phi_i[0] * phi_j[0] + phi_i[1] * phi_j[1] + phi_i[2] * phi_j[2];
            ret(j, i) = ret(i, j);
        }
    }
    
}

//elementSpaceIntegrand

void hpGemUIExtentions::anonymous3::elementIntegrand(Base::PhysicalElement<DIM>& el, LinearAlgebra::MiddleSizeVector& ret)
{
    
    const Base::Element* element = el.getElement();
    const Geometry::PointReference<DIM>& p = el.getPointReference();
    
    ret.resize(element->getNrOfBasisFunctions());
    PointPhysicalT pPhys;
    pPhys = element->referenceToPhysical(p);
    LinearAlgebra::SmallVector<DIM> val, phi;
    sourceTerm(pPhys, val);
    for (int i = 0; i < element->getNrOfBasisFunctions(); ++i)
    {
        el.basisFunction(i, phi);
        ret(i) = phi[0] * val[0] + phi[1] * val[1] + phi[2] * val[2];
    }
}

//initialConditionsIntegrand

void hpGemUIExtentions::anonymous4::elementIntegrand(Base::PhysicalElement<DIM>& el, LinearAlgebra::MiddleSizeVector& ret)
{
    const Base::Element* element = el.getElement();
    const Geometry::PointReference<DIM>& p = el.getPointReference();
    
    ret.resize(element->getNrOfBasisFunctions());
    PointPhysicalT pPhys;
    pPhys = element->referenceToPhysical(p);
    LinearAlgebra::SmallVector<DIM> val, phi;
    initialConditions(pPhys, val);

    for (int i = 0; i < element->getNrOfBasisFunctions(); ++i)
    {
        el.basisFunction(i, phi);
        ret(i) = phi[0] * val[0] + phi[1] * val[1] + phi[2] * val[2];
    }
}

//initialConditionsDerivIntegrand

void hpGemUIExtentions::anonymous5::elementIntegrand(Base::PhysicalElement<DIM>& el, LinearAlgebra::MiddleSizeVector& ret)
{
    const Base::Element* element = el.getElement();
    const Geometry::PointReference<DIM>& p = el.getPointReference();
    ret.resize(element->getNrOfBasisFunctions());
    PointPhysicalT pPhys;
    pPhys = element->referenceToPhysical(p);
    LinearAlgebra::SmallVector<DIM> val, phi;
    initialConditionsDeriv(pPhys, val);
    for (int i = 0; i < element->getNrOfBasisFunctions(); ++i)
    {
        el.basisFunction(i, phi);
        ret(i) = phi[0] * val[0] + phi[1] * val[1] + phi[2] * val[2];
    }
}

// faceStiffnessIntegrand

void hpGemUIExtentions::anonymous6::faceIntegrand(Base::PhysicalFace<DIM>& fa, LinearAlgebra::MiddleSizeMatrix& ret)
{
    const Base::Face* face = fa.getFace();
    LinearAlgebra::SmallVector<DIM> normal = fa.getUnitNormalVector();
    int M = face->getPtrElementLeft()->getNrOfBasisFunctions();
    LinearAlgebra::SmallVector<DIM> normedNormal;
    
    normedNormal[0] = (normal * (1 / Base::L2Norm(normal)))[0];
    normedNormal[1] = (normal * (1 / Base::L2Norm(normal)))[1];
    normedNormal[2] = (normal * (1 / Base::L2Norm(normal)))[2];

    if(face->isInternal())
    {
        M = face->getPtrElementLeft()->getNrOfBasisFunctions() + face->getPtrElementRight()->getNrOfBasisFunctions();
    }
    ret.resize(M, M);
    LinearAlgebra::SmallVector<DIM> phi_i_normal, phi_j_normal, phi_i_curl, phi_j_curl;
    for (int i = 0; i < M; ++i)
    {
        phi_i_curl = fa.basisFunctionCurl(i);
        fa.basisFunctionUnitNormal(i, phi_i_normal);
        //std::cout<<fa.getFace()->getPtrElementLeft()->getID()<<" "<<fa.getFace()->localFaceNumberLeft()<<" "<<i<<" "<<phi_i_normal<<std::endl;
        
        for (int j = i; j < M; ++j)
        {
            phi_j_curl = fa.basisFunctionCurl(j);
            fa.basisFunctionUnitNormal(j, phi_j_normal);
            
            ret(i, j) = -(face->isInternal() ? 0.5 : 1.) * (phi_i_normal[0] * phi_j_curl[0] + phi_i_normal[1] * phi_j_curl[1] + phi_i_normal[2] * phi_j_curl[2] + phi_j_normal[0] * phi_i_curl[0] + phi_j_normal[1] * phi_i_curl[1] + phi_j_normal[2] * phi_i_curl[2]);
            ret(j, i) = ret(i, j);
        }
    }
}

// faceStiffnessIntegrandIP

void hpGemUIExtentions::anonymous7::faceIntegrand(Base::PhysicalFace<DIM>& fa, LinearAlgebra::MiddleSizeMatrix& ret)
{
    const Base::Face* face = fa.getFace();
    LinearAlgebra::SmallVector<DIM> normal = fa.getNormalVector();
    LinearAlgebra::SmallVector<DIM> normedNormal;
    
    normedNormal[0] = (normal * (1 / Base::L2Norm(normal)))[0];
    normedNormal[1] = (normal * (1 / Base::L2Norm(normal)))[1];
    normedNormal[2] = (normal * (1 / Base::L2Norm(normal)))[2];

    
    int M = face->getPtrElementLeft()->getNrOfBasisFunctions();
    if(face->isInternal())
    {
         M = face->getPtrElementLeft()->getNrOfBasisFunctions() + face->getPtrElementRight()->getNrOfBasisFunctions();
    }
   
    ret.resize(M, M);
    LinearAlgebra::SmallVector<DIM> phi_i0, phi_j0, phi_i, phi_j;
    for (int i = 0; i < M; ++i)
    {
        fa.basisFunctionUnitNormal(i, phi_i);
        for (int j = i; j < M; ++j)
        {
            fa.basisFunctionUnitNormal(j, phi_j);
            
            ret(i, j) = MaxwellData::StabCoeff_ * (phi_i[0] * phi_j[0] + phi_i[1] * phi_j[1] + phi_i[2] * phi_j[2]);
            ret(j, i) = ret(i, j);
        }
    }
    
}
// faceSpaceIntegrandIP

void hpGemUIExtentions::anonymous8::faceIntegrand(Base::PhysicalFace<DIM>& fa, LinearAlgebra::MiddleSizeVector& ret)
{
    const Base::Face* face = fa.getFace();
    LinearAlgebra::SmallVector<DIM> normal = fa.getNormalVector();
    const Geometry::PointReference<DIM - 1>& p = fa.getPointReference();
    
    LinearAlgebra::SmallVector<DIM> normedNormal;
    
    normedNormal[0] = (normal * (1 / Base::L2Norm(normal)))[0];
    normedNormal[1] = (normal * (1 / Base::L2Norm(normal)))[1];
    normedNormal[2] = (normal * (1 / Base::L2Norm(normal)))[2];

    if(face->isInternal())
       {
           int M = face->getPtrElementLeft()->getNrOfBasisFunctions() + face->getPtrElementRight()->getNrOfBasisFunctions();
           ret.resize(M);
           for(int i = 0 ; i < M; ++i)
               ret(i) = 0;
       }
    else
       {
           ElementT* left = const_cast<ElementT*>(face->getPtrElementLeft());
           const PointFaceReferenceT& PLeft = face->mapRefFaceToRefElemL(p);
    
           PointPhysicalT PPhys;
           PPhys = left->referenceToPhysical(PLeft);
           LinearAlgebra::SmallVector<DIM> normedNormal;
    
           normedNormal[0] = (normal * (1 / Base::L2Norm(normal)))[0];
           normedNormal[1] = (normal * (1 / Base::L2Norm(normal)))[1];
           normedNormal[2] = (normal * (1 / Base::L2Norm(normal)))[2];
    
           LinearAlgebra::SmallVector<DIM> val, phi, phi_curl, dummy;
    
           boundaryConditions(PPhys, dummy); //assumes the initial conditions and the boundary conditions match
            
           OuterProduct(normedNormal, dummy, val);
           int n = face->getPtrElementLeft()->getNrOfBasisFunctions();
           ret.resize(n);
    
           for (int i = 0; i < n; ++i)
           {
               fa.basisFunctionUnitNormal(i, phi);
    
               phi_curl = fa.basisFunctionCurl(i);
        
               ret(i) = -(phi_curl[0] * val[0] + phi_curl[1] * val[1] + phi_curl[2] * val[2]) + MaxwellData::StabCoeff_ * (phi[0] * val[0] + phi[1] * val[1] + phi[2] * val[2]);
           }
       }
    
}
//faceStiffnessIntegrandBR

void hpGemUIExtentions::anonymous9::faceIntegrand(Base::PhysicalFace<DIM>& fa, LinearAlgebra::MiddleSizeMatrix& ret)
{
    const Base::Face* face = fa.getFace();
    LinearAlgebra::SmallVector<DIM> normal = fa.getNormalVector();
    ElementT* right;
    ElementT* left = const_cast<ElementT*>(face->getPtrElementLeft());
    ElementInfos* leftInfo = static_cast<ElementInfos*>(left->getUserData());
    ElementInfos* rightInfo;
    LinearAlgebra::SmallVector<DIM> normedNormal;
    
    normedNormal[0] = (normal * (1 / Base::L2Norm(normal)))[0];
    normedNormal[1] = (normal * (1 / Base::L2Norm(normal)))[1];
    normedNormal[2] = (normal * (1 / Base::L2Norm(normal)))[2];

    double localepsilon;
    
    if (face->isInternal())
    {
        right = const_cast<ElementT*>(face->getPtrElementRight());
        rightInfo = static_cast<ElementInfos*>(right->getUserData());
    }
    
    int n = face->getPtrElementLeft()->getNrOfBasisFunctions();
    int M = face->getNrOfBasisFunctions();
    ret.resize(M, M);
    
    LinearAlgebra::SmallVector<DIM> phi_i, phi_j;
    
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
        fa.basisFunction(i, phi_i);
        for (int j = 0; j < M; ++j)
        {
            fa.basisFunctionUnitNormal(j, phi_j);
            ret(j, i) = (face->isInternal() ? 1 : 2) * (phi_i[0] * phi_j[0] + phi_i[1] * phi_j[1] + phi_i[2] * phi_j[2]) * sqrt(localepsilon);
        }
    }
    
}

//faceSpaceIntegrandBR

void hpGemUIExtentions::anonymous10::faceIntegrand(Base::PhysicalFace<DIM>& fa, LinearAlgebra::MiddleSizeVector& ret)
{
    const Base::Face* face = fa.getFace();
    LinearAlgebra::SmallVector<DIM> normal = fa.getNormalVector();
    const Geometry::PointReference<DIM - 1>& p = fa.getPointReference();
    
    int n = face->getPtrElementLeft()->getNrOfBasisFunctions();
    
    ElementT* left = const_cast<ElementT*>(face->getPtrElementLeft());
    const PointElementReferenceT& PLeft = face->mapRefFaceToRefElemL(p);
    
    PointPhysicalT PPhys;
    PPhys = left->referenceToPhysical(PLeft);
    LinearAlgebra::SmallVector<DIM> normedNormal;
    
    normedNormal[0] = (normal * (1 / Base::L2Norm(normal)))[0];
    normedNormal[1] = (normal * (1 / Base::L2Norm(normal)))[1];
    normedNormal[2] = (normal * (1 / Base::L2Norm(normal)))[2];
    
    LinearAlgebra::SmallVector<DIM> val, dummy, phi;
    boundaryConditions(PPhys, dummy); //assumes the initial conditions and the boundary conditions match
    OuterProduct(normedNormal, dummy, val);
    ret.resize(n);
    for (int i = 0; i < n; ++i)
    {
       fa.basisFunction(i, phi);
       ret(i) = 2 * (phi[0] * val[0] + phi[1] * val[1] + phi[2] * val[2]);
    }
}

//faceSpaceIntegrand

void hpGemUIExtentions::anonymous11::faceIntegrand(Base::PhysicalFace<DIM>& fa, LinearAlgebra::MiddleSizeVector& ret)
{
    const Base::Face* face = fa.getFace();
    LinearAlgebra::SmallVector<DIM> normal = fa.getNormalVector();
    const Geometry::PointReference<DIM - 1>& p = fa.getPointReference();

    
    int n = face->getPtrElementLeft()->getNrOfBasisFunctions();
    ElementT* left = const_cast<ElementT*>(face->getPtrElementLeft());
    const PointElementReferenceT& PLeft = face->mapRefFaceToRefElemL(p);
    
    PointPhysicalT PPhys;
    PPhys = left->referenceToPhysical(PLeft);
    LinearAlgebra::SmallVector<DIM> normedNormal;
    
    normedNormal[0] = (normal * (1 / Base::L2Norm(normal)))[0];
    normedNormal[1] = (normal * (1 / Base::L2Norm(normal)))[1];
    normedNormal[2] = (normal * (1 / Base::L2Norm(normal)))[2];
    
    LinearAlgebra::SmallVector<DIM> val, dummy, phi_curl;
    boundaryConditions(PPhys, dummy); //assumes the initial conditions and the boundary conditions match
    OuterProduct(normedNormal, dummy, val);
    ret.resize(n);
    for (int i = 0; i < n; ++i)
    {
        phi_curl = fa.basisFunctionCurl(i);
        ret(i) = -(phi_curl[0] * val[0] + phi_curl[1] * val[1] + phi_curl[2] * val[2]);
    }
    
}

