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

#include "HpgemUISimplified.hpp"
#include "Integration/ReturnTrait1.hpp"
#include "Base/Element.hpp"
#include "Base/ShortTermStorageElementH1.hpp"
#include "Base/ShortTermStorageFaceH1.hpp"
#include "GlobalData.hpp"
#include "Integration/FaceIntegral.hpp"
#include "Integration/ElementIntegral.hpp"
#include "ConfigurationData.hpp"
#include "RectangularMeshDescriptor.hpp"
#include "ElementCacheData.hpp"
#include "FaceCacheData.hpp"
#include "cassert"
#include "CommandLineOptions.hpp"
#include "Output/TecplotDiscontinuousSolutionWriter.hpp"
#include "MpiContainer.hpp"

namespace Base
{
    class HpgemUISimplified;
    
        auto& numberOfSnapshots = Base::register_argument<std::size_t>(0, "nOutputFrames", "Number of frames to output", false,1);
        auto& endTime = Base::register_argument<double>(0, "endTime", "end time of the simulation", false,1);
        auto& startTime = Base::register_argument<double>(0, "startTime", "start time of the simulation", false,0);
        auto& dt = Base::register_argument<double>(0, "dt", "time step of the simulation", false);
        
    HpgemUISimplified::HpgemUISimplified(unsigned int DIM,int polynomialOrder)
      :  HpgemUI(new GlobalData,
                 new ConfigurationData(DIM,1,polynomialOrder,numberOfSnapshots.getValue())) 
    {
        endTime_=endTime.getValue();
        startTime_=startTime.getValue();
        if(!dt.isUsed()){
            ///TODO: compute CFL number
            dt_=1e-3;
        }else{
            dt_=dt.getValue();
        }
    }
    
    auto& outputName = Base::register_argument<std::string>(0, "outFile", "Name of the output file", false,"output.dat");
    
    bool
    HpgemUISimplified::solve()
    {
        initialise();
        assert(checkInitialisation());
        doAllElementIntegration();
        doAllFaceIntegration();
        beforeTimeIntegration();
        interpolate();
        
        std::string outFileName=outputName.getValue();
#ifdef HPGEM_USE_MPI
        outFileName=outFileName+"."+std::to_string(MPIContainer::Instance().getProcessorID());
#endif
        std::ofstream outFile(outFileName);
        std::string dimensions("012");
        Output::TecplotDiscontinuousSolutionWriter out(outFile, "solution of the problem", dimensions.substr(0,configData_->dimension_).c_str(), "u");
        double t=startTime_;
        double dtPlot;
        double origDt=dt_;
        LinearAlgebra::Matrix leftResidual,rightResidual,residual;
        if(numberOfSnapshots.getValue()>1L) {
            dtPlot=(endTime_-startTime_)/double(numberOfSnapshots.getValue()-1);
            out.write(meshes_[0],"t="+std::to_string(t),false,this);
        }else{
            dtPlot=endTime_-startTime_;
        }
        double tPlot=startTime_+dtPlot;
        while (t<endTime_) {
            t+=dt_;
            if(t>tPlot){
                t-=dt_;
                dt_=tPlot-t;
                t=tPlot;
            }
            
            computeLocalResidual();

            for(auto& it:meshes_[0]->getMesh().getSubmesh().getPullElements())
            {
                for(Base::Element* element:it.second)
                {
                    residual.resize(configData_->numberOfUnknowns_,element->getNrOfBasisFunctions());
                    element->setResidue(residual);
                }
            }
            
            //Now we need to perform the synchronisation between nodes.
            synchronize();
            
            computeFluxResidual();
            for(Base::Face* face:meshes_[0]->getFacesList()){
                if(face->isInternal()){
                    leftResidual=face->getPtrElementLeft()->getResidue();
                    rightResidual=face->getPtrElementRight()->getResidue();
                    residual=face->getResidue();
                    int n=face->getPtrElementLeft()->getNrOfBasisFunctions();
                    //can we get nicer matrix?
                    assert(n==leftResidual.getNCols());
                    assert(residual.getNCols()-n==rightResidual.getNCols());
                    for(std::size_t i=0;i<residual.getNRows();++i){
                        for(std::size_t j=0;j<n;++j){
                            leftResidual(i,j)+=residual(i,j);
                        }
                        for(std::size_t j=n;j<residual.getNCols();++j){
                            rightResidual(i,j-n)+=residual(i,j);
                        }
                    }
                    face->getPtrElementLeft()->setResidue(leftResidual);
                    face->getPtrElementRight()->setResidue(rightResidual);
                }else{
                    leftResidual=face->getPtrElementLeft()->getResidue();
                    residual=face->getResidue();
                    residual.axpy(1.,leftResidual);
                    face->getPtrElementLeft()->setResidue(residual);
                }
                
            }
            interpolate();
            
            if(t==tPlot){//yes, == for doubles, but see the start of the time loop
                tPlot+=dtPlot;
                out.write(meshes_[0],"t="+std::to_string(t),false,this);
                dt_=origDt;
                if(tPlot>endTime_){
                    tPlot=endTime_;
                }
            }
        }

        
        return true;
    }
    
    void
    HpgemUISimplified::doAllFaceIntegration(unsigned int meshID)
    {
        bool useCache   = false;
        
        LinearAlgebra::Matrix fMatrixData;
        LinearAlgebra::NumericalVector fVectorData;
        FaceIntegralT   faceIntegral(useCache);
        faceIntegral.setStorageWrapper(new Base::ShortTermStorageFaceH1(meshes_[meshID]->dimension()));
        
        for (MeshManipulator::FaceIterator citFe = Base::HpgemUI::faceColBegin(); citFe != Base::HpgemUI::faceColEnd(); ++citFe)
        {
        	int n=(*citFe)->getPtrElementLeft()->getNrOfUnknows()*(*citFe)->getPtrElementLeft()->getNrOfBasisFunctions();
        	if((*citFe)->isInternal())
        		n+=(*citFe)->getPtrElementRight()->getNrOfUnknows()*(*citFe)->getPtrElementRight()->getNrOfBasisFunctions();
        	fMatrixData.resize(n,n);
        	fVectorData.resize(n);
            faceIntegral.integrate<LinearAlgebra::Matrix>((*citFe), this, fMatrixData);
            (*citFe)->setFaceMatrix(fMatrixData);
            faceIntegral.integrate<LinearAlgebra::NumericalVector>((*citFe), this, fVectorData);
            (*citFe)->setFaceVector(fVectorData);
        }
    }
    
    void
    HpgemUISimplified::doAllElementIntegration(unsigned int meshID)
    {
        unsigned int ndof = HpgemUI::configData_->numberOfBasisFunctions_*HpgemUI::configData_->numberOfUnknowns_;
        LinearAlgebra::Matrix  	eMatrixData(ndof, ndof);
        LinearAlgebra::Matrix   initialData(configData_->numberOfUnknowns_,configData_->numberOfBasisFunctions_);
        LinearAlgebra::NumericalVector eVectorData(ndof);
        
        bool isUseCache(false);
        Integration::ElementIntegral 	elIntegral(isUseCache);
        elIntegral.setStorageWrapper(new ShortTermStorageElementH1(meshes_[meshID]->dimension()));
        
        for (ElementIterator it=HpgemUI::meshes_[meshID]->elementColBegin(); it!= HpgemUI::meshes_[meshID]->elementColEnd(); ++it)
        {
            elIntegral.integrate<LinearAlgebra::Matrix>((*it), this, eMatrixData);
            (*it)->setElementMatrix(eMatrixData);
            elIntegral.integrate<LinearAlgebra::NumericalVector>((*it), this, eVectorData);
            (*it)->setElementVector(eVectorData);
            if(configData_->numberOfUnknowns_==1){
                assert(initialData.getNCols()==eVectorData.size());
                for(std::size_t i=0;i<eVectorData.size();++i){
                    initialData(0,i)=eVectorData[i];
                }
                (*it)->setResidue(initialData);
            }
            

            //cout << result;
            //cout<< "#####################################END of ELEMENT######"<<endl;
        }
    }
    
    void HpgemUISimplified::synchronize(std::size_t meshID) {
#ifdef HPGEM_USE_MPI
        //Now, set it up.
        MeshManipulator * meshManipulator = meshes_[meshID];
        Submesh& mesh = meshManipulator->getMesh().getSubmesh();
        
        const auto& pushes = mesh.getPushElements();
        const auto& pulls = mesh.getPullElements();
        
        //recieve first for lower overhead
        for (const auto& it : pulls) {
            for (Element* el : it.second) {
                if(configData_->numberOfTimeLevels_>0)
                {
                    //std::cout<<"Receiving element "<<el->getID()<<" from process "<<it.first<<std::endl;
                    MPIContainer::Instance().receive(el->getTimeLevelData(0), it.first, el->getID() * 2 + 1);
                }
            }
        }
        for (const auto& it : pushes) {
            for (Element* el : it.second) {
                if(configData_->numberOfTimeLevels_>0)
                {
                    //std::cout<<"Sending element "<<el->getID()<<" to process "<<it.first<<std::endl;
                    MPIContainer::Instance().send(el->getTimeLevelData(0), it.first, el->getID() * 2 + 1);
                }
            }
        }
        MPIContainer::Instance().sync();
        
#endif
    }
    
    bool
    HpgemUISimplified::checkInitialisation()
    {
        if (HpgemUI::meshes_.size()==0)
        {
            std::cerr << "Error no mesh created : You need to create at least one mesh to solve a problem" << std::endl;
            return false;
        }
        return true;
    }
    
}
