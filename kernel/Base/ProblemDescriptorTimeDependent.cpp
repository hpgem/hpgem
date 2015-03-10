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

#include "ProblemDescriptorTimeDependent.h"
#include "Integration/ReturnTrait1.h"
#include "Base/Element.h"
#include "Base/ShortTermStorageElementH1.h"
#include "Base/ShortTermStorageFaceH1.h"
#include "GlobalData.h"
#include "Integration/FaceIntegral.h"
#include "Integration/ElementIntegral.h"
#include "ConfigurationData.h"
#include "RectangularMeshDescriptor.h"
#include "ElementCacheData.h"
#include "FaceCacheData.h"
#include "Logger.h"
#include "CommandLineOptions.h"
#include "Output/TecplotDiscontinuousSolutionWriter.h"
#include "MpiContainer.h"
#include "Output/VTKTimeDependentWriter.h"

namespace Base
{
    class ProblemDescriptorTimeDependent;
    
    extern CommandLineOption<std::size_t> &numberOfSnapshots;
    extern CommandLineOption<double> &endTime;
    extern CommandLineOption<double> &startTime;
    extern CommandLineOption<double> &dt;
    extern CommandLineOption<std::string> &outputName;
    
    ///Constructor, GlobalData and ConfigurationData is deleted in HpgemAPIBase
    
    ProblemDescriptorTimeDependent::ProblemDescriptorTimeDependent(std::size_t DIM, std::size_t polynomialOrder, const Base::ButcherTableau* integrator)
            : HpgemAPIBase(new GlobalData, new ConfigurationData(DIM, 1, polynomialOrder, integrator->getNumStages() + 1)), integrator_(integrator)
    {
        endTime_ = endTime.getValue();
        startTime_ = startTime.getValue();
        if (!dt.isUsed())
        {
            ///TODO: compute CFL number
            dt_ = 1e-3;
        }
        else
        {
            dt_ = dt.getValue();
        }
    }
    
    ///Function that executes one time step of an explicit Runge-Kutta like time
    ///integrator.
    ///\todo Find a way to nicely add the residues.
    void ProblemDescriptorTimeDependent::computeOneTimeStep()
    {
        //initialise the matrices for the right-hand side
        LinearAlgebra::NumericalVector leftResidual, rightResidual, residual;
        
        for (std::size_t level = 0; level < integrator_->getNumStages(); ++level)
        {
            //set the data we're going to work with, namely 
            //u_n + sum(a_{level,i} * k_{i+1}) 
            //with k_i the data of time level i.
            for (Base::Element *element : meshes_[0]->getElementsList())
            {
                LinearAlgebra::NumericalVector currentData = element->getTimeLevelData(0);
                for (std::size_t i = 0; i < level; ++i)
                {
                    currentData += dt_ * integrator_->getA(level, i) * element->getTimeLevelData(i + 1);
                }
                element->setCurrentData(currentData);
            }
            
            //Compute the local (element) part of the right hand side
            computeRhsLocal();
            
            //Now we need to perform the synchronisation between nodes.
            synchronize();
            
            //Compute the face (flux) part of the right hand side
            computeRhsFaces();
            
            //for every face, get the right hand side of both the local element and
            //the faces. Add the parts for both elements and save again.
            for (Base::Face *face : meshes_[0]->getFacesList())
            {
                if (face->isInternal())
                {
                    leftResidual = face->getPtrElementLeft()->getResidue();
                    
                    rightResidual = face->getPtrElementRight()->getResidue();
                    
                    residual = face->getResidue();
                    
                    std::size_t numBasisFuncsLeft = face->getPtrElementLeft()->getNrOfBasisFunctions();
                    //can we get nicer matrix?
                    leftResidual.resize(numBasisFuncsLeft);
                    rightResidual.resize(face->getNrOfBasisFunctions() - numBasisFuncsLeft);
                    
                    for (std::size_t j = 0; j < numBasisFuncsLeft; ++j)
                    {
                        leftResidual(j) += residual(j);
                    }
                    for (std::size_t j = numBasisFuncsLeft; j < residual.size(); ++j)
                    {
                        rightResidual(j - numBasisFuncsLeft) += residual(j);
                    }
                    
                    face->getPtrElementLeft()->setResidue(leftResidual);
                    face->getPtrElementRight()->setResidue(rightResidual);
                }
                else
                {
                    leftResidual = face->getPtrElementLeft()->getResidue();
                    residual = face->getResidue();
                    residual += leftResidual;
                    face->getPtrElementLeft()->setResidue(residual);
                }
            }
            
            //From the right-hand side computed earlier, compute M^{-1} * rhs
            //and save this as a time level.
            for (Base::Element *element : meshes_[0]->getElementsList())
            {
                LinearAlgebra::Matrix mass = element->getMassMatrix();
                LinearAlgebra::NumericalVector solution = element->getResidue();
                mass.solve(solution);
                element->setTimeLevelData(level + 1, solution);
            }
        }
        
        //combine all temporary solutions to the new solution:
        //u_{n+1} = u_n + sum(b_{level} k_{level}) with k being the timeLevelData
        for (Base::Element *element : meshes_[0]->getElementsList())
        {
            LinearAlgebra::NumericalVector nextTimeStep = element->getTimeLevelData(0);
            for (std::size_t level = 0; level < integrator_->getNumStages(); ++level)
            {
                nextTimeStep += dt_ * integrator_->getB(level) * element->getTimeLevelData(level + 1);
            }
            element->setTimeLevelData(0, nextTimeStep);
        }
    }
    
    ///Function that computes the grid and all matrices, executes the time stepping
    ///and writes the output files
    bool ProblemDescriptorTimeDependent::solve()
    {
        initialise(); //make the grid
        logger.assert(checkInitialisation(), "Something went wrong with initialisation of the problem."); //check if we made a grid
        doAllElementIntegration(); //compute all element integrals (except stiffness)
        doAllFaceIntegration(); //compute all face integrals (upwind flux)
        interpolate(); //compute the coefficients for the initial conditions
        
        //initialise the output files
        std::string outFileName = outputName.getValue();
#ifdef HPGEM_USE_MPI
        outFileName = outFileName + "." + std::to_string(MPIContainer::Instance().getProcessorID());
#endif
        std::ofstream outFile(outFileName + ".dat");
        std::string dimensions("012");
        Output::TecplotDiscontinuousSolutionWriter out(outFile, "solution of the problem", dimensions.substr(0, configData_->dimension_).c_str(), "u");
        if (outputName.isUsed())
        {
            outFileName = outputName.getValue();
        }
        else
        {
            outFileName = "VTK/output";
        }
        Output::VTKTimeDependentWriter VTKout(outFileName, meshes_[0]);
        
        //initialise and set the time-related variables
        double t = startTime_;
        double origDt = dt_;
        
        double dtPlot;
        if (numberOfSnapshots.getValue() > 1L)
        {
            dtPlot = (endTime_ - startTime_) / double(numberOfSnapshots.getValue() - 1);
            out.write(meshes_[0], "solution", false, this, t);
            VTKWrite(VTKout, t);
        }
        else
        {
            dtPlot = endTime_ - startTime_;
        }
        double tPlot = startTime_ + dtPlot;
        
        //start time stepping
        while (t < endTime_)
        {
            t += dt_;
            if (t > tPlot)
            {
                t -= dt_;
                dt_ = tPlot - t;
                t = tPlot;
            }
            
            computeOneTimeStep();
            
            if (t == tPlot)
            { //yes, == for doubles, but see the start of the time loop
                tPlot += dtPlot;
                out.write(meshes_[0], "solution", false, this, t);
                VTKWrite(VTKout, t);
                dt_ = origDt;
                if (tPlot > endTime_)
                {
                    tPlot = endTime_;
                }
            }
        }
        return true;
    }
    
    ///Compute the integrals on all faces. 
    void ProblemDescriptorTimeDependent::doAllFaceIntegration(std::size_t meshID)
    {
        bool useCache = false;
        
        LinearAlgebra::Matrix fMatrixData;
        LinearAlgebra::NumericalVector fVectorData;
        FaceIntegralT faceIntegral(useCache);
        faceIntegral.setStorageWrapper(new Base::ShortTermStorageFaceH1(meshes_[meshID]->dimension()));
        
        for (MeshManipulator::FaceIterator citFe = Base::HpgemAPIBase::faceColBegin(); citFe != Base::HpgemAPIBase::faceColEnd(); ++citFe)
        {
            std::size_t numBasisFuncs = (*citFe)->getNrOfBasisFunctions();
            fMatrixData = faceIntegral.integrate<LinearAlgebra::Matrix>((*citFe), this);
            (*citFe)->setFaceMatrix(fMatrixData);
            fVectorData = faceIntegral.integrate<LinearAlgebra::NumericalVector>((*citFe), this);
            (*citFe)->setFaceVector(fVectorData);
        }
    }
    
    ///Compute the integrals on all elements.
    void ProblemDescriptorTimeDependent::doAllElementIntegration(std::size_t meshID)
    {
        //numberOfUnknowns_ is the number of unknowns in the "real problem" you want a
        //solution for, this is automatically set to 1 in the contructor of hpgemUISimplified
        //ndof is now the size of your element matrix
        std::size_t ndof = HpgemAPIBase::configData_->numberOfBasisFunctions_;
        
        //initialise the element matrix and element vector.
        LinearAlgebra::Matrix eMatrixData(ndof, ndof);
        LinearAlgebra::NumericalVector eVectorData(ndof);
        
        bool isUseCache = false;
        Integration::ElementIntegral elIntegral(isUseCache);
        elIntegral.setStorageWrapper(new ShortTermStorageElementH1(meshes_[meshID]->dimension()));
        
        for (ElementIterator it = HpgemAPIBase::meshes_[meshID]->elementColBegin(); it != HpgemAPIBase::meshes_[meshID]->elementColEnd(); ++it)
        {
            eMatrixData = elIntegral.integrate<LinearAlgebra::Matrix>((*it), this);
            (*it)->setElementMatrix(eMatrixData);
            eVectorData = elIntegral.integrate<LinearAlgebra::NumericalVector>((*it), this);
            (*it)->setElementVector(eVectorData);
            (*it)->setResidue(eVectorData);
        }
    }
    
    ///Interpolates the new solution from the right hand side computed before.
    /// Default: solve Mass * x = rhs and set x as the new current data
    void ProblemDescriptorTimeDependent::interpolate()
    {
        for (Base::Element *element : meshes_[0]->getElementsList())
        {
            LinearAlgebra::Matrix mass = element->getMassMatrix();
            LinearAlgebra::NumericalVector solution = element->getResidue();
            mass.solve(solution);
            element->setTimeLevelData(0, solution);
        }
    }
    
    void ProblemDescriptorTimeDependent::synchronize(std::size_t meshID)
    {
#ifdef HPGEM_USE_MPI
        //Now, set it up.
        MeshManipulator * meshManipulator = meshes_[meshID];
        Submesh& mesh = meshManipulator->getMesh().getSubmesh();

        const auto& pushes = mesh.getPushElements();
        const auto& pulls = mesh.getPullElements();

        //recieve first for lower overhead
        for (const auto& it : pulls)
        {   
            for (Element* el : it.second)
            {   
                if (configData_->numberOfTimeLevels_ > 0)
                {   
                    logger.assert(el->getCurrentData().size() ==
                            configData_->numberOfBasisFunctions_ * configData_->numberOfUnknowns_,
                            "TimeLevelDataVector has the wrong size.");

                    //assert(el->getTimeLevelData(0).size()==configData_->numberOfBasisFunctions_);
                    //                    LinearAlgebra::NumericalVector timeLevelData(configData_->numberOfBasisFunctions_);
                    //std::cout<<"Receiving element "<<el->getID()<<" from process "<<it.first<<std::endl;
                    //                    MPIContainer::Instance().receive( timeLevelData , it.first, el->getID() * 2 + 1);
                    //                    el->setTimeLevelData( 0, timeLevelData );
                    MPIContainer::Instance().receive(el->getCurrentData(), it.first, el->getID() * 2 + 1);
                }
            }
        }
        for (const auto& it : pushes)
        {   
            for (Element* el : it.second)
            {   
                if (configData_->numberOfTimeLevels_ > 0)
                {   
                    logger.assert(el->getCurrentData().size() ==
                            configData_->numberOfBasisFunctions_ * configData_->numberOfUnknowns_,
                            "TimeLevelDataVector has the wrong size.");
                    MPIContainer::Instance().send(el->getCurrentData(), it.first, el->getID() * 2 + 1);
                    //std::cout<<"Sending element "<<el->getID()<<" to process "<<it.first<<std::endl;
                    //                    MPIContainer::Instance().send(el->getTimeLevelData(0), it.first, el->getID() * 2 + 1);
                }
            }
        }
        MPIContainer::Instance().sync();

#endif
    }
    
    bool ProblemDescriptorTimeDependent::checkInitialisation()
    {
        if (HpgemAPIBase::meshes_.size() == 0)
        {
            return false;
        }
        return true;
    }

}
