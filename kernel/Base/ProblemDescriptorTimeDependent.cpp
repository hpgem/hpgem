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

#include <cassert>

#include "ProblemDescriptorTimeDependent.hpp"
#include "Base/CommandLineOptions.hpp"
#include "Base/ConfigurationData.hpp"
#include "Base/Element.hpp"
#include "Base/ElementCacheData.hpp"
#include "Base/FaceCacheData.hpp"
#include "Base/GlobalData.hpp"
#include "Base/MpiContainer.hpp"
#include "Base/RectangularMeshDescriptor.hpp"
#include "Base/ShortTermStorageElementH1.hpp"
#include "Base/ShortTermStorageFaceH1.hpp"
#include "Base/TimeIntegration/AllTimeIntegrators.hpp"
#include "Integration/ElementIntegral.hpp"
#include "Integration/FaceIntegral.hpp"
#include "Integration/ReturnTrait1.hpp"
#include "Output/TecplotDiscontinuousSolutionWriter.hpp"
#include "Output/VTKTimeDependentWriter.hpp"

namespace Base
{
    class ProblemDescriptorTimeDependent;
    extern CommandLineOption<std::size_t> &numberOfSnapshots;
    extern CommandLineOption<double> &endTime;
    extern CommandLineOption<double> &startTime;
    extern CommandLineOption<double> &dt;
    extern CommandLineOption<std::string> &outputName;
    
    ///Constructor, GlobalData and ConfigurationData is deleted in HpgemUI
    ProblemDescriptorTimeDependent::ProblemDescriptorTimeDependent(std::size_t DIM, std::size_t polynomialOrder, const ButcherTableau* integrator)
    : HpgemUI(new GlobalData,
              new ConfigurationData(DIM, 1, polynomialOrder, integrator->numStages() + 1)),
              integrator_(integrator)
    {
        endTime_ = endTime.getValue();
        startTime_ = startTime.getValue();
        if (!dt.isUsed())
        {
            ///TODO: compute time step based on CFL number
            dt_ = 1e-3;
        }
        else
        {
            dt_ = dt.getValue();
        }
    }
    
    void ProblemDescriptorTimeDependent::interpolate()
    {
        for (Base::Element* element : meshes_[0]->getElementsList())
        {
            LinearAlgebra::Matrix mass = element->getMassMatrix();
            LinearAlgebra::NumericalVector solution = element->getResidue();
            mass.solve(solution);
            element->setTimeLevelData(0, solution);
        }
    }
    
    /// Compute the full right hand side of du/dt = f(u) by adding the contributions
    /// of the local element and the contributions of the faces
    void ProblemDescriptorTimeDependent::addRHSParts()
    {
        //initialise the matrices for the right-hand side
        LinearAlgebra::NumericalVector leftResidual, rightResidual, residual;

        //for every face, get the right hand side of both the local element and
        //the faces. Add the parts for both elements and save again.
        for (Base::Face* face : meshes_[0]->getFacesList())
        {
            leftResidual = face->getPtrElementLeft()->getResidue();
            residual = face->getResidue();
            if (face->isInternal())
            {
                rightResidual = face->getPtrElementRight()->getResidue();

                std::size_t numBasisFuncsLeft = face->getPtrElementLeft()->getNrOfBasisFunctions();

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
                residual += leftResidual;
                face->getPtrElementLeft()->setResidue(residual);
            }
        }
    }
    
    void ProblemDescriptorTimeDependent::computeOneTimeStep()
    {
        //iterate over the stages of the Runge Kutta method, compute temporary solutions
        for (std::size_t level = 0; level < integrator_->numStages(); ++level)
        {
            for (Base::Element* element : meshes_[0]->getElementsList())
            {
                LinearAlgebra::NumericalVector currentData = element->getTimeLevelData(0);
                for (std::size_t i = 0; i < level; ++i)
                {
                    currentData += dt_ * integrator_ -> a(level, i) * element->getTimeLevelData(i + 1);
                }
                element->setCurrentData(currentData);
            }

            //Compute the contribution of this element (local) to the right hand side
            computeRhsLocal();

            //Synchronize between nodes
            synchronize();

            //Compute the contribution of the faces (flux) to the right hand side
            computeRhsFaces();

            //Add both right hand side contributions
            addRHSParts();
            
            //Interpolate the new temporary solution from the right hand side computed before
            for (Base::Element* element : meshes_[0]->getElementsList())
            {
                LinearAlgebra::Matrix mass = element->getMassMatrix();
                LinearAlgebra::NumericalVector oldData = element->getCurrentData();
                LinearAlgebra::NumericalVector solution = element->getResidue();
                //M u_{n+1} = M u_n + dt_ * rhs
                solution = mass * oldData + dt_ * solution;
                mass.solve(solution);
                element->setTimeLevelData(level + 1, solution);
            }

        }
        
        //Combine all temporary solutions to the solution for the next time step
        for (Base::Element* element : meshes_[0]->getElementsList())
        {
            LinearAlgebra::NumericalVector newVals = element->getTimeLevelData(0);
            for (std::size_t i = 0; i < integrator_->numStages(); ++i)
            {
                newVals += dt_ * (integrator_->b(i)) * element->getTimeLevelData(i + 1);
            }
            element->setTimeLevelData(0,newVals);
        }
    }

    ///Function that computes the grid and all matrices, executes the time stepping
    ///and writes the output files
    ///
    /// \todo split this function in multiple parts.
    bool ProblemDescriptorTimeDependent::solve()
    {
        initialise(); //make the grid
        assert(checkInitialisation()); //check if we made a grid
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

            if (t == tPlot) //yes, == for doubles, but see the start of the time loop
            {
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

        for (MeshManipulator::FaceIterator citFe = Base::HpgemUI::faceColBegin(); citFe != Base::HpgemUI::faceColEnd(); ++citFe)
        {
            std::size_t numBasisFuncs = (*citFe)->getNrOfBasisFunctions();
            fMatrixData.resize(numBasisFuncs, numBasisFuncs);
            fVectorData.resize(numBasisFuncs);
            faceIntegral.integrate<LinearAlgebra::Matrix>((*citFe), this, fMatrixData);
            (*citFe)->setFaceMatrix(fMatrixData);
            faceIntegral.integrate<LinearAlgebra::NumericalVector>((*citFe), this, fVectorData);
            (*citFe)->setFaceVector(fVectorData);
        }
    }
    
    ///Compute the integrals on all elements.
    void ProblemDescriptorTimeDependent::doAllElementIntegration(std::size_t meshID)
    {
        //numberOfUnknowns_ is the number of unknowns in the "real problem" you want a
        //solution for, this is automatically set to 1 in the contructor of hpgemUISimplified
        //ndof is now the size of your element matrix
        std::size_t ndof = HpgemUI::configData_->numberOfBasisFunctions_;

        //initialise the element matrix and element vector.
        LinearAlgebra::Matrix eMatrixData(ndof, ndof);
        LinearAlgebra::NumericalVector eVectorData(ndof);

        bool isUseCache = false;
        Integration::ElementIntegral elIntegral(isUseCache);
        elIntegral.setStorageWrapper(new ShortTermStorageElementH1(meshes_[meshID]->dimension()));

        for (ElementIterator it = HpgemUI::meshes_[meshID]->elementColBegin(); it != HpgemUI::meshes_[meshID]->elementColEnd(); ++it)
        {
            elIntegral.integrate<LinearAlgebra::Matrix>((*it), this, eMatrixData);
            (*it)->setElementMatrix(eMatrixData);
            elIntegral.integrate<LinearAlgebra::NumericalVector>((*it), this, eVectorData);
            (*it)->setElementVector(eVectorData);
            (*it)->setResidue(eVectorData);
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
                    //@dducks: Shouldn't this be nobf * numberOfUnknowns?
                    assert(el->getTimeLevelDataMatrix(0).size()==configData_->numberOfBasisFunctions_);

                    //assert(el->getTimeLevelData(0).size()==configData_->numberOfBasisFunctions_);
//                    LinearAlgebra::NumericalVector timeLevelData(configData_->numberOfBasisFunctions_);
                    //std::cout<<"Receiving element "<<el->getID()<<" from process "<<it.first<<std::endl;
//                    MPIContainer::Instance().receive( timeLevelData , it.first, el->getID() * 2 + 1);
//                    el->setTimeLevelData( 0, timeLevelData );
                    MPIContainer::Instance().receive( el->getTimeLevelDataMatrix(0), it.first, el->getID() * 2 + 1);
                }
            }
        }
        for (const auto& it : pushes)
        {
            for (Element* el : it.second)
            {
                if (configData_->numberOfTimeLevels_ > 0)
                {
                    //@dducks: see note above with receive assert!
                    assert(el->getTimeLevelDataMatrix(0).size()==configData_->numberOfBasisFunctions_);
                    MPIContainer::Instance().send( el->getTimeLevelDataMatrix(0), it.first, el->getID() * 2 + 1);
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
        if (HpgemUI::meshes_.size() == 0)
        {
            return false;
        }
        return true;
    }

}


