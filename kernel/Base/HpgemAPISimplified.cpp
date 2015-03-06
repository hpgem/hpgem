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

#include "HpgemAPISimplified.hpp"

#include "Base/CommandLineOptions.hpp"
#include "Base/ConfigurationData.hpp"
#include "Base/Element.hpp"
#include "Base/Face.hpp"
#include "Base/MpiContainer.hpp"
#include "Base/RectangularMeshDescriptor.hpp"
#include "Base/TimeIntegration/AllTimeIntegrators.hpp"
#include "Geometry/PointReference.hpp"
#include "Integration/ElementIntegral.hpp"
#include "Integration/FaceIntegral.hpp"
#include "Output/TecplotDiscontinuousSolutionWriter.hpp"
#include "Output/TecplotSingleElementWriter.hpp"

//temporary; should be removed before the intesive week is over
#undef assert
#include "Logger.h"

namespace Base
{
    class HpgemAPISimplified;
    
    auto& numberOfSnapshots = Base::register_argument<std::size_t>(0, "nOutputFrames", "Number of frames to output", false, 1);
    auto& endTime = Base::register_argument<double>(0, "endTime", "end time of the simulation", false, 1);
    auto& startTime = Base::register_argument<double>(0, "startTime", "start time of the simulation", false, 0);
    auto& dt = Base::register_argument<double>(0, "dt", "time step of the simulation", false);
    auto& outputName = Base::register_argument<std::string>(0, "outFile", "Name of the output file (without extentions)", false, "output");
    
    /// \param[in] dimension Dimension of the domain
    /// \param[in] numberOfVariables Number of variables in the PDE
    /// \param[in] polynomialOrder Polynomial order of the basis functions
    /// \param[in] numberOfTimeLevels Number of time levels
    /// \param[in] useMatrixStorage Boolean to indicate if element and face matrices for the PDE should be stored
    HpgemAPISimplified::HpgemAPISimplified
    (
     const std::size_t dimension,
     const std::size_t numberOfVariables,
     const std::size_t polynomialOrder,
     const std::size_t numberOfTimeLevels,
     const bool useMatrixStorage
     ) :
    HpgemUI(new Base::GlobalData, new Base::ConfigurationData(dimension, numberOfVariables, polynomialOrder, numberOfTimeLevels)),
    startTime_(0.0),
    endTime_(0.0),
    useMatrixStorage_(useMatrixStorage),
    massMatrixID_(1),
    stiffnessElementMatrixID_(0),
    stiffnessFaceMatrixID_(0)
    {
    }
    
    void HpgemAPISimplified::createMassMatrices(const double time)
    {
        for(Base::Element *ptrElement : meshes_[0]->getElementsList())
        {
            LinearAlgebra::Matrix massMatrix(computeMassMatrixAtElement(ptrElement, time));
            ptrElement->setElementMatrix(massMatrix, massMatrixID_);
        }
    }
    
    /// \details Solve the equation \f$ Mu = r \f$ for \f$ u \f$, where \f$ r \f$ is the right-hand sid and \f$ M \f$ is the mass matrix.
    void HpgemAPISimplified::solveMassMatrixEquations(const std::size_t timeLevel, const double time)
    {
        for(Base::Element *ptrElement : meshes_[0]->getElementsList())
        {
            LinearAlgebra::NumericalVector &solutionCoefficients(ptrElement->getTimeLevelDataVector(timeLevel));
            
            if(useMatrixStorage_)
            {
                const LinearAlgebra::Matrix &massMatrix(ptrElement->getElementMatrix(massMatrixID_));
                massMatrix.solve(solutionCoefficients);
            }
            else
            {
                LinearAlgebra::Matrix massMatrix(computeMassMatrixAtElement(ptrElement, time));
                massMatrix.solve(solutionCoefficients);
            }
        }
    }
    
    void HpgemAPISimplified::integrateInitialSolution(const std::size_t timeLevelResult, const double startTime, const std::size_t orderTimeDerivative)
    {
        for(Base::Element *ptrElement : meshes_[0]->getElementsList())
        {
            LinearAlgebra::NumericalVector &solutionCoefficients = ptrElement->getTimeLevelDataVector(timeLevelResult);
            solutionCoefficients = integrateInitialSolutionAtElement(ptrElement, startTime, orderTimeDerivative);
        }
    }
    
    void HpgemAPISimplified::createStiffnessMatrices(const double time)
    {
        if(!useMatrixStorage_) return;
        
        std::cout << "- Creating stiffness matrices for the elements.\n";
        for(Base::Element *ptrElement : meshes_[0]->getElementsList())
        {
            LinearAlgebra::Matrix stiffnessMatrix( computeStiffnessMatrixAtElement(ptrElement, time) );
            ptrElement->setElementMatrix(stiffnessMatrix, stiffnessElementMatrixID_);
        }
        
        std::cout << "- Creating stiffness matrices for the faces.\n";
        for(Base::Face *ptrFace : meshes_[0]->getFacesList())
        {
            Base::FaceMatrix stiffnessFaceMatrix( computeStiffnessMatrixAtFace(ptrFace, time) );
            ptrFace->setFaceMatrix(stiffnessFaceMatrix, stiffnessFaceMatrixID_);
        }
    }
    
    /// \param[in] solutionTimeLevel Time level where the solution is stored.
    /// \param[in] time Time corresponding to the current solution.
    /// \details The energy-norm of the error is the square-root of the total energy of the error. The energy of the error is defined here as \f$ e^tMe\f$, where \f$ M\f$ is the mass matrix and \f$ e\f$ is the error. In most cases this is equivalent to the \f$ L^2 \f$ norm of the error.
    double HpgemAPISimplified::computeEnergyNormError(const std::size_t solutionTimeLevel, const double time)
    {
        LinearAlgebra::NumericalVector totalError(1);
        totalError(0) = 0;
        
        for(Base::Element *ptrElement : meshes_[0]->getElementsList())
        {
            LinearAlgebra::NumericalVector &solutionCoefficients = ptrElement->getTimeLevelDataVector(solutionTimeLevel);
            totalError += integrateEnergyErrorAtElement(ptrElement, time, solutionCoefficients);
        }
        
#ifdef HPGEM_USE_MPI
        int world_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        int world_size;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        
        double errorToAdd;
        for(std::size_t iRank = 1; iRank < world_size; iRank++)
        {
            if(world_rank == 0)
            {
                MPI_Recv(&errorToAdd, 1, MPI_DOUBLE, iRank, iRank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                totalError(0) += errorToAdd;
            }
            else if(world_rank == iRank)
            {
                errorToAdd = totalError(0);
                MPI_Send(&errorToAdd, 1, MPI_DOUBLE, 0, iRank, MPI_COMM_WORLD);
            }
        }
        
        if(world_rank == 0)
        {
            double error;
            if(totalError(0) >= 0)
            {
                error = std::sqrt(totalError(0));
            }
            else
            {
                logger(WARN,"Warning: the total energy of the error is negative.\n");
                error = std::sqrt(-totalError(0));
            }
            std::cout << "Energy norm of the error: " << error << ".\n";
            return error;
        }
        else
        {
            return -1.0;
        }
#endif
        double error;
        if(totalError(0) >= 0)
        {
            error = std::sqrt(totalError(0));
        }
        else
        {
            logger(WARN,"Warning: the total energy of the error is negative.\n");
            error = std::sqrt(-totalError(0));
        }
        std::cout << "Energy norm of the error: " << error << ".\n";
        return error;
    }
    
    /// \brief Compute the right hand side for the solution at time level 'timeLevelIn' and store the result at time level 'timeLevelResult'. Make sure timeLevelIn is different from timeLevelResult.
    void HpgemAPISimplified::computeRightHandSide(const std::size_t timeLevelIn, const std::size_t timeLevelResult, const double time)
    {
        // Apply the right hand side corresponding to integration on the elements.
        for(Base::Element *ptrElement : meshes_[0]->getElementsList())
        {
            LinearAlgebra::NumericalVector &solutionCoefficients(ptrElement->getTimeLevelDataVector(timeLevelIn));
            LinearAlgebra::NumericalVector &solutionCoefficientsNew(ptrElement->getTimeLevelDataVector(timeLevelResult));
            
            solutionCoefficientsNew = computeRightHandSideAtElement(ptrElement, time, solutionCoefficients);
        }
        
        // Apply the right hand side corresponding to integration on the faces.
        for(Base::Face *ptrFace : meshes_[0]->getFacesList())
        {
            LinearAlgebra::NumericalVector &solutionCoefficientsLeft( ptrFace->getPtrElementLeft()->getTimeLevelDataVector(timeLevelIn) );
            LinearAlgebra::NumericalVector &solutionCoefficientsRight( ptrFace->getPtrElementRight()->getTimeLevelDataVector(timeLevelIn) );
            LinearAlgebra::NumericalVector &solutionCoefficientsLeftNew( ptrFace->getPtrElementLeft()->getTimeLevelDataVector(timeLevelResult) );
            LinearAlgebra::NumericalVector &solutionCoefficientsRightNew( ptrFace->getPtrElementRight()->getTimeLevelDataVector(timeLevelResult) );
            
            solutionCoefficientsLeftNew += computeRightHandSideAtFace(ptrFace, time, Base::Side::LEFT, solutionCoefficientsLeft, solutionCoefficientsRight);
            solutionCoefficientsRightNew += computeRightHandSideAtFace(ptrFace, time, Base::Side::RIGHT, solutionCoefficientsLeft, solutionCoefficientsRight);
        }
    }
    
    LinearAlgebra::NumericalVector HpgemAPISimplified::getSolutionCoefficients(const Base::Element *ptrElement, const std::vector<std::size_t> timeLevelsIn, const std::vector<double> coefficientsTimeLevels)
    {
        logger.assert(timeLevelsIn.size() == coefficientsTimeLevels.size(),"Number of time levels and number of coefficients should be the same.");
        logger.assert(timeLevelsIn.size() > 0,"Number of time levels should be bigger than zero.");
        
        LinearAlgebra::NumericalVector solutionCoefficients(ptrElement->getTimeLevelDataVector(timeLevelsIn[0]) );
        solutionCoefficients *= coefficientsTimeLevels[0];
        for(std::size_t i = 1; i < timeLevelsIn.size(); i++)
        {
            solutionCoefficients.axpy(coefficientsTimeLevels[i], ptrElement->getTimeLevelDataVector(timeLevelsIn[i]) );
        }
        return solutionCoefficients;
    }
    
    /// \details Make sure timeLevelResult is different from the timeLevelsIn.
    void HpgemAPISimplified::computeRightHandSide(const std::vector<std::size_t> timeLevelsIn, const std::vector<double> coefficientsTimeLevels, const std::size_t timeLevelResult, const double time)
    {
        // Apply the right hand side corresponding to integration on the elements.
        for(Base::Element *ptrElement : meshes_[0]->getElementsList())
        {
            LinearAlgebra::NumericalVector solutionCoefficients(getSolutionCoefficients(ptrElement, timeLevelsIn, coefficientsTimeLevels) );
            LinearAlgebra::NumericalVector &solutionCoefficientsNew(ptrElement->getTimeLevelDataVector(timeLevelResult));
            
            solutionCoefficientsNew = computeRightHandSideAtElement(ptrElement, time, solutionCoefficients);
        }
        
        // Apply the right hand side corresponding to integration on the faces.
        for(Base::Face *ptrFace : meshes_[0]->getFacesList())
        {
            LinearAlgebra::NumericalVector solutionCoefficientsLeft( getSolutionCoefficients(ptrFace->getPtrElementLeft(), timeLevelsIn, coefficientsTimeLevels) );
            LinearAlgebra::NumericalVector solutionCoefficientsRight( getSolutionCoefficients(ptrFace->getPtrElementRight(), timeLevelsIn, coefficientsTimeLevels) );
            LinearAlgebra::NumericalVector &solutionCoefficientsLeftNew( ptrFace->getPtrElementLeft()->getTimeLevelDataVector(timeLevelResult) );
            LinearAlgebra::NumericalVector &solutionCoefficientsRightNew( ptrFace->getPtrElementRight()->getTimeLevelDataVector(timeLevelResult) );
            
            solutionCoefficientsLeftNew += computeRightHandSideAtFace(ptrFace, time, Base::Side::LEFT, solutionCoefficientsLeft, solutionCoefficientsRight);
            solutionCoefficientsRightNew += computeRightHandSideAtFace(ptrFace, time, Base::Side::RIGHT, solutionCoefficientsLeft, solutionCoefficientsRight);
        }
    }
    
    void HpgemAPISimplified::synchronize(const std::size_t timeLevel)
    {
#ifdef HPGEM_USE_MPI
        //Now, set it up.
        Base::MeshManipulator * meshManipulator = meshes_[0];
        Base::Submesh& mesh = meshManipulator->getMesh().getSubmesh();
        
        const auto& pushes = mesh.getPushElements();
        const auto& pulls = mesh.getPullElements();
        
        //recieve first for lower overhead
        for (const auto& it : pulls)
        {
            for (Base::Element* ptrElement : it.second)
            {
                logger.assert(ptrElement->getTimeLevelDataVector(timeLevel).size() == configData_->numberOfBasisFunctions_ * configData_->numberOfUnknowns_ , "Size of time level data vector is wrong.");
                
                Base::MPIContainer::Instance().receive(ptrElement->getTimeLevelDataVector(timeLevel), it.first, ptrElement->getID());
            }
        }
        for (const auto& it : pushes)
        {
            for (Base::Element* ptrElement : it.second)
            {
                logger.assert(ptrElement->getTimeLevelDataVector(timeLevel).size() == configData_->numberOfBasisFunctions_ * configData_->numberOfUnknowns_, "Size of time level data vector is wrong.");
                
                Base::MPIContainer::Instance().send(ptrElement->getTimeLevelDataVector(timeLevel), it.first, ptrElement->getID());
            }
        }
        Base::MPIContainer::Instance().sync();
#endif
    }
    
    void HpgemAPISimplified::scaleTimeLevel(const std::size_t timeLevel, const double scale)
    {
        for(Base::Element *ptrElement : meshes_[0]->getElementsList())
        {
            ptrElement->getTimeLevelDataVector(timeLevel) *= scale;
        }
        
        synchronize(timeLevel);
    }
    
    void HpgemAPISimplified::scaleAndAddTimeLevel(const std::size_t timeLevelToChange, const std::size_t timeLevelToAdd, const double scale)
    {
        for(Base::Element *ptrElement : meshes_[0]->getElementsList())
        {
            ptrElement->getTimeLevelDataVector(timeLevelToChange).axpy(scale, ptrElement->getTimeLevelDataVector(timeLevelToAdd));
        }
        
        synchronize(timeLevelToChange);
    }
    
    void HpgemAPISimplified::setInitialSolution(const std::size_t solutionTimeLevel, const double startTime)
    {
        integrateInitialSolution(solutionTimeLevel, startTime, 0);
        solveMassMatrixEquations(solutionTimeLevel, startTime);
        synchronize(solutionTimeLevel);
    }
    
    /// \details Computing the time derivative in this case means applying the right hand side for the solution at time level 'timeLevelIn' and then solving the mass matrix equations. The result is stored at time level 'timeLevelResult'.
    void HpgemAPISimplified::computeTimeDerivative(const std::size_t timeLevelIn, const std::size_t timeLevelResult, const double time)
    {
        computeRightHandSide(timeLevelIn, timeLevelResult, time);
        solveMassMatrixEquations(timeLevelResult, time);
        synchronize(timeLevelResult);
    }
    
    /// \details Computing the time derivative in this case means applying the right hand side for the linear combination of solutions at time levels 'timeLevelsIn' with coefficients given by coefficientsTimeLevels, and then solving the mass matrix equations. The result is stored at time level 'timeLevelResult'.
    void HpgemAPISimplified::computeTimeDerivative(const std::vector<std::size_t> timeLevelsIn, const std::vector<double> coefficientsTimeLevels, const std::size_t timeLevelResult, const double time)
    {
        computeRightHandSide(timeLevelsIn, coefficientsTimeLevels, timeLevelResult, time);
        solveMassMatrixEquations(timeLevelResult, time);
        synchronize(timeLevelResult);
    }
    
}

