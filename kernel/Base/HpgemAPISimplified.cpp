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

#include "HpgemAPISimplified.h"

#include "Base/CommandLineOptions.h"
#include "Base/ConfigurationData.h"
#include "Base/Element.h"
#include "Base/Face.h"
#include "Base/MpiContainer.h"
#include "Base/RectangularMeshDescriptor.h"
#include "Base/TimeIntegration/AllTimeIntegrators.h"
#include "Geometry/PointReference.h"
#include "Integration/ElementIntegral.h"
#include "Integration/FaceIntegral.h"
#include "Output/TecplotDiscontinuousSolutionWriter.h"
#include "Output/TecplotSingleElementWriter.h"
#include "Utilities/BasisFunctions1DH1ConformingLine.h"
#include "Utilities/BasisFunctions2DH1ConformingSquare.h"
#include "Utilities/BasisFunctions2DH1ConformingTriangle.h"
#include "Utilities/BasisFunctions3DH1ConformingCube.h"
#include "Utilities/BasisFunctions3DH1ConformingTetrahedron.h"

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
    /// \param[in] numOfVariables Number of variables in the PDE
    /// \param[in] polynomialOrder Polynomial order of the basis functions
    /// \param[in] butcherTableau A butcherTableau used to solve the PDE with a Runge-Kutta method.
    /// \param[in] numOfTimeLevels Number of time levels. If a butcherTableau is set and the number of time levels is too low, this will be corrected automatically.
    HpgemAPISimplified::HpgemAPISimplified
    (
     const std::size_t dimension,
     const std::size_t numOfVariables,
     const std::size_t polynomialOrder,
     const Base::ButcherTableau * const ptrButcherTableau,
     const std::size_t numOfTimeLevels
     ) :
    HpgemAPIBase(new Base::GlobalData, new Base::ConfigurationData(dimension, numOfVariables, polynomialOrder, (ptrButcherTableau->getNumStages() + 1 > numOfTimeLevels) ? ptrButcherTableau->getNumStages() + 1 : numOfTimeLevels)),
    ptrButcherTableau_(ptrButcherTableau)
    {
        solutionTimeLevel_ = 0;
        for (std::size_t i = 1; i < configData_->numberOfTimeLevels_; i++)
        {
            intermediateTimeLevels_.push_back(i);
        }
    }
    
    void HpgemAPISimplified::createMesh(const std::size_t numOfElementsPerDirection, const Base::MeshType meshType)
    {
        const Base::RectangularMeshDescriptor description = createMeshDescription(numOfElementsPerDirection);
        
        // Set the number of Element/Face Matrices/Vectors.
        std::size_t numOfElementMatrices = 0;
        std::size_t numOfElementVectors = 0;
        std::size_t numOfFaceMatrices = 0;
        std::size_t numOfFaceVectors = 0;
        
        // Create mesh and set basis functions.
        if (configData_->dimension_ == 2)
        {
            if (meshType == Base::MeshType::TRIANGULAR)
            {
                addMesh(description, Base::TRIANGULAR, numOfElementMatrices, numOfElementVectors, numOfFaceMatrices, numOfFaceVectors);
                meshes_[0]->setDefaultBasisFunctionSet(Utilities::createDGBasisFunctionSet2DH1Triangle(configData_->polynomialOrder_));
            }
            else if (meshType == Base::MeshType::RECTANGULAR)
            {
                addMesh(description, Base::RECTANGULAR, numOfElementMatrices, numOfElementVectors, numOfFaceMatrices, numOfFaceVectors);
                meshes_[0]->setDefaultBasisFunctionSet(Utilities::createDGBasisFunctionSet2DH1Square(configData_->polynomialOrder_));
            }
        }
        else if (configData_->dimension_ == 3)
        {
            if (meshType == Base::MeshType::TRIANGULAR)
            {
                addMesh(description, Base::TRIANGULAR, numOfElementMatrices, numOfElementVectors, numOfFaceMatrices, numOfFaceVectors);
                meshes_[0]->setDefaultBasisFunctionSet(Utilities::createDGBasisFunctionSet3DH1Tetrahedron(configData_->polynomialOrder_));
            }
            else if (meshType == Base::MeshType::RECTANGULAR)
            {
                addMesh(description, Base::RECTANGULAR, numOfElementMatrices, numOfElementVectors, numOfFaceMatrices, numOfFaceVectors);
                meshes_[0]->setDefaultBasisFunctionSet(Utilities::createDGBasisFunctionSet3DH1Cube(configData_->polynomialOrder_));
            }
        }
        
        std::size_t nElements = meshes_[0]->getNumberOfElements();
        logger(VERBOSE, "Total number of elements: %", nElements);
    }
    
    /// \details Solve the equation \f$ Mu = r \f$ for \f$ u \f$, where \f$ r \f$ is the right-hand sid and \f$ M \f$ is the mass matrix.
    void HpgemAPISimplified::solveMassMatrixEquations(const std::size_t timeLevel)
    {
        for (Base::Element *ptrElement : meshes_[0]->getElementsList())
        {
            LinearAlgebra::NumericalVector &solutionCoefficients(ptrElement->getTimeLevelDataVector(timeLevel));
            
            solveMassMatrixEquationsAtElement(ptrElement, solutionCoefficients);
        }
    }
    
    void HpgemAPISimplified::integrateInitialSolution(const std::size_t timeLevelResult, const double intialTime, const std::size_t orderTimeDerivative)
    {
        for (Base::Element *ptrElement : meshes_[0]->getElementsList())
        {
            LinearAlgebra::NumericalVector &solutionCoefficients = ptrElement->getTimeLevelDataVector(timeLevelResult);
            solutionCoefficients = integrateInitialSolutionAtElement(ptrElement, intialTime, orderTimeDerivative);
        }
    }
    
    /// \param[in] solutionTimeLevel Time level where the solution is stored.
    /// \param[in] time Time corresponding to the current solution.
    /// \details The square of the total error is defined as \f[ \int \|e\|^2 \,dV \f], where \f$\|e\|\f$ is some user-defined norm of the error.
    double HpgemAPISimplified::computeTotalError(const std::size_t solutionTimeLevel, const double time)
    {
        LinearAlgebra::NumericalVector totalError(1);
        totalError(0) = 0;
        
        for (Base::Element *ptrElement : meshes_[0]->getElementsList())
        {
            LinearAlgebra::NumericalVector &solutionCoefficients = ptrElement->getTimeLevelDataVector(solutionTimeLevel);
            totalError += integrateErrorAtElement(ptrElement, solutionCoefficients, time);
        }
        
        double error;
        
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
            if(totalError(0) >= 0)
            {   
                error = std::sqrt(totalError(0));
            }
            else
            {   
                logger(WARN,"Warning: the total energy of the error is negative.");
                error = std::sqrt(-totalError(0));
            }
        }
        
        
        for(std::size_t iRank = 1; iRank < world_size; iRank++)
        {
            if(world_rank == 0)
            {
                MPI_Send(&error, 1, MPI_DOUBLE, iRank, iRank, MPI_COMM_WORLD);
            }
            else if(world_rank == iRank)
            {
                MPI_Recv(&error, 1, MPI_DOUBLE, 0, iRank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        
        return error;
#endif
        
        if (totalError(0) >= 0)
        {
            error = std::sqrt(totalError(0));
        }
        else
        {
            logger(WARN, "Warning: the total energy of the error is negative.");
            error = std::sqrt(-totalError(0));
        }
        logger(INFO, "Energy norm of the error: %.", error);
        return error;
    }
    
    /// \brief Compute the right hand side for the solution at time level 'timeLevelIn' and store the result at time level 'timeLevelResult'. Make sure timeLevelIn is different from timeLevelResult.
    void HpgemAPISimplified::computeRightHandSide(const std::size_t timeLevelIn, const std::size_t timeLevelResult, const double time, const std::size_t orderTimeDerivative)
    {
        // Apply the right hand side corresponding to integration on the elements.
        for (Base::Element *ptrElement : meshes_[0]->getElementsList())
        {
            LinearAlgebra::NumericalVector &solutionCoefficients(ptrElement->getTimeLevelDataVector(timeLevelIn));
            LinearAlgebra::NumericalVector &solutionCoefficientsNew(ptrElement->getTimeLevelDataVector(timeLevelResult));
            
            solutionCoefficientsNew = computeRightHandSideAtElement(ptrElement,  solutionCoefficients, time, orderTimeDerivative);
        }
        
        // Apply the right hand side corresponding to integration on the faces.
        for (Base::Face *ptrFace : meshes_[0]->getFacesList())
        {
            LinearAlgebra::NumericalVector &solutionCoefficientsLeft(ptrFace->getPtrElementLeft()->getTimeLevelDataVector(timeLevelIn));
            LinearAlgebra::NumericalVector &solutionCoefficientsRight(ptrFace->getPtrElementRight()->getTimeLevelDataVector(timeLevelIn));
            LinearAlgebra::NumericalVector &solutionCoefficientsLeftNew(ptrFace->getPtrElementLeft()->getTimeLevelDataVector(timeLevelResult));
            LinearAlgebra::NumericalVector &solutionCoefficientsRightNew(ptrFace->getPtrElementRight()->getTimeLevelDataVector(timeLevelResult));
            
            solutionCoefficientsLeftNew += computeRightHandSideAtFace(ptrFace, Base::Side::LEFT, solutionCoefficientsLeft, solutionCoefficientsRight, time, orderTimeDerivative);
            solutionCoefficientsRightNew += computeRightHandSideAtFace(ptrFace, Base::Side::RIGHT, solutionCoefficientsLeft, solutionCoefficientsRight, time, orderTimeDerivative);
        }
    }
    
    LinearAlgebra::NumericalVector HpgemAPISimplified::getSolutionCoefficients(const Base::Element *ptrElement, const std::vector<std::size_t> timeLevelsIn, const std::vector<double> coefficientsTimeLevels)
    {
        logger.assert(timeLevelsIn.size() == coefficientsTimeLevels.size(), "Number of time levels and number of coefficients should be the same.");
        logger.assert(timeLevelsIn.size() > 0, "Number of time levels should be bigger than zero.");
        
        LinearAlgebra::NumericalVector solutionCoefficients(ptrElement->getTimeLevelDataVector(timeLevelsIn[0]));
        solutionCoefficients *= coefficientsTimeLevels[0];
        for (std::size_t i = 1; i < timeLevelsIn.size(); i++)
        {
            solutionCoefficients.axpy(coefficientsTimeLevels[i], ptrElement->getTimeLevelDataVector(timeLevelsIn[i]));
        }
        return solutionCoefficients;
    }
    
    /// \details Make sure timeLevelResult is different from the timeLevelsIn.
    void HpgemAPISimplified::computeRightHandSide(const std::vector<std::size_t> timeLevelsIn, const std::vector<double> coefficientsTimeLevels, const std::size_t timeLevelResult, const double time, const std::size_t orderTimeDerivative)
    {
        // Apply the right hand side corresponding to integration on the elements.
        for (Base::Element *ptrElement : meshes_[0]->getElementsList())
        {
            LinearAlgebra::NumericalVector solutionCoefficients(getSolutionCoefficients(ptrElement, timeLevelsIn, coefficientsTimeLevels));
            LinearAlgebra::NumericalVector &solutionCoefficientsNew(ptrElement->getTimeLevelDataVector(timeLevelResult));
            
            solutionCoefficientsNew = computeRightHandSideAtElement(ptrElement,  solutionCoefficients, time, orderTimeDerivative);
        }
        
        // Apply the right hand side corresponding to integration on the faces.
        for (Base::Face *ptrFace : meshes_[0]->getFacesList())
        {
            LinearAlgebra::NumericalVector solutionCoefficientsLeft(getSolutionCoefficients(ptrFace->getPtrElementLeft(), timeLevelsIn, coefficientsTimeLevels));
            LinearAlgebra::NumericalVector solutionCoefficientsRight(getSolutionCoefficients(ptrFace->getPtrElementRight(), timeLevelsIn, coefficientsTimeLevels));
            LinearAlgebra::NumericalVector &solutionCoefficientsLeftNew(ptrFace->getPtrElementLeft()->getTimeLevelDataVector(timeLevelResult));
            LinearAlgebra::NumericalVector &solutionCoefficientsRightNew(ptrFace->getPtrElementRight()->getTimeLevelDataVector(timeLevelResult));
            
            solutionCoefficientsLeftNew += computeRightHandSideAtFace(ptrFace, Base::Side::LEFT, solutionCoefficientsLeft, solutionCoefficientsRight, time, orderTimeDerivative);
            solutionCoefficientsRightNew += computeRightHandSideAtFace(ptrFace, Base::Side::RIGHT, solutionCoefficientsLeft, solutionCoefficientsRight, time, orderTimeDerivative);
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
            for (Base::Element *ptrElement : it.second)
            {   
                logger.assert(ptrElement->getTimeLevelDataVector(timeLevel).size() == configData_->numberOfBasisFunctions_ * configData_->numberOfUnknowns_ , "Size of time level data vector is wrong.");

                Base::MPIContainer::Instance().receive(ptrElement->getTimeLevelDataVector(timeLevel), it.first, ptrElement->getID());
            }
        }
        for (const auto& it : pushes)
        {   
            for (Base::Element *ptrElement : it.second)
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
        for (Base::Element *ptrElement : meshes_[0]->getElementsList())
        {
            ptrElement->getTimeLevelDataVector(timeLevel) *= scale;
        }
        
        synchronize(timeLevel);
    }
    
    void HpgemAPISimplified::scaleAndAddTimeLevel(const std::size_t timeLevelToChange, const std::size_t timeLevelToAdd, const double scale)
    {
        for (Base::Element *ptrElement : meshes_[0]->getElementsList())
        {
            ptrElement->getTimeLevelDataVector(timeLevelToChange).axpy(scale, ptrElement->getTimeLevelDataVector(timeLevelToAdd));
        }
        
        synchronize(timeLevelToChange);
    }
    
    void HpgemAPISimplified::setInitialSolution(const std::size_t solutionTimeLevel, const double intialTime, const std::size_t orderTimeDerivative)
    {
        integrateInitialSolution(solutionTimeLevel, intialTime, orderTimeDerivative);
        solveMassMatrixEquations(solutionTimeLevel);
        synchronize(solutionTimeLevel);
    }
    
    /// \details Computing the time derivative in this case means applying the right hand side for the solution at time level 'timeLevelIn' and then solving the mass matrix equations. The result is stored at time level 'timeLevelResult'.
    void HpgemAPISimplified::computeTimeDerivative(const std::size_t timeLevelIn, const std::size_t timeLevelResult, const double time, const std::size_t orderTimeDerivative)
    {
        computeRightHandSide(timeLevelIn, timeLevelResult, time, orderTimeDerivative);
        solveMassMatrixEquations(timeLevelResult);
        synchronize(timeLevelResult);
    }
    
    /// \details Computing the time derivative in this case means applying the right hand side for the linear combination of solutions at time levels 'timeLevelsIn' with coefficients given by coefficientsTimeLevels, and then solving the mass matrix equations. The result is stored at time level 'timeLevelResult'.
    void HpgemAPISimplified::computeTimeDerivative(const std::vector<std::size_t> timeLevelsIn, const std::vector<double> coefficientsTimeLevels, const std::size_t timeLevelResult, const double time, const std::size_t orderTimeDerivative)
    {
        computeRightHandSide(timeLevelsIn, coefficientsTimeLevels, timeLevelResult, time, orderTimeDerivative);
        solveMassMatrixEquations(timeLevelResult);
        synchronize(timeLevelResult);
    }
    
    void HpgemAPISimplified::computeOneTimeStep(double &time, const double dt)
    {
        std::size_t numOfStages = ptrButcherTableau_->getNumStages();
        
        // Compute intermediate Runge-Kutta stages
        for (std::size_t iStage = 0; iStage < numOfStages; iStage++)
        {
            double stageTime = time + ptrButcherTableau_->getC(iStage) * dt;
            
            std::vector<std::size_t> timeLevelsIn;
            std::vector<double> coefficientsTimeLevels;
            
            timeLevelsIn.push_back(solutionTimeLevel_);
            coefficientsTimeLevels.push_back(1);
            for (std::size_t jStage = 0; jStage < iStage; jStage++)
            {
                timeLevelsIn.push_back(intermediateTimeLevels_[jStage]);
                coefficientsTimeLevels.push_back(dt * ptrButcherTableau_->getA(iStage, jStage));
            }
            
            computeTimeDerivative(timeLevelsIn, coefficientsTimeLevels, intermediateTimeLevels_[iStage], stageTime, 0);
        }
        
        // Update the solution
        for (std::size_t jStage = 0; jStage < numOfStages; jStage++)
        {
            scaleAndAddTimeLevel(solutionTimeLevel_, intermediateTimeLevels_[jStage], dt * ptrButcherTableau_->getB(jStage));
        }
        
        // Update the time.
        time += dt;
    }
    
    /// \param[in] outputFileName Name of the output file (minus extensions like .dat).
    /// \param[in] internalFileTitle Title of the file as used by Tecplot internally.
    /// \param[in] solutionTitle Title of the solution.
    /// \param[in] variableNames String of variable names. The string should have the form "nameVar1,nameVar2,..,nameVarN".
    void HpgemAPISimplified::setOutputNames(std::string outputFileName, std::string internalFileTitle, std::string solutionTitle, std::string variableNames)
    {
        outputFileName_ = outputFileName;
#ifdef HPGEM_USE_MPI
        outputFileName_ = outputFileName_ + "." + std::to_string(Base::MPIContainer::Instance().getProcessorID());
#endif
        outputFileName_ = outputFileName_ + ".dat";
        
        internalFileTitle_ = internalFileTitle;
        solutionTitle_ = solutionTitle;
        variableNames_ = variableNames;
    }
    
    void HpgemAPISimplified::writeToTecplotFile(const ElementT *ptrElement, const PointReferenceT &pRef, std::ostream &out)
    {
        std::size_t numOfVariables = configData_->numberOfUnknowns_;
        
        LinearAlgebra::NumericalVector solution(numOfVariables);
        solution = ptrElement->getSolution(solutionTimeLevel_, pRef);
         
        std::size_t iV = 0; // Index for the variable
        out << solution(iV);
        for (iV = 1; iV < numOfVariables; iV++)
        {
            out << " " << solution(iV);
        }
    }
    
    /// \brief Solve the PDE over the time domain [intialTime, finalTime].
    /// \param[in] intialTime Initial time
    /// \param[in] finalTime End time
    /// \param[in] dt Size of the time step
    /// \param[in] numOutputFrames Number of times the solution is written to an output file.
    /// \param[in] doComputeError Boolean to indicate if the error should be computed.
    bool HpgemAPISimplified::solve(const double intialTime, const double finalTime, double dt, const std::size_t numOfOutputFrames, bool doComputeError)
    {
        // Create output file for Tecplot.
        std::ofstream outFile(outputFileName_);
        std::string dimensionsToWrite = "";
        for(std::size_t i = 0; i < configData_->dimension_; i++)
        {
            dimensionsToWrite = dimensionsToWrite + std::to_string(i);
        }
        Output::TecplotDiscontinuousSolutionWriter tecplotWriter(outFile, solutionTitle_, dimensionsToWrite, variableNames_);
        
        // Compute parameters for time integration
        double T = finalTime - intialTime;     // Time interval
        std::size_t numOfTimeSteps = T / dt;
        std::size_t numOfTimeStepsForOutput;
        if (numOfOutputFrames > 0)
        {
            // Round off to above such that the number of time steps is a multiple of the number of output frames.
            numOfTimeSteps += (numOfOutputFrames - (numOfTimeSteps % numOfOutputFrames)) % numOfOutputFrames;
            
            // Recompute dt.
            dt = T / numOfTimeSteps;
            
            // Compute the number of timesteps after which to create an output frame.
            numOfTimeStepsForOutput = (std::size_t) numOfTimeSteps / numOfOutputFrames;
        }
        
        // Set the initial time.
        double time = intialTime;
        
        // Create and Store things before solving the problem.
        tasksBeforeSolving();
        
        // Set the initial numerical solution.
        logger(INFO, "Computing and interpolating the initial solution.");
        setInitialSolution(solutionTimeLevel_, time, 0);
        tecplotWriter.write(meshes_[0], "discontinuous solution", false, this, time);
        
        // Solve the system of PDE's.
        logger(INFO,"Solving the system of PDE's.");
        logger(INFO, "dt: %.", dt);
        logger(INFO, "Total number of time steps: %.", numOfTimeSteps);
        logger(INFO, "Number of time steps for output: %.", numOfTimeStepsForOutput);
        for (std::size_t iT = 1; iT <= numOfTimeSteps; iT++)
        {
            computeOneTimeStep(time, dt);
            
            if (iT % numOfTimeStepsForOutput == 0)
            {
                tecplotWriter.write(meshes_[0], solutionTitle_, false, this, time);
            }
            showProgress(time, iT);
        }
        
        // Compute the energy norm of the error
        if(doComputeError)
        {
            double error = computeTotalError(solutionTimeLevel_, finalTime);
            logger(INFO, "Energy norm of the error: %.", error);
        }
        
        return true;
    }
}

