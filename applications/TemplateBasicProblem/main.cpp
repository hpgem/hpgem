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

#include <chrono>
#include <fstream>

#include "Base/HpgemAPISimplified.h"
#include "Base/TimeIntegration/AllTimeIntegrators.h"

#include "Logger.h"

//it is also possible to template ExampleProblem on the dimension if you want more flexibility, but if you do so you must prefix identifiers defined in the API with 'this->' due to the name look-up rules for templated classes
const std::size_t DIM = 2;

/// \brief Template for basis PDE problems that can be solved using HpgemAPISimplified.
/// \details This file is meant as a template for basic PDE problems that can be solved using HpgemAPISimplified.
/// When creating a solver for a certain PDE with HpgemAPISimplified one can copy this folder and fill in the blanks.
/// Please also change the name of the executable at the indicated locations in the CMakeLists.txt file.
/// You may need to re-run cmake to cause your new application to be detected.
/// Also see the documentation of HpgemAPISimplified for an overview of all the functionalities it offers.
class ExampleProblem : public Base::HpgemAPISimplified<DIM>
{
public:
    // constructor.
    /// \param[in] dimension Dimension of the domain
    /// \param[in] numOfUnknowns Number of variables in the PDE
    /// \param[in] polynomialOrder Polynomial order of the basis functions
    /// \param[in] ptrButcherTableau A butcherTableau used to solve the PDE with a Runge-Kutta method.
    /// \param[in] numOfTimeLevels Number of time levels. If a butcherTableau is set and the number of time levels is too low, this will be corrected automatically.
    ExampleProblem
    (
     const std::size_t numberOfUnknowns,
     const std::size_t polynomialOrder,
     const Base::ButcherTableau * const ptrButcherTableau = Base::AllTimeIntegrators::Instance().getRule(4, 4),
     const std::size_t numOfTimeLevels = 1
     ) :
    Base::HpgemAPISimplified<DIM>(numberOfUnknowns, polynomialOrder, ptrButcherTableau, numOfTimeLevels)
    {
        // Look at the constructor of HpgemAPISimplified to see what arguments are optional.
        logger(ERROR, "Remove this message. Make sure the constructor of this class is adapted to your purposes.");
    }
    
    /// \brief Create a rectangular mesh description
    Base::RectangularMeshDescriptor<DIM> createMeshDescription(const std::size_t numOfElementPerDirection) override final
    {
        //describes a rectangular domain
        Base::RectangularMeshDescriptor<DIM> description;
        
        for (std::size_t i = 0; i < configData_->dimension_; ++i)
        {
            description.bottomLeft_[i] = 0;
            description.topRight_[i] = 1;
            
            //Define elements in each direction.
            description.numElementsInDIM_[i] = numOfElementPerDirection;
            
            //Choose the type of boundary conditions (PERIODIC or SOLID_WALL) for each direction.
            description.boundaryConditions_[i] = Base::BoundaryType::PERIODIC;
        }
        
        return description;
    }
    
    /// \brief Compute the real solution at a given point in space and time.
    LinearAlgebra::MiddleSizeVector getExactSolution(const PointPhysicalT &pPhys, const double &time, const std::size_t orderTimeDerivative) override final
    {
        LinearAlgebra::MiddleSizeVector exactSolution(configData_->numberOfUnknowns_);
        
        // Implement here your exact solution if you have one. Otherwise delete this entire function.
        logger(ERROR, "No exact solution implemented.");
        
        return exactSolution;
    }
    
    /// \brief Compute the initial solution at a given point in space and time.
    LinearAlgebra::MiddleSizeVector getInitialSolution(const PointPhysicalT &pPhys, const double &startTime, const std::size_t orderTimeDerivative) override final
    {
        LinearAlgebra::MiddleSizeVector initialSolution(configData_->numberOfUnknowns_);
        
        // Compute the initial solution here
        logger(ERROR, "No initial solution implemented.");
        
        return initialSolution;
    }
    
    /// \brief Compute the right-hand side corresponding to an element
    LinearAlgebra::MiddleSizeVector computeRightHandSideAtElement
    (
     Base::Element *ptrElement,
     LinearAlgebra::MiddleSizeVector &solutionCoefficients,
     const double time
     ) override final
    {
        // Compute the number of basis functions
        std::size_t numOfBasisFunctions = ptrElement->getNumberOfBasisFunctions();
        
        // Declare the right-hand side at the element.
        LinearAlgebra::MiddleSizeVector rightHandSideAtElement(configData_->numberOfUnknowns_ * numOfBasisFunctions);
        
        // Compute the right hand side at the element here.
        logger(ERROR, "No function for computing the right-hand side at an element implemented.");
        
        return rightHandSideAtElement;
    }

    /// \brief Compute the right-hand side corresponding to a boundary face
    LinearAlgebra::MiddleSizeVector computeRightHandSideAtFace
    (
     Base::Face *ptrFace,
     LinearAlgebra::MiddleSizeVector &solutionCoefficients,
     const double time
     ) override final
    {
        // Compute the number of basis functions
        std::size_t numOfBasisFunctions = ptrFace->getPtrElementLeft()->getNumberOfBasisFunctions();
        
        // Declare the right-hand side at the boundary face.
        LinearAlgebra::MiddleSizeVector rightHandSideAtFace(configData_->numberOfUnknowns_ * numOfBasisFunctions);
        
        // Compute the right hand side at the boundary face here.
        logger(ERROR, "No function for computing the right-hand side at a boundary face implemented.");
        
        return rightHandSideAtFace;
    }
    
    /// \brief Compute the right-hand side corresponding to an internal face
    LinearAlgebra::MiddleSizeVector computeRightHandSideAtFace
    (
     Base::Face *ptrFace,
     const Base::Side side,
     LinearAlgebra::MiddleSizeVector &solutionCoefficientsLeft,
     LinearAlgebra::MiddleSizeVector &solutionCoefficientsRight,
     const double time
     ) override final
    {
        // Compute the number of basis functions corresponding to the element at the given side.
        std::size_t numOfBasisFunctions = ptrFace->getPtrElement(side)->getNumberOfBasisFunctions();
        
        // Declare the right-hand side at the internal face.
        LinearAlgebra::MiddleSizeVector rightHandSideAtFace(configData_->numberOfUnknowns_ * numOfBasisFunctions);
        
        // Compute the right hand side at the internal face here.
        logger(ERROR, "No function for computing the right-hand side at an internal face implemented.");
        
        return rightHandSideAtFace;
    }
    
    /// \brief Show the progress of the time integration.
    void showProgress(const double time, const std::size_t timeStepID) override final
    {
        // Choose how to show the progress.
        if (timeStepID % 10 == 0)
        {
            logger(INFO, "% time steps computed.", timeStepID);
        }
    }
    
};

auto& numOfElements = Base::register_argument<std::size_t>('n', "numElems", "number of elements per dimension", true);
auto& polynomialOrder = Base::register_argument<std::size_t>('p', "order", "polynomial order of the solution", true);

int main(int argc, char **argv)
{
    Base::parse_options(argc, argv);
    // Set parameters for the PDE. (dimension is set at the beginning of the file)
    const std::size_t numberOfVariables = 1;
    const Base::MeshType meshType = Base::MeshType::TRIANGULAR;    // Either TRIANGULAR or RECTANGULAR.
    const Base::ButcherTableau * const ptrButcherTableau = Base::AllTimeIntegrators::Instance().getRule(4, 4);
    const bool doComputeError = false;  // Set to true if you want to compute the error. (Requires exact solution)

    std::vector<std::string> variableNames(numberOfVariables);
    for(std::size_t i = 0; i < numberOfVariables; i++)
    {
        // Choose names for the variables
        logger(ERROR, "Remove this line. Make sure you have chosen names for the different variables.");
        variableNames[i] = "variable" + std::to_string(i);
    }

    // Create a problem solver, that can solve the implemented problem. Chooes your own name for the object.
    ExampleProblem problemSolver(numberOfVariables, polynomialOrder.getValue(), ptrButcherTableau);

    // Create the mesh
    problemSolver.createMesh(numOfElements.getValue(), meshType);

    // Set the names for the output file. Choose your own names.
    problemSolver.setOutputNames("output", "internalFileTitle", "solutionTitle", variableNames);

    // Solve the problem over time interval [startTime,endTime].
        // Start Measuring elapsed time
    std::chrono::time_point<std::chrono::system_clock> startClock, endClock;
    startClock = std::chrono::system_clock::now();

        // Solve the problem
    problemSolver.solve(Base::startTime.getValue(), Base::endTime.getValue(), Base::dt.getValue(), Base::numberOfSnapshots.getValue(), doComputeError);

        // Measure elapsed time
    endClock = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = endClock - startClock;
    std::cout << "Elapsed time for solving the PDE: " << elapsed_seconds.count() << "s\n";

    return 0;
}

