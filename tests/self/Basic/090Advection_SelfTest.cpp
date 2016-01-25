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

#include <cmath>
#include <functional>
#include <chrono>

#include "Base/CommandLineOptions.h"
#include "Base/ConfigurationData.h"
#include "Base/Element.h"
#include "Base/Face.h"
#include "Base/HpgemAPILinear.h"
#include "Base/RectangularMeshDescriptor.h"
#include "Integration/ElementIntegral.h"
#include "Integration/FaceIntegral.h"
#include "Output/TecplotDiscontinuousSolutionWriter.h"
#include "Output/TecplotSingleElementWriter.h"
#include "Utilities/BasisFunctions2DH1ConformingTriangle.h"
#include "Logger.h"

/// This class is used to test if the advection equation is solved correctly when using HpgemAPISimplified.
/// Linear advection equation: du/dt + a[0] du/dx + a[1] du/dy = 0, or equivalently: du/dt = - a[0] du/dx - a[1] du/dy = 0.

template<std::size_t DIM>
class Advection : public Base::HpgemAPISimplified<DIM>
{
public:
    
    using typename Base::HpgemAPIBase<DIM>::PointPhysicalT;
    using typename Base::HpgemAPIBase<DIM>::PointReferenceT;
    using typename Base::HpgemAPIBase<DIM>::PointReferenceOnFaceT;
    
    ///Constructor. Assign all private variables.
    Advection(const std::size_t n, const std::size_t p, const Base::MeshType meshType) :
    Base::HpgemAPISimplified<DIM>(1, p),
    n_(n),
    meshType_(meshType)
    {
        for (std::size_t i = 0; i < DIM; ++i)
        {
            a[i] = 0.1 + 0.1 * i;
        }
    }
    
    /// Create a mesh description
    Base::RectangularMeshDescriptor<DIM> createMeshDescription(const std::size_t numberOfElementsPerDirection) override final
    {
        //describes a rectangular domain
        Base::RectangularMeshDescriptor<DIM> description;
        
        //this demo will use the square [0,1]^2
        for (std::size_t i = 0; i < DIM; ++i)
        {
            description.bottomLeft_[i] = 0;
            description.topRight_[i] = 1;
            //Define elements in each direction.
            description.numberOfElementsInDIM_[i] = numberOfElementsPerDirection;
            
            //Choose whether you want periodic boundary conditions or other (solid wall)
            //boundary conditions.
            description.boundaryConditions_[i] = Base::BoundaryType::PERIODIC;
        }
        
        return description;
    }
    
    /// \brief Compute the integrand of the right-hand side associated with elements.
    LinearAlgebra::MiddleSizeVector computeIntegrandRightHandSideAtElement
    (
     Base::PhysicalElement<DIM> &element,
     const LinearAlgebra::MiddleSizeVector &inputFunctionCoefficients,
     const double time
     )
    {
        std::size_t numberOfBasisFunctions = element.getElement()->getNumberOfBasisFunctions();
        LinearAlgebra::MiddleSizeVector&  result = element.getResultVector();
        LinearAlgebra::MiddleSizeVector::type functionValue = 0;
        for(std::size_t j = 0; j < numberOfBasisFunctions; ++j)
        {
            functionValue += inputFunctionCoefficients(j) * element.basisFunction(j);
        }
        for (std::size_t i = 0; i < numberOfBasisFunctions; ++i)
        {
            result(i) =  functionValue * (a * element.basisFunctionDeriv(i));
        }
        
        return result;
    }
    
    /// \brief Compute the integrals of the right-hand side associated with elements.
    LinearAlgebra::MiddleSizeVector computeRightHandSideAtElement
    (
     Base::Element *ptrElement,
     LinearAlgebra::MiddleSizeVector &inputFunctionCoefficients,
     const double time
     ) override final
    {
        // Define a function for the integrand of the right hand side at the element.
        std::function<LinearAlgebra::MiddleSizeVector(Base::PhysicalElement<DIM>&)> integrandFunction = [=](Base::PhysicalElement<DIM> &element) -> LinearAlgebra::MiddleSizeVector {return this->computeIntegrandRightHandSideAtElement(element, inputFunctionCoefficients, time);};
        
        return this->elementIntegrator_.integrate(ptrElement, integrandFunction);
    }
    
    /// \brief Compute the integrals of the right-hand side associated with faces.
    LinearAlgebra::MiddleSizeVector computeIntegrandRightHandSideAtFace
    (
     Base::PhysicalFace<DIM> &face,
     const Base::Side iSide,
     const LinearAlgebra::MiddleSizeVector &inputFunctionCoefficientsLeft,
     const LinearAlgebra::MiddleSizeVector &inputFunctionCoefficientsRight,
     const double time
     )
    {
        //Get the number of basis functions of the elements at both sides.
        std::size_t numberOfTestFunctions = face.getFace()->getPtrElement(iSide)->getNumberOfBasisFunctions();
        std::size_t numberOfBasisFunctionsLeft = face.getFace()->getPtrElementLeft()->getNumberOfBasisFunctions();
        std::size_t numberOfBasisFunctionsRight = face.getFace()->getPtrElementRight()->getNumberOfBasisFunctions();
        
        //Resize the result to the correct size and set all elements to 0.
        LinearAlgebra::MiddleSizeVector integrandVal = face.getResultVector(iSide);
        integrandVal *= 0;
        
        // Check if the outward pointing normal vector of the left element is in the same direction as the advection term.
        const double A = a * face.getUnitNormalVector();
        
        // Compute the sign of the normal vector (1 if iSide is left, -1 if iSide is right)
        int iSign = 1;
        if(iSide == Base::Side::RIGHT)
        {
            iSign *= -1;
        }
        
        // Compute the value of the jump times the advection term
        LinearAlgebra::MiddleSizeVector::type jump = 0;
            //Advection in the same direction as outward normal of left element:
        if (A > 1e-12)
        {
            for(std::size_t j = 0; j < numberOfBasisFunctionsLeft; ++j)
            {
                jump += inputFunctionCoefficientsLeft(j) * face.basisFunction(Base::Side::LEFT, j);
            }
            jump *= A * iSign;
        }
            //Advection in the same direction as outward normal of right element:
        else if (A < -1e-12)
        {
            for(std::size_t j = 0; j < numberOfBasisFunctionsRight; ++j)
            {
                jump += inputFunctionCoefficientsRight(j) * face.basisFunction(Base::Side::RIGHT, j);
            }
            jump *= A * iSign;
        }
            //Advection orthogonal to normal:
        else if (std::abs(A) < 1e-12)
        {
            for(std::size_t j = 0; j < numberOfBasisFunctionsLeft; ++j)
            {
                jump += inputFunctionCoefficientsLeft(j) * face.basisFunction(Base::Side::LEFT, j);
            }
            for(std::size_t j = 0; j < numberOfBasisFunctionsRight; ++j)
            {
                jump += inputFunctionCoefficientsRight(j) * face.basisFunction(Base::Side::RIGHT, j);
            }
            jump *= A * iSign / 2;
        }
        
        //Compute all entries of the integrand at this point:
        for (std::size_t i = 0; i < numberOfTestFunctions; ++i)
        {
            integrandVal(i) = -jump * face.basisFunction(iSide,i);
        }
        
        return integrandVal;
    }
    
    /// \brief Compute the integrals of the right-hand side associated with faces.
    LinearAlgebra::MiddleSizeVector computeRightHandSideAtFace
    (
     Base::Face *ptrFace,
     const Base::Side iSide,
     LinearAlgebra::MiddleSizeVector &inputFunctionCoefficientsLeft,
     LinearAlgebra::MiddleSizeVector &inputFunctionCoefficientsRight,
     const double time
     ) override final
    {
        // Define a function for the integrand of the right hand side at the face.
        std::function<LinearAlgebra::MiddleSizeVector(Base::PhysicalFace<DIM>&)> integrandFunction = [=](Base::PhysicalFace<DIM> &face) -> LinearAlgebra::MiddleSizeVector {return this->computeIntegrandRightHandSideAtFace(face, iSide, inputFunctionCoefficientsLeft, inputFunctionCoefficientsRight, time);};
        
        return this->faceIntegrator_.integrate(ptrFace, integrandFunction);
    }
    
    
    
    /// Define a solution at time zero.
    double getSolutionAtTimeZero(const PointPhysicalT& point)
    {
        double solution;
        solution = std::sin(2 * M_PI * point[0]);
        for(std::size_t i=1; i<DIM; i++)
        {
            solution *= std::sin(2 * M_PI * point[i]);
        }
        return solution;
    }
    
    /// Define the exact solution. In this case that is \f$ u_0(\vec{x}-\vec{a}t) \f$, where \f$ u_0 \f$ is the solution at time zero.
    LinearAlgebra::MiddleSizeVector getExactSolution(const PointPhysicalT& point, const double &time, const std::size_t orderTimeDerivative) override final
    {
        LinearAlgebra::MiddleSizeVector exactSolution(1);
        if(orderTimeDerivative == 0)
        {
            PointPhysicalT displacement{a*time};
            exactSolution(0) = getSolutionAtTimeZero(point - displacement);
            return exactSolution;
        }
        else
        {
            logger(ERROR, "No exact solution for order time derivative % implemented", orderTimeDerivative);
            exactSolution(0) = 0;
            return exactSolution;
        }
    }
    
    /// Define the initial conditions. In this case it is just the exact solution at the start time.
    LinearAlgebra::MiddleSizeVector getInitialSolution(const PointPhysicalT& point, const double &startTime, const std::size_t orderTimeDerivative) override final
    {
        return getExactSolution(point, startTime, orderTimeDerivative);
    }
    
    
    /// \brief Create a mesh, solve the problem and return the total error.
    LinearAlgebra::MiddleSizeVector::type createAndSolve
    (
     const double T, // final time
     const std::size_t nT // number of time steps
    )
    {
        this->createMesh(n_, meshType_);
        this->solve(0, T,  (T / nT), 0, false);
        return this->computeTotalError(this->solutionVectorId_, T);
    }
    
private:
    /// Number of elements per direction
    std::size_t n_;
    
    /// Mesh type
    Base::MeshType meshType_;
    
    ///Advective vector
    LinearAlgebra::SmallVector<DIM> a;
};


int main(int argc, char **argv)
{
    Base::parse_options(argc, argv);
    
    // Define test parameters
    const std::size_t numberOfTests = 14;
    std::array<std::size_t, numberOfTests> dim = {1,1,1,1,1, 2,2,2,2,2, 3,3,3,3};
    std::array<std::size_t, numberOfTests> p = {1,2,3,4,5, 1,2,3,4,5, 1,2,3,4};
    std::array<std::size_t, numberOfTests> n = {16,8,4,4,4, 16,8,4,4,4, 16,8,4,4};
    std::array<double, numberOfTests> T = {0.1,0.1,0.1,0.1,0.1, 0.1,0.1,0.1,0.1,0.1, 0.02,0.02,0.02,0.01};
    std::array<std::size_t, numberOfTests> nT = {10,10,10,10,10, 10,10,10,10,10, 2,2,2,1};
    std::array<std::size_t, numberOfTests> shapeId = {1,0,1,0,1, 0,1,0,1,0, 1,0,1,0};
    std::array<double, numberOfTests> errors =
    {
        0.00322137,
        0.000642002,
        0.000417087,
        3.90839e-05,
        2.83226e-06,
        0.00361477,
        0.000847142,
        0.00236532,
        4.79497e-05,
        7.70023e-05,
        0.00148109,
        0.00113544,
        0.000181175,
        0.00032545
    };
    
    // Define clocks for measuring simulation time.
    std::chrono::time_point<std::chrono::system_clock> startClock, endClock;
    std::chrono::duration<double> elapsed_seconds;
    
    // Define error
    double error;
    
    // Define mesh type
    Base::MeshType meshType;
    
    
    startClock = std::chrono::system_clock::now();
    
    // Perform tests
    for(std::size_t i=0; i<numberOfTests; i++)
    {
        
        if(shapeId[i] == 0)
        {
            meshType = Base::MeshType::TRIANGULAR;
        }
        else
        {
            meshType = Base::MeshType::RECTANGULAR;
        }
        
        if(dim[i] == 1)
        {
            Advection<1> test(n[i], p[i], meshType);
            error = test.createAndSolve(T[i], nT[i]);
            std::cout << "Error: " << error << "\n";
            logger.assert_always((std::abs(error - errors[i]) < 1e-8), "comparison to old results");
        }
        else if(dim[i] == 2)
        {
            Advection<2> test(n[i], p[i], meshType);
            error = test.createAndSolve(T[i], nT[i]);
            std::cout << "Error: " << error << "\n";
            logger.assert_always((std::abs(error - errors[i]) < 1e-8), "comparison to old results");
        }
        else
        {
            Advection<3> test(n[i], p[i], meshType);
            error = test.createAndSolve(T[i], nT[i]);
            std::cout << "Error: " << error << "\n";
            logger.assert_always((std::abs(error - errors[i]) < 1e-8), "comparison to old results");
        }
        
    }
    endClock = std::chrono::system_clock::now();
    elapsed_seconds = endClock - startClock;
    std::cout << "Elapsed time for solving the PDE: " << elapsed_seconds.count() << "s\n";
    
    
    return 0;
}
