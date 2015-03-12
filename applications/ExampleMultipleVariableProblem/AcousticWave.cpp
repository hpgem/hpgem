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

#include "AcousticWave.h"

/// \param[in] dimension Dimension of the domain
/// \param[in] numberOfVariables Number of variables in the PDE
/// \param[in] polynomialOrder Polynomial order of the basis functions
/// \param[in] useMatrixStorage Boolean to indicate if element and face matrices for the PDE should be stored
/// \param[in] ptrButcherTableau Pointer to a Butcher Tableau used to do the time integration with a Runge-Kutta scheme. By default this is a RK4 scheme.
AcousticWave::AcousticWave
(
 const std::size_t dimension,
 const std::size_t numOfVariables,
 const std::size_t polynomialOrder,
 const Base::ButcherTableau * const ptrButcherTableau
 ) :
HpgemAPISimplified(dimension, numOfVariables, polynomialOrder, ptrButcherTableau),
DIM_(dimension),
numOfVariables_(numOfVariables),
elementIntegrator_(),
faceIntegrator_(),
cInv_(1.0)
{
}

Base::RectangularMeshDescriptor AcousticWave::createMeshDescription(const std::size_t numOfElementPerDirection)
{
    // Create the domain. In this case the domain is the square [0,1]^DIM and periodic.
    Base::RectangularMeshDescriptor description(DIM_);
    for (int i = 0; i < DIM_; ++i)
    {
        description.bottomLeft_[i] = 0;
        description.topRight_[i] = 1;
        description.numElementsInDIM_[i] = numOfElementPerDirection;
        description.boundaryConditions_[i] = Base::RectangularMeshDescriptor::PERIODIC;
    }
    
    return description;
}

LinearAlgebra::NumericalVector AcousticWave::getExactSolution(const PointPhysicalT &pPhys, const double &time, const std::size_t orderTimeDerivative)
{
    LinearAlgebra::NumericalVector realSolution(numOfVariables_);
    double c = std::sqrt(1.0 / cInv_); // Wave velocity.
    
    double xWave = 0; // The inner product of the physical point and some direction vector. In this case the direction (1,..,1).
    for (std::size_t iD = 0; iD < DIM_; iD++) // Index for the dimension.
    {
        xWave += pPhys[iD];
    }
    
    for (std::size_t iV = 0; iV < numOfVariables_; iV++) // iV is the index for the variable.
    {
        realSolution(iV) = 2 * M_PI * std::cos(2 * M_PI * (xWave - std::sqrt(DIM_) * c * time));
        if (iV == 0)
        {
            realSolution(iV) *= -std::sqrt(DIM_) * c;
        }
        else
        {
            realSolution(iV) /= cInv_;
        }
    }
    
    return realSolution;
}

/// \brief Compute the initial solution at a given point in space and time.
LinearAlgebra::NumericalVector AcousticWave::getInitialSolution(const PointPhysicalT &pPhys, const double &startTime, const std::size_t orderTimeDerivative)
{
    return getExactSolution(pPhys, startTime, orderTimeDerivative);
}

/// \details The integrand for the reference element is the same as the physical element, but scaled with the reference-to-physical element scale, which is the determinant of the jacobian of the reference-to-physical element mapping.
LinearAlgebra::Matrix AcousticWave::integrandMassMatrixOnRefElement(const Base::Element *ptrElement, const Geometry::PointReference &pRef)
{
    std::size_t numOfBasisFunctions = ptrElement->getNrOfBasisFunctions();
    LinearAlgebra::Matrix integrand(numOfVariables_ * numOfBasisFunctions, numOfVariables_ * numOfBasisFunctions);
    Geometry::PointPhysical pPhys = ptrElement->referenceToPhysical(pRef);
    
    std::size_t iVB, jVB; // indices for both variable and basis function.
    for (std::size_t iV = 0; iV < numOfVariables_; iV++)
    {
        for (std::size_t iB = 0; iB < numOfBasisFunctions; iB++)
        {
            for (std::size_t jB = 0; jB < numOfBasisFunctions; jB++)
            {
                iVB = ptrElement->convertToSingleIndex(iB, iV);
                jVB = ptrElement->convertToSingleIndex(jB, iV);
                integrand(iVB, jVB) = ptrElement->basisFunction(iB, pRef) * ptrElement->basisFunction(jB, pRef);
                if (iV > 0)
                {
                    integrand(iVB, jVB) *= getCInv(pPhys);
                }
            }
        }
    }
    
    // Scale with the reference-to-physical element ratio.
    Geometry::Jacobian jac = ptrElement->calcJacobian(pRef);
    integrand *= jac.determinant();
    
    return integrand;
}

/// \details The integrand for the initial solution is the exact solution at time 0 multiplied by a test function. The integrand is then scaled by the reference-to-physical element scale, since we compute the integral on a reference element.
LinearAlgebra::NumericalVector AcousticWave::integrandInitialSolutionOnRefElement
(const Base::Element *ptrElement, const double &startTime, const Geometry::PointReference &pRef)
{
    std::size_t numOfBasisFunctions = ptrElement->getNrOfBasisFunctions();
    
    LinearAlgebra::NumericalVector integrand(numOfVariables_ * numOfBasisFunctions);
    
    Geometry::PointPhysical pPhys = ptrElement->referenceToPhysical(pRef);
    
    LinearAlgebra::NumericalVector initialSolution(getInitialSolution(pPhys, startTime));
    
    std::size_t iVB; // Index for both variable and basis function.
    for (std::size_t iV = 0; iV < numOfVariables_; iV++)
    {
        for (std::size_t iB = 0; iB < numOfBasisFunctions; iB++)
        {
            iVB = ptrElement->convertToSingleIndex(iB, iV);
            integrand(iVB) = ptrElement->basisFunction(iB, pRef) * initialSolution(iV);
            if (iV > 0)
            {
                integrand(iVB) *= getCInv(pPhys);
            }
        }
    }
    
    // Scale with the reference-to-physical element ratio.
    Geometry::Jacobian jac = ptrElement->calcJacobian(pRef);
    integrand *= jac.determinant();
    
    return integrand;
}

/// \details The integrand for the reference element is the same as the physical element, but scaled with the reference-to-physical element scale, which is the determinant of the jacobian of the reference-to-physical element mapping.
LinearAlgebra::NumericalVector AcousticWave::integrandRightHandSideOnRefElement(const Base::Element *ptrElement, const double &time, const Geometry::PointReference &pRef, const LinearAlgebra::NumericalVector &solutionCoefficients)
{
    std::size_t numOfBasisFunctions = ptrElement->getNrOfBasisFunctions();
    
    LinearAlgebra::NumericalVector integrand(numOfVariables_ * numOfBasisFunctions);
    
    LinearAlgebra::NumericalVector gradientBasisFunction(DIM_);
    LinearAlgebra::NumericalVector gradientScalarFunction(DIM_);
    double divergenceVectorFunction = 0;
    
    // Compute the gradient of the scalar function and the divergence of the vector function.
    std::size_t jVB; // Index for both basis function and variable
    for (std::size_t jD = 0; jD < DIM_; jD++) // Index for derivatives
    {
        gradientScalarFunction(jD) = 0;
    }
    for (std::size_t jB = 0; jB < numOfBasisFunctions; jB++) // Index for basis functions
    {
        gradientBasisFunction = ptrElement->basisFunctionDeriv(jB, pRef);
        
        for (std::size_t jD = 0; jD < DIM_; jD++) // Index for derivatives
        {
            jVB = ptrElement->convertToSingleIndex(jB, 0);
            gradientScalarFunction(jD) += gradientBasisFunction(jD) * solutionCoefficients(jVB);
            
            jVB = ptrElement->convertToSingleIndex(jB, jD + 1);
            divergenceVectorFunction += gradientBasisFunction(jD) * solutionCoefficients(jVB);
        }
    }
    
    // Compute integrand on the physical element.
    std::size_t iVB; // Index for both basis function and variable
    for (std::size_t iB = 0; iB < numOfBasisFunctions; iB++) // Index for basis function
    {
        iVB = ptrElement->convertToSingleIndex(iB, 0);
        integrand(iVB) = ptrElement->basisFunction(iB, pRef) * divergenceVectorFunction;
        
        for (std::size_t iD = 0; iD < DIM_; iD++) // Index for derivative
        {
            iVB = ptrElement->convertToSingleIndex(iB, iD + 1);
            integrand(iVB) = ptrElement->basisFunction(iB, pRef) * gradientScalarFunction(iD);
        }
    }
    
    // Scale with the reference-to-physical element ratio.
    Geometry::Jacobian jac = ptrElement->calcJacobian(pRef);
    integrand *= jac.determinant();
    
    return integrand;
}

/// \details The integrand for the reference face is the same as the physical face, but scaled with the reference-to-physical face scale. This face scale is absorbed in the normal vector, since it is relatively cheap to compute the normal vector with a length (L2-norm) equal to the reference-to-physical face scale.
LinearAlgebra::NumericalVector AcousticWave::integrandRightHandSideOnRefFace(const Base::Face *ptrFace, const double &time, const Geometry::PointReference &pRef, const Base::Side &iSide, const LinearAlgebra::NumericalVector &solutionCoefficientsLeft, const LinearAlgebra::NumericalVector &solutionCoefficientsRight)
{
    std::size_t numOfTestBasisFunctions = ptrFace->getPtrElement(iSide)->getNrOfBasisFunctions();
    std::size_t numOfSolutionBasisFunctionsLeft = ptrFace->getPtrElementLeft()->getNrOfBasisFunctions();
    std::size_t numOfSolutionBasisFunctionsRight = ptrFace->getPtrElementRight()->getNrOfBasisFunctions();
    
    LinearAlgebra::NumericalVector integrand(numOfVariables_ * numOfTestBasisFunctions);
    
    LinearAlgebra::NumericalVector numericalSolutionLeft(numOfVariables_);
    LinearAlgebra::NumericalVector numericalSolutionRight(numOfVariables_);
    
    // Compute the numerical solution at the given point at the left and right side.
    std::size_t jVB; // Index for both variable and basis function.
    for (std::size_t jV = 0; jV < numOfVariables_; jV++)
    {
        numericalSolutionLeft(jV) = 0;
        numericalSolutionRight(jV) = 0;
        for (std::size_t jB = 0; jB < numOfSolutionBasisFunctionsLeft; jB++)
        {
            jVB = ptrFace->getPtrElementLeft()->convertToSingleIndex(jB, jV);
            numericalSolutionLeft(jV) += ptrFace->basisFunction(Base::Side::LEFT, jB, pRef) * solutionCoefficientsLeft(jVB);
        }
        for (std::size_t jB = 0; jB < numOfSolutionBasisFunctionsRight; jB++)
        {
            jVB = ptrFace->getPtrElementRight()->convertToSingleIndex(jB, jV);
            numericalSolutionRight(jV) += ptrFace->basisFunction(Base::Side::RIGHT, jB, pRef) * solutionCoefficientsRight(jVB);
        }
    }
    
    // Compute normal vector, with size of the ref-to-phys face scale, pointing outward of the left element.
    LinearAlgebra::NumericalVector normal = ptrFace->getNormalVector(pRef);
    
    // Compute the jump of the scalar function and the vector function.
    LinearAlgebra::NumericalVector jumpScalarFunction(DIM_);
    double jumpVectorFunction = 0;
    for (std::size_t jD = 0; jD < DIM_; jD++)
    {
        jumpScalarFunction(jD) = normal(jD) * (numericalSolutionLeft(0) - numericalSolutionRight(0));
        jumpVectorFunction += normal(jD) * (numericalSolutionLeft(jD + 1) - numericalSolutionRight(jD + 1));
    }
    
    // Compute integrand on the reference element.
    std::size_t iVB; // Index for both variable and basis function.
    for (std::size_t iB = 0; iB < numOfTestBasisFunctions; iB++)
    {
        iVB = ptrFace->getPtrElement(iSide)->convertToSingleIndex(iB, 0);
        integrand(iVB) = -0.5 * ptrFace->basisFunction(iSide, iB, pRef) * jumpVectorFunction;
        
        for (std::size_t iD = 0; iD < DIM_; iD++) // Index for direction
        {
            iVB = ptrFace->getPtrElement(iSide)->convertToSingleIndex(iB, iD + 1);
            integrand(iVB) = -0.5 * ptrFace->basisFunction(iSide, iB, pRef) * jumpScalarFunction(iD);
        }
    }
    return integrand;
}


/// \details The integrand for the reference element is the same as the physical element, but scaled with the reference-to-physical element scale, which is the determinant of the jacobian of the reference-to-physical element mapping.
LinearAlgebra::NumericalVector AcousticWave::integrandErrorOnRefElement
(
 const Base::Element *ptrElement,
 const double &time,
 const Geometry::PointReference &pRef,
 const LinearAlgebra::NumericalVector &solutionCoefficients
 )
{
    std::size_t numOfBasisFunctions = ptrElement->getNrOfBasisFunctions();
    
    LinearAlgebra::NumericalVector integrand(1);
    integrand(0) = 0;
    
    Geometry::PointPhysical pPhys = ptrElement->referenceToPhysical(pRef);
    
    LinearAlgebra::NumericalVector realSolution(getExactSolution(pPhys, time));
    LinearAlgebra::NumericalVector numericalSolution(numOfVariables_);
    
    std::size_t jVB; // Index for both variable and basis function.
    for (std::size_t jV = 0; jV < numOfVariables_; jV++)
    {
        numericalSolution(jV) = 0;
        for (std::size_t jB = 0; jB < numOfBasisFunctions; jB++)
        {
            jVB = ptrElement->convertToSingleIndex(jB, jV);
            numericalSolution(jV) += ptrElement->basisFunction(jB, pRef) * solutionCoefficients(jVB);
        }
        if (jV > 0)
        {
            integrand(0) += std::pow(numericalSolution(jV) - realSolution(jV), 2);
        }
        else
        {
            integrand(0) += getCInv(pPhys) * std::pow(numericalSolution(jV) - realSolution(jV), 2);
        }
    }
    
    // Scale with the reference-to-physical element ratio.
    Geometry::Jacobian jac = ptrElement->calcJacobian(pRef);
    integrand *= jac.determinant();
    
    return integrand;
}

LinearAlgebra::Matrix AcousticWave::computeMassMatrixAtElement(const Base::Element *ptrElement)
{
    std::function<LinearAlgebra::Matrix(const Geometry::PointReference &)> integrandFunction = [=](const Geometry::PointReference & pRef) -> LinearAlgebra::Matrix{ return this -> integrandMassMatrixOnRefElement(ptrElement, pRef);};
    
    return elementIntegrator_.referenceElementIntegral(ptrElement->getGaussQuadratureRule(), integrandFunction);
}

/// \details Solve the equation \f$ Mu = r \f$ for \f$ u \f$ for a single element, where \f$ r \f$ is the right-hand sid and \f$ M \f$ is the mass matrix. The input is the right hand side here called 'solutionCoefficients' and the result is returned in this same vector.
void AcousticWave::solveMassMatrixEquationsAtElement(const Base::Element *ptrElement, LinearAlgebra::NumericalVector &solutionCoefficients)
{
    LinearAlgebra::Matrix massMatrix(computeMassMatrixAtElement(ptrElement));
    massMatrix.solve(solutionCoefficients);
}

LinearAlgebra::NumericalVector AcousticWave::integrateInitialSolutionAtElement(const Base::Element * ptrElement, const double startTime, const std::size_t orderTimeDerivative)
{
    // Define the integrand function for the the initial solution integral.
    std::function<LinearAlgebra::NumericalVector(const Geometry::PointReference &)> integrandFunction = [=](const Geometry::PointReference & pRef) -> LinearAlgebra::NumericalVector { return this -> integrandInitialSolutionOnRefElement(ptrElement, startTime, pRef);};
    
    return elementIntegrator_.referenceElementIntegral(ptrElement->getGaussQuadratureRule(), integrandFunction);
}

/// \details The error is defined as error = realSolution - numericalSolution. The energy of the vector (u, s0, s1) is defined as u^2 + c^{-1} * |s|^2.
LinearAlgebra::NumericalVector AcousticWave::integrateErrorAtElement(const Base::Element *ptrElement, LinearAlgebra::NumericalVector &solutionCoefficients, double time)
{
    // Define the integrand function for the error energy.
    std::function<LinearAlgebra::NumericalVector(const Geometry::PointReference &)> integrandFunction = [=](const Geometry::PointReference & pRef) -> LinearAlgebra::NumericalVector
    {
        return this->integrandErrorOnRefElement(ptrElement, time, pRef, solutionCoefficients);
    };
    
    return elementIntegrator_.referenceElementIntegral(ptrElement->getGaussQuadratureRule(), integrandFunction);
}

LinearAlgebra::NumericalVector AcousticWave::computeRightHandSideAtElement(const Base::Element *ptrElement, LinearAlgebra::NumericalVector &solutionCoefficients, const double time, const std::size_t orderTimeDerivative)
{
    // Define the integrand function for the right hand side for the reference element.
    std::function<LinearAlgebra::NumericalVector(const Geometry::PointReference &)> integrandFunction = [=](const Geometry::PointReference & pRef) -> LinearAlgebra::NumericalVector
    {   return this->integrandRightHandSideOnRefElement(ptrElement, time, pRef, solutionCoefficients);};
    
    return elementIntegrator_.referenceElementIntegral(ptrElement->getGaussQuadratureRule(), integrandFunction);
}

LinearAlgebra::NumericalVector AcousticWave::computeRightHandSideAtFace
(
 const Base::Face *ptrFace,
 const Base::Side side,
 LinearAlgebra::NumericalVector &solutionCoefficientsLeft,
 LinearAlgebra::NumericalVector &solutionCoefficientsRight,
 const double time,
 const std::size_t orderTimeDerivative
 )
{
    // Define the integrand function for the right hand side for the reference face.
    std::function<LinearAlgebra::NumericalVector(const Geometry::PointReference &)> integrandFunction = [=](const Geometry::PointReference &pRef) -> LinearAlgebra::NumericalVector
    {   return this->integrandRightHandSideOnRefFace(ptrFace, time, pRef, side, solutionCoefficientsLeft, solutionCoefficientsRight);};
    
    return faceIntegrator_.referenceFaceIntegral(ptrFace->getGaussQuadratureRule(), integrandFunction);
}





