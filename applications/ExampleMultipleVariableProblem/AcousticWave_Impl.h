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
template<std::size_t DIM>
AcousticWave<DIM>::AcousticWave
(
 const std::size_t dimension,
 const std::size_t numOfVariables,
 const std::size_t polynomialOrder,
 const Base::ButcherTableau * const ptrButcherTableau
 ) :
Base::HpgemAPISimplified<DIM>(dimension, numOfVariables, polynomialOrder, ptrButcherTableau),
numOfVariables_(numOfVariables),
cInv_(1.0)
{
}

template<std::size_t DIM>
Base::RectangularMeshDescriptor<DIM> AcousticWave<DIM>::createMeshDescription(const std::size_t numOfElementPerDirection)
{
    // Create the domain. In this case the domain is the square [0,1]^DIM and periodic.
    Base::RectangularMeshDescriptor<DIM> description;
    for (std::size_t i = 0; i < DIM; ++i)
    {
        description.bottomLeft_[i] = 0;
        description.topRight_[i] = 1;
        description.numElementsInDIM_[i] = numOfElementPerDirection;
        if(i == 0)
        {
            description.boundaryConditions_[i] = Base::BoundaryType::SOLID_WALL;
        }
        else
        {
            description.boundaryConditions_[i] = Base::BoundaryType::PERIODIC;
        }
    }
    
    return description;
}

template<std::size_t DIM>
LinearAlgebra::MiddleSizeVector AcousticWave<DIM>::getExactSolution(const PointPhysicalT &pPhys, const double &time, const std::size_t orderTimeDerivative)
{
    LinearAlgebra::MiddleSizeVector realSolution(numOfVariables_);
    double c = std::sqrt(1.0 / cInv_); // Wave velocity.
    
    double x0 = pPhys[0];
    double x1 = 0;
    for (std::size_t iD = 1; iD < DIM; iD++) // Index for the dimension.
    {
        x1 += pPhys[iD];
    }
    
    realSolution(0) = - std::sqrt(DIM) * c * (2 * M_PI) * std::cos(2 * M_PI * x0) * std::cos(2 * M_PI * (x1 - std::sqrt(DIM) * c * time));
    realSolution(1) = - (2 * M_PI) * std::sin(2 * M_PI * x0) * std::sin(2 * M_PI * (x1 - std::sqrt(DIM) * c * time)) / cInv_;
    for (std::size_t iV = 2; iV < numOfVariables_; iV++) // iV is the index for the variable.
    {
        realSolution(iV) = (2 * M_PI) * std::cos(2 * M_PI * x0) * std::cos(2 * M_PI * (x1 - std::sqrt(DIM) * c * time)) / cInv_;
    }
    
    return realSolution;
}

/// \brief Compute the initial solution at a given point in space and time.
template<std::size_t DIM>
LinearAlgebra::MiddleSizeVector AcousticWave<DIM>::getInitialSolution(const PointPhysicalT &pPhys, const double &startTime, const std::size_t orderTimeDerivative)
{
    return getExactSolution(pPhys, startTime, orderTimeDerivative);
}

/// \details The integrand for the reference element is the same as the physical element, but scaled with the reference-to-physical element scale, which is the determinant of the jacobian of the reference-to-physical element mapping.
template<std::size_t DIM>
LinearAlgebra::MiddleSizeMatrix AcousticWave<DIM>::integrandMassMatrixOnRefElement(const Base::Element *ptrElement, const PointReferenceT &pRef)
{
    std::size_t numOfBasisFunctions = ptrElement->getNrOfBasisFunctions();
    LinearAlgebra::MiddleSizeMatrix integrand(numOfVariables_ * numOfBasisFunctions, numOfVariables_ * numOfBasisFunctions);
    PointPhysicalT pPhys = ptrElement->referenceToPhysical(pRef);
    
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
    Geometry::Jacobian<DIM, DIM> jac = ptrElement->calcJacobian(pRef);
    integrand *= std::abs(jac.determinant());
    
    return integrand;
}

/// \details The integrand for the initial solution is the exact solution at time 0 multiplied by a test function. The integrand is then scaled by the reference-to-physical element scale, since we compute the integral on a reference element.
template<std::size_t DIM>
LinearAlgebra::MiddleSizeVector AcousticWave<DIM>::integrandInitialSolutionOnRefElement
(const Base::Element *ptrElement, const double &startTime, const PointReferenceT &pRef)
{
    std::size_t numOfBasisFunctions = ptrElement->getNrOfBasisFunctions();
    
    LinearAlgebra::MiddleSizeVector integrand(numOfVariables_ * numOfBasisFunctions);
    
    PointPhysicalT pPhys = ptrElement->referenceToPhysical(pRef);
    
    LinearAlgebra::MiddleSizeVector initialSolution(getInitialSolution(pPhys, startTime));
    
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
    Geometry::Jacobian<DIM, DIM> jac = ptrElement->calcJacobian(pRef);
    integrand *= std::abs(jac.determinant());
    
    return integrand;
}

/// \details The integrand for the reference element is the same as the physical element, but scaled with the reference-to-physical element scale, which is the determinant of the jacobian of the reference-to-physical element mapping.
template<std::size_t DIM>
LinearAlgebra::MiddleSizeVector AcousticWave<DIM>::integrandRightHandSideOnRefElement(const Base::Element *ptrElement, const double &time, const PointReferenceT &pRef, const LinearAlgebra::MiddleSizeVector &solutionCoefficients)
{
    std::size_t numOfBasisFunctions = ptrElement->getNrOfBasisFunctions();
    
    LinearAlgebra::MiddleSizeVector integrand(numOfVariables_ * numOfBasisFunctions);
    
    LinearAlgebra::SmallVector<DIM> gradientBasisFunction;
    LinearAlgebra::SmallVector<DIM> gradientScalarFunction;
    double divergenceVectorFunction = 0;
    
    // Compute the gradient of the scalar function and the divergence of the vector function.
    std::size_t jVB; // Index for both basis function and variable
    for (std::size_t jD = 0; jD < DIM; jD++) // Index for derivatives
    {
        gradientScalarFunction(jD) = 0;
    }
    for (std::size_t jB = 0; jB < numOfBasisFunctions; jB++) // Index for basis functions
    {
        gradientBasisFunction = ptrElement->basisFunctionDeriv(jB, pRef);
        
        for (std::size_t jD = 0; jD < DIM; jD++) // Index for derivatives
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
        
        for (std::size_t iD = 0; iD < DIM; iD++) // Index for derivative
        {
            iVB = ptrElement->convertToSingleIndex(iB, iD + 1);
            integrand(iVB) = ptrElement->basisFunction(iB, pRef) * gradientScalarFunction(iD);
        }
    }
    
    // Scale with the reference-to-physical element ratio.
    Geometry::Jacobian<DIM, DIM> jac = ptrElement->calcJacobian(pRef);
    integrand *= std::abs(jac.determinant());
    
    return integrand;
}

/// \details The integrand for the reference face is the same as the physical face, but scaled with the reference-to-physical face scale. This face scale is absorbed in the normal vector, since it is relatively cheap to compute the normal vector with a length (L2-norm) equal to the reference-to-physical face scale.
template<std::size_t DIM>
LinearAlgebra::MiddleSizeVector AcousticWave<DIM>::integrandRightHandSideOnRefFace
(
 const Base::Face *ptrFace,
 const double &time,
 const PointReferenceOnFaceT &pRef,
 const LinearAlgebra::MiddleSizeVector &solutionCoefficients
 )
{
    std::size_t numOfBasisFunctions = ptrFace->getPtrElementLeft()->getNrOfBasisFunctions();
    
    LinearAlgebra::MiddleSizeVector integrand(numOfVariables_ * numOfBasisFunctions);
    
    LinearAlgebra::MiddleSizeVector numericalSolution(numOfVariables_);
    
    // Compute the numerical solution at the given point.
    std::size_t jVB; // Index for both variable and basis function.
    for (std::size_t jV = 0; jV < numOfVariables_; jV++)
    {
        numericalSolution(jV) = 0;
        for (std::size_t jB = 0; jB < numOfBasisFunctions; jB++)
        {
            jVB = ptrFace->getPtrElementLeft()->convertToSingleIndex(jB, jV);
            numericalSolution(jV) += ptrFace->basisFunction(Base::Side::LEFT, jB, pRef) * solutionCoefficients(jVB);
        }
    }
    
    // Compute normal vector, with size of the ref-to-phys face scale, pointing outward of the left element.
    LinearAlgebra::SmallVector<DIM> normal = ptrFace->getNormalVector(pRef);
    
    // Compute the jump of the vector function.
    double jumpVectorFunction = 0;
    for (std::size_t jD = 0; jD < DIM; jD++)
    {
        jumpVectorFunction += normal(jD) * (numericalSolution(jD + 1));
    }
    
    // Compute integrand on the reference element.
    std::size_t iVB; // Index for both variable and basis function.
    for (std::size_t iB = 0; iB < numOfBasisFunctions; iB++)
    {
        iVB = ptrFace->getPtrElementLeft()->convertToSingleIndex(iB, 0);
        integrand(iVB) = - ptrFace->basisFunction(Base::Side::LEFT, iB, pRef) * jumpVectorFunction;
        
        for (std::size_t iD = 0; iD < DIM; iD++) // Index for direction
        {
            iVB = ptrFace->getPtrElementLeft()->convertToSingleIndex(iB, iD + 1);
            integrand(iVB) = 0;
        }
    }
    return integrand;
}

/// \details The integrand for the reference face is the same as the physical face, but scaled with the reference-to-physical face scale. This face scale is absorbed in the normal vector, since it is relatively cheap to compute the normal vector with a length (L2-norm) equal to the reference-to-physical face scale.
template<std::size_t DIM>
LinearAlgebra::MiddleSizeVector AcousticWave<DIM>::integrandRightHandSideOnRefFace
(
 const Base::Face *ptrFace,
 const double &time,
 const PointReferenceOnFaceT &pRef,
 const Base::Side &iSide,
 const LinearAlgebra::MiddleSizeVector &solutionCoefficientsLeft,
 const LinearAlgebra::MiddleSizeVector &solutionCoefficientsRight
 )
{
    std::size_t numOfTestBasisFunctions = ptrFace->getPtrElement(iSide)->getNrOfBasisFunctions(); // Get the number of test basis functions on a given side, iSide
    std::size_t numOfSolutionBasisFunctionsLeft = ptrFace->getPtrElementLeft()->getNrOfBasisFunctions(); //Get the number of basis functions on the left
    std::size_t numOfSolutionBasisFunctionsRight = ptrFace->getPtrElementRight()->getNrOfBasisFunctions(); //Get the number of basis functions on the right side
    
    LinearAlgebra::MiddleSizeVector integrand(numOfVariables_ * numOfTestBasisFunctions); // Integrand value based on n number of testbasisfunctions from element corresponding to side iSide
    
    LinearAlgebra::MiddleSizeVector numericalSolutionLeft(numOfVariables_);
    LinearAlgebra::MiddleSizeVector numericalSolutionRight(numOfVariables_);
    
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
    LinearAlgebra::SmallVector<DIM> normal = ptrFace->getNormalVector(pRef);
    
    // Compute the jump of the scalar function and the vector function.
    LinearAlgebra::SmallVector<DIM> jumpScalarFunction;
    double jumpVectorFunction = 0;
    for (std::size_t jD = 0; jD < DIM; jD++)
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
        
        for (std::size_t iD = 0; iD < DIM; iD++) // Index for direction
        {
            iVB = ptrFace->getPtrElement(iSide)->convertToSingleIndex(iB, iD + 1);
            integrand(iVB) = -0.5 * ptrFace->basisFunction(iSide, iB, pRef) * jumpScalarFunction(iD);
        }
    }
    return integrand;
}


/// \details The integrand for the reference element is the same as the physical element, but scaled with the reference-to-physical element scale, which is the determinant of the jacobian of the reference-to-physical element mapping.
template<std::size_t DIM>
double AcousticWave<DIM>::integrandErrorOnRefElement
(
 const Base::Element *ptrElement,
 const double &time,
 const PointReferenceT &pRef,
 const LinearAlgebra::MiddleSizeVector &solutionCoefficients
 )
{
    std::size_t numOfBasisFunctions = ptrElement->getNrOfBasisFunctions();
    
    double integrand = 0.;
    
    PointPhysicalT pPhys = ptrElement->referenceToPhysical(pRef);
    
    LinearAlgebra::MiddleSizeVector realSolution(getExactSolution(pPhys, time));
    LinearAlgebra::MiddleSizeVector numericalSolution(numOfVariables_);
    
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
            integrand += std::pow(numericalSolution(jV) - realSolution(jV), 2);
        }
        else
        {
            integrand += getCInv(pPhys) * std::pow(numericalSolution(jV) - realSolution(jV), 2);
        }
    }
    
    // Scale with the reference-to-physical element ratio.
    Geometry::Jacobian<DIM, DIM> jac = ptrElement->calcJacobian(pRef);
    integrand *= std::abs(jac.determinant());
    
    return integrand;
}

template<std::size_t DIM>
LinearAlgebra::MiddleSizeMatrix AcousticWave<DIM>::computeMassMatrixAtElement(Base::Element *ptrElement)
{
    std::function<LinearAlgebra::MiddleSizeMatrix(const PointReferenceT &)> integrandFunction = [=](const PointReferenceT & pRef) -> LinearAlgebra::MiddleSizeMatrix{ return this -> integrandMassMatrixOnRefElement(ptrElement, pRef);};
    
    return this->elementIntegrator_.referenceElementIntegral(ptrElement->getGaussQuadratureRule(), integrandFunction);
}

template<std::size_t DIM>
LinearAlgebra::MiddleSizeVector AcousticWave<DIM>::integrateInitialSolutionAtElement(Base::Element * ptrElement, const double startTime, const std::size_t orderTimeDerivative)
{
    // Define the integrand function for the the initial solution integral.
    std::function<LinearAlgebra::MiddleSizeVector(const PointReferenceT &)> integrandFunction = [=](const PointReferenceT & pRef) -> LinearAlgebra::MiddleSizeVector { return this -> integrandInitialSolutionOnRefElement(ptrElement, startTime, pRef);};
    
    return this->elementIntegrator_.referenceElementIntegral(ptrElement->getGaussQuadratureRule(), integrandFunction);
}

/// \details The error is defined as error = realSolution - numericalSolution. The energy of the vector (u, s0, s1) is defined as u^2 + c^{-1} * |s|^2.
template<std::size_t DIM>
double AcousticWave<DIM>::integrateErrorAtElement(Base::Element *ptrElement, LinearAlgebra::MiddleSizeVector &solutionCoefficients, double time)
{
    // Define the integrand function for the error energy.
    std::function<double(const PointReferenceT &)> integrandFunction = [=](const PointReferenceT & pRef) -> double
    {
        return this->integrandErrorOnRefElement(ptrElement, time, pRef, solutionCoefficients);
    };
    
    return this->elementIntegrator_.referenceElementIntegral(ptrElement->getGaussQuadratureRule(), integrandFunction);
}

template<std::size_t DIM>
LinearAlgebra::MiddleSizeVector AcousticWave<DIM>::computeRightHandSideAtElement(Base::Element *ptrElement, LinearAlgebra::MiddleSizeVector &solutionCoefficients, const double time)
{
    // Define the integrand function for the right hand side for the reference element.
    std::function<LinearAlgebra::MiddleSizeVector(const PointReferenceT &)> integrandFunction = [=](const PointReferenceT & pRef) -> LinearAlgebra::MiddleSizeVector
    {   return this->integrandRightHandSideOnRefElement(ptrElement, time, pRef, solutionCoefficients);};
    
    return this->elementIntegrator_.referenceElementIntegral(ptrElement->getGaussQuadratureRule(), integrandFunction);
}

template<std::size_t DIM>
LinearAlgebra::MiddleSizeVector AcousticWave<DIM>::computeRightHandSideAtFace
(
 Base::Face *ptrFace,
 LinearAlgebra::MiddleSizeVector &solutionCoefficients,
 const double time
 )
{
    // Define the integrand function for the right hand side for the reference face.
    std::function<LinearAlgebra::MiddleSizeVector(const PointReferenceOnFaceT &)> integrandFunction = [=](const PointReferenceOnFaceT &pRef) -> LinearAlgebra::MiddleSizeVector
    {   return this->integrandRightHandSideOnRefFace(ptrFace, time, pRef, solutionCoefficients);};
    
    return this->faceIntegrator_.referenceFaceIntegral(ptrFace->getGaussQuadratureRule(), integrandFunction);
}

template<std::size_t DIM>
LinearAlgebra::MiddleSizeVector AcousticWave<DIM>::computeRightHandSideAtFace
(
 Base::Face *ptrFace,
 const Base::Side side,
 LinearAlgebra::MiddleSizeVector &solutionCoefficientsLeft,
 LinearAlgebra::MiddleSizeVector &solutionCoefficientsRight,
 const double time
 )
{
    // Define the integrand function for the right hand side for the reference face.
    std::function<LinearAlgebra::MiddleSizeVector(const PointReferenceOnFaceT &)> integrandFunction = [=](const PointReferenceOnFaceT &pRef) -> LinearAlgebra::MiddleSizeVector
    {   return this->integrandRightHandSideOnRefFace(ptrFace, time, pRef, side, solutionCoefficientsLeft, solutionCoefficientsRight);};
    
    return this->faceIntegrator_.referenceFaceIntegral(ptrFace->getGaussQuadratureRule(), integrandFunction);
}





