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

#include "AcousticWaveLinear.h"

/// \param[in] dimension Dimension of the domain
/// \param[in] numberOfVariables Number of variables in the PDE
/// \param[in] polynomialOrder Polynomial order of the basis functions
/// \param[in] ptrButcherTableau Pointer to a Butcher Tableau used to do the time integration with a Runge-Kutta scheme. By default this is a RK4 scheme.
/// \param[in] useSourceTerm Boolean to indicate if a source term is used.
template<std::size_t DIM>
AcousticWaveLinear<DIM>::AcousticWaveLinear
(
 const std::size_t numberOfVariables,
 const std::size_t polynomialOrder,
 const TimeIntegration::ButcherTableau * const ptrButcherTableau,
 const bool useSourceTerm
 ) :
Base::HpgemAPILinear<DIM>(numberOfVariables, polynomialOrder, ptrButcherTableau, 1, useSourceTerm),
numberOfVariables_(numberOfVariables),
cInv_(1.0)
{
}

template<std::size_t DIM>
Base::RectangularMeshDescriptor<DIM> AcousticWaveLinear<DIM>::createMeshDescription(const std::size_t numberOfElementsPerDirection)
{
    // Create the domain. In this case the domain is the square [0,1]^DIM and periodic.
    Base::RectangularMeshDescriptor<DIM> description;
    for (std::size_t i = 0; i < DIM; ++i)
    {
        description.bottomLeft_[i] = 0;
        description.topRight_[i] = 1;
        description.numElementsInDIM_[i] = numberOfElementsPerDirection;
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
LinearAlgebra::MiddleSizeVector AcousticWaveLinear<DIM>::getExactSolution(const PointPhysicalT &pPhys, const double &time, const std::size_t orderTimeDerivative)
{
    LinearAlgebra::MiddleSizeVector realSolution(numberOfVariables_);
    double c = std::sqrt(1.0 / cInv_); // Wave velocity.
    
    double x0 = pPhys[0];
    double x1 = 0;
    for (std::size_t iD = 1; iD < DIM; iD++) // Index for the dimension.
    {
        x1 += pPhys[iD];
    }
    
    realSolution(0) = - std::sqrt(DIM) * c * (2 * M_PI) * std::cos(2 * M_PI * x0) * std::cos(2 * M_PI * (x1 - std::sqrt(DIM) * c * time));
    realSolution(1) = - (2 * M_PI) * std::sin(2 * M_PI * x0) * std::sin(2 * M_PI * (x1 - std::sqrt(DIM) * c * time)) / cInv_;
    for (std::size_t iV = 2; iV < numberOfVariables_; iV++) // iV is the index for the variable.
    {
        realSolution(iV) = (2 * M_PI) * std::cos(2 * M_PI * x0) * std::cos(2 * M_PI * (x1 - std::sqrt(DIM) * c * time)) / cInv_;
    }
    
    return realSolution;
}

/// \brief Compute the initial solution at a given point in space and time.
template<std::size_t DIM>
LinearAlgebra::MiddleSizeVector AcousticWaveLinear<DIM>::getInitialSolution(const PointPhysicalT &pPhys, const double &startTime, const std::size_t orderTimeDerivative)
{
    return getExactSolution(pPhys, startTime, orderTimeDerivative);
}

/// \details The integrand for the reference element is the same as the physical element, but scaled with the reference-to-physical element scale, which is the determinant of the jacobian of the reference-to-physical element mapping.
template<std::size_t DIM>
LinearAlgebra::MiddleSizeMatrix AcousticWaveLinear<DIM>::integrandMassMatrixOnRefElement(Base::PhysicalElement<DIM>& element)
{
    std::size_t numberOfBasisFunctions = element.getNumberOfBasisFunctions();
    LinearAlgebra::MiddleSizeMatrix& integrand = element.getResultMatrix();
    const PointPhysicalT& pPhys = element.getPointPhysical();
    
    std::size_t iVB, jVB; // indices for both variable and basis function.
    for (std::size_t iV = 0; iV < numberOfVariables_; iV++)
    {
        for (std::size_t iB = 0; iB < numberOfBasisFunctions; iB++)
        {
            for (std::size_t jB = 0; jB < numberOfBasisFunctions; jB++)
            {
                iVB = element.convertToSingleIndex(iB, iV);
                jVB = element.convertToSingleIndex(jB, iV);
                integrand(iVB, jVB) = element.basisFunction(iB) * element.basisFunction(jB);
                if (iV > 0)
                {
                    integrand(iVB, jVB) *= getCInv(pPhys);
                }
            }
        }
    }
    
    // Scale with the reference-to-physical element ratio.
    integrand *= element.getJacobianAbsDet();
    
    return integrand;
}

/// \details The integrand for the initial solution is the exact solution at time 0 multiplied by a test function. The integrand is then scaled by the reference-to-physical element scale, since we compute the integral on a reference element.
template<std::size_t DIM>
LinearAlgebra::MiddleSizeVector AcousticWaveLinear<DIM>::integrandInitialSolutionOnRefElement
(Base::PhysicalElement<DIM>& element, const double &startTime)
{
    std::size_t numberOfBasisFunctions = element.getNumberOfBasisFunctions();
    
    LinearAlgebra::MiddleSizeVector integrand = element.getResultVector();
    
    PointPhysicalT pPhys = element.getPointPhysical();
    
    LinearAlgebra::MiddleSizeVector initialSolution(getInitialSolution(pPhys, startTime));
    
    std::size_t iVB; // Index for both variable and basis function.
    for (std::size_t iV = 0; iV < numberOfVariables_; iV++)
    {
        for (std::size_t iB = 0; iB < numberOfBasisFunctions; iB++)
        {
            iVB = element.convertToSingleIndex(iB, iV);
            integrand(iVB) = element.basisFunction(iB) * initialSolution(iV);
            if (iV > 0)
            {
                integrand(iVB) *= getCInv(pPhys);
            }
        }
    }
    
    // Scale with the reference-to-physical element ratio.
    integrand *= element.getJacobianAbsDet();
    
    return integrand;
}

/// \details The integrand for the reference element is the same as the physical element, but scaled with the reference-to-physical element scale, which is the determinant of the jacobian of the reference-to-physical element mapping.
template<std::size_t DIM>
LinearAlgebra::MiddleSizeMatrix AcousticWaveLinear<DIM>::integrandStiffnessMatrixOnRefElement(Base::PhysicalElement<DIM>& element)
{
    std::size_t numberOfBasisFunctions = element.getNumberOfBasisFunctions();
    LinearAlgebra::MiddleSizeMatrix& integrand = element.getResultMatrix();
    LinearAlgebra::SmallVector<DIM> gradientBasisFunction;
    double valueTestFunction;
    
    std::size_t iVB, jVB; // Indices for both variable and basisfunction
    for (std::size_t jB = 0; jB < numberOfBasisFunctions; jB++)
    {
        gradientBasisFunction = element.basisFunctionDeriv(jB);
        for (std::size_t iB = 0; iB < numberOfBasisFunctions; iB++)
        {
            valueTestFunction = element.basisFunction(iB);
            
            iVB = element.convertToSingleIndex(iB, 0);
            for (std::size_t jD = 0; jD < DIM; jD++) // Index for the derivatives
            {
                jVB = element.convertToSingleIndex(jB, jD + 1);
                integrand(iVB, jVB) = gradientBasisFunction(jD) * valueTestFunction;
                
            }
            
            jVB = element.convertToSingleIndex(jB, 0);
            for (std::size_t iD = 0; iD < DIM; iD++) // Index for the derivatives
            {
                iVB = element.convertToSingleIndex(iB, iD + 1);
                integrand(iVB, jVB) = gradientBasisFunction(iD) * valueTestFunction;
            }
        }
    }
    
    // Scale with the reference-to-physical element ratio.
    integrand *= element.getJacobianAbsDet();
    
    return integrand;
}

/// \details The integrand for the reference face is the same as the physical face, but scaled with the reference-to-physical face scale. This face scale is absorbed in the normal vector, since it is relatively cheap to compute the normal vector with a length (L2-norm) equal to the reference-to-physical face scale.
template<std::size_t DIM>
LinearAlgebra::MiddleSizeMatrix AcousticWaveLinear<DIM>::integrandStiffnessMatrixOnRefFace(Base::PhysicalFace<DIM>& face, const Base::Side &iSide, const Base::Side &jSide)
{
    std::size_t numberOfSolutionBasisFunctions = face.getFace()->getPtrElement(jSide)->getNumberOfBasisFunctions();
    std::size_t numberOfTestBasisFunctions = face.getFace()->getPtrElement(iSide)->getNumberOfBasisFunctions();
    
    LinearAlgebra::MiddleSizeMatrix integrand = face.getResultMatrix(iSide, jSide);
    
    double valueTestFunction;
    double valueBasisFunction;
    
    // Compute normal vector, with size of the ref-to-phys face scale, pointing outward of the left element.
    LinearAlgebra::SmallVector<DIM> normal = face.getNormalVector();
    if (jSide == Base::Side::RIGHT)
    {
        normal *= -1;
    };
    
    // Compute the integrand
    std::size_t iVB, jVB; // Indices for both variable and basisfunction.
    for (std::size_t jB = 0; jB < numberOfSolutionBasisFunctions; jB++) // iB and jB are indices for the basis function.
    {
        valueBasisFunction = face.basisFunction(jSide, jB);
        for (std::size_t iB = 0; iB < numberOfTestBasisFunctions; iB++)
        {
            valueTestFunction = face.basisFunction(iSide, iB);
            
            iVB = face.getFace()->getPtrElement(iSide)->convertToSingleIndex(iB, 0);
            for (std::size_t jD = 0; jD < DIM; jD++) // index for the direction
            {
                jVB = face.getFace()->getPtrElement(jSide)->convertToSingleIndex(jB, jD + 1);
                if(face.isInternal())
                {
                    integrand(iVB, jVB) = -0.5 * normal(jD) * valueBasisFunction * valueTestFunction;
                }
                else
                {
                    integrand(iVB, jVB) = - normal(jD) * valueBasisFunction * valueTestFunction;
                }
                
            }
            
            jVB = face.getFace()->getPtrElement(jSide)->convertToSingleIndex(jB, 0);
            for (std::size_t iD = 0; iD < DIM; iD++) // index for the direction
            {
                iVB = face.getFace()->getPtrElement(iSide)->convertToSingleIndex(iB, iD + 1);
                if(face.isInternal())
                {
                    integrand(iVB, jVB) = -0.5 * normal(iD) * valueBasisFunction * valueTestFunction;
                }
                else
                {
                    integrand(iVB, jVB) = 0;
                }
            }
        }
    }
    
    return integrand;
}

/// \details The integrand for the reference element is the same as the physical element, but scaled with the reference-to-physical element scale, which is the determinant of the jacobian of the reference-to-physical element mapping.
template<std::size_t DIM>
double AcousticWaveLinear<DIM>::integrandErrorOnRefElement
(
 Base::PhysicalElement<DIM>& element,
 const double &time,
 const LinearAlgebra::MiddleSizeVector &solutionCoefficients
 )
{
    double integrand = 0.;
    
    const PointPhysicalT& pPhys = element.getPointPhysical();
    
    LinearAlgebra::MiddleSizeVector realSolution(getExactSolution(pPhys, time));
    const LinearAlgebra::MiddleSizeVector& numericalSolution = element.getSolution();
    
    for (std::size_t jV = 0; jV < numberOfVariables_; jV++)
    {
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
    integrand *= element.getJacobianAbsDet();
    
    return integrand;
}

template<std::size_t DIM>
LinearAlgebra::MiddleSizeMatrix AcousticWaveLinear<DIM>::computeMassMatrixAtElement(Base::Element *ptrElement)
{
    std::function<LinearAlgebra::MiddleSizeMatrix(Base::PhysicalElement<DIM>&)> integrandFunction =
        [=](Base::PhysicalElement<DIM>& element) -> LinearAlgebra::MiddleSizeMatrix
        { return this -> integrandMassMatrixOnRefElement(element);};
    
    return this->elementIntegrator_.integrate(ptrElement, integrandFunction);
}

template<std::size_t DIM>
LinearAlgebra::MiddleSizeVector AcousticWaveLinear<DIM>::integrateInitialSolutionAtElement(Base::Element * ptrElement, const double startTime, const std::size_t orderTimeDerivative)
{
    // Define the integrand function for the the initial solution integral.
    std::function<LinearAlgebra::MiddleSizeVector(Base::PhysicalElement<DIM>&)> integrandFunction = [=](Base::PhysicalElement<DIM>& element) -> LinearAlgebra::MiddleSizeVector { return this -> integrandInitialSolutionOnRefElement(element, startTime);};
    
    return this->elementIntegrator_.integrate(ptrElement, integrandFunction);
}

/// \details The error is defined as error = realSolution - numericalSolution. The energy of the vector (u, s0, s1) is defined as u^2 + c^{-1} * |s|^2.
template<std::size_t DIM>
double AcousticWaveLinear<DIM>::integrateErrorAtElement(Base::Element *ptrElement, LinearAlgebra::MiddleSizeVector &solutionCoefficients, double time)
{
    // Define the integrand function for the error energy.
    std::function<double(Base::PhysicalElement<DIM>&)> integrandFunction = [=](Base::PhysicalElement<DIM>& element) -> double
    {
        return this->integrandErrorOnRefElement(element, time, solutionCoefficients);
    };
    
    return this->elementIntegrator_.integrate(ptrElement, integrandFunction);
}

template<std::size_t DIM>
LinearAlgebra::MiddleSizeMatrix AcousticWaveLinear<DIM>::computeStiffnessMatrixAtElement(Base::Element *ptrElement)
{
    // Define the integrand function for the stiffness matrix for the element.
    std::function<LinearAlgebra::MiddleSizeMatrix(Base::PhysicalElement<DIM>&)> integrandFunction = [=](Base::PhysicalElement<DIM>& element) -> LinearAlgebra::MiddleSizeMatrix
    {return this->integrandStiffnessMatrixOnRefElement(element);};
    
    return this->elementIntegrator_.integrate(ptrElement, integrandFunction);
}

template<std::size_t DIM>
Base::FaceMatrix AcousticWaveLinear<DIM>::computeStiffnessMatrixAtFace(Base::Face *ptrFace)
{
    std::size_t numberOfBasisFunctionsLeft = 0;
    std::size_t numberOfBasisFunctionsRight = 0;
    
    numberOfBasisFunctionsLeft = ptrFace->getPtrElementLeft()->getNumberOfBasisFunctions();
    if (ptrFace->isInternal())
    {
        numberOfBasisFunctionsRight = ptrFace->getPtrElementRight()->getNumberOfBasisFunctions();
    }
    
    std::vector<Base::Side> allSides; // Vector with all sides of the face.
    allSides.push_back(Base::Side::LEFT);
    if (ptrFace->isInternal())
    {
        allSides.push_back(Base::Side::RIGHT);
    }
    
    Base::FaceMatrix stiffnessFaceMatrix(numberOfBasisFunctionsLeft * numberOfVariables_, numberOfBasisFunctionsRight * numberOfVariables_);
    
    for (Base::Side iSide : allSides)
    {
        for (Base::Side jSide : allSides)
        {
            // Define the integrand function for the stiffness matrix for the face.
            std::function<LinearAlgebra::MiddleSizeMatrix(Base::PhysicalFace<DIM>&)> integrandFunction = [=](Base::PhysicalFace<DIM>& face) -> LinearAlgebra::MiddleSizeMatrix
            {   return this->integrandStiffnessMatrixOnRefFace(face, iSide, jSide);};
            
            LinearAlgebra::MiddleSizeMatrix stiffnessMatrix(this->faceIntegrator_.integrate(ptrFace, integrandFunction));
            
            stiffnessFaceMatrix.setElementMatrix(stiffnessMatrix, iSide, jSide);
        }
    }
    return stiffnessFaceMatrix;
}
