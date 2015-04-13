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

#include "SavageHutterRightHandSideComputer.h"
#include "SavageHutter.h"
#include "Logger.h"

using LinearAlgebra::NumericalVector;

/// \details The integrand for the reference element is the same as the physical element, but scaled with the reference-to-physical element scale, which is the determinant of the jacobian of the reference-to-physical element mapping.
NumericalVector SavageHutterRightHandSideComputer::integrandRightHandSideOnRefElement
(const Base::Element *ptrElement, const double &time, const Geometry::PointReference &pRef, const NumericalVector &solutionCoefficients)
{
    const std::size_t numBasisFuncs = ptrElement->getNrOfBasisFunctions();
    
    NumericalVector integrand(numOfVariables_ * numBasisFuncs);
    const NumericalVector numericalSolution = computeNumericalSolution(ptrElement, pRef, solutionCoefficients);
    const NumericalVector physicalFlux = computePhysicalFlux(numericalSolution);
    const NumericalVector source = computeSourceTerm(numericalSolution);
    
    // Compute integrand on the physical element.
    std::size_t iVB; // Index for both basis function and variable
    for (std::size_t iB = 0; iB < numBasisFuncs; iB++) // Index for basis function
    {
        iVB = ptrElement->convertToSingleIndex(iB, 0);
        integrand(iVB) = physicalFlux(0) * ptrElement->basisFunctionDeriv(iB, pRef)(0);
        integrand(iVB) += source(0) * ptrElement->basisFunction(iB, pRef);

        iVB = ptrElement->convertToSingleIndex(iB, 1);
        integrand(iVB) = physicalFlux(1) * ptrElement->basisFunctionDeriv(iB, pRef)(0);
        integrand(iVB) += source(1) * ptrElement->basisFunction(iB, pRef);
    }
    
    logger(DEBUG, "Integrand on element: %", integrand);
    return integrand;
}

/// \details The integrand for the reference face is the same as the physical face, but scaled with the reference-to-physical face scale. This face scale is absorbed in the normal vector, since it is relatively cheap to compute the normal vector with a length (L2-norm) equal to the reference-to-physical face scale.
NumericalVector SavageHutterRightHandSideComputer::integrandRightHandSideOnRefFace
( const Base::Face *ptrFace, const Base::Side &iSide, const NumericalVector &normalVec, const Geometry::PointReference &pRef, const NumericalVector &solutionCoefficientsLeft, const NumericalVector &solutionCoefficientsRight)
{
    logger.assert(ptrFace != nullptr, "gave an empty face");
    logger.assert(ptrFace->getPtrElement(iSide) != nullptr, "acquired an empty element");
    double normal = normalVec(0);
    const std::size_t numOfTestBasisFunctions = ptrFace->getPtrElement(iSide)->getNrOfBasisFunctions();
    const std::size_t numOfBasisFuncsLeft = ptrFace->getPtrElement(Base::Side::LEFT)->getNrOfBasisFunctions();
    const std::size_t numOfBasisFuncsRight = ptrFace->getPtrElement(Base::Side::RIGHT)->getNrOfBasisFunctions();
    
    NumericalVector solutionLeft(2);
    for (std::size_t i = 0; i < numOfBasisFuncsLeft; ++i)    
    {
        std::size_t iH = ptrFace->getPtrElement(Base::Side::LEFT)->convertToSingleIndex(i, 0);
        solutionLeft(0) += solutionCoefficientsLeft(iH) * ptrFace->basisFunction(Base::Side::LEFT, i, pRef);
        std::size_t iHu = ptrFace->getPtrElement(Base::Side::LEFT)->convertToSingleIndex(i, 1);
        solutionLeft(1) += solutionCoefficientsLeft(iHu) * ptrFace->basisFunction(Base::Side::LEFT, i, pRef);
    }
    
    NumericalVector solutionRight(2);
    for (std::size_t i = 0; i < numOfBasisFuncsRight; ++i)    
    {
        std::size_t iH = ptrFace->getPtrElement(Base::Side::RIGHT)->convertToSingleIndex(i, 0);
        solutionRight(0) += solutionCoefficientsRight(iH) * ptrFace->basisFunction(Base::Side::RIGHT, i, pRef);
        std::size_t iHu = ptrFace->getPtrElement(Base::Side::RIGHT)->convertToSingleIndex(i, 1);
        solutionRight(1) += solutionCoefficientsRight(iHu) * ptrFace->basisFunction(Base::Side::RIGHT, i, pRef);
    }
    NumericalVector flux(2);
    if (normal > 0)
    {
        flux = localLaxFriedrichsFlux(solutionLeft, solutionRight);
    }
    else
    {
        flux = localLaxFriedrichsFlux(solutionRight, solutionLeft);
    }
    
    if (iSide == Base::Side::RIGHT) //the normal is defined for the left element
    {
        normal *= -1;
    }
    NumericalVector integrand(numOfVariables_ * numOfTestBasisFunctions); // Integrand value based on n number of testbasisfunctions from element corresponding to side iSide

    for (std::size_t iFun = 0; iFun < numOfTestBasisFunctions; ++iFun)
    {
        for (std::size_t iVar = 0; iVar < numOfVariables_; ++iVar)
        {
            std::size_t iVarFun = ptrFace->getPtrElement(iSide)->convertToSingleIndex(iFun, iVar);
            integrand(iVarFun) = -flux(iVar) * ptrFace->basisFunction(iSide, iFun, Geometry::PointReference(0)) * normal;            
        }
    }
    
    return integrand;
}

NumericalVector SavageHutterRightHandSideComputer::integrandRightHandSideOnRefFace
    (
     const Base::Face *ptrFace,
     const Geometry::PointReference &pRef,
     const NumericalVector &solutionCoefficients
     )
{
    double normal = ptrFace->getNormalVector(Geometry::PointReference(0))(0);
    const std::size_t numBasisFuncs = ptrFace->getNrOfBasisFunctions();
    NumericalVector solution(2);
    for (std::size_t i = 0; i < numBasisFuncs; ++i)    
    {
        std::size_t iH = ptrFace->getPtrElement(Base::Side::LEFT)->convertToSingleIndex(i, 0);
        solution(0) += solutionCoefficients(iH) * ptrFace->basisFunction(i, pRef);
        std::size_t iHu = ptrFace->getPtrElement(Base::Side::LEFT)->convertToSingleIndex(i, 1);
        solution(1) += solutionCoefficients(iHu) * ptrFace->basisFunction(i, pRef);
    }
    NumericalVector flux(2);
    if (normal > 0) //outflow
    {
        flux = localLaxFriedrichsFlux(solution, solution);
    }
    else //inflow
    {
        flux = localLaxFriedrichsFlux(NumericalVector({2, 5}), solution);
    }
    
    NumericalVector integrand(numOfVariables_ * numBasisFuncs);

    for (std::size_t iFun = 0; iFun < numBasisFuncs; ++iFun)
    {
        for (std::size_t iVar = 0; iVar < numOfVariables_; ++iVar)
        {
            std::size_t iVarFun = ptrFace->getPtrElementLeft()->convertToSingleIndex(iFun, iVar);
            integrand(iVarFun) = -flux(iVar) * ptrFace->basisFunction(iFun, Geometry::PointReference(0)) * normal;            
        }
    }
    
    return integrand;
}

NumericalVector SavageHutterRightHandSideComputer::computePhysicalFlux(const NumericalVector &numericalSolution)
{    
    const double h = numericalSolution(0);
    const double hu = numericalSolution(1);
    double u = 0;
    if (h > 1e-10)
    {
        u = hu/h;
    }
    NumericalVector flux(2);
    flux(0) = hu;
    flux(1) = hu * u + epsilon_/2 * std::cos(theta_) * h * h;
    logger(DEBUG, "flux values: %, %", flux(0), flux(1));
    return flux;
}

NumericalVector SavageHutterRightHandSideComputer::computeSourceTerm(const NumericalVector& numericalSolution)
{
    logger.assert(theta_ < 3.1415 / 2, "Angle must be in radians, not degrees!");
    const double h = numericalSolution(0);
    const double hu = numericalSolution(1);
    double u = 0;
    if (h > 1e-10)
    {
        u = hu/h;
    }
    double mu = computeFriction(numericalSolution);
    const int signU = (numericalSolution(1) > 0) ? 1 : -1;
    double sourceX = h * std::sin(theta_) - h * mu * signU * std::cos(theta_);
    logger(DEBUG, "Source: %", sourceX);
    return NumericalVector({0, sourceX});
}

NumericalVector SavageHutterRightHandSideComputer::computeNumericalSolution(const Base::Element *ptrElement, const Geometry::PointReference &pRef, const NumericalVector& solutionCoefficients)
{    
    logger.assert(1 == pRef.size(), "Empty reference point given.");
    const std::size_t numBasisFuns = ptrElement->getNrOfBasisFunctions();
    double h = 0;
    double hu = 0;
    for (std::size_t i = 0; i < numBasisFuns; ++i)
    {
        std::size_t iH = ptrElement->convertToSingleIndex(i, 0);
        h += solutionCoefficients(iH) * ptrElement->basisFunction(i, pRef);
        std::size_t iHu = ptrElement->convertToSingleIndex(i, 1);
        hu += solutionCoefficients(iHu) * ptrElement->basisFunction(i, pRef);
    }
    logger(DEBUG, "h: %, hu: %", h, hu);
    return NumericalVector({h,hu});
}

NumericalVector SavageHutterRightHandSideComputer::localLaxFriedrichsFlux(const NumericalVector& numericalSolutionLeft, const NumericalVector& numericalSolutionRight)
{
    double uLeft = 0;
    if (numericalSolutionLeft(0) > 1e-10)
    {
        uLeft = numericalSolutionLeft(1) / numericalSolutionLeft(0);
    }
    
    double uRight = 0;
    if (numericalSolutionRight(0) > 1e-10)
    {
        uRight = numericalSolutionRight(1) / numericalSolutionRight(0);
    }
    
    const double alpha = std::max(std::abs(uLeft) + std::sqrt(epsilon_ * numericalSolutionLeft(0)), 
                      std::abs(uRight) + std::sqrt(epsilon_ * numericalSolutionRight(0)));
        
    NumericalVector diffSolutions = numericalSolutionRight - numericalSolutionLeft;
    
    const NumericalVector numericalFlux = 0.5 * 
        (computePhysicalFlux(numericalSolutionLeft) + computePhysicalFlux(numericalSolutionRight)
            - alpha * (diffSolutions));
    
    return numericalFlux;
}

double SavageHutterRightHandSideComputer::computeFriction(const NumericalVector& numericalSolution)
{
    return std::tan(theta_);
}

