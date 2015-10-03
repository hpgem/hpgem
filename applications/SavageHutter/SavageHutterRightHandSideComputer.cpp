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

#include "SavageHutter.h"
#include "SavageHutterRightHandSideComputer.h"
#include "Logger.h"
#include "Base/L2Norm.h"
#include "Geometry/Mappings/MappingReferenceToPhysical.h"

using LinearAlgebra::MiddleSizeVector;

/// \details The integrand for the reference element is the same as the physical element, but scaled with the reference-to-physical element scale, which is the determinant of the jacobian of the reference-to-physical element mapping.
MiddleSizeVector SavageHutterRightHandSideComputer::integrandRightHandSideOnElement
(Base::PhysicalElement<DIM>& element, const double &time, const MiddleSizeVector &solutionCoefficients)
{
    const std::size_t numBasisFuncs = element.getElement()->getNrOfBasisFunctions();
    
    MiddleSizeVector& integrand = element.getResultVector(); //just to have the correct length
    const MiddleSizeVector numericalSolution = element.getSolution();
    const MiddleSizeVector physicalFlux = computePhysicalFlux(numericalSolution);
    const PointPhysicalT& pPhys = element.getPointPhysical();
    const MiddleSizeVector source = computeSourceTerm(numericalSolution, pPhys, time);
    logger.assert(Base::L2Norm(source) < 1e-10, "Source non-zero: %", source);
    
    // Compute integrand on the physical element.
    std::size_t iVB; // Index for both basis function and variable
    for (std::size_t iB = 0; iB < numBasisFuncs; iB++) // Index for basis function
    {
        for (std::size_t iV = 0; iV < numOfVariables_; ++iV)
        {
            iVB = element.getElement()->convertToSingleIndex(iB, iV);
            integrand(iVB) = physicalFlux(iV) * element.basisFunctionDeriv(iB)(0);
            integrand(iVB) += source(iV) * element.basisFunction(iB);
        }
    }
    
    logger(DEBUG, "Integrand on element: %", integrand);
    return integrand;
}

/// \details The integrand for the reference face is the same as the physical face, but scaled with the reference-to-physical face scale. This face scale is absorbed in the normal vector, since it is relatively cheap to compute the normal vector with a length (L2-norm) equal to the reference-to-physical face scale.
MiddleSizeVector SavageHutterRightHandSideComputer::integrandRightHandSideOnRefFace
( Base::PhysicalFace<DIM>& face, const Base::Side &iSide, const MiddleSizeVector &solutionCoefficientsLeft, const MiddleSizeVector &solutionCoefficientsRight)
{
    double normal = face.getNormalVector()[0];
    const std::size_t numTestBasisFuncs = face.getFace()->getPtrElement(iSide)->getNrOfBasisFunctions();
    const std::size_t numBasisFuncsLeft = face.getFace()->getPtrElement(Base::Side::LEFT)->getNrOfBasisFunctions();
    const std::size_t numBasisFuncsRight = face.getFace()->getPtrElement(Base::Side::RIGHT)->getNrOfBasisFunctions();

    MiddleSizeVector solutionLeft = face.getSolution(Base::Side::LEFT);
    MiddleSizeVector solutionRight = face.getSolution(Base::Side::RIGHT);
    
    logger(DEBUG, "face: %, uL: %, uR:%", face.getFace()->getID(), solutionLeft, solutionRight);
    MiddleSizeVector flux(2);
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
    MiddleSizeVector& integrand = face.getResultVector(iSide); // Integrand value based on n number of testbasisfunctions from element corresponding to side iSide

    for (std::size_t iFun = 0; iFun < numTestBasisFuncs; ++iFun)
    {
        for (std::size_t iVar = 0; iVar < numOfVariables_; ++iVar)
        {
            std::size_t iVarFun = face.getFace()->getPtrElement(iSide)->convertToSingleIndex(iFun, iVar);
            const PointReferenceOnFaceT& center = face.getFace()->getReferenceGeometry()->getCenter();
            integrand(iVarFun) = -flux(iVar) * face.getFace()->basisFunction(iSide, iFun, center) * normal;
        }
    }
    
    return integrand;
}

MiddleSizeVector SavageHutterRightHandSideComputer::integrandRightHandSideOnRefFace
    (
     Base::PhysicalFace<DIM>& face,
     const MiddleSizeVector &solutionCoefficients
     )
{
    double normal = face.getNormalVector()[0];
    const std::size_t numBasisFuncs = face.getFace()->getNrOfBasisFunctions();
    MiddleSizeVector solution = face.getSolution(Base::Side::LEFT); //note that at the boundary, the element is the left element by definition
    
    MiddleSizeVector flux(2);
    double u = 0;
    if (solution(0) > 1e-10)
    {
        u = solution(1)/solution(0);
    }
    
    //outflow
    if (normal > 0)
    {
        flux = localLaxFriedrichsFlux(solution, solution);
    }
    else //inflow
    {
        flux = localLaxFriedrichsFlux(getInflowBC(), solution);
    }
    
    MiddleSizeVector integrand(numOfVariables_ * numBasisFuncs);

    for (std::size_t iFun = 0; iFun < numBasisFuncs; ++iFun)
    {
        for (std::size_t iVar = 0; iVar < numOfVariables_; ++iVar)
        {
            std::size_t iVarFun = face.getFace()->getPtrElementLeft()->convertToSingleIndex(iFun, iVar);
            integrand(iVarFun) = -flux(iVar) * face.basisFunction(iFun) * normal;
        }
    }
    
    return integrand;
}

MiddleSizeVector SavageHutterRightHandSideComputer::computePhysicalFlux(const MiddleSizeVector &numericalSolution)
{    
    const double h = numericalSolution(0);
    logger.assert(h > -1e-16, "Negative height (%)", h);
    const double hu = numericalSolution(1);
    double u = 0;
    if (h > 1e-10)
    {
        u = hu/h;
    }
    MiddleSizeVector flux(2);
    flux(0) = hu;
    flux(1) = hu * u + epsilon_/2 * std::cos(theta_) * h * h;
    logger(DEBUG, "flux values: %, %", flux(0), flux(1));
    return flux;
}

MiddleSizeVector SavageHutterRightHandSideComputer::computeSourceTerm(const MiddleSizeVector& numericalSolution, const PointPhysicalT& pPhys, const double time)
{
    logger.assert(theta_ < M_PI / 2, "Angle must be in radians, not degrees!");
    const double h = numericalSolution(0);
    const double hu = numericalSolution(1);
    double u = 0;
    if (h > 1e-10)
    {
        u = hu/h;
    }
    double mu = computeFriction(numericalSolution);
    const int signU = (u > -1e-16) ? 1 : -1;
    double sourceX = h * std::sin(theta_) - h * mu * signU * std::cos(theta_);
    logger(DEBUG, "Source: %, h: %", sourceX, h);
    return MiddleSizeVector({0, sourceX});
}

MiddleSizeVector SavageHutterRightHandSideComputer::localLaxFriedrichsFlux(const MiddleSizeVector& numericalSolutionLeft, const MiddleSizeVector& numericalSolutionRight)
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
    
    const double alpha = std::max(std::abs(uLeft) + std::sqrt(epsilon_ * std::max(0.,numericalSolutionLeft(0))), 
                      std::abs(uRight) + std::sqrt(epsilon_ * std::max(0.,numericalSolutionRight(0))));
    
    logger(DEBUG, "alpha: %", alpha);
        
    MiddleSizeVector diffSolutions = numericalSolutionRight - numericalSolutionLeft;
    
    const MiddleSizeVector numericalFlux = 0.5 * 
        (computePhysicalFlux(numericalSolutionLeft) + computePhysicalFlux(numericalSolutionRight)
            - alpha * (diffSolutions));
    
    return numericalFlux;
}

double SavageHutterRightHandSideComputer::computeFriction(const MiddleSizeVector& numericalSolution)
{
    return std::tan(theta_);
}

LinearAlgebra::MiddleSizeVector SavageHutterRightHandSideComputer::getInflowBC()
{
    return LinearAlgebra::MiddleSizeVector({.1, 0});
}
