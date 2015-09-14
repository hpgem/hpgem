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
#include "Logger.h"
#include "Base/L2Norm.h"
#include "Geometry/Mappings/MappingReferenceToPhysical.h"
#include "HelperFunctions.h"

using LinearAlgebra::MiddleSizeVector;

/// \details The integrand for the reference element is the same as the physical element, but scaled with the reference-to-physical element scale, which is the determinant of the jacobian of the reference-to-physical element mapping.
MiddleSizeVector SavageHutterRightHandSideComputer::integrandRightHandSideOnElement
(Base::PhysicalElement<DIM>& element, const double &time, const MiddleSizeVector &solutionCoefficients)
{
    const std::size_t numBasisFuncs = element.getElement()->getNumberOfBasisFunctions();
    
    MiddleSizeVector& integrand = element.getResultVector(); //just to have the correct length    
    const PointPhysicalT& pPhys = element.getPointPhysical();
    const PointReferenceT& pRef = element.getPointReference();
    const MiddleSizeVector numericalSolution = Helpers::getSolution<DIM>(element.getElement(), solutionCoefficients, pRef, numOfVariables_);
    const MiddleSizeVector physicalFlux = computePhysicalFlux(numericalSolution);
    const MiddleSizeVector source = computeSourceTerm(numericalSolution, pPhys, time);
    
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
    
    return integrand;
}

/// \details The integrand for the reference face is the same as the physical face, but scaled with the reference-to-physical face scale. This face scale is absorbed in the normal vector, since it is relatively cheap to compute the normal vector with a length (L2-norm) equal to the reference-to-physical face scale.
MiddleSizeVector SavageHutterRightHandSideComputer::integrandRightHandSideOnRefFace
( Base::PhysicalFace<DIM>& face, const Base::Side &iSide, const MiddleSizeVector &solutionCoefficientsLeft, const MiddleSizeVector &solutionCoefficientsRight)
{
    double normal = face.getNormalVector()[0];
    const std::size_t numTestBasisFuncs = face.getFace()->getPtrElement(iSide)->getNumberOfBasisFunctions();

    //compute numerical solution at the left side and right side of this face
    const PointReferenceOnFaceT& pRef = face.getPointReference();
    const PointReferenceT& pRefL = face.getFace()->mapRefFaceToRefElemL(pRef);
    const PointReferenceT& pRefR = face.getFace()->mapRefFaceToRefElemR(pRef);
    MiddleSizeVector solutionLeft = Helpers::getSolution<DIM>(face.getFace()->getPtrElementLeft(), solutionCoefficientsLeft, pRefL, numOfVariables_);    
    MiddleSizeVector solutionRight = Helpers::getSolution<DIM>(face.getFace()->getPtrElementRight(), solutionCoefficientsRight, pRefR, numOfVariables_);
    
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
            integrand(iVarFun) = -flux(iVar) * face.basisFunction(iSide, iFun) * normal;
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
    const std::size_t numBasisFuncs = face.getFace()->getNumberOfBasisFunctions();
    
    const PointReferenceOnFaceT& pRef = face.getPointReference();
    //note that at the boundary, the element is the left element by definition
    const PointReferenceT& pRefL = face.getFace()->mapRefFaceToRefElemL(pRef);
    MiddleSizeVector solution = Helpers::getSolution<DIM>(face.getFace()->getPtrElementLeft(), solutionCoefficients, pRefL, numOfVariables_);
    
    MiddleSizeVector flux(2);
    double u = 0;
    if (solution(0) > minH_)
    {
        u = solution(1)/solution(0);
    }
    
    //outflow
    if (normal*solution(1) > 0)
    {
        const Base::Face *otherFace = face.getPhysicalElement(Base::Side::LEFT).getElement()->getFace(0);
        logger(DEBUG, "face before last: %, position %", otherFace->getID(), otherFace->referenceToPhysical(pRef));
        MiddleSizeVector otherSideSolution = 0.5*(otherFace->getPtrElementRight()->getSolution(0, otherFace->mapRefFaceToRefElemR(pRef)) + otherFace->getPtrElementLeft()->getSolution(0, otherFace->mapRefFaceToRefElemL(pRef)));
        const Base::Element *otherElement = otherFace->getPtrElementLeft();
        logger(DEBUG, "one back: %", otherElement->getID());
        otherFace = otherElement->getFace(0);
        MiddleSizeVector twoBackSolution = 0.5*(otherFace->getPtrElementRight()->getSolution(0, otherFace->mapRefFaceToRefElemR(pRef)) + otherFace->getPtrElementLeft()->getSolution(0, otherFace->mapRefFaceToRefElemL(pRef)));
        otherElement = otherFace->getPtrElementLeft();
        logger(DEBUG, "two back: %", otherElement->getID());
        otherFace = otherElement->getFace(0);
        MiddleSizeVector threeBackSolution = 0.5*(otherFace->getPtrElementRight()->getSolution(0, otherFace->mapRefFaceToRefElemR(pRef)) + otherFace->getPtrElementLeft()->getSolution(0, otherFace->mapRefFaceToRefElemL(pRef)));
        otherElement = otherFace->getPtrElementLeft();
        logger(DEBUG, "three back: %", otherElement->getID());
        otherFace = otherElement->getFace(0);
        MiddleSizeVector fourBackSolution = 0.5*(otherFace->getPtrElementRight()->getSolution(0, otherFace->mapRefFaceToRefElemR(pRef)) + otherFace->getPtrElementLeft()->getSolution(0, otherFace->mapRefFaceToRefElemL(pRef)));
        otherElement = otherFace->getPtrElementLeft();
        logger(DEBUG, "four back: %", otherElement->getID());
        otherFace = otherElement->getFace(0);
        MiddleSizeVector fiveBackSolution = 0.5*(otherFace->getPtrElementRight()->getSolution(0, otherFace->mapRefFaceToRefElemR(pRef)) + otherFace->getPtrElementLeft()->getSolution(0, otherFace->mapRefFaceToRefElemL(pRef)));
        otherElement = otherFace->getPtrElementLeft();
        logger(DEBUG, "five back: %", otherElement->getID());
        otherFace = otherElement->getFace(0);
        MiddleSizeVector sixBackSolution = 0.5*(otherFace->getPtrElementRight()->getSolution(0, otherFace->mapRefFaceToRefElemR(pRef)) + otherFace->getPtrElementLeft()->getSolution(0, otherFace->mapRefFaceToRefElemL(pRef)));
        
        MiddleSizeVector ghostSolution =  solution;
        
        //logger(INFO, "new flux %", flux);
        logger(DEBUG, "solution: %, otherSideSolution: % \n", solution, otherSideSolution);
        //subcritical
        if (solution[1]/solution[0] < std::sqrt(solution[0]*std::cos(chuteAngle_)*epsilon_))
        {
            double huNew = 1;
            double hNew = 1.14;//bisectionHFinder(huNew, solution);
            double uNew = solution[1] / solution[0] + 2*std::sqrt(epsilon_*std::cos(chuteAngle_))*(std::sqrt(solution[0]) - std::sqrt(hNew));
            auto stateNew = MiddleSizeVector({hNew, hNew * uNew});
            flux = localLaxFriedrichsFlux(solution, stateNew);
            
            logger(INFO, "reached subcritical! old state: % new state: %",solution, stateNew);
        }
        else
        {            
            flux = localLaxFriedrichsFlux(solution, solution);
        }
        //flux = computePhysicalFlux(ghostSolution);
        //logger(INFO, "old flux %", flux);
    }
    else //inflow
    {
        //subcritical
        if (solution[1]/solution[0] < std::sqrt(solution[0]*std::cos(chuteAngle_)*epsilon_))
        {
            double hNew = inflowBC_[0];
            double uNew = solution[1] / solution[0] - 2*std::sqrt(epsilon_*std::cos(chuteAngle_))*(std::sqrt(solution[0]) - std::sqrt(hNew));
            auto stateNew = MiddleSizeVector({hNew, hNew * uNew});
            flux = localLaxFriedrichsFlux(solution, stateNew);
        }
        else
        {            
            flux = localLaxFriedrichsFlux(inflowBC_, inflowBC_);
        }
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
    if (h > minH_)
    {
        u = hu/h;
    }
    MiddleSizeVector flux(2);
    flux(0) = hu;
    flux(1) = hu * u * alpha_ + epsilon_/2 * std::cos(chuteAngle_) * h * h;
    logger(DEBUG, "flux values: %, %", flux(0), flux(1));
    return flux;
}

MiddleSizeVector SavageHutterRightHandSideComputer::computeSourceTerm(const MiddleSizeVector& numericalSolution, const PointPhysicalT& pPhys, const double time)
{
    logger.assert(chuteAngle_ < M_PI, "Angle must be in radians, not degrees!");
    const double h = numericalSolution(0);
    const double hu = numericalSolution(1);
    double u = 0;
    if (h > minH_)
    {
        u = hu/h;
    }
    double mu = computeFrictionCoulomb(numericalSolution);
    const int signU = Helpers::sign(u);
    double sourceX = h * std::sin(chuteAngle_) - h * mu * signU * std::cos(chuteAngle_);
    logger(DEBUG, "Source: %, h: %", sourceX, h);
    return MiddleSizeVector({0, sourceX});
}

MiddleSizeVector SavageHutterRightHandSideComputer::localLaxFriedrichsFlux(const MiddleSizeVector& numericalSolutionLeft, const MiddleSizeVector& numericalSolutionRight)
{
    double uLeft = 0;
    if (numericalSolutionLeft(0) > minH_)
    {
        uLeft = numericalSolutionLeft(1) / numericalSolutionLeft(0);
    }
    
    double uRight = 0;
    if (numericalSolutionRight(0) > minH_)
    {
        uRight = numericalSolutionRight(1) / numericalSolutionRight(0);
    }
    
    const double alpha = std::max(std::abs(uLeft) + std::sqrt(epsilon_ * std::max(0.,numericalSolutionLeft(0))), 
                      std::abs(uRight) + std::sqrt(epsilon_ * std::max(0.,numericalSolutionRight(0))));
    
    logger(DEBUG, "alpha: %", alpha);
        
    logger(DEBUG, "physical fluxes: %, %", computePhysicalFlux(numericalSolutionLeft), computePhysicalFlux(numericalSolutionRight));
    MiddleSizeVector diffSolutions = numericalSolutionRight - numericalSolutionLeft;
    const MiddleSizeVector numericalFlux = 0.5 * 
        (computePhysicalFlux(numericalSolutionLeft) + computePhysicalFlux(numericalSolutionRight)
            - alpha * (diffSolutions));
    
    return numericalFlux;
}

double SavageHutterRightHandSideComputer::computeFrictionCoulomb(const MiddleSizeVector& numericalSolution)
{
    return std::tan(60./180*M_PI);
}

/// Compute friction as described in Weinhart (2012), eq (50) with lambda = 1, d = 1
/// Notice that the minus sign for gamma in eq (50) is wrong.
double SavageHutterRightHandSideComputer::computeFriction(const MiddleSizeVector& numericalSolution)
{
    const double delta1 = 17.561 / 180 * M_PI ;
    const double delta2 = 32.257 / 180 * M_PI;
    const double A = 3.836;
    const double beta = 0.191;
    const double gamma = -.045;
    const double d = 1.;
    const double h = numericalSolution[0];
    if (h < minH_)
        return std::tan(delta1);
    const double u  = numericalSolution[1] / h;
    const double F = u /std::sqrt(epsilon_*std::cos(chuteAngle_) * h);
    return std::tan(delta1) + (std::tan(delta2) - std::tan(delta1))/(beta*h/(A*d*(F - gamma)) + 1);
}

/// Compute friction as in Anthony's draft, eq (3.9)
double SavageHutterRightHandSideComputer::computeFrictionExponential(const MiddleSizeVector& numericalSolution)
{
    const double delta1 = 27. / 180 * M_PI ;
    const double delta2 = 37. / 180 * M_PI;
    const double beta = 0.136;
    const double L = 2.;
    const double h = numericalSolution[0];
    if (h < minH_)
        return std::tan(delta1);
    const double u  = numericalSolution[1] / h;
    if (std::abs(u) < 1e-16)
        return std::tan(delta1);
    return std::tan(delta1) + (std::tan(delta2) - std::tan(delta1))*std::exp(-beta*std::pow(epsilon_*h, 1.5)/(L * std::abs(u)));
}
