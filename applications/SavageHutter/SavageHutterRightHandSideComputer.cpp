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
NumericalVector SavageHutterRightHandSideComputer::integrandRightHandSideOnRefElement(const Base::Element *ptrElement, const double &time, const Geometry::PointReference &pRef, const NumericalVector &solutionCoefficients)
{
    const std::size_t numOfBasisFunctions = ptrElement->getNrOfBasisFunctions();
    
    NumericalVector integrand(numOfVariables_ * numOfBasisFunctions);
    const NumericalVector numericalSolution = computeNumericalSolution(ptrElement, pRef, solutionCoefficients);
    const NumericalVector physicalFlux = computePhysicalFlux(numericalSolution);
    
    // Compute integrand on the physical element.
    std::size_t iVB; // Index for both basis function and variable
    for (std::size_t iB = 0; iB < numOfBasisFunctions; iB++) // Index for basis function
    {
        iVB = ptrElement->convertToSingleIndex(iB, 0);
        integrand(iVB) = physicalFlux(0) * ptrElement->basisFunctionDeriv(iB, 0, pRef);

        iVB = ptrElement->convertToSingleIndex(iB, 1);
        integrand(iVB) = physicalFlux(1) * ptrElement->basisFunctionDeriv(iB, 0, pRef);
    }
    Geometry::Jacobian jac = ptrElement->calcJacobian(pRef);
    integrand *= jac.determinant();
    
    logger(DEBUG, "Integrand on element: %", integrand);
    return integrand;
}

/// \details The integrand for the reference face is the same as the physical face, but scaled with the reference-to-physical face scale. This face scale is absorbed in the normal vector, since it is relatively cheap to compute the normal vector with a length (L2-norm) equal to the reference-to-physical face scale.
NumericalVector SavageHutterRightHandSideComputer::computeRightHandSideOnRefFace
(
 const Base::Face *ptrFace,
 const Base::Side &iSide,
 const NumericalVector &solutionCoefficientsLeft,
 const NumericalVector &solutionCoefficientsRight
 )
{
    Geometry::PointReference pRef(1);
    double normal = ptrFace->getNormalVector(Geometry::PointReference(0))(0);
    if(iSide == Base::Side::LEFT) //right face of element
    {
        pRef.setCoordinate(0, 1);
    }
    else
    {
        pRef.setCoordinate(0, -1);
    }
    const std::size_t numOfTestBasisFunctions = ptrFace->getPtrElement(iSide)->getNrOfBasisFunctions(); // Get the number of test basis functions on a given side, iSide
    NumericalVector solutionLeft = computeNumericalSolution(ptrFace->getPtrElementLeft(), pRef, solutionCoefficientsLeft);
    NumericalVector solutionRight = computeNumericalSolution(ptrFace->getPtrElementRight(), pRef, solutionCoefficientsRight);
    NumericalVector flux = localLaxFriedrichsFlux(solutionLeft, solutionRight, normal);
    
    if (iSide == Base::Side::LEFT)
    {
        normal *= -1;
    }
    NumericalVector integrand(numOfVariables_ * numOfTestBasisFunctions); // Integrand value based on n number of testbasisfunctions from element corresponding to side iSide

    for (std::size_t iFun = 0; iFun < numOfTestBasisFunctions; ++iFun)
    {
        for (std::size_t iVar = 0; iVar < numOfVariables_; ++iVar)
        {
            std::size_t iVarFun = ptrFace->getPtrElement(iSide)->convertToSingleIndex(iFun, iVar);
            integrand(iVarFun) = flux(iVar) * ptrFace->basisFunction(iSide, iFun, Geometry::PointReference(0)) * normal;            
        }
    }
    logger(DEBUG, "Values of rhs on face % : %, %", ptrFace->getID(), integrand(0), integrand(1));
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
    flux(1) = hu * u + epsilon_/2 * h * h;
    logger(DEBUG, "flux values: %, %", flux(0), flux(1));
    return flux;
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
    return NumericalVector({h,hu});
}

NumericalVector SavageHutterRightHandSideComputer::localLaxFriedrichsFlux(const NumericalVector& numericalSolutionLeft, const NumericalVector& numericalSolutionRight, double normal)
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
    
    double alpha = std::max(std::abs(uLeft) + std::sqrt(epsilon_ * numericalSolutionLeft(0)), 
                      std::abs(uRight) + std::sqrt(epsilon_ * numericalSolutionRight(0)));
        
    NumericalVector diffSolutions(2);
    if (normal == 1)
    {
        diffSolutions = numericalSolutionRight - numericalSolutionLeft;
    }
    else
    {
        diffSolutions = numericalSolutionLeft - numericalSolutionRight;
    }
    NumericalVector numericalFlux = 0.5 * 
        (computePhysicalFlux(numericalSolutionLeft) + computePhysicalFlux(numericalSolutionRight)
            - alpha * (diffSolutions));
    logger(DEBUG, "%, %, %, %, %", computePhysicalFlux(numericalSolutionLeft), computePhysicalFlux(numericalSolutionRight), alpha, numericalSolutionRight, numericalSolutionLeft);
    return numericalFlux;
}