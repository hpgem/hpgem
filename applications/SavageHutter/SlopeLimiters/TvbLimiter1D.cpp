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

#include "TvbLimiter1D.h"
#include "../GlobalConstants.h"

void TvbLimiter1D::limitSlope(Base::Element* element) 
{
    logger.assert(1 == DIM, "Slope limiter for the wrong dimension");
    for (std::size_t iVar = 0; iVar < numOfVariables_; ++iVar)
    {
        if (!hasSmallSlope(element, iVar))
        {
            limitWithMinMod(element, iVar);
        }
    }
}


void TvbLimiter1D::limitWithMinMod(Base::Element* element, const std::size_t iVar)
{    
    
    const LinearAlgebra::MiddleSizeVector &solutionCoefficients = element->getTimeLevelDataVector(0);
    const PointReferenceT &pRefL = element->getReferenceGeometry()->getReferenceNodeCoordinate(0); 
    const PointReferenceT &pRefR = element->getReferenceGeometry()->getReferenceNodeCoordinate(1);
    logger.assert(std::abs(-1  - pRefL[0]) < 1e-10, "xi_L != -1");    
    logger.assert(std::abs(1 - pRefR[0]) < 1e-10, "xi_R != 1");
    
    //this does not work for first boundary, probably also not for triangular mesh.
    //maybe write getNeighbour(face)?
    const Base::Element *elemL;
    const Base::Element *elemR;
    if (element->getFace(0)->isInternal())
        elemL  = element->getFace(0)->getPtrElementLeft();
    else
        return; //for now, just don't use the limiter here...
    if (element->getFace(1)->isInternal())
        elemR = element->getFace(1)->getPtrElementRight();
    else
        return; //for now, just don't use the limiter here...
    
    logger.assert(elemR->getID() == element->getID() + 1 && element->getID() == elemL->getID() + 1, "elements not in correct order");

    const double uPlus = element->getSolution(0, pRefR)[iVar];
    const double uMinus = element->getSolution(0, pRefL)[iVar];
    
    const double u0 = Helpers::computeAverageOfSolution<1>(element, solutionCoefficients, elementIntegrator_)(iVar);
    const double uElemR = Helpers::computeAverageOfSolution<1>(const_cast<Base::Element*>(elemR), solutionCoefficients, elementIntegrator_)(iVar);
    const double uElemL = Helpers::computeAverageOfSolution<1>(const_cast<Base::Element*>(elemL), solutionCoefficients, elementIntegrator_)(iVar);
    logger(DEBUG, "coefficients: %", element->getTimeLevelData(0, iVar));
    logger(DEBUG, "uPlus: %, basis function vals: %, %", uPlus , element->basisFunction(0,pRefR), element->basisFunction(1,pRefR));
    logger(DEBUG, "uMinus: %, basis function vals: %, %", uMinus, element->basisFunction(0,pRefL), element->basisFunction(1,pRefL));
    logger(INFO, "Element %: %, %, %",element->getID(), (uPlus - uMinus), uElemR - u0, u0 - uElemL);
    double slope =  0;
    if (Helpers::sign(uElemR - u0) == Helpers::sign(u0 - uElemL) && Helpers::sign(uElemR - u0) == Helpers::sign(uPlus - uMinus))
    {
        slope = Helpers::sign(uElemR - u0) * std::min(std::abs(uPlus - uMinus), std::min(std::abs(uElemR - u0), std::abs(u0 - uElemL)));
    }
    
    if (std::abs(std::abs(slope) - std::abs(uPlus - uMinus)) > 1e-16 )
    {
        //replace coefficients with "u0 + slope/2 * xi" coefficients
        std::function<double(const PointReferenceT&)> newFun = [=] (const PointReferenceT& pRef) {return u0 + slope*pRef[0];};
        LinearAlgebra::MiddleSizeVector newCoeffs = Helpers::projectOnBasisFuns<1>(element, newFun, elementIntegrator_);
        element->setTimeLevelData(0, iVar, newCoeffs);
        logger(INFO, "element % is being limited.", element->getID());
    }
}


bool TvbLimiter1D::hasSmallSlope(const Base::Element* element, const std::size_t iVar)
{
    const PointReferenceT &pRefL = element->getReferenceGeometry()->getReferenceNodeCoordinate(0); 
    const PointReferenceT &pRefR = element->getReferenceGeometry()->getReferenceNodeCoordinate(1);
    
    const double uPlus = element->getSolution(0, pRefR)[iVar];
    const double uMinus = element->getSolution(0, pRefL)[iVar];
    const double M = .1;
    return (std::abs(uPlus - uMinus) < M * (2*element->calcJacobian(pRefL).determinant())* (2*element->calcJacobian(pRefL).determinant()));
}