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

#include "AverageValuesNonNegativeLimiter.h"
#include "../HelperFunctions.h"
#include "../GlobalConstants.h"
void AverageValuesNonNegativeLimiter::limit(Base::Element *element, LinearAlgebra::MiddleSizeVector &solutionCoefficients)
{
    const double minimumHeight = getMinimumHeight(element);
    if (minimumHeight >= minH_)
        return;

    logger(DEBUG, "solution coefficients before limiting element %: %", element->getID(), solutionCoefficients);
    const double averageHeight = Helpers::computeAverageOfSolution<DIM>(element, solutionCoefficients, elementIntegrator_)(0);
    const double averageDischarge = Helpers::computeAverageOfSolution<DIM>(element, solutionCoefficients, elementIntegrator_)(1);
    logger(DEBUG, "average height: %", averageHeight);
    LinearAlgebra::MiddleSizeVector heightCoefficients;
    LinearAlgebra::MiddleSizeVector dischargeCoefficients;

    heightCoefficients =
        Helpers::projectOnBasisFuns<DIM>(element, [ = ](const PointReferenceT & pRef){return averageHeight;}, elementIntegrator_);
    if (false && averageHeight < minH_)
    {
        dischargeCoefficients = Helpers::projectOnBasisFuns<DIM>(element, [ = ](const PointReferenceT & pRef){return 0;}, elementIntegrator_);
    }
    else
    {
        dischargeCoefficients = Helpers::projectOnBasisFuns<DIM>(element, [ = ](const PointReferenceT & pRef){return averageDischarge;}, elementIntegrator_);
    }

    for (std::size_t iFun = 0; iFun < element->getNrOfBasisFunctions(); ++iFun)
    {
        std::size_t iVF = element->convertToSingleIndex(iFun, 0);
        solutionCoefficients[iVF] = heightCoefficients[iFun];
        iVF = element->convertToSingleIndex(iFun, 1);
        solutionCoefficients[iVF] = dischargeCoefficients[iFun];
    }
    logger(DEBUG, "solution coefficients after limiting element %: %", element->getID(), solutionCoefficients);
}
double AverageValuesNonNegativeLimiter::getMinimumHeight(const Base::Element* element)
{
    const PointReferenceT &pRefL = element->getReferenceGeometry()->getReferenceNodeCoordinate(0);
    const PointReferenceT &pRefR = element->getReferenceGeometry()->getReferenceNodeCoordinate(1);
    const LinearAlgebra::MiddleSizeVector &solutionCoefficients = element->getTimeLevelDataVector(0);
    const std::size_t numOfVariables = element->getNrOfUnknowns();
    const double solutionLeft = Helpers::getSolution<DIM>(element, solutionCoefficients, pRefL, numOfVariables)(0);
    const double solutionRight = Helpers::getSolution<DIM>(element, solutionCoefficients, pRefR, numOfVariables)(0);
    double minimum = std::min(solutionLeft, solutionRight);
    for (std::size_t p = 0; p < element->getGaussQuadratureRule()->nrOfPoints(); ++p)
    {
        const PointReferenceT& pRef = element->getGaussQuadratureRule()->getPoint(p);
        minimum = std::min(minimum, Helpers::getSolution<DIM>(element, solutionCoefficients, pRef, numOfVariables)(0));
    }
    logger(DEBUG, "Minimum in element %: %", element->getID(), minimum);
    return minimum;
}

