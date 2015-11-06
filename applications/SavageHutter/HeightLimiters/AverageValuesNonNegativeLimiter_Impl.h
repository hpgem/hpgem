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

#include "../HelperFunctions.h"

///\details Limit the solution by changing the solution coefficients that are given to this function.
template <std::size_t DIM>
void AverageValuesNonNegativeLimiter<DIM>::limit(Base::Element *element, LinearAlgebra::MiddleSizeVector &solutionCoefficients)
{
    const double minimumHeight = getMinimumHeight(element);
    if (minimumHeight >= minH_)
        return;

    logger(DEBUG, "solution coefficients before limiting element %: %", element->getID(), solutionCoefficients);
    const LinearAlgebra::MiddleSizeVector averages = Helpers::computeAverageOfSolution(element, solutionCoefficients, elementIntegrator_);
    const double averageHeight = averages(0);
    const double averageDischargeX = averages(1);
    double averageDischargeY;
    if (DIM == 2)    
        averageDischargeY = averages(2);
    logger(DEBUG, "average height: %", averageHeight);
    LinearAlgebra::MiddleSizeVector heightCoefficients;
    LinearAlgebra::MiddleSizeVector dischargeCoefficientsX;
    LinearAlgebra::MiddleSizeVector dischargeCoefficientsY;

    //compute the new solution coefficients
    heightCoefficients =
        Helpers::projectOnBasisFuns<DIM>(element, [ = ](const PointReferenceT & pRef){return averageHeight;}, elementIntegrator_);
    dischargeCoefficientsX = Helpers::projectOnBasisFuns<DIM>(element, [ = ](const PointReferenceT & pRef){return averageDischargeX;}, elementIntegrator_);
    if (DIM == 2) 
        dischargeCoefficientsY = Helpers::projectOnBasisFuns<DIM>(element, [ = ](const PointReferenceT & pRef){return averageDischargeY;}, elementIntegrator_);

    // set the solution coefficients
    for (std::size_t iFun = 0; iFun < element->getNumberOfBasisFunctions(); ++iFun)
    {
        std::size_t iVF = element->convertToSingleIndex(iFun, 0);
        solutionCoefficients[iVF] = heightCoefficients[iFun];
        iVF = element->convertToSingleIndex(iFun, 1);
        solutionCoefficients[iVF] = dischargeCoefficientsX[iFun];
        if (DIM == 2)
        {
            iVF = element->convertToSingleIndex(iFun, 2);
            solutionCoefficients[iVF] = dischargeCoefficientsY[iFun];
        }
    }
    logger(DEBUG, "solution coefficients after limiting element %: %", element->getID(), solutionCoefficients);
}

template <std::size_t DIM>
double AverageValuesNonNegativeLimiter<DIM>::getMinimumHeight(const Base::Element* element)
{
    const PointReferenceT &pRefL = element->getReferenceGeometry()->getReferenceNodeCoordinate(0);
    const LinearAlgebra::MiddleSizeVector &solutionCoefficients = element->getTimeIntegrationVector(0);
    const double solutionLeft = element->getSolution(0, pRefL)(0);//Helpers::getSolution<DIM>(element, solutionCoefficients, pRefL, numOfVariables)(0);    
    double minimum = solutionLeft;
    for (std::size_t iPoint = 1; iPoint < element->getReferenceGeometry()->getNumberOfNodes(); ++iPoint)
    {
        const PointReferenceT &pRef = element->getReferenceGeometry()->getReferenceNodeCoordinate(iPoint);
        minimum = std::min(minimum, element->getSolution(0, pRef)(0));
    }
    for (std::size_t p = 0; p < element->getGaussQuadratureRule()->getNumberOfPoints(); ++p)
    {
        const PointReferenceT& pRef = element->getGaussQuadratureRule()->getPoint(p);
        minimum = std::min(minimum, element->getSolution(0, pRef)(0));
    }
    logger(DEBUG, "Minimum in element %: %", element->getID(), minimum);
    return minimum;
}

