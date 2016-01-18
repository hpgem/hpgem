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

#include "BidispersedLimiter1D.h"
#include "../HelperFunctions.h"
#include "PositiveLayerLimiter.h"

void BidispersedLimiter::limit(Base::Element* element, LinearAlgebra::MiddleSizeVector& solutionCoefficients)
{
    limitHeight(element, solutionCoefficients);
    //limitSmallHeight(element, solutionCoefficients);
    
    const LinearAlgebra::MiddleSizeVector averages = Helpers::computeAverageOfSolution(element, solutionCoefficients, elementIntegrator_);
    if (averages(0) + minH_ < averages(2))
    {
        setSmallHeigthToHeight(element, solutionCoefficients);
    }
    else if (averages(2) < -minH_)
    {
        setSmallHeightToZero(element, solutionCoefficients);
    }
    else if (smallHeightDoesNotFitBetweenHeightAndZero(element, solutionCoefficients))
    {
        squeezeSmallHeight(element, solutionCoefficients);
    }
    
}

void BidispersedLimiter::limitHeight(Base::Element* element, LinearAlgebra::MiddleSizeVector& solutionCoefficients)
{
    PositiveLayerLimiter<1> heightLimiter(minH_);
    heightLimiter.limit(element, solutionCoefficients);
}

void BidispersedLimiter::setSmallHeigthToHeight(Base::Element* element, LinearAlgebra::MiddleSizeVector& solutionCoefficients)
{
    for (std::size_t i = 0; i < element->getNumberOfBasisFunctions(); ++i)
    {
        const std::size_t iVFHeight = element->convertToSingleIndex(i, 0);
        const std::size_t iVFSmallHeight = element->convertToSingleIndex(i,2);
        solutionCoefficients[iVFSmallHeight] = solutionCoefficients[iVFHeight];
    }
}

void BidispersedLimiter::setSmallHeightToZero(Base::Element* element, LinearAlgebra::MiddleSizeVector& solutionCoefficients)
{
    for (std::size_t i = 0; i < element->getNumberOfBasisFunctions(); ++i)
    {
        const std::size_t iVFSmallHeight = element->convertToSingleIndex(i,2);
        solutionCoefficients[iVFSmallHeight] = 0;
    }
}

void BidispersedLimiter::squeezeSmallHeight(Base::Element* element, LinearAlgebra::MiddleSizeVector& solutionCoefficients)
{
    const double squeezeFactor = findSqueezeFactorSmallHeight(element, solutionCoefficients);    
    const LinearAlgebra::MiddleSizeVector averages = Helpers::computeAverageOfSolution(element, solutionCoefficients, elementIntegrator_);
    std::function<double(const PointReferenceT&) > smallHeightFun = [ = ] (const PointReferenceT & pRef)
        {
            return averages(2);
        };
    const LinearAlgebra::MiddleSizeVector smallHeightCoefficients = Helpers::projectOnBasisFuns<1>(element, smallHeightFun, elementIntegrator_);
    for (std::size_t iFun = 0; iFun < element->getNumberOfBasisFunctions(); ++iFun)
    {
        std::size_t iVF = element->convertToSingleIndex(iFun, 2);
        solutionCoefficients[iVF] = smallHeightCoefficients[iFun];
    }
}

bool BidispersedLimiter::smallHeightDoesNotFitBetweenHeightAndZero(Base::Element *element, LinearAlgebra::MiddleSizeVector &solutionCoefficients)
{
    //for both end points of the element, check if it is between 0 and h
    for (std::size_t iPoint = 0; iPoint < element->getReferenceGeometry()->getNumberOfNodes(); ++iPoint)
    {
        const PointReferenceT &pRef = element->getReferenceGeometry()->getReferenceNodeCoordinate(iPoint);
        const LinearAlgebra::MiddleSizeVector solution = Helpers::getSolution(element, solutionCoefficients, pRef);
        if (solution(0) + minH_ < solution(2) || solution(2) < -minH_)
            return true;
    }
    //for all quadrature points, check if it is between 0 and h
    for (std::size_t p = 0; p < element->getGaussQuadratureRule()->getNumberOfPoints(); ++p)
    {
        const PointReferenceT& pRef = element->getGaussQuadratureRule()->getPoint(p);
        const LinearAlgebra::MiddleSizeVector solution = Helpers::getSolution(element, solutionCoefficients, pRef);
        if (solution(0) + minH_ < solution(2) || solution(2) < -minH_)
            return true;
    }
    return false;
}

double BidispersedLimiter::findSqueezeFactorSmallHeight(Base::Element* element, LinearAlgebra::MiddleSizeVector& solutionCoefficients)
{
    //for now, just set the smallHeight to its average. This has to be solved in the future.
    //we want to fit the smallHeight between 0 and height, but I'm currently not sure how to do that or if it's even possible in the general case.
    return 0;
}

void BidispersedLimiter::limitSmallHeight(Base::Element* element, LinearAlgebra::MiddleSizeVector& solutionCoefficients)
{
    const double averageSmallHeight = Helpers::computeAverageOfSolution<1>(element, solutionCoefficients, elementIntegrator_)(2);
    LinearAlgebra::MiddleSizeVector heightCoefficients;
    if (averageSmallHeight < -minH_)
    {
        setSmallHeightToZero(element, solutionCoefficients);
    }
    else
    {
        const PointReferenceT &pRefL = element->getReferenceGeometry()->getReferenceNodeCoordinate(0);
        const double solutionLeft = element->getSolution(0, pRefL)(2);
        const PointReferenceT &pRefR = element->getReferenceGeometry()->getReferenceNodeCoordinate(1);
        const double solutionRight = element->getSolution(0, pRefR)(2);

        double slopeSmallHeight = std::min((solutionLeft + solutionRight) / 2, averageSmallHeight);

        if (solutionRight < solutionLeft)
        {
            slopeSmallHeight *= -1;
        }
        std::function<double(const PointReferenceT&) > heightFun = [ = ] (const PointReferenceT & pRef){return averageSmallHeight + slopeSmallHeight * pRef[0];};
        heightCoefficients = Helpers::projectOnBasisFuns<1>(element, heightFun, elementIntegrator_);
        logger(DEBUG, "averages: %", Helpers::computeAverageOfSolution<1>(element, solutionCoefficients, elementIntegrator_));
        for (std::size_t iFun = 0; iFun < element->getNumberOfBasisFunctions(); ++iFun)
        {
            std::size_t iVF = element->convertToSingleIndex(iFun, 2);
            solutionCoefficients[iVF] = heightCoefficients[iFun];
        }
    }
}