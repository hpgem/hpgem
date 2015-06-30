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

#include "PositiveLayerLimiter.h"
#include "HelperFunctions.h"

double PositiveLayerLimiter::getMinimumHeight(const Base::Element* element)
{
    const PointReferenceT &pRefL = element->getReferenceGeometry()->getReferenceNodeCoordinate(0); 
    const PointReferenceT &pRefR = element->getReferenceGeometry()->getReferenceNodeCoordinate(1);
    
    double minimum = std::min(element->getSolution(0,pRefL)(0), element->getSolution(0,pRefR)(0));
    for (std::size_t p = 0; p < element->getGaussQuadratureRule()->nrOfPoints(); ++p)
    {
        const PointReferenceT& pRef = element->getGaussQuadratureRule()->getPoint(p);
        minimum = std::min(minimum, element->getSolution(0,pRef)(0));
    }
    logger(DEBUG, "Minimum in element %: %", element->getID(), minimum);
    return minimum;
}

void PositiveLayerLimiter::limitHeight(Base::Element* element)
{
    double minimum = getMinimumHeight(element);
    logger(DEBUG, "minimum before adaption in element %: %",element->getID(), minimum);
    const double averageH = Helpers::computeAverageOfSolution<1>(element, elementIntegrator_)(0);
    const double averageHU = Helpers::computeAverageOfSolution<1>(element, elementIntegrator_)(1);
    if (averageH < minH_)
    {
        //solution is constant with value average
        LinearAlgebra::MiddleSizeVector newCoeffsH = Helpers::projectOnBasisFuns<1>(element, [=](const PointReferenceT& pRef){return averageH;}, elementIntegrator_);
        element->setTimeLevelData(0,0,newCoeffsH);
        return;
    }
    //squeeze around average, such that the minimum height is at least minH_
    const double theta = 0;//std::min(1.0, (averageH - minH_)/(averageH - minimum));
    logger(DEBUG, "Theta: %", theta);
    std::function<double(const PointReferenceT&)> newFun = [=] (const PointReferenceT& pRef){return theta*(element->getSolution(0,pRef)(0) - averageH) + averageH;};
    LinearAlgebra::MiddleSizeVector newCoeffsH = Helpers::projectOnBasisFuns<1>(element, newFun, elementIntegrator_);
    element->setTimeLevelData(0, 0, newCoeffsH);
    logger(DEBUG, "minimum after adaption: % (average %)", getMinimumHeight(element), averageH);
}

void PositiveLayerLimiter::limitDischarge(Base::Element* element)
{
    const PointReferenceT &pRefL = element->getReferenceGeometry()->getReferenceNodeCoordinate(0); 
    const PointReferenceT &pRefR = element->getReferenceGeometry()->getReferenceNodeCoordinate(1);
    if (element->getSolution(0,pRefL)(0) <= minH_ && element->getSolution(0,pRefR)(0) <= minH_ )
    {
        LinearAlgebra::MiddleSizeVector newCoeffsHU = Helpers::projectOnBasisFuns<1>(element, [=](const PointReferenceT& pRef){return 0;}, elementIntegrator_);
        element->setTimeLevelData(0,1,newCoeffsHU);
        return;
    }
    if (element->getSolution(0,pRefL)(0) > minH_ && element->getSolution(0,pRefR)(0) > minH_ )
    {
        return;
    }
    double huL = element->getSolution(0,pRefL)(1);
    double huR = element->getSolution(0,pRefR)(1);
    double average = (huL + huR) / 2;
    double slope = 0;//(huL + huR)/2;

    if (element->getSolution(0, pRefR)(0) <= minH_)
    {
        slope *= -1;
    }
    std::function<double(const PointReferenceT&) > newFun = [ = ] (const PointReferenceT & pRef){return average + slope * pRef[0];};
    LinearAlgebra::MiddleSizeVector newCoeffs = Helpers::projectOnBasisFuns<1>(element, newFun, elementIntegrator_);
    element->setTimeLevelData(0, 1, newCoeffs);
    logger(DEBUG, "Values of hu: [%,%]", element->getSolution(0,pRefL)(1), element->getSolution(0,pRefR)(1));
}

//Note that this only holds for a flat bottom!
void PositiveLayerLimiter::isDryElement(Base::Element* elt)
{
    bool flag = static_cast<Helpers::DryFlag*>(elt->getUserData())->isDry;
    logger(DEBUG, "isDry: %", flag);
    double averageH = Helpers::computeAverageOfSolution<1>(elt, elementIntegrator_)(0);
    flag = (averageH > minH_);
}