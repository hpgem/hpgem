/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2014, University of Twente
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.

 3. Neither the name of the copyright holder nor the names of its contributors
 may be used to endorse or promote products derived from this software without
 specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "../HelperFunctions.h"

using namespace hpgem;

void SqueezeLimiterWithLayer::limit(
    Base::Element *element,
    LinearAlgebra::MiddleSizeVector &solutionCoefficients) {
    const double minimumHeight = getMinimumHeight(element);
    if (minimumHeight >= minH_) return;

    const double averageHeight = Helpers::computeAverageOfSolution<1>(
        element, solutionCoefficients, elementIntegrator_)(0);
    const double averageDischarge = Helpers::computeAverageOfSolution<1>(
        element, solutionCoefficients, elementIntegrator_)(1);
    LinearAlgebra::MiddleSizeVector heightCoefficients;
    LinearAlgebra::MiddleSizeVector dischargeCoefficients;
    if (averageHeight < minH_) {
        heightCoefficients = Helpers::projectOnBasisFuns<1>(
            element, [=](const PointReferenceT &pRef) { return averageHeight; },
            elementIntegrator_);
        dischargeCoefficients = Helpers::projectOnBasisFuns<1>(
            element, [=](const PointReferenceT &pRef) { return 0; },
            elementIntegrator_);
    } else {
        const PointReferenceT &pRefL =
            element->getReferenceGeometry()->getReferenceNodeCoordinate(0);
        const double solutionLeft = element->getSolution(0, pRefL)(0);
        const PointReferenceT &pRefR =
            element->getReferenceGeometry()->getReferenceNodeCoordinate(1);
        const double solutionRight = element->getSolution(0, pRefR)(0);

        double slopeDischarge = averageDischarge;
        double slopeHeight =
            std::min((solutionLeft + solutionRight - 2 * minH_) / 2,
                     averageHeight - minH_);
        double squeezeFactorH =
            (averageHeight - minH_) / (averageHeight - minimumHeight);
        const PointReferenceT &pRef = *getMinimumHeightPoint(element);
        double squeezeFactorHu =
            std::min(1., averageDischarge / (averageDischarge -
                                             element->getSolution(0, pRef)(1)));

        std::function<double(const PointReferenceT &)> heightFun =
            [=](const PointReferenceT &pRef) {
                return averageHeight +
                       squeezeFactorH *
                           (element->getSolution(0, pRef)(0) - averageHeight);
            };
        heightCoefficients = Helpers::projectOnBasisFuns<1>(element, heightFun,
                                                            elementIntegrator_);
        std::function<double(const PointReferenceT &)> dischargeFun =
            [=](const PointReferenceT &pRef) {
                return averageDischarge +
                       squeezeFactorHu * (element->getSolution(0, pRef)[1] -
                                          averageDischarge);
            };
        dischargeCoefficients = Helpers::projectOnBasisFuns<1>(
            element, dischargeFun, elementIntegrator_);
        logger(DEBUG, "averages: %",
               Helpers::computeAverageOfSolution<1>(
                   element, solutionCoefficients, elementIntegrator_));
        logger(DEBUG, "slopes: [% %]", slopeHeight, slopeDischarge);
        logger(DEBUG, "heightCoefficients: %", heightCoefficients);
        logger(DEBUG, "discharge coefficients: %", dischargeCoefficients);
    }

    for (std::size_t iFun = 0; iFun < element->getNumberOfBasisFunctions();
         ++iFun) {
        std::size_t iVF = element->convertToSingleIndex(iFun, 0);
        solutionCoefficients[iVF] = heightCoefficients[iFun];
        iVF = element->convertToSingleIndex(iFun, 1);
        solutionCoefficients[iVF] = dischargeCoefficients[iFun];
    }
}

double SqueezeLimiterWithLayer::getMinimumHeight(const Base::Element *element) {
    const PointReferenceT &pRefL =
        element->getReferenceGeometry()->getReferenceNodeCoordinate(0);
    const LinearAlgebra::MiddleSizeVector &solutionCoefficients =
        element->getTimeIntegrationVector(0);
    const std::size_t numOfVariables = element->getNumberOfUnknowns();
    const double solutionLeft = element->getSolution(
        0, pRefL)(0);  // Helpers::getSolution<DIM>(element,
                       // solutionCoefficients, pRefL, numOfVariables)(0);
    double minimum = solutionLeft;
    for (std::size_t iPoint = 1;
         iPoint < element->getReferenceGeometry()->getNumberOfNodes();
         ++iPoint) {
        const PointReferenceT &pRef =
            element->getReferenceGeometry()->getReferenceNodeCoordinate(iPoint);
        minimum = std::min(minimum, element->getSolution(0, pRef)(0));
    }
    for (std::size_t p = 0;
         p < element->getGaussQuadratureRule()->getNumberOfPoints(); ++p) {
        const PointReferenceT &pRef =
            element->getGaussQuadratureRule()->getPoint(p);
        minimum = std::min(minimum, element->getSolution(0, pRef)(0));
    }
    logger(DEBUG, "Minimum in element %: %", element->getID(), minimum);
    return minimum;
}

///\return a POINTER to the reference point where the height is smallest. This
/// cannot be a reference, since the reference is not allowed to be copied.
const Geometry::PointReference<1>
    *SqueezeLimiterWithLayer::getMinimumHeightPoint(
        const Base::Element *element) {
    const PointReferenceT &pRefL =
        element->getReferenceGeometry()->getReferenceNodeCoordinate(0);
    const PointReferenceT &pRefR =
        element->getReferenceGeometry()->getReferenceNodeCoordinate(1);
    const LinearAlgebra::MiddleSizeVector &solutionCoefficients =
        element->getTimeIntegrationVector(0);
    const std::size_t numOfVariables = element->getNumberOfUnknowns();
    const double solutionLeft = element->getSolution(0, pRefL)(0);
    const double solutionRight = element->getSolution(0, pRefL)(1);
    double minimum = solutionLeft;
    PointReferenceT const *minPoint = &pRefL;
    if (solutionRight < solutionLeft) {
        minimum = solutionRight;
        minPoint = &pRefR;
    }
    for (std::size_t p = 0;
         p < element->getGaussQuadratureRule()->getNumberOfPoints(); ++p) {
        const PointReferenceT &pRef =
            element->getGaussQuadratureRule()->getPoint(p);
        if (minimum > element->getSolution(0, pRef)(0)) {
            minimum = element->getSolution(0, pRef)(0);
            minPoint = &pRef;
        }
    }
    return minPoint;
}