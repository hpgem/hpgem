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

#include "TvbLimiter1D.h"

using namespace hpgem;

TvbLimiter1D::TvbLimiter1D(std::size_t numberOfVariables)
    : SlopeLimiter(numberOfVariables) {}

void TvbLimiter1D::limitSlope(Base::Element *element) {
    for (std::size_t iVar = 0; iVar < numberOfVariables_; ++iVar) {
        if (!hasSmallSlope(element, iVar)) {
            limitWithMinMod(element, iVar);
        }
    }
}

///\details Doing it the DG way for p=1 does not lead to conservation of mass,
/// so for now, do it the FVM way. Don't look at the slope in this element at
/// all, just at the averages of this element and the adjacent elements. We also
/// assume that we limit only once per time step. Lastly, it is assumed that the
/// first two basis functions are the nodal basis functions
void TvbLimiter1D::limitWithMinMod(Base::Element *element,
                                   const std::size_t iVar) {
    // this does not work for boundaries, probably also not for triangular mesh.
    // maybe write getNeighbour(face), which deals with nullpointers in case of
    // a boundary?
    const Base::Element *elemL;
    const Base::Element *elemR;
    if (element->getFace(0)->isInternal())
        elemL = element->getFace(0)->getPtrElementLeft();
    else
        return;  // for now, just don't use the limiter here...
    if (element->getFace(1)->isInternal())
        elemR = element->getFace(1)->getPtrElementRight();
    else
        return;  // for now, just don't use the limiter here...

    logger.assert_always(elemR->getID() == element->getID() + 1 &&
                             element->getID() == elemL->getID() + 1,
                         "elements not in correct order");

    const double u0 = Helpers::computeAverageOfSolution<1>(
        element, element->getTimeIntegrationVector(0),
        elementIntegrator_)(iVar);
    const double uElemR = Helpers::computeAverageOfSolution<1>(
        const_cast<Base::Element *>(elemR), elemR->getTimeIntegrationVector(0),
        elementIntegrator_)(iVar);
    const double uElemL = Helpers::computeAverageOfSolution<1>(
        const_cast<Base::Element *>(elemL), elemL->getTimeIntegrationVector(0),
        elementIntegrator_)(iVar);

    logger(INFO, "uLeft: %\nu0: %\nuRight: %", uElemL, u0, uElemR);

    LinearAlgebra::MiddleSizeVector newCoeffs =
        LinearAlgebra::MiddleSizeVector(element->getNumberOfBasisFunctions());
    if (Helpers::sign(uElemR - u0) != Helpers::sign(u0 - uElemL))  // phi(r) = 0
    {
        newCoeffs[0] = u0;
        newCoeffs[1] = u0;
        logger(INFO, "phi = 0");
    } else {
        if ((u0 - uElemL) / (uElemR - u0) < 1)  // phi(r) = r
        {
            newCoeffs[0] = .5 * u0 + .5 * uElemL;
            newCoeffs[1] = 1.5 * u0 - .5 * uElemL;
            logger(INFO, "phi = r");
        } else  // phi(r) = 1
        {
            newCoeffs[0] = 1.5 * u0 - .5 * uElemR;
            newCoeffs[1] = .5 * u0 + .5 * uElemR;
            logger(INFO, "phi = 1");
        }
    }
    logger(INFO, "new coefficients: % \n \n", newCoeffs);
    element->setTimeIntegrationSubvector(0, iVar, newCoeffs);
}

bool TvbLimiter1D::hasSmallSlope(const Base::Element *element,
                                 const std::size_t iVar) {
    const PointReferenceT &pRefL =
        element->getReferenceGeometry()->getReferenceNodeCoordinate(0);
    const PointReferenceT &pRefR =
        element->getReferenceGeometry()->getReferenceNodeCoordinate(1);

    const double uPlus = element->getSolution(0, pRefR)[iVar];
    const double uMinus = element->getSolution(0, pRefL)[iVar];
    const double M = 1;
    return (std::abs(uPlus - uMinus) <
            M * (2 * element->calcJacobian(pRefL).determinant()) *
                (2 * element->calcJacobian(pRefL).determinant()));
}