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

#include "TvbLimiterWithDetector1D.h"
#include "HelperFunctions.h"


void TvbLimiterWithDetector1D::limitSlope(Base::Element *element)
{
    for (std::size_t i = 0; i < numOfVariables_; ++i)
    {
        if (detectDiscontinuity(element)[i] && !hasSmallSlope(element, i))
        {
            limitWithMinMod(element, i);
        }
    } 
}

///\details In here, the discontinuity detector of Krivodonova et. al (2004) is
///implemented to determine whether or not the given element needs limiting. If
///limiting is needed, this function calls limitWithMinmod
std::vector<bool> TvbLimiterWithDetector1D::detectDiscontinuity(Base::Element* element)
{
    LinearAlgebra::MiddleSizeVector totalIntegral(numOfVariables_);
    LinearAlgebra::MiddleSizeVector numericalSolution(numOfVariables_);
    std::vector<bool> isDiscontinuous(numOfVariables_, false);
    const LinearAlgebra::MiddleSizeVector &solutionCoefficients = element->getTimeIntegrationVector(0);

    //For every face of this element, check the size of the jump
    for (const Base::Face *face : element->getFacesList())
    {
        //Compute the numerical solution and the (outward) normal of the face
        const Geometry::PointReference<0>& pRefForInflowTest = face->getReferenceGeometry()->getCenter();
        Base::Side sideOfElement = (face->getPtrElementLeft() == element) ? Base::Side::LEFT : Base::Side::RIGHT;
        LinearAlgebra::SmallVector<1> normal = face->getNormalVector(pRefForInflowTest);
        if (sideOfElement == Base::Side::LEFT)
        {
            numericalSolution = Helpers::getSolution<1>(face->getPtrElement(sideOfElement), solutionCoefficients, face->mapRefFaceToRefElemL(pRefForInflowTest), numOfVariables_);
        }
        else
        {
            numericalSolution = Helpers::getSolution<1>(face->getPtrElement(sideOfElement), solutionCoefficients, face->mapRefFaceToRefElemR(pRefForInflowTest), numOfVariables_);
            normal *= -1;
        }

        //determine the direction of the flow across this face
        LinearAlgebra::SmallVector<1> velocity = computeVelocity(numericalSolution);
        //if this is an inflow face, integrate the difference in {h,hu} over the face
        if (velocity * normal < -1e-14)
        {
            logger(DEBUG, "Face % is an inflow face for element %.", face->getID(), element->getID());
            LinearAlgebra::MiddleSizeVector numericalSolutionOther(numOfVariables_);
            if (face->isInternal())
            {
                if (sideOfElement == Base::Side::LEFT)
                {
                    LinearAlgebra::MiddleSizeVector solutionCoefficientsOther = face->getPtrElement(Base::Side::RIGHT)->getTimeIntegrationVector(0);
                    numericalSolutionOther = Helpers::getSolution<1>(face->getPtrElement(Base::Side::RIGHT), solutionCoefficientsOther, face->mapRefFaceToRefElemR(pRefForInflowTest), numOfVariables_);
                }
                else
                {
                    LinearAlgebra::MiddleSizeVector solutionCoefficientsOther = face->getPtrElement(Base::Side::LEFT)->getTimeIntegrationVector(0);
                    numericalSolutionOther = Helpers::getSolution<1>(face->getPtrElement(Base::Side::LEFT), solutionCoefficientsOther, face->mapRefFaceToRefElemL(pRefForInflowTest), numOfVariables_);
                }
            }
            else
            {
                numericalSolutionOther = inflowBC_;
            }

            ///\todo make this multi-dimensional by integrating over the face.
            totalIntegral += numericalSolution - numericalSolutionOther;
        }
    }

    logger(DEBUG, "Integral over all inflow boundaries of difference for element %: %", element->getID(), totalIntegral);

    //divide the integral by dx^{(p+1)/2}, the norm of {h,hu} and the size of the face
    ///\todo make dx the maximum of the edge-lengths
    ///\todo check if dx is computed correctly
    const Geometry::PointPhysical<1>& pPhys0 = element->getPhysicalGeometry()->getLocalNodeCoordinates(0);
    const Geometry::PointPhysical<1>& pPhys1 = element->getPhysicalGeometry()->getLocalNodeCoordinates(1);
    const double dx = Base::L2Norm(pPhys0 - pPhys1);
    logger(DEBUG, "grid size: %", dx);
    totalIntegral /= std::pow(dx, (polynomialOrder_ + 1.) / 2);

    LinearAlgebra::MiddleSizeVector average = Helpers::computeAverageOfSolution<1>(element, solutionCoefficients, elementIntegrator_);
    for (std::size_t i = 0; i < numOfVariables_; ++i)
    {
        if (average(i) > 1e-10)
            totalIntegral(i) /= std::abs(average(i));
    }
    ///\todo divide by the size of the faces for DIM > 1

    //compare this with 1: if >1, discontinuous, if <1 smooth
    for (std::size_t i = 0; i < numOfVariables_; ++i)
    {
        if (std::abs(totalIntegral(i)) > 1)
        {
            logger(DEBUG, "Element % with variable % will be limited.", element->getID(), i);
            isDiscontinuous[i] = true;
        }
    }
 
        return isDiscontinuous;
}

LinearAlgebra::SmallVector<1> TvbLimiterWithDetector1D::computeVelocity(LinearAlgebra::MiddleSizeVector numericalSolution)
{
    if (numericalSolution(0) < 1e-5)
    {
        return LinearAlgebra::SmallVector<1>({0.});
    }
    
    return LinearAlgebra::SmallVector<1>({numericalSolution(1) / numericalSolution(0)});
}

bool TvbLimiterWithDetector1D::hasSmallSlope(Base::Element* element, std::size_t iVar)
{
    const PointReferenceT &pRefL = element->getReferenceGeometry()->getReferenceNodeCoordinate(0); 
    const PointReferenceT &pRefR = element->getReferenceGeometry()->getReferenceNodeCoordinate(1);
    
    const double uPlus = element->getSolution(0, pRefR)[iVar];
    const double uMinus = element->getSolution(0, pRefL)[iVar];
    const double M = 500;
    return (std::abs(uPlus - uMinus) < M * (2*element->calcJacobian(pRefL).determinant())* (2*element->calcJacobian(pRefL).determinant()));
}

void TvbLimiterWithDetector1D::limitWithMinMod(Base::Element* element, const std::size_t iVar)
{    
    
    const LinearAlgebra::MiddleSizeVector &solutionCoefficients = element->getTimeIntegrationVector(0);
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

    //check if left is indeed of the left side of right //TODO: this may be false for general meshes (or may even examine the same node twice)
    logger.assert((PointPhysicalT(elemL->getPhysicalGeometry()->getLocalNodeCoordinates(0)))[0] < PointPhysicalT(elemR->getPhysicalGeometry()->getLocalNodeCoordinates(0))[0], "elements left/right in wrong order");

    const double uPlus = element->getSolution(0, pRefR)[iVar];
    const double uMinus = element->getSolution(0, pRefL)[iVar];
    
    const double u0 = Helpers::computeAverageOfSolution<1>(element, solutionCoefficients, elementIntegrator_)(iVar);
    const double uElemR = Helpers::computeAverageOfSolution<1>(const_cast<Base::Element*>(elemR), solutionCoefficients, elementIntegrator_)(iVar);
    const double uElemL = Helpers::computeAverageOfSolution<1>(const_cast<Base::Element*>(elemL), solutionCoefficients, elementIntegrator_)(iVar);
    logger(DEBUG, "coefficients: %", element->getTimeIntegrationSubvector(0, iVar));
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
        element->setTimeIntegrationSubvector(0, iVar, newCoeffs);
        logger(INFO, "Element % is limited, slope %!", element->getID(), slope);
    }
}
