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
#include "SavageHutter1DBase.h"
#include "HelperFunctions.h"
#include <cmath>

Base::RectangularMeshDescriptor<1> SavageHutter1DBase::createMeshDescription(const std::size_t numberOfElements, const double endOfDomain, const Base::BoundaryType boundary)
{
    // Make a description of the domain [0,endOfDomain] with the correct boundary conditions and number of elements
    Base::RectangularMeshDescriptor<1> description;
    description.boundaryConditions_[0] = boundary;
    description.numElementsInDIM_[0] = numberOfElements;
    description.bottomLeft_[0] = 0;
    description.topRight_[0] = endOfDomain;
    return description;
}

/// \details The integrand for the reference element is the same as the physical element, but scaled with the reference-to-physical element scale, which is the determinant of the jacobian of the reference-to-physical element mapping.
const LinearAlgebra::MiddleSizeVector SavageHutter1DBase::integrandRightHandSideOnElement
(Base::PhysicalElement<1>& element, const double &time, const LinearAlgebra::MiddleSizeVector &solutionCoefficients)
{
    const std::size_t numBasisFuncs = element.getElement()->getNumberOfBasisFunctions();

    MiddleSizeVector& integrand = element.getResultVector(); //just to have the correct length    
    const PointPhysicalT& pPhys = element.getPointPhysical();
    const LinearAlgebra::MiddleSizeVector numericalSolution = Helpers::getSolution<1>(element, solutionCoefficients, numberOfVariables_);
    const LinearAlgebra::MiddleSizeVector physicalFlux = computePhysicalFlux(numericalSolution, element.getPointPhysical());
    const LinearAlgebra::MiddleSizeVector source = computeSourceTerm(numericalSolution, pPhys, time);

    // Compute integrand on the physical element.
    std::size_t iVB; // Index for both basis function and variable
    for (std::size_t iB = 0; iB < numBasisFuncs; iB++) // Index for basis function
    {
        for (std::size_t iV = 0; iV < numberOfVariables_; ++iV)
        {
            iVB = element.getElement()->convertToSingleIndex(iB, iV);
            integrand(iVB) = physicalFlux(iV) * element.basisFunctionDeriv(iB)(0);
            integrand(iVB) += source(iV) * element.basisFunction(iB);
        }
    }

    return integrand;
}

/// \details The integrand for the reference face is the same as the physical face, but scaled with the reference-to-physical face scale. This face scale is absorbed in the normal vector, since it is relatively cheap to compute the normal vector with a length (L2-norm) equal to the reference-to-physical face scale.
const LinearAlgebra::MiddleSizeVector SavageHutter1DBase::integrandRightHandSideOnRefFace
(Base::PhysicalFace<1>& face, const Base::Side &iSide, const LinearAlgebra::MiddleSizeVector &solutionCoefficientsLeft, const LinearAlgebra::MiddleSizeVector &solutionCoefficientsRight)
{
    double normal = face.getNormalVector()[0];
    const std::size_t numTestBasisFuncs = face.getFace()->getPtrElement(iSide)->getNumberOfBasisFunctions();

    //compute numerical solution at the left side and right side of this face
    LinearAlgebra::MiddleSizeVector solutionLeft = Helpers::getSolution<1>(face.getPhysicalElement(Base::Side::LEFT), solutionCoefficientsLeft, numberOfVariables_);
    LinearAlgebra::MiddleSizeVector solutionRight = Helpers::getSolution<1>(face.getPhysicalElement(Base::Side::RIGHT), solutionCoefficientsRight, numberOfVariables_);

    LinearAlgebra::MiddleSizeVector flux = hllcFlux(solutionLeft, solutionRight, face.getUnitNormalVector()[0], face);

    if (iSide == Base::Side::RIGHT) //the normal is defined for the left element
    {
        normal *= -1;
    }
    LinearAlgebra::MiddleSizeVector& integrand = face.getResultVector(iSide); // Integrand value based on n number of testbasisfunctions from element corresponding to side iSide

    for (std::size_t iFun = 0; iFun < numTestBasisFuncs; ++iFun)
    {
        for (std::size_t iVar = 0; iVar < numberOfVariables_; ++iVar)
        {
            std::size_t iVarFun = face.getFace()->getPtrElement(iSide)->convertToSingleIndex(iFun, iVar);
            integrand(iVarFun) = -flux(iVar) * face.basisFunction(iSide, iFun) * normal;
        }
    }

    return integrand;
}

///\details Compute the flux on the boundary. Note that there is a lot of dead code
///in here, this can be deleted when everything proves to work fine or in April 2016, whichever comes first.
const LinearAlgebra::MiddleSizeVector SavageHutter1DBase::integrandRightHandSideOnRefFace
(
 Base::PhysicalFace<1>& face,
 const LinearAlgebra::MiddleSizeVector &solutionCoefficients,
 const double time
 )
{
    
    /* Dead code that is a back-up of the boundary conditions.
    MiddleSizeVector flux(numberOfVariables_);
    double u = 0;
    if (solution(0) > dryLimit_)
    {
        u = solution(1) / solution(0);
    }
    
    //outflow
    if (normal > 0)
    {
        //solid wall on the end
        if (time > 10 && time < 7)
        {
            LinearAlgebra::MiddleSizeVector reflection;
            if (time < 10)
            {
                reflection = LinearAlgebra::MiddleSizeVector({solution[0], -solution[1]});
            }
            else
            {
                reflection = LinearAlgebra::MiddleSizeVector({solution[0], 20 * (time - 0.15) * solution[1]});
            }
            flux = hllcFlux(solution, reflection, 1);
        }
        else
        {
            //reference for how to do interpolation with finite differences of faces:
            const Base::Face *otherFace = face.getPhysicalElement(Base::Side::LEFT).getElement()->getFace(0);
            logger(DEBUG, "face before last: %, position %", otherFace->getID(), otherFace->referenceToPhysical(pRef));
            MiddleSizeVector otherSideSolution = 0.5 * (otherFace->getPtrElementRight()->getSolution(0, otherFace->mapRefFaceToRefElemR(pRef)) + otherFace->getPtrElementLeft()->getSolution(0, otherFace->mapRefFaceToRefElemL(pRef)));
            const Base::Element *otherElement = otherFace->getPtrElementLeft();
            logger(DEBUG, "one back: %", otherElement->getID());
            otherFace = otherElement->getFace(0);
            MiddleSizeVector twoBackSolution = 0.5 * (otherFace->getPtrElementRight()->getSolution(0, otherFace->mapRefFaceToRefElemR(pRef)) + otherFace->getPtrElementLeft()->getSolution(0, otherFace->mapRefFaceToRefElemL(pRef)));
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
            
            MiddleSizeVector ghostSolution = 2 * solution - 1 * otherSideSolution + 0 * twoBackSolution;

            //subcritical outflow:
            if ((solution[1] / solution[0] < std::sqrt(solution[0] * std::cos(chuteAngle_) * epsilon_)))
            {
                double uIn = solution[1] / solution[0];
                double invariantIn = uIn + 2 * std::sqrt(epsilon_ * std::cos(chuteAngle_) * solution[0]);
                double froudeIn = solution[1] / solution[0] / std::sqrt(epsilon_ * std::cos(chuteAngle_) * solution[0]);
                double froudePrescribed = 1;
                double dischargePrescribed = .5;
                //double hOut = bisectionHFinder(dischargePrescribed, solution);
                double hOut = (invariantIn / (2 + froudePrescribed)) * (invariantIn / (2 + froudePrescribed)) / (epsilon_ * std::cos(chuteAngle_));
                //double hOut = solution[0];
                double uOut = invariantIn - 2 * std::sqrt(epsilon_ * std::cos(chuteAngle_) * hOut);
                logger(DEBUG, "new h and u and F: %, %, %", hOut, uOut, uOut / std::sqrt(epsilon_ * std::cos(chuteAngle_) * hOut));
                auto stateNew = MiddleSizeVector({hOut, hOut * uOut});
                flux = hllcFlux(solution, stateNew, 1);
            }
            else
            {
                //supercritical outflow:
                flux = hllcFlux(solution, solution, face.getUnitNormalVector()[0]);
            }
        }
    }
    else //inflow
    {
        //subcritical inflow:
        if (solution[1] / solution[0] < std::sqrt(solution[0] * std::cos(chuteAngle_) * epsilon_))
        {
            double hNew = inflowBC_[0];
            double uNew = solution[1] / solution[0] - 2 * std::sqrt(epsilon_ * std::cos(chuteAngle_))*(std::sqrt(solution[0]) - std::sqrt(hNew));
            auto stateNew = MiddleSizeVector({hNew, hNew * uNew});
            flux = hllcFlux(stateNew, solution, 1);
        }
        else
        {
            flux = hllcFlux(inflowBC_, solution, 1);
        }
    }*/
    
    double normal = face.getNormalVector()[0];
    const std::size_t numberOfBasisFuncs = face.getFace()->getNumberOfBasisFunctions();

    //note that at the boundary, the element is the left element by definition
    LinearAlgebra::MiddleSizeVector solution = Helpers::getSolution<1>(face.getPhysicalElement(Base::Side::LEFT), solutionCoefficients, numberOfVariables_);
    LinearAlgebra::MiddleSizeVector flux = hllcFlux(solution, computeGhostSolution(solution, normal, time), 1, face);
    //enforce inflow bc strongly
    if (std::abs(computeGhostSolution(solution, normal, time)[1] - inflowBC_[1]) < 1e-10)
    {
        flux = computePhysicalFlux(inflowBC_, face.getPointPhysical());
    }
    LinearAlgebra::MiddleSizeVector integrand(numberOfVariables_ * numberOfBasisFuncs);

    for (std::size_t iFun = 0; iFun < numberOfBasisFuncs; ++iFun)
    {
        for (std::size_t iVar = 0; iVar < numberOfVariables_; ++iVar)
        {
            std::size_t iVarFun = face.getFace()->getPtrElementLeft()->convertToSingleIndex(iFun, iVar);
            integrand(iVarFun) = -flux(iVar) * face.basisFunction(iFun) * normal;
        }
    }
    
    return integrand;
}



LinearAlgebra::MiddleSizeVector SavageHutter1DBase::localLaxFriedrichsFlux(const LinearAlgebra::MiddleSizeVector& numericalSolutionLeft, const LinearAlgebra::MiddleSizeVector& numericalSolutionRight, Base::PhysicalFace<1> &face)
{
    double uLeft = 0;
    if (numericalSolutionLeft(0) > dryLimit_)
    {
        uLeft = numericalSolutionLeft(1) / numericalSolutionLeft(0);
    }

    double uRight = 0;
    if (numericalSolutionRight(0) > dryLimit_)
    {
        uRight = numericalSolutionRight(1) / numericalSolutionRight(0);
    }

    const double alpha = std::max(std::abs(uLeft) + std::sqrt(epsilon_ * std::max(0., numericalSolutionLeft(0))),
                                  std::abs(uRight) + std::sqrt(epsilon_ * std::max(0., numericalSolutionRight(0))));

    logger(DEBUG, "alpha: %", alpha);

    logger(DEBUG, "physical fluxes: %, %", computePhysicalFlux(numericalSolutionLeft, face.getPointPhysical()), computePhysicalFlux(numericalSolutionRight, face.getPointPhysical()));
    LinearAlgebra::MiddleSizeVector diffSolutions = numericalSolutionRight - numericalSolutionLeft;
    const LinearAlgebra::MiddleSizeVector numericalFlux = 0.5 *
        (computePhysicalFlux(numericalSolutionLeft, face.getPointPhysical()) + computePhysicalFlux(numericalSolutionRight, face.getPointPhysical())
         - alpha * (diffSolutions));

    return numericalFlux;
}

LinearAlgebra::MiddleSizeVector SavageHutter1DBase::hllcFlux(const LinearAlgebra::MiddleSizeVector& numericalSolutionLeft, const LinearAlgebra::MiddleSizeVector& numericalSolutionRight, const double normal, Base::PhysicalFace<1> &face)
{
    const double hLeft = numericalSolutionLeft[0];
    const double hRight = numericalSolutionRight[0];
    double normalSpeedLeft = 0;
    if (numericalSolutionLeft(0) > dryLimit_)
    {
        normalSpeedLeft = normal * numericalSolutionLeft(1) / numericalSolutionLeft(0);
    }

    double normalSpeedRight = 0;
    if (numericalSolutionRight(0) > dryLimit_)
    {
        normalSpeedRight = normal * numericalSolutionRight(1) / numericalSolutionRight(0);
    }
    LinearAlgebra::MiddleSizeVector fluxNormalLeft = normal * computePhysicalFlux(numericalSolutionLeft, face.getPointPhysical());
    LinearAlgebra::MiddleSizeVector fluxNormalRight = normal * computePhysicalFlux(numericalSolutionRight, face.getPointPhysical());
    double phaseSpeedLeft = std::sqrt(epsilon_ * std::cos(chuteAngle_) * hLeft);
    double phaseSpeedRight = std::sqrt(epsilon_ * std::cos(chuteAngle_) * hRight);
    double sl = std::min(normalSpeedLeft - phaseSpeedLeft, normalSpeedRight - phaseSpeedRight);
    double sr = std::max(normalSpeedLeft + phaseSpeedLeft, normalSpeedRight + phaseSpeedRight);
    if (sl >= 0)
    {
        return fluxNormalLeft;
    }
    if (sr <= 0)
    {
        return fluxNormalRight;
    }
    return ((sr * fluxNormalLeft - sl * fluxNormalRight + sl * sr * (numericalSolutionRight - numericalSolutionLeft)) / (sr - sl));
}


double SavageHutter1DBase::computeFrictionCoulomb(const LinearAlgebra::MiddleSizeVector& numericalSolution, const double frictionAngle)
{
    return std::tan(frictionAngle);
}

/// Compute friction as described in Weinhart et al (2012), eq (50)
/// Notice that the minus sign for gamma in eq (50) is wrong.
double SavageHutter1DBase::computeFriction(const LinearAlgebra::MiddleSizeVector& numericalSolution)
{
    const double delta1 = 17.518 / 180 * M_PI;
    const double delta2 = 29.712 / 180 * M_PI;
    const double A = 5.29;
    const double beta = 0.189;
    const double gamma = -.080;
    const double d = 2;
    const double h = numericalSolution[0];
    if (h < dryLimit_)
        return std::tan(delta1);
    const double u = numericalSolution[1] / h;
    const double F = u / std::sqrt(epsilon_ * std::cos(chuteAngle_) * h);
    return std::tan(delta1) + (std::tan(delta2) - std::tan(delta1)) / (beta * h / (A * d * (F - gamma)) + 1);
}

/// Compute friction as in Anthony's draft, eq (3.9)
double SavageHutter1DBase::computeFrictionExponential(const LinearAlgebra::MiddleSizeVector& numericalSolution)
{
    const double delta1 = 27. / 180 * M_PI;
    const double delta2 = 37. / 180 * M_PI;
    const double beta = 0.136;
    const double L = 2.;
    const double h = numericalSolution[0];
    if (h < dryLimit_)
        return std::tan(delta1);
    const double u = numericalSolution[1] / h;
    if (std::abs(u) < 1e-16)
        return std::tan(delta1);
    return std::tan(delta1) + (std::tan(delta2) - std::tan(delta1)) * std::exp(-beta * std::pow(epsilon_*h, 1.5) / (L * std::abs(u)));
}

