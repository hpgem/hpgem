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
#include "SavageHutter1DBase.h"
#include "HelperFunctions.h"
#include <cmath>

/// \details The integrand for the reference element is the same as the physical
/// element, but scaled with the reference-to-physical element scale, which is
/// the determinant of the jacobian of the reference-to-physical element
/// mapping.
const LinearAlgebra::MiddleSizeVector
    SavageHutter1DBase::integrandRightHandSideOnElement(
        Base::PhysicalElement<1> &element, const double &time,
        const LinearAlgebra::MiddleSizeVector &solutionCoefficients) {
    const std::size_t numberOfBasisFunctions =
        element.getNumberOfBasisFunctions();

    MiddleSizeVector &integrand =
        element.getResultVector();  // just to have the correct length
    const PointPhysicalT &pPhys = element.getPointPhysical();
    const LinearAlgebra::MiddleSizeVector numericalSolution =
        Helpers::getSolution<1>(element, solutionCoefficients,
                                numberOfVariables_);
    const LinearAlgebra::MiddleSizeVector physicalFlux =
        computePhysicalFlux(numericalSolution, element.getPointPhysical());
    const LinearAlgebra::MiddleSizeVector source =
        computeSourceTerm(numericalSolution, pPhys, time);

    // Compute integrand on the physical element.
    for (std::size_t iFun = 0; iFun < numberOfBasisFunctions; iFun++) {
        for (std::size_t iVar = 0; iVar < numberOfVariables_; ++iVar) {
            const std::size_t iVarFun =
                element.convertToSingleIndex(iFun, iVar);
            integrand(iVarFun) =
                physicalFlux(iVar) * element.basisFunctionDeriv(iFun)(0);
            integrand(iVarFun) += source(iVar) * element.basisFunction(iFun);
        }
    }

    return integrand;
}

/// \details The integrand for the reference face is the same as the physical
/// face, but scaled with the reference-to-physical face scale. This face scale
/// is absorbed in the normal vector, since it is relatively cheap to compute
/// the normal vector with a length (L2-norm) equal to the reference-to-physical
/// face scale.
const LinearAlgebra::MiddleSizeVector
    SavageHutter1DBase::integrandRightHandSideOnRefFace(
        Base::PhysicalFace<1> &face, const Base::Side &iSide,
        const LinearAlgebra::MiddleSizeVector &solutionCoefficientsLeft,
        const LinearAlgebra::MiddleSizeVector &solutionCoefficientsRight) {
    const std::size_t numberOfTestBasisFunctions =
        face.getPhysicalElement(iSide).getNumberOfBasisFunctions();

    // compute numerical solution at the left side and right side of this face
    const LinearAlgebra::MiddleSizeVector solutionLeft =
        Helpers::getSolution<1>(face.getPhysicalElement(Base::Side::LEFT),
                                solutionCoefficientsLeft, numberOfVariables_);
    const LinearAlgebra::MiddleSizeVector solutionRight =
        Helpers::getSolution<1>(face.getPhysicalElement(Base::Side::RIGHT),
                                solutionCoefficientsRight, numberOfVariables_);

    const LinearAlgebra::MiddleSizeVector flux = hllcFlux(
        solutionLeft, solutionRight, face.getUnitNormalVector()[0], face);

    LinearAlgebra::MiddleSizeVector &integrand = face.getResultVector(
        iSide);  // Integrand value based on n number of testbasisfunctions from
                 // element corresponding to side iSide

    for (std::size_t iFun = 0; iFun < numberOfTestBasisFunctions; ++iFun) {
        for (std::size_t iVar = 0; iVar < numberOfVariables_; ++iVar) {
            std::size_t iVarFun = face.convertToSingleIndex(iSide, iFun, iVar);
            integrand(iVarFun) =
                -flux(iVar) * face.basisFunctionNormal(iSide, iFun)(0);
        }
    }

    return integrand;
}

std::pair<LinearAlgebra::MiddleSizeVector, LinearAlgebra::MiddleSizeVector>
    SavageHutter1DBase::integrandsAtFace(
        Base::PhysicalFace<1> &face, const double &time,
        const LinearAlgebra::MiddleSizeVector &solutionCoefficientsLeft,
        const LinearAlgebra::MiddleSizeVector &solutionCoefficientsRight) {
    // compute numerical solution at the left side and right side of this face
    const LinearAlgebra::MiddleSizeVector solutionLeft =
        Helpers::getSolution<1>(face.getPhysicalElement(Base::Side::LEFT),
                                solutionCoefficientsLeft, numberOfVariables_);
    const LinearAlgebra::MiddleSizeVector solutionRight =
        Helpers::getSolution<1>(face.getPhysicalElement(Base::Side::RIGHT),
                                solutionCoefficientsRight, numberOfVariables_);

    const LinearAlgebra::MiddleSizeVector flux = hllcFlux(
        solutionLeft, solutionRight, face.getUnitNormalVector()[0], face);

    const std::size_t numOfTestBasisFunctionsLeft =
        face.getPhysicalElement(Base::Side::LEFT).getNumberOfBasisFunctions();
    const std::size_t numOfTestBasisFunctionsRight =
        face.getPhysicalElement(Base::Side::RIGHT).getNumberOfBasisFunctions();
    std::pair<LinearAlgebra::MiddleSizeVector, LinearAlgebra::MiddleSizeVector>
        integrands(std::piecewise_construct,
                   std::forward_as_tuple(numberOfVariables_ *
                                         numOfTestBasisFunctionsLeft),
                   std::forward_as_tuple(numberOfVariables_ *
                                         numOfTestBasisFunctionsRight));

    for (std::size_t iVar = 0; iVar < numberOfVariables_; ++iVar) {
        for (std::size_t iFun = 0; iFun < numOfTestBasisFunctionsLeft; ++iFun) {
            const std::size_t iVarFunLeft =
                face.getPhysicalElement(Base::Side::LEFT)
                    .convertToSingleIndex(iFun, iVar);
            integrands.first(iVarFunLeft) =
                -flux(iVar) *
                face.basisFunctionNormal(Base::Side::LEFT, iFun)(0);
        }
        for (std::size_t iFun = 0; iFun < numOfTestBasisFunctionsRight;
             ++iFun) {
            const std::size_t iVarFunRight =
                face.getPhysicalElement(Base::Side::RIGHT)
                    .convertToSingleIndex(iFun, iVar);
            integrands.second(iVarFunRight) =
                -flux(iVar) *
                face.basisFunctionNormal(Base::Side::RIGHT, iFun)(0);
        }
    }
    return integrands;
}

///\details Compute the flux on the boundary.
///\todo split in multiple functions
const LinearAlgebra::MiddleSizeVector
    SavageHutter1DBase::integrandRightHandSideOnRefFace(
        Base::PhysicalFace<1> &face,
        const LinearAlgebra::MiddleSizeVector &solutionCoefficients,
        const double &time) {
    MiddleSizeVector flux(numberOfVariables_);

    const LinearAlgebra::MiddleSizeVector solution =
        Helpers::getSolution<1>(face.getPhysicalElement(Base::Side::LEFT),
                                solutionCoefficients, numberOfVariables_);
    const double u = (solution(0) > dryLimit_) ? solution(1) / solution(0) : 0;

    // outflow
    if (face.getNormalVector()(0) > 0) {
        // solid wall on the end, can be put in as a function of time
        if (false) {
            const LinearAlgebra::MiddleSizeVector reflection =
                LinearAlgebra::MiddleSizeVector({solution[0], -solution[1]});
            flux = hllcFlux(solution, reflection, 1, face);
        } else {
            // subcritical outflow:
            if ((solution[1] / solution[0] <
                 std::sqrt(solution[0] * std::cos(chuteAngle_) * epsilon_))) {
                const double uIn = solution[1] / solution[0];
                const double invariantIn =
                    uIn + 2 * std::sqrt(epsilon_ * std::cos(chuteAngle_) *
                                        solution[0]);
                // const double froudePrescribed = .75;
                // const double hOut = (invariantIn / (2 + froudePrescribed)) *
                // (invariantIn / (2 + froudePrescribed)) / (epsilon_ *
                // std::cos(chuteAngle_));
                double hOut = 1.35;
                const double uOut =
                    invariantIn -
                    2 * std::sqrt(epsilon_ * std::cos(chuteAngle_) * hOut);
                logger(
                    DEBUG, "new h and u and F: %, %, %", hOut, uOut,
                    uOut / std::sqrt(epsilon_ * std::cos(chuteAngle_) * hOut));
                const auto stateNew = MiddleSizeVector({hOut, hOut * uOut});
                flux = hllcFlux(solution, stateNew, 1, face);
            } else {
                // supercritical outflow:
                flux = hllcFlux(solution, solution,
                                face.getUnitNormalVector()[0], face);
            }
        }
    } else  // inflow
    {
        // subcritical inflow:
        if (false &&
            solution[1] / solution[0] <
                std::sqrt(solution[0] * std::cos(chuteAngle_) * epsilon_)) {
            const double hNew = inflowBC_[0];
            const double uNew =
                solution[1] / solution[0] -
                2 * std::sqrt(epsilon_ * std::cos(chuteAngle_)) *
                    (std::sqrt(solution[0]) - std::sqrt(hNew));
            const auto stateNew = MiddleSizeVector({hNew, hNew * uNew});
            flux = hllcFlux(stateNew, solution, 1, face);
        } else {
            flux = hllcFlux(inflowBC_, solution, 1, face);
        }
    }

    const std::size_t numberOfBasisFunctions = face.getNumberOfBasisFunctions();
    LinearAlgebra::MiddleSizeVector integrand(numberOfVariables_ *
                                              numberOfBasisFunctions);

    for (std::size_t iFun = 0; iFun < numberOfBasisFunctions; ++iFun) {
        for (std::size_t iVar = 0; iVar < numberOfVariables_; ++iVar) {
            const std::size_t iVarFun =
                face.convertToSingleIndex(Base::Side::LEFT, iFun, iVar);
            integrand(iVarFun) =
                -flux(iVar) * face.basisFunctionNormal(iFun)[0];
        }
    }

    return integrand;
}

LinearAlgebra::MiddleSizeVector SavageHutter1DBase::localLaxFriedrichsFlux(
    const LinearAlgebra::MiddleSizeVector &numericalSolutionLeft,
    const LinearAlgebra::MiddleSizeVector &numericalSolutionRight,
    Base::PhysicalFace<1> &face) {
    const double hLeft = numericalSolutionLeft(0);
    const double uLeft =
        (hLeft > dryLimit_) ? numericalSolutionLeft(1) / hLeft : 0;

    const double hRight = numericalSolutionRight(0);
    const double uRight =
        (hRight > dryLimit_) ? numericalSolutionRight(1) / hRight : 0;

    const double alpha = std::max(
        std::abs(uLeft) +
            std::sqrt(epsilon_ * std::max(0., numericalSolutionLeft(0))),
        std::abs(uRight) +
            std::sqrt(epsilon_ * std::max(0., numericalSolutionRight(0))));

    logger(DEBUG, "alpha: %", alpha);

    logger(
        DEBUG, "physical fluxes: %, %",
        computePhysicalFlux(numericalSolutionLeft, face.getPointPhysical()),
        computePhysicalFlux(numericalSolutionRight, face.getPointPhysical()));
    const LinearAlgebra::MiddleSizeVector diffSolutions =
        numericalSolutionRight - numericalSolutionLeft;
    const LinearAlgebra::MiddleSizeVector numericalFlux =
        0.5 *
        (computePhysicalFlux(numericalSolutionLeft, face.getPointPhysical()) +
         computePhysicalFlux(numericalSolutionRight, face.getPointPhysical()) -
         alpha * (diffSolutions));

    return numericalFlux;
}

LinearAlgebra::MiddleSizeVector SavageHutter1DBase::hllcFlux(
    const LinearAlgebra::MiddleSizeVector &numericalSolutionLeft,
    const LinearAlgebra::MiddleSizeVector &numericalSolutionRight,
    const double normal, Base::PhysicalFace<1> &face) {
    const double hLeft = numericalSolutionLeft(0);
    const double hRight = numericalSolutionRight(0);

    const double normalSpeedLeft =
        (hLeft > dryLimit_) ? normal * numericalSolutionLeft(1) / hLeft : 0;
    const double normalSpeedRight =
        (hRight > dryLimit_) ? normal * numericalSolutionRight(1) / hRight : 0;

    const LinearAlgebra::MiddleSizeVector fluxNormalLeft =
        normal *
        computePhysicalFlux(numericalSolutionLeft, face.getPointPhysical());
    const LinearAlgebra::MiddleSizeVector fluxNormalRight =
        normal *
        computePhysicalFlux(numericalSolutionRight, face.getPointPhysical());

    const double phaseSpeedLeft =
        std::sqrt(epsilon_ * std::cos(chuteAngle_) * hLeft);
    const double phaseSpeedRight =
        std::sqrt(epsilon_ * std::cos(chuteAngle_) * hRight);

    const double sl = std::min(normalSpeedLeft - phaseSpeedLeft,
                               normalSpeedRight - phaseSpeedRight);
    const double sr = std::max(normalSpeedLeft + phaseSpeedLeft,
                               normalSpeedRight + phaseSpeedRight);

    if (sl >= 0) {
        return fluxNormalLeft;
    }
    if (sr <= 0) {
        return fluxNormalRight;
    }
    return ((sr * fluxNormalLeft - sl * fluxNormalRight +
             sl * sr * (numericalSolutionRight - numericalSolutionLeft)) /
            (sr - sl));
}

double SavageHutter1DBase::computeFrictionCoulomb(
    const LinearAlgebra::MiddleSizeVector &numericalSolution,
    const double &frictionAngle) {
    return std::tan(frictionAngle);
}

/// Compute friction as described in Weinhart et al (2012), eq (50)
/// Notice that the minus sign for gamma in eq (50) is wrong.
double SavageHutter1DBase::computeFriction(
    const LinearAlgebra::MiddleSizeVector &numericalSolution) {
    const double delta1 = 17.518 / 180 * M_PI;
    const double delta2 = 29.712 / 180 * M_PI;
    const double A = 5.29;
    const double beta = 0.189;
    const double gamma = -.080;
    const double d = 2;
    const double h = numericalSolution[0];
    if (h < dryLimit_) return std::tan(delta1);
    const double u = numericalSolution[1] / h;
    const double F = u / std::sqrt(epsilon_ * std::cos(chuteAngle_) * h);
    return std::tan(delta1) + (std::tan(delta2) - std::tan(delta1)) /
                                  (beta * h / (A * d * (F - gamma)) + 1);
}

/// Compute friction as in Anthony's draft, eq (3.9)
double SavageHutter1DBase::computeFrictionExponential(
    const LinearAlgebra::MiddleSizeVector &numericalSolution) {
    const double delta1 = 27. / 180 * M_PI;
    const double delta2 = 37. / 180 * M_PI;
    const double beta = 0.136;
    const double L = 2.;
    const double h = numericalSolution[0];
    if (h < dryLimit_) return std::tan(delta1);
    const double u = numericalSolution[1] / h;
    if (std::abs(u) < 1e-16) return std::tan(delta1);
    return std::tan(delta1) + (std::tan(delta2) - std::tan(delta1)) *
                                  std::exp(-beta * std::pow(epsilon_ * h, 1.5) /
                                           (L * std::abs(u)));
}
