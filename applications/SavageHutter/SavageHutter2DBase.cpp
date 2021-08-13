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

#include "SavageHutter2DBase.h"
#include "HelperFunctions.h"

using namespace hpgem;

const LinearAlgebra::MiddleSizeVector
    SavageHutter2DBase::integrandRightHandSideOnElement(
        Base::PhysicalElement<2>& element, const double& time,
        const LinearAlgebra::MiddleSizeVector& solutionCoefficients) {
    const std::size_t numBasisFuncs =
        element.getElement()->getNumberOfBasisFunctions();

    MiddleSizeVector& integrand =
        element.getResultVector();  // just to have the correct length
    const PointPhysicalT& pPhys = element.getPointPhysical();
    const LinearAlgebra::MiddleSizeVector numericalSolution =
        Helpers::getSolution<2>(element, solutionCoefficients,
                                numberOfVariables_);
    const LinearAlgebra::MiddleSizeVector physicalFlux =
        computePhysicalFlux(numericalSolution);
    const LinearAlgebra::MiddleSizeVector source =
        computeSourceTerm(numericalSolution, pPhys, time);

    // Compute integrand on the physical element.
    std::size_t iVB;  // Index for both basis function and variable
    for (std::size_t iB = 0; iB < numBasisFuncs;
         iB++)  // Index for basis function
    {
        for (std::size_t iV = 0; iV < numberOfVariables_; ++iV) {
            iVB = element.getElement()->convertToSingleIndex(iB, iV);
            integrand(iVB) =
                physicalFlux(2 * iV) * element.basisFunctionDeriv(iB)(0);
            integrand(iVB) +=
                physicalFlux(2 * iV + 1) * element.basisFunctionDeriv(iB)(1);
            integrand(iVB) += source(iV) * element.basisFunction(iB);
        }
    }
    return integrand;
}

const LinearAlgebra::MiddleSizeVector
    SavageHutter2DBase::integrandRightHandSideOnRefFace(
        Base::PhysicalFace<2>& face, const Base::Side& iSide,
        const MiddleSizeVector& solutionCoefficientsLeft,
        const MiddleSizeVector& solutionCoefficientsRight) {
    const std::size_t numTestBasisFuncs =
        face.getFace()->getPtrElement(iSide)->getNumberOfBasisFunctions();

    // compute numerical solution at the left side and right side of this face
    const PointReferenceOnFaceT& pRef = face.getPointReference();
    LinearAlgebra::MiddleSizeVector solutionLeft =
        Helpers::getSolution<2>(face.getPhysicalElement(Base::Side::LEFT),
                                solutionCoefficientsLeft, numberOfVariables_);
    LinearAlgebra::MiddleSizeVector solutionRight =
        Helpers::getSolution<2>(face.getPhysicalElement(Base::Side::RIGHT),
                                solutionCoefficientsRight, numberOfVariables_);

    LinearAlgebra::MiddleSizeVector flux =
        hllcFlux(solutionLeft, solutionRight, face.getUnitNormalVector());

    if (iSide == Base::Side::RIGHT)  // the normal is defined for the left
                                     // element
    {
        flux *= -1;
    }

    LinearAlgebra::MiddleSizeVector& integrand = face.getResultVector(
        iSide);  // Integrand value based on n number of testbasisfunctions from
                 // element corresponding to side iSide

    for (std::size_t iFun = 0; iFun < numTestBasisFuncs; ++iFun) {
        for (std::size_t iVar = 0; iVar < numberOfVariables_; ++iVar) {
            std::size_t iVarFun =
                face.getFace()->getPtrElement(iSide)->convertToSingleIndex(
                    iFun, iVar);
            integrand(iVarFun) = -flux(iVar) * face.basisFunction(iSide, iFun);
        }
    }
    return integrand;
}

std::pair<LinearAlgebra::MiddleSizeVector, LinearAlgebra::MiddleSizeVector>
    SavageHutter2DBase::integrandsAtFace(
        Base::PhysicalFace<2>& face, const double& time,
        const LinearAlgebra::MiddleSizeVector& solutionCoefficientsLeft,
        const LinearAlgebra::MiddleSizeVector& solutionCoefficientsRight) {

    // compute numerical solution at the left side and right side of this face
    LinearAlgebra::MiddleSizeVector solutionLeft =
        Helpers::getSolution<2>(face.getPhysicalElement(Base::Side::LEFT),
                                solutionCoefficientsLeft, numberOfVariables_);
    LinearAlgebra::MiddleSizeVector solutionRight =
        Helpers::getSolution<2>(face.getPhysicalElement(Base::Side::RIGHT),
                                solutionCoefficientsRight, numberOfVariables_);

    LinearAlgebra::MiddleSizeVector flux =
        hllcFlux(solutionLeft, solutionRight, face.getUnitNormalVector());

    std::size_t numOfTestBasisFunctionsLeft =
        face.getPhysicalElement(Base::Side::LEFT).getNumberOfBasisFunctions();
    std::size_t numOfTestBasisFunctionsRight =
        face.getPhysicalElement(Base::Side::RIGHT).getNumberOfBasisFunctions();
    std::pair<LinearAlgebra::MiddleSizeVector, LinearAlgebra::MiddleSizeVector>
        integrands(std::piecewise_construct,
                   std::forward_as_tuple(numberOfVariables_ *
                                         numOfTestBasisFunctionsLeft),
                   std::forward_as_tuple(numberOfVariables_ *
                                         numOfTestBasisFunctionsRight));

    for (std::size_t iVar = 0; iVar < numberOfVariables_; ++iVar) {
        for (std::size_t iFun = 0; iFun < numOfTestBasisFunctionsLeft; ++iFun) {
            std::size_t iVarFun = face.getFace()
                                      ->getPtrElement(Base::Side::LEFT)
                                      ->convertToSingleIndex(iFun, iVar);
            integrands.first(iVarFun) =
                -flux(iVar) * face.basisFunction(Base::Side::LEFT, iFun);
        }
        for (std::size_t iFun = 0; iFun < numOfTestBasisFunctionsRight;
             ++iFun) {
            std::size_t iVarFun = face.getFace()
                                      ->getPtrElement(Base::Side::RIGHT)
                                      ->convertToSingleIndex(iFun, iVar);
            integrands.second(iVarFun) =
                flux(iVar) * face.basisFunction(Base::Side::RIGHT, iFun);
        }
    }
    return integrands;
}

const LinearAlgebra::MiddleSizeVector
    SavageHutter2DBase::integrandRightHandSideOnRefFace(
        Base::PhysicalFace<2>& face,
        const LinearAlgebra::MiddleSizeVector& solutionCoefficients,
        const double& time) {
    double normalX = face.getUnitNormalVector()[0];
    double normalY = face.getUnitNormalVector()[1];
    const std::size_t numberOfBasisFunctions =
        face.getFace()->getNumberOfBasisFunctions();

    const PointReferenceOnFaceT& pRef = face.getPointReference();
    // note that at the boundary, the element is the left element by definition
    LinearAlgebra::MiddleSizeVector solution =
        Helpers::getSolution<2>(face.getPhysicalElement(Base::Side::LEFT),
                                solutionCoefficients, numberOfVariables_);

    LinearAlgebra::MiddleSizeVector flux(3);

    // outflow
    if (normalX > .5 && std::abs(normalY) < 1e-5) {
        if (solution[1] / solution[0] <
            std::sqrt(epsilon_ * std::cos(chuteAngle_) *
                      solution[0]))  // subcritical
        {
            logger(DEBUG, "subcritical outflow");
            double uIn = solution[1] / solution[0];
            double invariantIn =
                uIn +
                2 * std::sqrt(epsilon_ * std::cos(chuteAngle_) * solution[0]);
            double pseudoFroude = 1;
            // double hOut = (invariantIn / (2+pseudoFroude)) * (invariantIn /
            // (2+pseudoFroude)) / (epsilon_ * std::cos(chuteAngle_)); double
            // dischargeOut = 1; double hTest =
            // solveCubicEquation(4*epsilon_*std::cos(chuteAngle_), -invariantIn
            // * invariantIn, 2*invariantIn*dischargeOut,
            // -dischargeOut*dischargeOut);
            double hOut = 3.772;
            double uOut =
                invariantIn -
                2 * std::sqrt(epsilon_ * std::cos(chuteAngle_) * hOut);
            auto stateNew = MiddleSizeVector(
                {hOut, hOut * uOut, solution[2]});  // keep hv continuous
            flux = hllcFlux(solution, stateNew, face.getUnitNormalVector());
        } else {
            flux = hllcFlux(solution, solution, face.getUnitNormalVector());
        }
    } else if (std::abs(normalY) < 1e-5)  // inflow
    {
        if (solution[1] / solution[0] <
            std::sqrt(epsilon_ * std::cos(chuteAngle_) *
                      solution[0]))  // subcritical
        {
            double hNew = inflowBC_[0];
            double uNew = solution[1] / solution[0] -
                          2 * std::sqrt(epsilon_ * std::cos(chuteAngle_)) *
                              (std::sqrt(solution[0]) - std::sqrt(hNew));
            auto stateNew = MiddleSizeVector({hNew, hNew * uNew, solution[2]});
            flux = hllcFlux(solution, stateNew, face.getUnitNormalVector());
        } else {
            flux = hllcFlux(inflowBC_, inflowBC_, face.getUnitNormalVector());
        }
    } else  // solid wall
    {
        LinearAlgebra::SmallVector<2> velocity(
            {solution[1] / solution[0], solution[2] / solution[0]});
        LinearAlgebra::SmallVector<2> velocityReflected =
            velocity - 2 * (velocity * face.getUnitNormalVector()) *
                           face.getUnitNormalVector();
        const LinearAlgebra::MiddleSizeVector& reflection =
            MiddleSizeVector({solution[0], velocityReflected[0] * solution[0],
                              velocityReflected[1] * solution[0]});
        flux = hllcFlux(solution, reflection, face.getUnitNormalVector());
    }

    LinearAlgebra::MiddleSizeVector integrand(numberOfVariables_ *
                                              numberOfBasisFunctions);

    for (std::size_t iFun = 0; iFun < numberOfBasisFunctions; ++iFun) {
        for (std::size_t iVar = 0; iVar < numberOfVariables_; ++iVar) {
            std::size_t iVarFun =
                face.getFace()->getPtrElementLeft()->convertToSingleIndex(iFun,
                                                                          iVar);
            integrand(iVarFun) = -flux(iVar) * face.basisFunction(iFun);
        }
    }

    return integrand;
}

LinearAlgebra::MiddleSizeVector SavageHutter2DBase::localLaxFriedrichsFlux(
    const LinearAlgebra::MiddleSizeVector& numericalSolutionLeft,
    const LinearAlgebra::MiddleSizeVector& numericalSolutionRight,
    const LinearAlgebra::SmallVector<2>& normal) {
    logger.assert_debug(std::abs(normal.l2Norm() - 1) < 1e-16,
                        "LLF flux needs a unit normal vector");
    const double hLeft = numericalSolutionLeft[0];
    const double hRight = numericalSolutionRight[0];
    double uLeft = 0;
    double vLeft = 0;
    if (hLeft > dryLimit_) {
        uLeft = numericalSolutionLeft(1) / hLeft;
        vLeft = numericalSolutionLeft(2) / hLeft;
    }

    double uRight = 0;
    double vRight = 0;
    if (hRight > dryLimit_) {
        uRight = numericalSolutionRight(1) / hRight;
        vRight = numericalSolutionRight(2) / hRight;
    }

    // take the maximum of |u+-sqrt(epsilon cos(theta) h)| and |v+-sqrt(epsilon
    // cos(theta) h)| for left and right side
    const double eigenSpeed1 = std::max(
        std::abs(uLeft + std::sqrt(epsilon_ * std::cos(chuteAngle_) * hLeft)),
        std::abs(uLeft - std::sqrt(epsilon_ * std::cos(chuteAngle_) * hLeft)));
    const double eigenSpeed2 = std::max(
        std::abs(vLeft + std::sqrt(epsilon_ * std::cos(chuteAngle_) * hLeft)),
        std::abs(vLeft - std::sqrt(epsilon_ * std::cos(chuteAngle_) * hLeft)));
    const double eigenSpeed3 = std::max(
        std::abs(uRight + std::sqrt(epsilon_ * std::cos(chuteAngle_) * hRight)),
        std::abs(uRight -
                 std::sqrt(epsilon_ * std::cos(chuteAngle_) * hRight)));
    const double eigenSpeed4 = std::max(
        std::abs(vRight + std::sqrt(epsilon_ * std::cos(chuteAngle_) * hRight)),
        std::abs(vRight -
                 std::sqrt(epsilon_ * std::cos(chuteAngle_) * hRight)));

    // alpha is the maximum eigenvalue of the system
    double alpha = std::max(eigenSpeed1, eigenSpeed2);
    alpha = std::max(alpha, eigenSpeed3);
    alpha = std::max(alpha, eigenSpeed4);

    LinearAlgebra::MiddleSizeVector fluxLeft =
        computePhysicalFlux(numericalSolutionLeft);
    LinearAlgebra::MiddleSizeVector fluxRight =
        computePhysicalFlux(numericalSolutionRight);
    LinearAlgebra::MiddleSizeVector fluxNormalLeft(numberOfVariables_);
    LinearAlgebra::MiddleSizeVector fluxNormalRight(numberOfVariables_);
    for (std::size_t i = 0; i < numberOfVariables_; ++i) {
        fluxNormalLeft(i) =
            fluxLeft(2 * i) * normal(0) + fluxLeft(2 * i + 1) * normal(1);
        fluxNormalRight(i) =
            fluxRight(2 * i) * normal(0) + fluxRight(2 * i + 1) * normal(1);
    }

    const LinearAlgebra::MiddleSizeVector numericalFlux =
        0.5 * (fluxNormalLeft + fluxNormalRight -
               alpha * (numericalSolutionRight - numericalSolutionLeft));

    return numericalFlux;
}

/// The HLLC flux is an approximate Riemann solver
LinearAlgebra::MiddleSizeVector SavageHutter2DBase::hllcFlux(
    const LinearAlgebra::MiddleSizeVector& numericalSolutionLeft,
    const LinearAlgebra::MiddleSizeVector& numericalSolutionRight,
    const LinearAlgebra::SmallVector<2>& normal) {
    logger.assert_debug(std::abs(normal.l2Norm() - 1) < 1e-14,
                        "hllc flux needs unit normal vector");
    const double nx = normal[0];
    const double ny = normal[1];
    const double hLeft = numericalSolutionLeft[0];
    const double hRight = numericalSolutionRight[0];
    double uLeft = 0;
    double vLeft = 0;
    if (hLeft > dryLimit_) {
        uLeft = numericalSolutionLeft(1) / hLeft;
        vLeft = numericalSolutionLeft(2) / hLeft;
    }
    double uRight = 0;
    double vRight = 0;
    if (hRight > dryLimit_) {
        uRight = numericalSolutionRight(1) / hRight;
        vRight = numericalSolutionRight(2) / hRight;
    }
    LinearAlgebra::MiddleSizeVector fluxLeft =
        computePhysicalFlux(numericalSolutionLeft);
    LinearAlgebra::MiddleSizeVector fluxRight =
        computePhysicalFlux(numericalSolutionRight);
    LinearAlgebra::MiddleSizeVector fluxNormalLeft(numberOfVariables_);
    LinearAlgebra::MiddleSizeVector fluxNormalRight(numberOfVariables_);
    for (std::size_t i = 0; i < numberOfVariables_; ++i) {
        fluxNormalLeft(i) = fluxLeft(2 * i) * nx + fluxLeft(2 * i + 1) * ny;
        fluxNormalRight(i) = fluxRight(2 * i) * nx + fluxRight(2 * i + 1) * ny;
    }
    double normalSpeedLeft = uLeft * nx + vLeft * ny;
    double normalSpeedRight = uRight * nx + vRight * ny;
    double phaseSpeedLeft = std::sqrt(epsilon_ * std::cos(chuteAngle_) * hLeft);
    double phaseSpeedRight =
        std::sqrt(epsilon_ * std::cos(chuteAngle_) * hRight);
    double sl = std::min(normalSpeedLeft - phaseSpeedLeft,
                         normalSpeedRight - phaseSpeedRight);
    double sr = std::max(normalSpeedLeft + phaseSpeedLeft,
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

double SavageHutter2DBase::computeFriction(
    const LinearAlgebra::MiddleSizeVector& numericalSolution) {
    const double delta1 = 17.518 / 180 * M_PI;
    const double delta2 = 29.712 / 180 * M_PI;
    const double h = numericalSolution[0];
    if (h < dryLimit_) return std::tan(delta1);
    const double u = numericalSolution[1] / h;
    const double v = numericalSolution[2] / h;
    const double froude = std::sqrt(u * u + v * v) /
                          std::sqrt(epsilon_ * std::cos(chuteAngle_) * h);
    const double A = 5.29;
    const double beta = 0.189;
    const double gamma = -0.080;
    const double d = 2;

    // return std::tan(45./180*M_PI);
    return std::tan(delta1) + (std::tan(delta2) - std::tan(delta1)) /
                                  (beta * h / (A * d * (froude - gamma)) + 1);
}

///\details function that computes the width-average of the solution by simply
/// adding
/// the values of all elements and then divide by the number of nodes that have
/// been added. Finally, each point is multiplied by the width of the chute at
/// that point.
// Sorry for the ugliness...
std::vector<std::pair<double, LinearAlgebra::MiddleSizeVector>>
    SavageHutter2DBase::widthAverage() {

    // Since we're using a rectangular grid, we can get the number of nodes in x
    // direction and y direction The number of elements in x direction is given
    // by the user in the commandline. extern
    // Base::CommandLineOption<std::size_t>& numberOfElements; const std::size_t
    // nodesInXDirection = numberOfElements.getValue() + 1;
    const std::size_t nodesInXDirection = 21;
    const std::size_t elementsInYDirection =
        this->meshes_[0]->getNumberOfElements(Base::IteratorType::GLOBAL) /
        (nodesInXDirection - 1);
    logger(DEBUG, "elements in y direction: % ", elementsInYDirection);

    // make xs
    ///\todo insert length of the domain here automatically instead of hardcoded
    const double dx = 11. / (nodesInXDirection - 1);
    std::vector<std::pair<double, LinearAlgebra::MiddleSizeVector>> totals;
    for (std::size_t i = 0; i < nodesInXDirection; ++i) {
        totals.push_back(std::make_pair(
            i * dx, LinearAlgebra::MiddleSizeVector(numberOfVariables_)));
    }
    logger(DEBUG, "size of totals: %", totals.size());

    // add all values at a certain x-coordinate. To do that, first check if this
    // value of x is already in the vector. If not, make a pair of this x-value
    // and the value of the variables if there was already an entry for this x,
    // add the value of the current point
    for (Base::Element* element : this->meshes_[0]->getElementsList()) {
        const Geometry::ReferenceGeometry* referenceElement =
            element->getReferenceGeometry();
        for (std::size_t i = 0; i < referenceElement->getNumberOfNodes(); ++i) {
            const Geometry::PointReference<2>& nodeReference =
                referenceElement->getReferenceNodeCoordinate(i);
            const LinearAlgebra::MiddleSizeVector& value =
                element->getSolution(0, nodeReference);
            logger(DEBUG, "value: %", value);
            const Geometry::PointPhysical<2>& node =
                element->referenceToPhysical(nodeReference);

            const auto xPosInVector = std::find_if(
                totals.begin(), totals.end(),
                [=](const std::pair<double, LinearAlgebra::MiddleSizeVector>
                        current) {
                    return std::abs(current.first - (node)[0]) < 1e-10;
                });
            if (xPosInVector == totals.end()) {
                totals.push_back(std::make_pair(node[0], value));
                logger(WARN, "x = % not found", node[0]);
            } else {
                (*xPosInVector).second += value;
            }
        }
    }
#ifdef HPGEM_USE_MPI
    int world_rank = Base::MPIContainer::Instance().getProcessorID();
    auto& comm = Base::MPIContainer::Instance().getComm();

    // split the pairs, since MPI can't send over a vector of pairs
    std::vector<LinearAlgebra::MiddleSizeVector> solutions;
    solutions.reserve(totals.size());
    for (std::pair<double, LinearAlgebra::MiddleSizeVector> p : totals) {
        solutions.push_back(p.second);
    }
    logger(DEBUG, "solutions size: %", solutions.size());

    std::vector<LinearAlgebra::MiddleSizeVector> globalSolutions(
        solutions.size(), LinearAlgebra::MiddleSizeVector(numberOfVariables_));
    logger(DEBUG, "global solutions size %", globalSolutions.size());
    for (std::size_t i = 0; i < solutions.size(); ++i) {
        LinearAlgebra::MiddleSizeVector v = solutions[i];
        logger(DEBUG, "v pointer: %, v: %", v.data(), v);
        logger(DEBUG, "pointer to global solutions: %",
               globalSolutions[i].data());
        MPI_Reduce(v.data(), globalSolutions[i].data(), v.size(),
                   Base::Detail::toMPIType(*v.data()), MPI_SUM, 0, comm);
    }
    logger(DEBUG, "I'm still here!");

    if (world_rank == 0) {
        logger.assert_debug(totals.size() == globalSolutions.size(),
                            "wrong size");
        for (std::size_t i = 0; i < globalSolutions.size(); ++i) {
            totals[i].second = globalSolutions[i];
        }
#endif
        // divide by the number of times a value for the given x-point is added,
        // which is 2 times the number of elements in y-direction for boundary
        // nodes and 4 times the number of elements in y-direction otherwise
        (*totals.begin()).second *= 2;
        (totals.back()).second *= 2;
        for (std::pair<double, LinearAlgebra::MiddleSizeVector>& val : totals) {
            val.second /= 4 * elementsInYDirection;
            logger(INFO, "x: %, average: %", val.first, val.second);
        }
        // just make sure we did not forget any points or that points that are
        // the same have been put in different rows
        logger.assert_debug(
            totals.size() == nodesInXDirection,
            "wrong number of points in vector of width-averaging");
#ifdef HPGEM_USE_MPI
    }
#endif
    return totals;
}
