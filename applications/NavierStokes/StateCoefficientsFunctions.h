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

#ifndef HPGEM_APP_STATECOEFFICIENTSFUNCTIONS_H
#define HPGEM_APP_STATECOEFFICIENTSFUNCTIONS_H

#include "Base/PhysicalElement.h"
#include "Base/PhysicalFace.h"
#include "LinearAlgebra/MiddleSizeVector.h"
#include "LinearAlgebra/MiddleSizeMatrix.h"

/// *********************************************************
/// ***      General Solution Coefficient Functions       ***
/// *********************************************************
/// The general solution coefficient functions are required for every integrand
/// calculation \brief Compute state at an element

template <std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
LinearAlgebra::MiddleSizeVector computeStateOnElement(
    Base::PhysicalElement<DIM> &element,
    const LinearAlgebra::MiddleSizeVector &stateCoefficients) {
    std::size_t numberOfBasisFunctions = element.getNumOfBasisFunctions();
    LinearAlgebra::MiddleSizeVector elementState(NUMBER_OF_VARIABLES);
    std::size_t iVB;  // Index in state coefficients for variable i and
                      // basisfunction j
    // todo: check if this function can be made more efficient
    for (std::size_t iV = 0; iV < NUMBER_OF_VARIABLES; iV++) {
        elementState(iV) = 0.0;
        for (std::size_t iB = 0; iB < numberOfBasisFunctions; iB++) {
            iVB = element.convertToSingleIndex(iB, iV);
            elementState(iV) +=
                stateCoefficients(iVB) * element.basisFunction(iB);
        }
    }

    return elementState;
}

/// \brief Compute state derivatives at an element.
template <std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
LinearAlgebra::MiddleSizeMatrix computeStateJacobianAtElement(
    Base::PhysicalElement<DIM> &element,
    const LinearAlgebra::MiddleSizeVector &stateCoefficients) {
    std::size_t numberOfBasisFunctions = element.getNumOfBasisFunctions();
    std::size_t iVB;  // Index in state coefficients for variable i and
                      // basisfunction j

    LinearAlgebra::MiddleSizeMatrix stateJacobian(NUMBER_OF_VARIABLES, DIM);
    LinearAlgebra::SmallVector<DIM> gradientBasisFunction;

    for (std::size_t iB = 0; iB < numberOfBasisFunctions; iB++) {
        gradientBasisFunction = element.basisFunctionDeriv(iB);
        for (std::size_t iV = 0; iV < NUMBER_OF_VARIABLES; iV++) {
            iVB = element.convertToSingleIndex(iB, iV);
            for (std::size_t iD = 0; iD < DIM; iD++) {
                stateJacobian(iV, iD) +=
                    stateCoefficients(iVB) * gradientBasisFunction(iD);
            }
        }
    }

    return stateJacobian;
}

/// \brief Compute state at a face
template <std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
LinearAlgebra::MiddleSizeVector computeStateOnFace(
    Base::PhysicalFace<DIM> &face, const Base::Side &iSide,
    const LinearAlgebra::MiddleSizeVector &stateCoefficients) {
    std::size_t numOfBasisFunctions =
        face.getPhysicalElement(iSide).getNumOfBasisFunctions();
    LinearAlgebra::MiddleSizeVector state(NUMBER_OF_VARIABLES);
    std::size_t iVB;  // Index in solution coefficients for variable i and
                      // basisfunction j

    for (std::size_t iV = 0; iV < NUMBER_OF_VARIABLES; iV++) {
        for (std::size_t iB = 0; iB < numOfBasisFunctions; iB++) {
            iVB = face.getPhysicalElement(iSide).convertToSingleIndex(iB, iV);
            state(iV) +=
                stateCoefficients(iVB) *
                face.basisFunction(iSide,
                                   iB);  // basisFunction returns physical value
        }
    }

    return state;
}

/// \brief Compute derivatives of the state with respect to the coordinates, at
/// a face
template <std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
LinearAlgebra::MiddleSizeMatrix computeStateJacobianAtFace(
    Base::PhysicalFace<DIM> &face, const Base::Side &iSide,
    const LinearAlgebra::MiddleSizeVector &stateCoefficients) {
    std::size_t numberOfBasisFunctions =
        face.getPhysicalElement(iSide).getNumOfBasisFunctions();
    std::size_t iVB;  // Index in state coefficients for variable i and
                      // basisfunction j

    LinearAlgebra::MiddleSizeMatrix stateJacobian(NUMBER_OF_VARIABLES, DIM);
    LinearAlgebra::SmallVector<DIM> gradientBasisFunction;

    for (std::size_t iB = 0; iB < numberOfBasisFunctions; iB++) {
        gradientBasisFunction = face.basisFunctionDeriv(iSide, iB);
        for (std::size_t iV = 0; iV < NUMBER_OF_VARIABLES; iV++) {
            iVB = face.getPhysicalElement(iSide).convertToSingleIndex(iB, iV);
            for (std::size_t iD = 0; iD < DIM; iD++) {
                stateJacobian(iV, iD) +=
                    stateCoefficients(iVB) * gradientBasisFunction(iD);
            }
        }
    }

    return stateJacobian;
}

/// \brief Computes the partial states: all states except the density are
/// divided by the density
template <std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
LinearAlgebra::MiddleSizeVector computePartialState(
    const LinearAlgebra::MiddleSizeVector &state) {
    LinearAlgebra::MiddleSizeVector partialState(NUMBER_OF_VARIABLES);
    double q1Inverse = 1.0 / state(0);

    partialState(0) = state(0);

    for (std::size_t iD = 0; iD < NUMBER_OF_VARIABLES - 1; iD++) {
        partialState(iD + 1) = state(iD + 1) * q1Inverse;
    }
    return partialState;
}

/// ************************************************************
/// ***      Additional Solution Coefficient Functions       ***
/// ************************************************************
/// These functions generally do not have to be computed everytime (or do they?)

/// \brief Computes the Jacobian of the velocity matrix
template <std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
LinearAlgebra::MiddleSizeMatrix computeVelocityJacobian(
    const LinearAlgebra::MiddleSizeVector &partialState,
    const LinearAlgebra::MiddleSizeMatrix &stateJacobian) {

    LinearAlgebra::MiddleSizeMatrix velocityJacobian(DIM, DIM);
    // First put in the velocity Jacobian matrix L
    double inverseDensity = 1.0 / partialState(0);
    for (std::size_t iD1 = 0; iD1 < DIM; iD1++) {
        for (std::size_t iD2 = 0; iD2 < DIM; iD2++) {
            velocityJacobian(iD1, iD2) =
                inverseDensity *
                (stateJacobian(iD1 + 1, iD2) -
                 partialState(iD1 + 1) * stateJacobian(0, iD2));
        }
    }
    return velocityJacobian;
}

/// \brief This function computes the state on a face given a matrix of
/// stateCoefficients. Implemented for efficiency reasons when calculating
/// stabilityparameters
template <std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
LinearAlgebra::MiddleSizeMatrix computeStateOnFace(
    Base::PhysicalFace<DIM> &face, const Base::Side &iSide,
    const LinearAlgebra::MiddleSizeMatrix &stateCoefficients) {
    std::size_t numOfBasisFunctions =
        face.getPhysicalElement(iSide).getNumOfBasisFunctions();
    LinearAlgebra::MiddleSizeMatrix state(NUMBER_OF_VARIABLES, DIM);
    std::size_t iVB;  // Index in solution coefficients for variable i and
                      // basisfunction j

    for (std::size_t iV = 0; iV < NUMBER_OF_VARIABLES; iV++) {
        for (std::size_t iD = 0; iD < DIM; iD++) {
            for (std::size_t iB = 0; iB < numOfBasisFunctions; iB++) {
                iVB =
                    face.getPhysicalElement(iSide).convertToSingleIndex(iB, iV);
                state(iV, iD) +=
                    stateCoefficients(iVB, iD) *
                    face.basisFunction(
                        iSide, iB);  // basisFunction returns physical value
                /*						std::cout <<
                   "stateCoefficient[ " << iVB << "," << iD << "]: " <<
                   stateCoefficients(iVB,iD) << std::endl; std::cout <<
                   "basisFunction[" << iB << "]: " << face.basisFunction(iSide,
                   iB) << std::endl;*/
            }
        }
    }
    return state;
}

// A_ikrs = A_(iV)(iD)(iVm)(iDm)
// This is a naive implementation, but can be used for debugging fast
// implementations
template <std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
LinearAlgebra::MiddleSizeMatrix computeEllipticTensorMatrixContraction(
    const std::vector<LinearAlgebra::MiddleSizeMatrix> &ellipticTensor,
    const LinearAlgebra::MiddleSizeMatrix &matrix) {
    LinearAlgebra::MiddleSizeMatrix result(NUMBER_OF_VARIABLES, DIM);
    double pos;

    for (std::size_t iD = 0; iD < DIM; iD++) {
        for (std::size_t iV = 0; iV < NUMBER_OF_VARIABLES; iV++) {
            for (std::size_t iDm = 0; iDm < DIM; iDm++) {
                for (std::size_t iVm = 0; iVm < NUMBER_OF_VARIABLES; iVm++) {
                    pos = (DIM)*iD + iDm;
                    result(iV, iD) +=
                        ellipticTensor[pos](iV, iVm) * matrix(iVm, iDm);
                }
            }
        }
    }

    return result;
}

#endif // HPGEM_APP_STATECOEFFICIENTSFUNCTIONS_H
