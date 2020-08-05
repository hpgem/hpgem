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

#include "HpgemAPINonLinearSteadyState.h"
#include "Utilities/GlobalVector.h"

#include "Logger.h"
#include <iostream>  //remove this after this API is official
#include <fstream>   //remove this after this API is official
namespace hpgem {
namespace Base {
// todo: Possibly add restart from data file to create a p-multigrid solution
// strategy todo: Add function that can add setup functions to the solve
// function

// note: compute error is always on, compute both faces is always on
template <std::size_t DIM>
HpgemAPINonLinearSteadyState<DIM>::HpgemAPINonLinearSteadyState(
    const std::size_t numberOfVariables, const std::size_t polynomialOrder,
    const bool computeBothFaces)
    : HpgemAPISimplified<DIM>(numberOfVariables, polynomialOrder, 4, 0,
                              computeBothFaces),
      sourceElementVectorID_(0),
      sourceFaceVectorID_(0) {}

template <std::size_t DIM>
void HpgemAPINonLinearSteadyState<DIM>::computeJacobian() {
    // Compute new Jacobian element matrices for J(u)
    for (Base::Element *ptrElement : this->meshes_[0]->getElementsList()) {
        // Get the solution coefficients required for the calculation
        LinearAlgebra::MiddleSizeVector &solutionCoefficients(
            ptrElement->getTimeIntegrationVector(this->solutionVectorId_));
        // Calculate and store the matrix
        ptrElement->setElementMatrix(
            computeJacobianAtElement(ptrElement, solutionCoefficients, 0),
            jacobianElementMatrixID_);
    }

    // Compute the local Jacobian Matrix: Face integrals
    for (Base::Face *ptrFace : this->meshes_[0]->getFacesList()) {
        // grab the coefficients
        LinearAlgebra::MiddleSizeVector solutionCoefficientsLeft(
            ptrFace->getPtrElementLeft()->getTimeIntegrationVector(
                this->solutionVectorId_));
        LinearAlgebra::MiddleSizeVector solutionCoefficientsRight(
            ptrFace->getPtrElementRight()->getTimeIntegrationVector(
                this->solutionVectorId_));

        // Calculate face Matrix
        int numberOfDOFLeft =
            ptrFace->getPtrElement(Base::Side::LEFT)
                ->getNumberOfBasisFunctions() *
            ptrFace->getPtrElement(Base::Side::LEFT)->getNumberOfUnknowns();
        int numberOfDOFRight =
            ptrFace->getPtrElement(Base::Side::RIGHT)
                ->getNumberOfBasisFunctions() *
            ptrFace->getPtrElement(Base::Side::RIGHT)->getNumberOfUnknowns();
        Base::FaceMatrix faceMatrix(numberOfDOFLeft, numberOfDOFRight);
        // Left element. Side indexes are: elementSide, derivativeSide
        faceMatrix.setElementMatrix(
            computeJacobianAtFace(ptrFace, solutionCoefficientsLeft,
                                  solutionCoefficientsRight, 0,
                                  Base::Side::LEFT, Base::Side::LEFT),
            Base::Side::LEFT, Base::Side::LEFT);
        faceMatrix.setElementMatrix(
            computeJacobianAtFace(ptrFace, solutionCoefficientsLeft,
                                  solutionCoefficientsRight, 0,
                                  Base::Side::LEFT, Base::Side::RIGHT),
            Base::Side::LEFT, Base::Side::RIGHT);
        // Right element
        faceMatrix.setElementMatrix(
            computeJacobianAtFace(ptrFace, solutionCoefficientsLeft,
                                  solutionCoefficientsRight, 0,
                                  Base::Side::RIGHT, Base::Side::LEFT),
            Base::Side::RIGHT, Base::Side::LEFT);
        faceMatrix.setElementMatrix(
            computeJacobianAtFace(ptrFace, solutionCoefficientsLeft,
                                  solutionCoefficientsRight, 0,
                                  Base::Side::RIGHT, Base::Side::RIGHT),
            Base::Side::RIGHT, Base::Side::RIGHT);

        // Store face Matrix
        ptrFace->setFaceMatrix(faceMatrix, jacobianFaceMatrixID_);  // time=0
    }
}

template <std::size_t DIM>
bool HpgemAPINonLinearSteadyState<DIM>::solve(bool doComputeInitialCondition,
                                              bool doComputeError,
                                              bool doUseJacobian) {

    logger(ERROR,
           "Sundials is required to solve the non linear steady state problem "
           "using this function.");
    return false;
}
}  // namespace Base
}  // namespace hpgem