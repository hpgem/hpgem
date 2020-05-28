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

#include "AirParameters.h"
#include "../StateCoefficientsFunctions.h"
#include "StateCoefficientsStructAir.h"
#include <cmath>

StateCoefficientsStructAir::StateCoefficientsStructAir(
    Base::PhysicalElement<DIM> &element,
    const LinearAlgebra::MiddleSizeVector &stateCoefficients, const double time)
    : StateCoefficientsStruct<DIM, NUMBER_OF_VARIABLES>(
          element, stateCoefficients, time) {
    pressure_ = this->computePressure(state_);
    speedOfSound_ = this->computeSpeedOfSound(state_, pressure_);
    viscosity_ = this->computeViscosity(state_, partialState_, stateJacobian_,
                                        pressure_);
    hyperbolicMatrix_ = this->computeHyperbolicMatrix(state_, pressure_);
    ellipticTensor_ = this->computeEllipticTensor(partialState_, viscosity_);
}

StateCoefficientsStructAir::StateCoefficientsStructAir(
    Base::PhysicalFace<DIM> &face,
    const LinearAlgebra::MiddleSizeVector &stateCoefficients,
    const Base::Side side, const double time)
    : StateCoefficientsStruct<DIM, NUMBER_OF_VARIABLES>(face, stateCoefficients,
                                                        side, time) {
    pressure_ = this->computePressure(state_);
    speedOfSound_ = this->computeSpeedOfSound(state_, pressure_);
    viscosity_ = this->computeViscosity(state_, partialState_, stateJacobian_,
                                        pressure_);
    hyperbolicMatrix_ = this->computeHyperbolicMatrix(state_, pressure_);
    ellipticTensor_ = this->computeEllipticTensor(partialState_, viscosity_);
}

// todo: NOTE viscosity is based on the stateJacobian, but it is not implemented
// for this case(!) Not sure how to resolve this issue yet
StateCoefficientsStructAir::StateCoefficientsStructAir(
    const LinearAlgebra::MiddleSizeVector &stateBoundary, const double time)
    : StateCoefficientsStruct<DIM, NUMBER_OF_VARIABLES>(stateBoundary, time)

{
    pressure_ = this->computePressure(state_);
    speedOfSound_ = this->computeSpeedOfSound(state_, pressure_);
    viscosity_ = this->computeViscosity(state_, partialState_, stateJacobian_,
                                        pressure_);
    hyperbolicMatrix_ = this->computeHyperbolicMatrix(state_, pressure_);
    ellipticTensor_ = this->computeEllipticTensor(partialState_, viscosity_);
}

StateCoefficientsStructAir::~StateCoefficientsStructAir() {}

/// ***************************************
/// ***      Constutive relations       ***
/// ***************************************

double StateCoefficientsStructAir::computePressure(
    const LinearAlgebra::MiddleSizeVector &state) const {

    // Compute the velocity components
    double velocitySquared = 0;
    for (std::size_t iD = 0; iD < DIM; iD++) {
        velocitySquared += state(iD + 1) * state(iD + 1);
    }

    double pressure =
        (GAMMA - 1) * (state(DIM + 1) - velocitySquared / (2 * state(0)));

    // std::cout << "======Computing pressure=====" << std::endl;
    // std::cout << "state: " << state << std::endl;
    // std::cout << "velocity squared: " << velocitySquared << std::endl;

    if (pressure <= 0) {
        std::cout << "Non-physical behaviour: pressure is zero" << std::endl;
        std::exit(-1);
    }
    // std::cout << "real pressure: " << pressure << std::endl;
    return pressure;
}

double StateCoefficientsStructAir::computeSpeedOfSound(
    const LinearAlgebra::MiddleSizeVector &state, const double pressure) const {
    return std::sqrt(GAMMA * pressure / state(0));
}

LinearAlgebra::MiddleSizeMatrix
    StateCoefficientsStructAir::computeHyperbolicMatrix(
        const LinearAlgebra::MiddleSizeVector &state, const double pressure) {
    LinearAlgebra::MiddleSizeMatrix fluxMatrix(NUMBER_OF_VARIABLES, DIM);
    double densityInverse = 1.0 / state(0);

    for (std::size_t iD = 0; iD < DIM; iD++) {
        // Mass flux
        fluxMatrix(0, iD) = state(iD + 1);  // rho*u, rho*v, rho*w

        // Momentum flux: convection part
        for (std::size_t iD2 = 0; iD2 < DIM; iD2++) {
            fluxMatrix(iD2 + 1, iD) =
                state(iD2 + 1) * state(iD + 1) *
                densityInverse;  // rho*u*u, rho*u*v, rho*u*w for iD2 = 0;
        }

        // Momentum flux: Add pressure part
        fluxMatrix(iD + 1, iD) += pressure;

        // Energy
        fluxMatrix(DIM + 1, iD) =
            (state(DIM + 1) + pressure) * state(iD + 1) * densityInverse;
    }

    return fluxMatrix;
}

double StateCoefficientsStructAir::computeViscosity(
    const LinearAlgebra::MiddleSizeVector &state,
    const LinearAlgebra::MiddleSizeVector &partialState,
    const LinearAlgebra::MiddleSizeMatrix stateJacobian,
    const double pressure) {
    // todo: this temperature definition is bullshit, see FlowBetweenPlates
    double temperature = MACH * MACH * GAMMA * pressure / state(0);

    /*
     * http://ocw.mit.edu/courses/mathematics/18-102-introduction-to-functional-analysis-spring-2009/lecture-notes/
            std::cout << "====Computing Viscosity====" << std::endl;
            std::cout << "MACH: " << MACH << std::endl;
            std::cout << "GAMMA: " << GAMMA << std::endl;
            std::cout << "pressure: " << pressure << std::endl;
            std::cout << "inverse density: " << 1.0/state(0) << std::endl;
            std::cout << "ViscosityScaling: " << VISCOSITY_SCALING << std::endl;
            std::cout << "temperature: " << temperature << std::endl;
            std::cout << "viscosity: " <<
     VISCOSITY_SCALING*(1+THETA_S)/(temperature +
     THETA_S)*std::pow(temperature,3.0/20) << std::endl;
    */
    // todo: VISCOSITY_SCALING has been altered for a test case with high
    // viscosity (!!)
    return VISCOSITY_SCALING;  //*(1+THETA_S)/(temperature +
                               //THETA_S)*std::pow(temperature,3.0/20);
}

std::vector<LinearAlgebra::MiddleSizeMatrix>
    StateCoefficientsStructAir::computeEllipticTensor(
        const LinearAlgebra::MiddleSizeVector partialState,
        const double viscosity) {
    // todo: Check if this also works correctly in 3D

    // std::cout << "viscosity: " << viscosity << std::endl;

    std::vector<LinearAlgebra::MiddleSizeMatrix> ellipticTensor_(DIM * DIM);
    double velocityNormSquared = 0.0;
    double thermalFactor = GAMMA / (REYNOLDS * PRANDTL);

    double pos1;
    double pos2;

    double factor43 = 4.0 / 3.0;
    double vis23 = 2.0 / 3.0 * viscosity;
    double inverseRho = 1.0 / partialState(0);

    for (std::size_t iD = 0; iD < DIM; iD++) {
        velocityNormSquared += partialState(iD + 1) * partialState(iD + 1);
    }

    LinearAlgebra::MiddleSizeMatrix APartial1(DIM + 2, DIM + 2);
    LinearAlgebra::MiddleSizeMatrix APartial2(DIM + 2, DIM + 2);

    // A11 A22 en A33: For documentation see the full matrix in Klaij et al.
    // 2006
    for (std::size_t iD = 0; iD < DIM; iD++) {
        // Reset matrix
        APartial1 *= 0.0;

        // Tensor index
        pos1 = (DIM)*iD + iD;

        // (1) viscosity contributions
        for (std::size_t iD2 = 0; iD2 < DIM; iD2++) {
            APartial1(iD2 + 1, iD2 + 1) = viscosity;
            APartial1(iD2 + 1, 0) = -viscosity * partialState(iD2 + 1);
            APartial1(DIM + 1, iD2 + 1) = viscosity;
        }

        // Multiply by correct value of the dominant component
        APartial1(iD + 1, 0) *= factor43;
        APartial1(iD + 1, iD + 1) *= factor43;
        APartial1(DIM + 1, iD + 1) *= factor43;

        // (2) temperature contribution
        for (std::size_t iD2 = 0; iD2 < DIM; iD2++) {
            APartial1(DIM + 1, iD2 + 1) += -thermalFactor;
        }
        APartial1(DIM + 1, DIM + 1) = thermalFactor;

        // Complete energy part by multiplying by velocity.
        for (std::size_t iD2 = 0; iD2 < DIM; iD2++) {
            APartial1(DIM + 1, iD2 + 1) *= partialState(iD2 + 1);
        }

        APartial1(DIM + 1, 0) =
            -(1.0 / 3.0) * viscosity * partialState(iD + 1) *
                partialState(iD + 1) -
            viscosity * velocityNormSquared -
            thermalFactor * (partialState(DIM + 1) - velocityNormSquared);

        // Divide by rho
        APartial1 *= inverseRho;
        ellipticTensor_[pos1] = APartial1;
    }

    // A12 A13 A23 ? Note: A12 --> A(iD1)(iD2)
    for (std::size_t iD1 = 0; iD1 < DIM - 1; iD1++) {
        for (std::size_t iD2 = iD1 + 1; iD2 < DIM; iD2++) {
            // Reset matrix
            APartial1 *= 0.0;
            APartial2 *= 0.0;

            // Tensor index
            pos1 = (DIM)*iD1 + iD2;
            pos2 = (DIM)*iD2 + iD1;

            // viscosity contributions for A(iD1)(iD2)
            APartial1(iD1 + 1, 0) = vis23 * partialState(iD2 + 1);
            APartial1(iD2 + 1, 0) = -viscosity * partialState(iD1 + 1);

            APartial1(iD1 + 1, iD2 + 1) = -vis23;
            APartial1(iD2 + 1, iD1 + 1) = viscosity;

            APartial1(DIM + 1, 0) = -1.0 / 3.0 * viscosity *
                                    partialState(iD1 + 1) *
                                    partialState(iD2 + 1);
            APartial1(DIM + 1, iD2 + 1) = -vis23 * partialState(iD1 + 1);
            APartial1(DIM + 1, iD1 + 1) = viscosity * partialState(iD2 + 1);

            // viscosity contributions for A(iD2)(iD1)
            APartial2(iD1 + 1, 0) = -viscosity * partialState(iD2 + 1);
            APartial2(iD2 + 1, 0) = vis23 * partialState(iD1 + 1);

            APartial2(iD1 + 1, iD2 + 1) = viscosity;
            APartial2(iD2 + 1, iD1 + 1) = -vis23;

            APartial2(DIM + 1, 0) = -1.0 / 3.0 * viscosity *
                                    partialState(iD1 + 1) *
                                    partialState(iD2 + 1);
            APartial2(DIM + 1, iD2 + 1) = viscosity * partialState(iD1 + 1);
            APartial2(DIM + 1, iD1 + 1) = -vis23 * partialState(iD2 + 1);

            // Divide by rho
            APartial1 *= inverseRho;
            APartial2 *= inverseRho;

            // Assign matrices to the tensor vector
            ellipticTensor_[pos1] = APartial1;
            ellipticTensor_[pos2] = APartial2;
        }
    }
    return ellipticTensor_;
}

/// ***************************************
/// ***      Additional Functions       ***
/// ***************************************

LinearAlgebra::MiddleSizeMatrix
    StateCoefficientsStructAir::computeEllipticTensorMatrixContractionFast(
        const LinearAlgebra::MiddleSizeMatrix &matrix) const {
    LinearAlgebra::MiddleSizeMatrix result(NUMBER_OF_VARIABLES, DIM);
    double pos;

    // A11 A22 and A33
    for (std::size_t iD = 0; iD < DIM; iD++) {
        std::size_t iDm = iD;
        pos = (DIM)*iD + iDm;
        // Velocity part
        for (std::size_t it1 = 0; it1 < DIM; it1++) {
            result(it1 + 1, iD) +=
                ellipticTensor_[pos](it1 + 1, 0) * matrix(0, iDm);
            result(it1 + 1, iD) +=
                ellipticTensor_[pos](it1 + 1, it1 + 1) * matrix(it1 + 1, iDm);
        }
        // Energy part
        for (std::size_t iVm = 0; iVm < NUMBER_OF_VARIABLES; iVm++) {
            result(DIM + 1, iD) +=
                ellipticTensor_[pos](DIM + 1, iVm) * matrix(iVm, iDm);
        }
    }

    // Other
    for (std::size_t iD = 0; iD < DIM; iD++) {
        for (std::size_t iDm = iD + 1; iDm < DIM; iDm++) {

            // A12 A13 and A23
            pos = (DIM)*iD + iDm;
            // Velocity part
            result(iD + 1, iD) +=
                ellipticTensor_[pos](iD + 1, 0) * matrix(0, iDm);
            result(iD + 1, iD) +=
                ellipticTensor_[pos](iD + 1, iDm + 1) * matrix(iDm + 1, iDm);
            result(iDm + 1, iD) +=
                ellipticTensor_[pos](iDm + 1, 0) * matrix(0, iDm);
            result(iDm + 1, iD) +=
                ellipticTensor_[pos](iDm + 1, iD + 1) * matrix(iD + 1, iDm);
            // Energy part
            for (std::size_t iVm = 0; iVm < NUMBER_OF_VARIABLES; iVm++) {
                result(DIM + 1, iD) +=
                    ellipticTensor_[pos](DIM + 1, iVm) * matrix(iVm, iDm);
            }

            // A21 A31 and A32
            // Note: Same structure as above, however iD and iDm are switched
            pos = (DIM)*iDm + iD;
            // Velocity part
            result(iDm + 1, iDm) +=
                ellipticTensor_[pos](iDm + 1, 0) * matrix(0, iD);
            result(iDm + 1, iDm) +=
                ellipticTensor_[pos](iDm + 1, iD + 1) * matrix(iD + 1, iD);
            result(iD + 1, iDm) +=
                ellipticTensor_[pos](iD + 1, 0) * matrix(0, iD);
            result(iD + 1, iDm) +=
                ellipticTensor_[pos](iD + 1, iDm + 1) * matrix(iDm + 1, iD);
            // Energy part
            for (std::size_t iVm = 0; iVm < NUMBER_OF_VARIABLES; iVm++) {
                result(DIM + 1, iDm) +=
                    ellipticTensor_[pos](DIM + 1, iVm) * matrix(iVm, iD);
            }
        }
    }

    return result;
}
