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
#include "Air.h"
#include "../UnsteadyNavierStokesAPI.h"
#include "StateCoefficientsStructAir.h"

Air::Air(const std::size_t numOfVariables, const double endTime,
         const std::size_t polynomialOrder,
         const TimeIntegration::ButcherTableau *const ptrButcherTableau,
         const bool computeBothFaces)
    : UnsteadyNavierStokesAPI<DIM, NUMBER_OF_VARIABLES>(
          numOfVariables, endTime, polynomialOrder, ptrButcherTableau,
          computeBothFaces) {}

/// \brief This function couples the constitutive functions to the NavierStokes
/// API
StateCoefficientsStruct<DIM, NUMBER_OF_VARIABLES>
    Air::constutiveRelationsStruct() {
    StateCoefficientsStructAir result;
    return result;
}

/// \brief This function couples the constitutive relations to the NavierStokes
/// API
StateCoefficientsStruct<DIM, NUMBER_OF_VARIABLES>
    Air::computeElementStateStruct(
        Base::PhysicalElement<DIM> &element,
        const LinearAlgebra::MiddleSizeVector &stateCoefficients,
        const double time) {
    StateCoefficientsStructAir result(element, stateCoefficients, time);
    return result;
}

StateCoefficientsStruct<DIM, NUMBER_OF_VARIABLES> Air::computeFaceStateStruct(
    Base::PhysicalFace<DIM> &face,
    const LinearAlgebra::MiddleSizeVector &stateCoefficients,
    const Base::Side side, const double time) {
    StateCoefficientsStructAir result(face, stateCoefficients, side, time);

    return result;
}

StateCoefficientsStruct<DIM, NUMBER_OF_VARIABLES>
    Air::computeBoundaryFaceStateStruct(
        const LinearAlgebra::MiddleSizeVector &stateBoundary,
        const double time) {
    StateCoefficientsStructAir result(stateBoundary, time);

    return result;
}

/// \brief computes the source at an element
LinearAlgebra::MiddleSizeVector Air::integrandSourceAtElement(
    Base::PhysicalElement<DIM> &element,
    const StateCoefficientsStruct<DIM, NUMBER_OF_VARIABLES> &elementStateStruct,
    const double &time) {
    std::size_t numOfBasisFunctions = element.getNumberOfBasisFunctions();

    LinearAlgebra::MiddleSizeVector integrandSource(NUMBER_OF_VARIABLES *
                                                    numOfBasisFunctions);

    // Convert pRef to pPhys
    Geometry::PointPhysical<DIM> pPhys = element.getPointPhysical();

    // create datastructures
    double A = DENSITY_REF / DENSITY_REF;
    double C = 0.1;
    double D = 0.5;
    double B = 10.0;
    double C_t = 0;
    double pi = M_PI;
    LinearAlgebra::MiddleSizeVector zVal(4);
    LinearAlgebra::MiddleSizeVector exactState(4);

    double zT = time;
    double zT_t = 1.0;

    double k1 = 2 * pi;
    double k2 = 2 * pi;

    zVal(0) = std::cos(k1 * pPhys[0]) * std::cos(k2 * pPhys[1]);
    zVal(1) = std::sin(k1 * pPhys[0]) * std::cos(k2 * pPhys[1]);
    zVal(2) = std::cos(k1 * pPhys[0]) * std::sin(k2 * pPhys[1]);
    zVal(3) = std::sin(k1 * pPhys[0]) * std::sin(k2 * pPhys[1]);

    // compute state
    exactState(0) = A * (1 + zT * C * zVal(0));
    exactState(1) = D * (1 + zT * C * zVal(0));
    exactState(2) = D * (1 + zT * C * zVal(0));
    exactState(3) = B * (1 + zT * C * zVal(0));

    LinearAlgebra::MiddleSizeVector state_t(4);
    LinearAlgebra::MiddleSizeVector state_x(4);
    LinearAlgebra::MiddleSizeVector state_y(4);
    LinearAlgebra::MiddleSizeVector state_xx(4);
    LinearAlgebra::MiddleSizeVector state_yy(4);
    LinearAlgebra::MiddleSizeVector state_xy(4);

    LinearAlgebra::MiddleSizeVector partialState(4);
    LinearAlgebra::MiddleSizeVector partialState_x(4);
    LinearAlgebra::MiddleSizeVector partialState_y(4);
    LinearAlgebra::MiddleSizeVector partialState_xx(4);
    LinearAlgebra::MiddleSizeVector partialState_yy(4);
    LinearAlgebra::MiddleSizeVector partialState_xy(4);

    double inverseDensity = 1.0 / exactState(0);

    /*	//Compute state_t
            state_t(0) 	= zT_t*A*C*zVal(0);
            state_t(1) 	= zT_t*D*C*zVal(1);
            state_t(2) 	= zT_t*D*C*zVal(2);
            state_t(3)  = zT_t*B*C*zVal(3);

            //compute state_x
            state_x(0) = -k1*zT*C*A*zVal(1);
            state_x(1) = k1*zT*C*D*zVal(0);
            state_x(2) = -k1*zT*C*D*zVal(3);
            state_x(3) = k1*zT*C*B*zVal(2);

            //compute state_y
            state_y(0) = -k2*zT*C*A*zVal(2);
            state_y(1) = -k2*zT*C*D*zVal(3);
            state_y(2) = k2*zT*C*D*zVal(0);
            state_y(3) = k2*zT*C*B*zVal(1);

            //compute state_xx
            state_xx(0) = -k1*k1*zT*C*A*zVal(0);
            state_xx(1) = -k1*k1*zT*C*D*zVal(1);
            state_xx(2) = -k1*k1*zT*C*D*zVal(2);
            state_xx(3) = -k1*k1*zT*C*B*zVal(3);

            //compute state_yy
            state_yy(0) = -k2*k2*zT*C*A*zVal(0);
            state_yy(1) = -k2*k2*zT*C*D*zVal(1);
            state_yy(2) = -k2*k2*zT*C*D*zVal(2);
            state_yy(3) = -k2*k2*zT*C*B*zVal(3);

            //compute state_xy
            state_xy(0) = k1*k2*C*A*zVal(3);
            state_xy(1) = -k1*k2*C*D*zVal(2);
            state_xy(2) = -k1*k2*C*D*zVal(1);
            state_xy(3) = k1*k2*C*B*zVal(0);*/

    // Compute state_t
    state_t(0) = zT_t * A * C * zVal(0);
    state_t(1) = zT_t * D * C * zVal(0);
    state_t(2) = zT_t * D * C * zVal(0);
    state_t(3) = zT_t * B * C * zVal(0);

    // compute state_x
    state_x(0) = -k1 * zT * C * A * zVal(1);
    state_x(1) = -k1 * zT * C * D * zVal(1);
    state_x(2) = -k1 * zT * C * D * zVal(1);
    state_x(3) = -k1 * zT * C * B * zVal(1);

    // compute state_y
    state_y(0) = -k2 * zT * C * A * zVal(2);
    state_y(1) = -k2 * zT * C * D * zVal(2);
    state_y(2) = -k2 * zT * C * D * zVal(2);
    state_y(3) = -k2 * zT * C * B * zVal(2);

    // compute state_xx
    state_xx(0) = -k1 * k1 * zT * C * A * zVal(0);
    state_xx(1) = -k1 * k1 * zT * C * D * zVal(0);
    state_xx(2) = -k1 * k1 * zT * C * D * zVal(0);
    state_xx(3) = -k1 * k1 * zT * C * B * zVal(0);

    // compute state_yy
    state_yy(0) = -k2 * k2 * zT * C * A * zVal(0);
    state_yy(1) = -k2 * k2 * zT * C * D * zVal(0);
    state_yy(2) = -k2 * k2 * zT * C * D * zVal(0);
    state_yy(3) = -k2 * k2 * zT * C * B * zVal(0);

    // compute state_xy
    state_xy(0) = k1 * k2 * C * A * zVal(3);
    state_xy(1) = k1 * k2 * C * D * zVal(3);
    state_xy(2) = k1 * k2 * C * D * zVal(3);
    state_xy(3) = k1 * k2 * C * B * zVal(3);

    // compute partialState
    partialState(0) = exactState(0);
    partialState(1) = exactState(1) * inverseDensity;
    partialState(2) = exactState(2) * inverseDensity;
    partialState(3) = exactState(3) * inverseDensity;

    // compute partialState_x
    partialState_x(0) = state_x(0);
    partialState_x(1) =
        inverseDensity * (state_x(1) - partialState(1) * state_x(0));
    partialState_x(2) =
        inverseDensity * (state_x(2) - partialState(2) * state_x(0));
    partialState_x(3) =
        inverseDensity * (state_x(3) - partialState(3) * state_x(0));

    // compute partialState_y
    partialState_y(0) = state_y(0);
    partialState_y(1) =
        inverseDensity * (state_y(1) - partialState(1) * state_y(0));
    partialState_y(2) =
        inverseDensity * (state_y(2) - partialState(2) * state_y(0));
    partialState_y(3) =
        inverseDensity * (state_y(3) - partialState(3) * state_y(0));

    // compute partialState_xx
    partialState_xx(0) = state_xx(0);
    partialState_yy(0) = state_yy(0);
    partialState_xy(0) = state_xy(0);

    for (std::size_t iV = 1; iV < NUMBER_OF_VARIABLES; iV++) {
        partialState_xx(iV) =
            inverseDensity *
            (state_xx(iV) - 2.0 * (partialState_x(iV) * state_x(0)) -
             partialState(iV) * state_xx(0));
        partialState_yy(iV) =
            inverseDensity *
            (state_yy(iV) - 2.0 * (partialState_y(iV) * state_y(0)) -
             partialState(iV) * state_yy(0));
        partialState_xy(iV) =
            inverseDensity * (state_xy(iV) - (partialState_x(iV) * state_y(0)) -
                              (partialState_y(iV) * state_x(0)) -
                              partialState(iV) * state_xy(0));
    }

    //*************************************
    //***	Compute closure variables	***
    //*************************************

    double pressure =
        (GAMMA - 1) * (exactState(3) - 0.5 *
                                           (partialState(1) * partialState(1) +
                                            partialState(2) * partialState(2)) *
                                           (exactState(0)));

    double e = partialState(3) - 0.5 * (partialState(1) * partialState(1) +
                                        partialState(2) * partialState(2));
    double e_x = partialState_x(3) - partialState(1) * partialState_x(1) -
                 partialState(2) * partialState_x(2);
    double e_y = partialState_y(3) - partialState(1) * partialState_y(1) -
                 partialState(2) * partialState_y(2);
    double e_xx = partialState_xx(3) - partialState_x(1) * partialState_x(1) -
                  partialState(1) * partialState_xx(1) -
                  partialState_x(2) * partialState_x(2) -
                  partialState(2) * partialState_xx(2);
    double e_yy = partialState_yy(3) - partialState_y(1) * partialState_y(1) -
                  partialState(1) * partialState_yy(1) -
                  partialState_y(2) * partialState_y(2) -
                  partialState(2) * partialState_yy(2);

    double temperatureScaling = GAMMA * (GAMMA - 1) * MACH * MACH;
    double T = e * temperatureScaling;
    double T_x = e_x * temperatureScaling;
    double T_y = e_y * temperatureScaling;
    double T_xx = e_xx * temperatureScaling;
    double T_yy = e_yy * temperatureScaling;

    double lambda = -2.0 / 3.0;

    double mu = VISCOSITY_SCALING;  // VISCOSITY_SCALING*(1+THETA_S)/(T +
                                    // THETA_S)*std::pow(T,3.0/20);
    double mu_x =
        0.0;  // VISCOSITY_SCALING*(1+THETA_S)/(T + THETA_S)*( -(1/(T +
              // THETA_S))*std::pow(T,3.0/2.0) + 3.0/2.0*std::sqrt(T))*T_x;
    double mu_y =
        0.0;  // VISCOSITY_SCALING*(1+THETA_S)/(T + THETA_S)*( -(1/(T +
              // THETA_S))*std::pow(T,3.0/2.0) + 3.0/2.0*std::sqrt(T))*T_y;

    // todo: this is a more elegant way, replace it after the code is verified
    double mu_xnew = (-mu / (T + THETA_S) + 3.0 * mu / (2.0 * T)) * T_x;

    /*	std::cout << "mu_x old" << mu_x << std::endl;
            std::cout << "mu_x new" << mu_xnew << std::endl;
            std::exit(-1);*/
    // std::cout << "E: " << partialState_x(3) << std::endl;
    // std::cout << "u^2: " << partialState(1)*partialState(1) << std::endl;
    // std::cout << "v^2: " << partialState(2)*partialState(2) << std::endl;
    // std::cout << "e_x: " << e_x << std::endl;
    // std::cout << "e_y: " << e_y << std::endl;

    // std::cout << "mu_x: " << mu_x << std::endl;
    // std::cout << "mu_y: " << mu_y << std::endl;
    // std::cout << "partialState_x" << partialState_x << std::endl;
    // std::cout << "partialState_y" << partialState_y << std::endl;

    /*	std::cout << "partialState_xx" << partialState_xx << std::endl;
            std::cout << "partialState_yy" << partialState_yy << std::endl;
            std::cout << "partialState_xy" << partialState_xy << std::endl;*/

    // std::exit(-1);

    double Txx =
        (2 + lambda) * mu * partialState_x(1) + lambda * mu * partialState_y(2);
    double Txy = mu * (partialState_y(1) + partialState_x(2));
    double Tyy =
        (2 + lambda) * mu * partialState_y(2) + lambda * mu * partialState_x(1);

    //*****************************
    //*** 	Convection terms	***
    //*****************************

    double conv_uu_x =
        exactState(1) * partialState_x(1) + partialState(1) * state_x(1);
    double conv_uu_y =
        exactState(1) * partialState_y(1) + partialState(1) * state_y(1);
    double conv_vu_y =
        exactState(1) * partialState_y(2) + partialState(2) * state_y(1);
    double conv_uv_x =
        exactState(1) * partialState_x(2) + partialState(2) * state_x(1);
    double conv_vv_x =
        exactState(2) * partialState_x(2) + partialState(2) * state_x(2);
    double conv_vv_y =
        exactState(2) * partialState_y(2) + partialState(2) * state_y(2);

    double p_x = (GAMMA - 1) * (state_x(3) - 0.5 * conv_uu_x - 0.5 * conv_vv_x);
    double p_y = (GAMMA - 1) * (state_y(3) - 0.5 * conv_uu_y - 0.5 * conv_vv_y);

    //**************************
    //***	Viscous terms	 ***
    //**************************

    double Txx_x = (2 + lambda) * mu_x * partialState_x(1) +
                   (2 + lambda) * mu * partialState_xx(1) +
                   lambda * mu_x * partialState_y(2) +
                   lambda * mu * partialState_xy(2);
    double Txy_x = mu_x * (partialState_y(1) + partialState_x(2)) +
                   mu * (partialState_xy(1) + partialState_xx(2));
    double Tyx_y = mu_y * (partialState_y(1) + partialState_x(2)) +
                   mu * (partialState_yy(1) + partialState_xy(2));
    double Tyy_y = (2 + lambda) * mu_y * partialState_y(2) +
                   (2 + lambda) * mu * partialState_yy(2) +
                   lambda * mu_y * partialState_x(1) +
                   lambda * mu * partialState_xy(1);

    double H = partialState(3) + pressure / exactState(0);
    double H_x = partialState_x(3) + p_x / exactState(0) -
                 pressure / (exactState(0) * exactState(0)) * state_x(0);
    double H_y = partialState_y(3) + p_y / exactState(0) -
                 pressure / (exactState(0) * exactState(0)) * state_y(0);
    double conv_Hu_x = exactState(1) * H_x + H * state_x(1);
    double conv_Hv_y = exactState(2) * H_y + H * state_y(2);

    double Ax_x = partialState(1) * Txx_x + Txx * partialState_x(1) +
                  partialState(2) * Txy_x + Txy * partialState_x(2) +
                  1.0 / (PRANDTL * REYNOLDS * MACH * MACH * (GAMMA - 1)) * T_xx;

    double Ay_y = partialState(1) * Tyx_y + Txy * partialState_y(1) +
                  partialState(2) * Tyy_y + Tyy * partialState_y(2) +
                  1.0 / (PRANDTL * REYNOLDS * MACH * MACH * (GAMMA - 1)) * T_yy;

    /*	std::cout << "============" << std::endl;
            std::cout << "partialState" << partialState << std::endl;
            std::cout << "partialState_x" << partialState_x << std::endl;
        std::cout << "partialState_y" << partialState_y << std::endl;
            std::cout << "state_x " << state_x << std::endl;
            std::cout << "state_y " << state_y << std::endl;
            std::cout << "state_xx " << state_xx << std::endl;
            std::cout << "state_yy " << state_yy << std::endl;
            std::cout << "state_xy " << state_xy << std::endl;
            std::cout << "--------" << std::endl;
            std::cout << "Txx = " << Txx << std::endl;
            std::cout << "Tyy = " << Tyy << std::endl;
            std::cout << "Txt = " << Txy << std::endl;
            std::cout << "Txx_x = " << Txx_x << std::endl;
            std::cout << "Txy_x = " << Txy_x << std::endl;
            std::cout << "Tyx_y = " << Tyx_y << std::endl;
            std::cout << "Tyy_y = " << Tyy_y << std::endl;
            std::cout << "Ax_x" << Ax_x << std::endl;
            std::cout << "Ay_y" << Ay_y << std::endl;*/

    //*********************************
    //***	Compute Source term		***
    //*********************************

    double sDensity = state_t(0) + state_x(1) + state_y(2);
    double sMomentumX =
        state_t(1) + conv_uu_x + conv_vu_y + p_x;  // - (Txx_x + Tyx_y);
    double sMomentumY =
        state_t(2) + conv_uv_x + conv_vv_y + p_y;         // - (Txy_x + Tyy_y);
    double sEnergy = state_t(3) + conv_Hu_x + conv_Hv_y;  // - (Ax_x + Ay_y );

    // std::cout << "---" << std::endl;
    // std::cout << "pressure: " << pressure << std::endl;

    /*	std::cout << "Ax_x: " << Ax_x << std::endl;
            std::cout << "conv_Hu_x: " << conv_Hu_x << std::endl;
            std::cout << "Txx_x: " << Txx_x << std::endl;*/

    //*************************************
    //***	Compute source integrand	***
    //*************************************

    //	std::cout << "time " << time << std::endl;

    std::size_t iVB;
    for (std::size_t iB = 0; iB < numOfBasisFunctions; iB++) {
        // Density
        iVB = element.convertToSingleIndex(iB, 0);
        integrandSource(iVB) = sDensity * element.basisFunction(iB);

        // Momentum x
        iVB = element.convertToSingleIndex(iB, 1);
        integrandSource(iVB) = sMomentumX * element.basisFunction(iB);

        // Momentum Y
        iVB = element.convertToSingleIndex(iB, 2);
        integrandSource(iVB) = sMomentumY * element.basisFunction(iB);

        // Energy
        iVB = element.convertToSingleIndex(iB, 3);
        integrandSource(iVB) = sEnergy * element.basisFunction(iB);
    }

    return integrandSource;
}

LinearAlgebra::MiddleSizeVector Air::getInitialSolution(
    const Geometry::PointPhysical<DIM> &pPhys, const double &startTime,
    const std::size_t orderTimeDerivative) {
    // create datastructures
    double A = DENSITY_REF / DENSITY_REF;
    double C = 0.1;
    double D = 0.5;
    double B = 10.0;
    double pi = M_PI;
    LinearAlgebra::MiddleSizeVector zVal(4);
    LinearAlgebra::MiddleSizeVector initialState(4);

    /*    double zT = 1 + 0.5*std::cos(2*pi*startTime);
        double zT_t = -pi*std::sin(2*pi*startTime);*/

    double zT = startTime;
    double zT_t = 1.0;

    double k1 = 2 * pi;
    double k2 = 2 * pi;

    zVal(0) = std::cos(k1 * pPhys[0]) * std::cos(k2 * pPhys[1]);
    zVal(1) = std::sin(k1 * pPhys[0]) * std::cos(k2 * pPhys[1]);
    zVal(2) = std::cos(k1 * pPhys[0]) * std::sin(k2 * pPhys[1]);
    zVal(3) = std::sin(k1 * pPhys[0]) * std::sin(k2 * pPhys[1]);

    // compute state
    initialState(0) = A * (1 + zT * C * zVal(0));
    initialState(1) = D * (1 + zT * C * zVal(0));
    initialState(2) = D * (1 + zT * C * zVal(0));
    initialState(3) = B * (1 + zT * C * zVal(0));

    return initialState;
}

LinearAlgebra::MiddleSizeVector Air::getExactSolution(
    const Geometry::PointPhysical<DIM> &pPhys, const double &time,
    const std::size_t orderTimeDerivative) {

    double A = DENSITY_REF / DENSITY_REF;
    double C = 0.1;
    double D = 0.5;
    double B = 10.0;
    double pi = M_PI;
    LinearAlgebra::MiddleSizeVector zVal(4);
    LinearAlgebra::MiddleSizeVector exactState(4);

    double zT = time;
    double zT_t = 1.0;

    double k1 = 2 * pi;
    double k2 = 2 * pi;

    zVal(0) = std::cos(k1 * pPhys[0]) * std::cos(k2 * pPhys[1]);
    zVal(1) = std::sin(k1 * pPhys[0]) * std::cos(k2 * pPhys[1]);
    zVal(2) = std::cos(k1 * pPhys[0]) * std::sin(k2 * pPhys[1]);
    zVal(3) = std::sin(k1 * pPhys[0]) * std::sin(k2 * pPhys[1]);

    // compute state
    exactState(0) = A * (1 + zT * C * zVal(0));
    exactState(1) = D * (1 + zT * C * zVal(0));
    exactState(2) = D * (1 + zT * C * zVal(0));
    exactState(3) = B * (1 + zT * C * zVal(0));

    return exactState;
}

/// \brief Compute the dirichlet boundary condition on the face
// todo: add time
LinearAlgebra::MiddleSizeVector Air::computeBoundaryState(
    Base::PhysicalFace<DIM> &face,
    const StateCoefficientsStruct<DIM, NUMBER_OF_VARIABLES> faceStateStuctLeft,
    const double time) {
    const Geometry::PointPhysical<DIM> pPhys = face.getPointPhysical();

    return getExactSolution(pPhys, time, 0);

    /*	//Check which boundary it is
            if(pPhys[1] == 1)
            {
                    //std::cout << "top BC!" << std::endl;

            }
            else if(pPhys[1] == 0)
            {
                    //std::cout << "bottom BC!" << std::endl;
            }
            else if(pPhys[0] == 0)
            {
                    //std::cout << "left BC!" << std::endl;
            }
            else if(pPhys[0] == 1)
            {
                    //std::cout << "right BC!" << std::endl;
            }
            else
            {
                    logger(ERROR,"Incorrect Boundary.");
                    LinearAlgebra::MiddleSizeVector stateBoundary;
                    return stateBoundary;
            }*/
}

/// \brief Output function. This fucntion gives the state as output, and
/// additionally also the error if specified.
void Air::registerVTKWriteFunctions() {
    HpgemAPISimplified::registerVTKWriteFunctions();

    /*    registerVTKWriteFunction([ = ](Base::Element* element, const
       Geometry::PointReference<1>& pRef, std::size_t timeLevel) -> double
        {
            if (element->getSolution(timeLevel, pRef)[0] > 1e-5)
                return std::real(element->getSolution(timeLevel, pRef)[1] /
       element->getSolution(timeLevel, pRef)[0]); return 0;
        }, "u");*/

    registerVTKWriteFunction(
        [=](Base::Element *element, const Geometry::PointReference<DIM> &pRef,
            std::size_t timeLevel) -> double {
            LinearAlgebra::MiddleSizeVector exactSolution =
                getExactSolution(element->referenceToPhysical(pRef), time_, 0);

            return std::real(element->getSolution(timeLevel, pRef)[0] -
                             exactSolution(0));
        },
        "density error");

    registerVTKWriteFunction(
        [=](Base::Element *element, const Geometry::PointReference<DIM> &pRef,
            std::size_t timeLevel) -> double {
            LinearAlgebra::MiddleSizeVector exactSolution =
                getExactSolution(element->referenceToPhysical(pRef), time_, 0);

            return std::real(element->getSolution(timeLevel, pRef)[1] -
                             exactSolution(1));
        },
        "momentum x error");

    registerVTKWriteFunction(
        [=](Base::Element *element, const Geometry::PointReference<DIM> &pRef,
            std::size_t timeLevel) -> double {
            LinearAlgebra::MiddleSizeVector exactSolution =
                getExactSolution(element->referenceToPhysical(pRef), time_, 0);

            return std::real(element->getSolution(timeLevel, pRef)[2] -
                             exactSolution(2));
        },
        "momentum y error");

    registerVTKWriteFunction(
        [=](Base::Element *element, const Geometry::PointReference<DIM> &pRef,
            std::size_t timeLevel) -> double {
            LinearAlgebra::MiddleSizeVector exactSolution =
                getExactSolution(element->referenceToPhysical(pRef), time_, 0);

            return std::real(element->getSolution(timeLevel, pRef)[3] -
                             exactSolution(3));
        },
        "energy error");
}
