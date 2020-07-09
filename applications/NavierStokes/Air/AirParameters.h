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

#ifndef HPGEM_APP_AIRPARAMETERS_H
#define HPGEM_APP_AIRPARAMETERS_H

#include <cstddef>

/// Simulation parameters
const static std::size_t DIM = 2;
const static std::size_t NUMBER_OF_VARIABLES = DIM + 2;

/// reference input parameters
const static double DENSITY_REF = 1.0;      // [kg/m^3] at T = 308.16 [K]
const static double VELOCITY_REF = 1.0;     // [m/s]
const static double LENGTH_REF = 1.0;       // [m]
const static double TEMPERATURE_REF = 1.0;  // [K]

/// Dimensionless input parameters
const static double MACH = 0.5;  // Rati o of wall velocity with speed of sound
const static double REYNOLDS = 500;  // Ratio of inertia and viscous forces

/// Dimensionless groups
const static double GAMMA = 1.4;    // Ratio of specific heats cp/cv
const static double THETA_S = 0.4;  // Ratio in Sutherlands law
const static double PRANDTL =
    0.72;  // Ratio of viscous diffusion and thermal diffusion

/// Derived from the input:
const static double SOUND_VELOCITY = VELOCITY_REF / MACH;
const static double PRESSURE_WALL =
    DENSITY_REF * SOUND_VELOCITY * SOUND_VELOCITY / GAMMA;
const static double GAS_CONSTANT =
    PRESSURE_WALL / (DENSITY_REF * TEMPERATURE_REF);
const static double SPECIFIC_HEAT_PRESSURE = GAMMA * GAS_CONSTANT / (GAMMA - 1);
const static double MU_WALL =
    DENSITY_REF * VELOCITY_REF * LENGTH_REF / REYNOLDS;
const static double KAPPA_WALL = SPECIFIC_HEAT_PRESSURE * MU_WALL / PRANDTL;
const static double SUTHERLAND_TEMPERATURE = THETA_S * TEMPERATURE_REF;

/// Often used parameters:
const static double VISCOSITY_SCALING = 0.01;  // 1.0/REYNOLDS;

#endif  // HPGEM_APP_AIRPARAMETERS_H
