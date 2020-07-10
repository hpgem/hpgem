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

#include "CompressibleNavierStokes.h"
#include "Logger.h"
#include <iostream>
#include <chrono>

using namespace hpgem;

// todo: decide if choosing the dimension during run time is worth the design
// overhead
auto& name = Base::register_argument<std::string>(
    'n', "meshName", "name of the mesh file", true);
auto& polynomialOrder = Base::register_argument<std::size_t>(
    'p', "order", "polynomial order of the solution", true);

auto& numberOfOutputFrames = Base::register_argument<std::size_t>(
    'O', "numOfOutputFrames", "Number of frames to output", false, 100);
auto& startTime = Base::register_argument<double>(
    'S', "startTime", "start time of the simulation", false, 0.0);
auto& endTime = Base::register_argument<double>(
    'T', "endTime", "end time of the simulation", false, 0.001);
auto& dt = Base::register_argument<double>(
    'd', "timeStepSize", "time step of the simulation", false, 0.00001);

int main(int argc, char** argv) {

    // Start measuring elapsed time
    std::chrono::time_point<std::chrono::system_clock> startClock, endClock;
    startClock = std::chrono::system_clock::now();

    // todo: Reduce the number of numericalvector creations
    // todo: optimize the code, placement of variables
    // todo: make timestep CFL dependent
    Base::parse_options(argc, argv);

    // logger(WARN,"WARNING: Timestep is determined a priori. Stability Criteria
    // might not be satisfied!");
    const TimeIntegration::ButcherTableau* const ptrButcherTableau =
        TimeIntegration::AllTimeIntegrators::Instance().getRule(3, 3, true);

    // Set variable names and number of parameters
    std::vector<std::string> variableNames;
    variableNames.push_back("q1");

    for (std::size_t i = 1; i <= DIM + 1; i++)  // +1 is for the energy equation
    {
        std::string variableName = "q" + std::to_string(i + 1);
        variableNames.push_back(variableName);
    }

    // Calculate number of variables
    const std::size_t numOfVariables = 2 + DIM;

    // Create problem solver 'test', that can solve the compressible
    // Navier-Stokes equations.
    CompressibleNavierStokes test(numOfVariables, endTime.getValue(),
                                  polynomialOrder.getValue(), ptrButcherTableau,
                                  true);

    // Create the mesh, a simple square domain
    test.readMesh(name.getValue());

    // Sets the mass matrix required for the stability parameters in the viscous
    // class. Slight hack. todo: improve this
    test.setStabilityMassMatrix();

    // Solve the problem over time interval [startTime,endTime].
    test.solve(startTime.getValue(), endTime.getValue(), dt.getValue(),
               numberOfOutputFrames.getValue(), true);

    // Compute errors at the end of the simulation
    LinearAlgebra::MiddleSizeVector maxError = test.Error(endTime.getValue());
    std::cout << maxError << std::endl;

    // Measure elapsed time
    endClock = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = endClock - startClock;
    logger(INFO, "Elapsed time for solving the PDE: % s",
           elapsed_seconds.count());

    return 0;
}
