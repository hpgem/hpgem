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

#include "Air.h"
#include "Logger.h"
#include <iostream>
#include <chrono>
#include <fstream>

using namespace hpgem;

auto& name = Base::register_argument<std::string>(
    'n', "meshName", "name of the mesh file", true);
auto& polynomialOrder = Base::register_argument<std::size_t>(
    'p', "order", "polynomial order of the solution", true);

auto& numOfOutputFrames = Base::register_argument<std::size_t>(
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

    // Read input parameters
    Base::parse_options(argc, argv);

    // Select explicit integration method
    const TimeIntegration::ButcherTableau* const ptrButcherTableau =
        TimeIntegration::AllTimeIntegrators::Instance().getRule(3, 3, true);

    // Set variable names and number of parameters
    std::vector<std::string> variableNames;
    for (std::size_t i = 0; i < NUMBER_OF_VARIABLES;
         i++)  // +1 is for the energy equation
    {
        std::string variableName = "q" + std::to_string(i + 1);
        variableNames.push_back(variableName);
    }

    // Calculate number of variables
    const std::size_t numOfVariables = NUMBER_OF_VARIABLES;

    // Create problem solver 'test', that can solve the compressible
    // Navier-Stokes equations.
    Air problem(numOfVariables, endTime.getValue(), polynomialOrder.getValue(),
                ptrButcherTableau, true);

    // Create the mesh, a simple square domain
    problem.readMesh(name.getValue());

    // Solve the problem over time interval [startTime,endTime].
    problem.solve(startTime.getValue(), endTime.getValue(), dt.getValue(),
                  numOfOutputFrames.getValue(), true);

    // Measure elapsed time
    endClock = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = endClock - startClock;
    logger(INFO, "Elapsed time for solving the PDE: % s",
           elapsed_seconds.count());

    // Output dimensionfull variables:
    std::cout << "TEMPERATURE_WALL: " << TEMPERATURE_REF << std::endl;
    std::cout << "SOUND_VELOCITY: " << SOUND_VELOCITY << std::endl;
    std::cout << "Pressure_wall: " << PRESSURE_WALL << std::endl;
    std::cout << "Viscosity_wall: " << MU_WALL << std::endl;
    std::cout << "gas constant: " << GAS_CONSTANT << std::endl;

    // Output error to file for testing purposes
    LinearAlgebra::MiddleSizeVector errors =
        problem.computeMaxError(0, endTime.getValue());
    std::fstream myFile;
    myFile.open("errorResult", std::fstream::out | std::fstream::app);
    myFile << "=========START=========" << std::endl;
    myFile << "Density: " << errors(0) << std::endl;
    myFile << "Momentum x: " << errors(1) << std::endl;
    myFile << "Momentum y: " << errors(2) << std::endl;
    myFile << "Total Energy: " << errors(3) << std::endl;
    myFile << "p = " << polynomialOrder.getValue() << std::endl;
    myFile << "Mesh name = " << name.getValue() << std::endl;
    myFile << "=========END===========" << std::endl;
    myFile.close();

    return 0;
}
