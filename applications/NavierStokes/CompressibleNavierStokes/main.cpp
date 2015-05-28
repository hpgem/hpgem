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

#include "CompressibleNavierStokes.h"
#include "Logger.h"
#include <iostream>

auto& dimension = Base::register_argument<std::size_t>('D', "dim", "number of dimensions in the problem");
auto& numOfElements = Base::register_argument<std::size_t>('n', "numElems", "number of elements per dimension", true);
auto& polynomialOrder = Base::register_argument<std::size_t>('p', "order", "polynomial order of the solution", true);

auto& numOfOutputFrames = Base::register_argument<std::size_t>('O', "numOfOutputFrames", "Number of frames to output", false, 200);
auto& startTime = Base::register_argument<double>('S', "startTime", "start time of the simulation", false, 0.0);
auto& endTime = Base::register_argument<double>('T', "endTime", "end time of the simulation", false, 5.0);
auto& dt = Base::register_argument<double>('d', "timeStepSize", "time step of the simulation", false, 0.001);

int main (int argc, char **argv){

	// todo: Reduce the number of numericalvector creations
	// todo: optimize the use of qinverse
	// todo: optimize the code, placement of variables
	// todo: make timestep CFL dependent
	// todo: reorganise the code into more intuitief constructions
	Base::parse_options(argc, argv);

	logger(WARN,"WARNING: Timestep is determined a priori. Stability Criteria might not be satisfied!");
    // Set parameters for the PDE.
    const Base::MeshType meshType = Base::MeshType::TRIANGULAR;
    const Base::ButcherTableau * const ptrButcherTableau = Base::AllTimeIntegrators::Instance().getRule(3,3,true);

    //Set variable names and number of parameters
    std::vector<std::string> variableNames;
    variableNames.push_back("q1");

    for(std::size_t i = 1; i <= dimension.getValue()+1; i++) // +1 is for the energy equation
    {
        std::string variableName = "q" + std::to_string(i+1);
        variableNames.push_back(variableName);
    }

    // Calculate number of variables
    const std::size_t numOfVariables = 2+dimension.getValue();

    // Create problem solver 'test', that can solve the compressible Navier-Stokes equations.
    CompressibleNavierStokes test(dimension.getValue(), numOfVariables, endTime.getValue(), polynomialOrder.getValue(), ptrButcherTableau);

    // Create the mesh, a simple square domain
    test.createMesh(numOfElements.getValue(), meshType);

    // Solve the problem over time interval [startTime,endTime].
    test.solve(startTime.getValue(), endTime.getValue(), dt.getValue(), numOfOutputFrames.getValue(), true);

    //Compute errors at the end of the simulation
    LinearAlgebra::NumericalVector maxError = test.Error(endTime.getValue());
    std::cout << maxError << std::endl;

	std::cout << "Hi, I am Navier-Stokes, but I am still on holiday, next month I will do my job!" << std::endl;

    return 0;
}
