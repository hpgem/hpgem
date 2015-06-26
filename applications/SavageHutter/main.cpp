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

#include <chrono>
#include <fstream>

#include "SavageHutter.h"
#include "Base/TimeIntegration/AllTimeIntegrators.h"

#include "Logger.h"

auto& numOfElements = Base::register_argument<std::size_t>('n', "numElems", "number of elements per dimension", false, 10);
auto& polynomialOrder = Base::register_argument<std::size_t>('p', "order", "polynomial order of the solution", false, 2);
auto& numOfOutputFrames = Base::register_argument<std::size_t>('O', "numOfOutputFrames", "Number of frames to output", false, 1);
auto& startTime = Base::register_argument<double>('S', "startTime", "start time of the simulation", false, 0.0);
auto& endTime = Base::register_argument<double>('T', "endTime", "end time of the simulation", false, 0.001);
auto& dt = Base::register_argument<double>('d', "timeStepSize", "time step of the simulation", false, 0.001);

//This code does not work correctly yet!
int main(int argc, char **argv)
{
    Base::parse_options(argc, argv);
    
    // Set parameters for the PDE.
    SHConstructorStruct inputVals;
    //DIM is declared in SavageHutterRightHandSideComputer.h
    inputVals.numOfVariables = 2;    
    inputVals.polyOrder = polynomialOrder.getValue();
    inputVals.numElements = numOfElements.getValue();
    inputVals.meshType = Base::MeshType::RECTANGULAR; // Either TRIANGULAR or RECTANGULAR.
    inputVals.ptrButcherTableau = Base::AllTimeIntegrators::Instance().getRule(1,1);
    
    //Construct the problem and output generator
    SavageHutter test(inputVals);    
    std::vector<std::string> variableNames = {"h", "hu"};
    test.setOutputNames("output", "SavageHutter", "SavageHutter", variableNames);

    // Start measuring elapsed time
    std::chrono::time_point<std::chrono::system_clock> startClock, endClock;
    startClock = std::chrono::system_clock::now();

    // Solve the problem over time interval [startTime,endTime].
    test.solve(startTime.getValue(), endTime.getValue(), dt.getValue(), numOfOutputFrames.getValue(), false);

    // Measure elapsed time
    endClock = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = endClock - startClock;
    logger(INFO, "Elapsed time for solving the PDE: % s", elapsed_seconds.count());

    return 0;
}
