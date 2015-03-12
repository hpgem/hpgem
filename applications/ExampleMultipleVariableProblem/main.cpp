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

#include <fstream>

#include "AcousticWave.h"
#include "AcousticWaveLinear.h"
#include "Base/TimeIntegration/AllTimeIntegrators.h"

#include "Logger.h"


auto& numOfElements = Base::register_argument<std::size_t>('n', "numElems", "number of elements per dimension", true);
auto& polynomialOrder = Base::register_argument<std::size_t>('p', "order", "polynomial order of the solution", true);

// Parse options seem to do weird things. Flag 'h' also communicates with Petsc. Several other flags are also already defined in the library, because they are defined in class HpgemUISimplified.
auto& numOfOutputFrames = Base::register_argument<std::size_t>('O', "numOfOutputFrames", "Number of frames to output", false, 1);
auto& startTime = Base::register_argument<double>('S', "startTime", "start time of the simulation", false, 0.0);
auto& endTime = Base::register_argument<double>('T', "endTime", "end time of the simulation", false, 1.0);
auto& dt = Base::register_argument<double>('d', "timeStepSize", "time step of the simulation", false, 0.01);

int main(int argc, char **argv)
{
    Base::parse_options(argc, argv);
    try
    {
        const bool useLinearSolver = false;
        
        // Set parameters for the PDE.
        const std::size_t dimension = 2;
        const Base::MeshType meshType = Base::MeshType::TRIANGULAR;
        const Base::ButcherTableau * const ptrButcherTableau = Base::AllTimeIntegrators::Instance().getRule(4, 4);
        const double c = 1.0;
        
        std::string variableString;
        if (dimension == 2)
        {
            variableString = "v,s0,s1";
        }
        else if (dimension == 3)
        {
            variableString = "v,s0,s1,s2";
        }
        
        // Compute parameters for PDE
        const std::size_t numOfVariables = dimension + 1;
        
        if(useLinearSolver)
        {
            // Create problem solver 'test', that can solve the acoustic wave equations.
            AcousticWaveLinear test(dimension, numOfVariables, polynomialOrder.getValue(), ptrButcherTableau);
            
            // Create the mesh
            test.createMesh(numOfElements.getValue(), meshType);
            
            // Set the material parameter
            test.setMaterialParameter(c);
            
            // Set the names for the output file
            test.setOutputNames("output","acousticWaveLinear","acousticWaveLinear",variableString);
            
            // Solve the problem over time interval [startTime,endTime].
            test.solve(startTime.getValue(), endTime.getValue(), dt.getValue(), numOfOutputFrames.getValue(), true);
        }
        else
        {
            // Create problem solver 'test', that can solve the acoustic wave equations.
            AcousticWave test(dimension, numOfVariables, polynomialOrder.getValue(), ptrButcherTableau);
            
            // Create the mesh
            test.createMesh(numOfElements.getValue(), meshType);
            
            // Set the material parameter
            test.setMaterialParameter(c);
            
            // Set the names for the output file
            test.setOutputNames("output","acousticWave","acousticWave",variableString);
            
            // Solve the problem over time interval [startTime,endTime].
            test.solve(startTime.getValue(), endTime.getValue(), dt.getValue(), numOfOutputFrames.getValue(), true);
            
        }
        
        return 0;
    }
    catch (const char* e)
    {
        std::cout << e;
    }
}
