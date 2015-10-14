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

#ifndef SavageHutterH
#define SavageHutterH

#include "SavageHutterBase.h"

class SavageHutter : public SavageHutterBase
{
public:
    
    ///\brief Constructor that takes an object specially designed to contain all values needed for construction of this kind of problem.
    SavageHutter(const SHConstructorStruct& inputValues);
    
private:       
    ///\brief Create the slope limiter that will be used in this simulation.
    SlopeLimiter * createSlopeLimiter(const SHConstructorStruct &inputValues) override final;
    
    ///\brief Create the non-negativity limiter that will be used in this simulation.
    HeightLimiter * createHeightLimiter(const SHConstructorStruct &inputValues) override final;
    
    ///\brief Create the object that can compute the right hand side of the differential equation for this simulation.
    RightHandSideComputer * createRightHandSideComputer(const SHConstructorStruct &inputValues) override final;

    ///\brief Compute the initial solution at a given point in space and time.
    LinearAlgebra::MiddleSizeVector getInitialSolution(const PointPhysicalT &pPhys, const double &startTime, const std::size_t orderTimeDerivative = 0) override final;
    
    ///\brief Show the progress of the time integration.
    void showProgress(const double time, const std::size_t timeStepID) override final;
    
    LinearAlgebra::MiddleSizeVector getExactSolution(const PointPhysicalT &pPhys, const double &time, const std::size_t orderTimeDerivative = 0) override final;
    
    void registerVTKWriteFunctions() override final;
    
    void setInflowBC(double time) override final
    {
        
    }

};

#endif
