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

#include "SavageHutter.h"
#include "SlopeLimiters/TvbLimiterWithDetector1D.h"
#include "HeightLimiters/PositiveLayerLimiter.h"
#include "SavageHutterRightHandSideComputer.h"
#include "HeightLimiters/EmptyHeightLimiter.h"
#include "SlopeLimiters/EmptySlopeLimiter.h"
#include "SavageHutterRHS2D.h"
#include "SlopeLimiters/TvbLimiter1D.h"
#include "HeightLimiters/SqueezeLimiterWithLayer.h"
#include "HeightLimiters/AverageValuesNonNegativeLimiter.h"

using LinearAlgebra::MiddleSizeVector;

///\details Construct the simulation by constructing the basic functionality of SavageHutterBase and constructing the right hand side computer, slope limiter and non-negativity limiter.
SavageHutter::SavageHutter(const SHConstructorStruct& inputValues) :
SavageHutterBase(inputValues)
{
    rhsComputer_ = createRightHandSideComputer(inputValues);
    slopeLimiter_ = createSlopeLimiter(inputValues);
    heightLimiter_ = createHeightLimiter(inputValues);
}

///\details Actual creation of the slope limiter. The delete is called in the class SavageHutterBase, since that's also where the slope limiter resides.
SlopeLimiter * SavageHutter::createSlopeLimiter(const SHConstructorStruct &inputValues)
{
    const PointPhysicalT &pPhys = createMeshDescription(1).bottomLeft_;
    LinearAlgebra::MiddleSizeVector inflowBC = getInitialSolution(pPhys, 0.);
    return (new EmptySlopeLimiter());
    //return (new TvbLimiter1D(numOfVariables_));
    //return new TvbLimiterWithDetector1D(numOfVariables_, inflowBC, inputValues.polyOrder);
}

///\details Actual creation of the non-negativity limiter. The delete is called in the class SavageHutterBase, since that's also where the non-negativity limiter resides.
HeightLimiter * SavageHutter::createHeightLimiter(const SHConstructorStruct& inputValues)
{
    if (inputValues.polyOrder == 0)
        return new EmptyHeightLimiter();
    return new PositiveLayerLimiter(1e-5);
}

///\details Actual creation of the right hand side computer. The delete is called in the class SavageHutterBase, since that's also where the right hand side computer resides.
RightHandSideComputer * SavageHutter::createRightHandSideComputer(const SHConstructorStruct& inputValues)
{
    const PointPhysicalT &pPhys = createMeshDescription(1).bottomLeft_;
    LinearAlgebra::MiddleSizeVector inflowBC = getInitialSolution(pPhys, 0);
    //magic numbers: epsilon and chute angle (in radians)
    if (DIM == 1)
        return new SavageHutterRightHandSideComputer(inputValues.numOfVariables, 0.1, 32./180*M_PI, inflowBC);
    
    return new SavageHutterRHS2D(inputValues.numOfVariables, 1, 0., inflowBC);
}



///\details Show the number of time steps that have been computed on the console.
void SavageHutter::showProgress(const double time, const std::size_t timeStepID)
{
    if (timeStepID % 100 == 0)
    {
        logger(INFO, "% time steps computed.", timeStepID);
    }
}

void SavageHutter::applyInflowConditions()
{
    if (false && DIM == 1)
    {
        Base::Element *element0 = meshes_[0]->getElementsList()[0];
        const LinearAlgebra::MiddleSizeVector hCoeffs = Helpers::projectOnBasisFuns<DIM>(element0, [ = ](const PointReferenceT & pRef){return getInitialSolution(createMeshDescription(1).bottomLeft_, 0)(0);}, elementIntegrator_);
        element0->setTimeLevelData(0, 0, hCoeffs);
        const LinearAlgebra::MiddleSizeVector huCoeffs = Helpers::projectOnBasisFuns<DIM>(element0, [ = ](const PointReferenceT & pRef){return getInitialSolution(createMeshDescription(1).bottomLeft_, 0)(1);}, elementIntegrator_);
        element0->setTimeLevelData(0, 1, huCoeffs);
    }
}

LinearAlgebra::MiddleSizeVector SavageHutter::getInitialSolution(const PointPhysicalT& pPhys, const double& startTime, const std::size_t orderTimeDerivative)
{
    double h = 0;
    const double x = pPhys[0];
    if (x > 0.5 && x < 2.5)
        h = 1.-(x-1.5) * (x-1.5);
    LinearAlgebra::MiddleSizeVector initialCondition(numOfVariables_);
    initialCondition(0) = h;
    initialCondition(1) = 0;
    return initialCondition;
}

/*
//Dam break
///\details Given the physical point in the domain and the start time of the simulation, compute the initial conditions for the simulation
LinearAlgebra::MiddleSizeVector SavageHutter::getInitialSolution(const PointPhysicalT &pPhys, const double &startTime, const std::size_t orderTimeDerivative)
{    
    return getExactSolution(pPhys, startTime);
}

LinearAlgebra::MiddleSizeVector SavageHutter::getExactSolution(const PointPhysicalT& pPhys, const double& time, const std::size_t orderTimeDerivative)
{
    const double c0 = 0.7;
    const double x = pPhys[0];
    double h;
    double u;
    if (x <= 0.5 - c0 * time)
    {
        h = c0*c0;
        u = 0;
    }
    else
    {
        if (x >= 0.5 + 2 * c0 * time)
        {
            h = 0;
            u = 0;
        }
        else
        {
            h = 1./9. * ((x-0.5)/time - 2.*c0) * ((x-0.5)/time - 2.*c0);
            
            u = 2./3 * ((x-0.5)/time + c0);
            logger(DEBUG, "analytical solution on x = %: %, %", x, h, h*u);
        }
    }
    return LinearAlgebra::MiddleSizeVector({h, h*u, 0});
}

void SavageHutter::registerVTKWriteFunctions()
{
    for (std::size_t iV = 0; iV < configData_->numberOfUnknowns_; iV++)
    {
        registerVTKWriteFunction([ = ](Base::Element* element, const Geometry::PointReference<DIM>& pRef, std::size_t timeLevel) -> double
        {
            return std::real(element->getSolution(timeLevel, pRef)[iV]);
        }, variableNames_[iV]);

    }
    registerVTKWriteFunction([ = ](Base::Element* element, const Geometry::PointReference<DIM>& pRef, std::size_t timeLevel) -> double
    {
        const PointPhysicalT &pPhys = element->referenceToPhysical(pRef);
                             return std::real(getExactSolution(pPhys, time_)[0]);
    }, "analytical height");
    registerVTKWriteFunction([ = ](Base::Element* element, const Geometry::PointReference<DIM>& pRef, std::size_t timeLevel) -> double
    {
        const PointPhysicalT &pPhys = element->referenceToPhysical(pRef);
                             return std::real(getExactSolution(pPhys, time_)[1]);
    }, "analytical discharge");
}
*/