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

#include "SavageHutter1DBidispersed.h"
#include "HelperFunctions.h"
#include "HeightLimiters/PositiveLayerLimiter.h"

///\details In this constructor, some of the parameters for the problem are set.
///Most of these parameters are declared in SavageHutterBase, but since they are protected
///they can also be used here.
SavageHutter1DBidispersed::SavageHutter1DBidispersed(std::size_t polyOrder, std::size_t numberOfElements)
: SavageHutter1DBase(3, polyOrder)
{
    alpha_ = .1;
    chuteAngle_ = M_PI / 180 * 29.6484;
    epsilon_ = .1;
    const PointPhysicalT &pPhys = createMeshDescription(1).bottomLeft_;
    inflowBC_ = getInitialSolution(pPhys, 0);    
    
    std::vector<std::string> variableNames = {"h", "hu", "eta"};
    setOutputNames("output1D", "SavageHutter", "SavageHutter", variableNames);
    
    createMesh(numberOfElements, Base::MeshType::RECTANGULAR);
}

///\details In this function, the mesh gets described. 
///The domain is given by [0, endOfDomain], it consists of numberOfElementsPerDirection
///elements (in this case lines) and has either periodic or non-periodic boundary
///conditions. If the boundary conditions are not periodic, set the BoundaryType
///in this function to Base::BoundaryType::SOLID_WALL and make sure that the function
///SavageHutter1DBase::integrandRightHandSideOnRefFace is set correctly for boundary
///faces. 
Base::RectangularMeshDescriptor<1> SavageHutter1DBidispersed::createMeshDescription(const std::size_t numOfElementsPerDirection)
{
    const double endOfDomain = 8;
    const Base::BoundaryType boundary = Base::BoundaryType::SOLID_WALL;
    return SavageHutter1DBase::createMeshDescription(numOfElementsPerDirection, endOfDomain, boundary);
}

///\details Gives the initial solution for the problem. One could also call getExactSolution
///here to get the analytical solution at the start time.
 LinearAlgebra::MiddleSizeVector SavageHutter1DBidispersed::getInitialSolution(const PointPhysicalT &pPhys, const double &startTime, const std::size_t orderTimeDerivative)
{
    const double h = 0;
    const double hu = 0;
    const double eta = 0;
    return LinearAlgebra::MiddleSizeVector({h, hu, eta});
}

void SavageHutter1DBidispersed::setInflowBC(double time)
{
    const double inflowTurnOnRate = 2;
    const double hIn = 1 - std::exp(-inflowTurnOnRate * time);
    const double uIn = .884404;
    const double phiIn = 1;
    inflowBC_ = LinearAlgebra::MiddleSizeVector({hIn, hIn * uIn, hIn * phiIn});
}

 ///\details Gives the exact (analytical) solution for the test problem. This function
 ///is not necessary, unless the last flag in main::solve has been set to true.
LinearAlgebra::MiddleSizeVector SavageHutter1DBidispersed::getExactSolution(const PointPhysicalT& pPhys, const double& time, const std::size_t orderTimeDerivative)
{
    const double h = 1;
    const double hu = 1;
    const double eta = 1;
    return LinearAlgebra::MiddleSizeVector({h, hu, eta});
}

///\details Constructs the slope limiter, available slope limiters can be found 
///in the folder SlopeLimiters.
SlopeLimiter* SavageHutter1DBidispersed::createSlopeLimiter()
{
    //Little hack: the polynomial order is the number of basis functions minus one.
    return new EmptySlopeLimiter;
}

///\details Constructs the non-negativity limiter, available non-negativity limiters
///can be found in the folder HeightLimiters.
HeightLimiter* SavageHutter1DBidispersed::createHeightLimiter()
{
    return new PositiveLayerLimiter<1>(1e-5);
}

///\details Write the values of the variables to the VTK file. The method of 
///HpgemAPISimplified writes all the unknowns of the system (here h and hu), here
///the concentration phi is also written to the VTK file.
void SavageHutter1DBidispersed::registerVTKWriteFunctions()
{
    HpgemAPISimplified::registerVTKWriteFunctions();

    registerVTKWriteFunction([ = ](Base::Element* element, const Geometry::PointReference<1>& pRef, std::size_t timeLevel) -> double
    {
        if (element->getSolution(timeLevel, pRef)[0] > 1e-5)
            return std::real(element->getSolution(timeLevel, pRef)[2] / element->getSolution(timeLevel, pRef)[0]);
        return 0;
    }, "phi");
}

///\details Compute the source term of the 1D shallow granular flow system, namely
///h(\sin\theta - \mu\sign(u)\cos\theta). It is important to note here that one can
///choose between a Coulomb-type friction law, a Pouliquen friction law as can be found
///in Weinhart et. al. (2012) or an exponential form for the friction. The respective functions
///are defined in the class SavageHutter1DBase.
LinearAlgebra::MiddleSizeVector SavageHutter1DBidispersed::computeSourceTerm(const LinearAlgebra::MiddleSizeVector& numericalSolution, const PointPhysicalT& pPhys, const double time)
{
    logger.assert(chuteAngle_ < M_PI, "Angle must be in radians, not degrees!");
    const double h = numericalSolution(0);
    double sourceX;
    if (h > dryLimit_)
    {
        const double hu = numericalSolution(1);
        const double u = hu / h;
        const double mu = computeFriction(numericalSolution);
        const int signU = Helpers::sign(u);
        sourceX = h * std::sin(chuteAngle_) - h * mu * signU * std::cos(chuteAngle_);
    }
    else
    {
        sourceX = 0;
    }
    
    logger(DEBUG, "Source: %, h: %", sourceX, h);
    return MiddleSizeVector({0, sourceX, 0});
}

///\details Compute the function f(h,hu) = {hu,  hu^2 + h^2/2 \epsilon \cos \theta, eta*u - (1-alpha)*eta*u*(1-eta/h)}
LinearAlgebra::MiddleSizeVector SavageHutter1DBidispersed::computePhysicalFlux(const MiddleSizeVector &numericalSolution, const PointPhysicalT& pPhys)
{
    double h = numericalSolution(0);
    logger.assert(h > -1e-16, "Negative height (%)", h);
    const double hu = numericalSolution(1);
    const double smallHeight = numericalSolution(2);
    MiddleSizeVector flux(3);
    if (h > dryLimit_)
    {
        double u = hu / h;
        double phi = smallHeight / h;

        flux(0) = hu;
        flux(1) = hu * u + epsilon_ / 2 * std::cos(chuteAngle_) * h * h;
        flux(2) = hu * phi * (alpha_ + (1 - alpha_) * phi);
    }
    else
    {
        flux = LinearAlgebra::MiddleSizeVector({hu, 0, 0});
    }
    return flux;
}
