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

#ifndef HPGEM_APP_SAVAGEHUTTER1DBIDISPERSED_H
#define HPGEM_APP_SAVAGEHUTTER1DBIDISPERSED_H

#include "SavageHutter1DBase.h"

///\brief Class to solve the bidispersed Savage-Hutter equations with
/// segregation equation. \details The idea is that almost everything you want
/// to change in your application can be changed in this class: the domain,
/// initial
/// solution, parameter values and output functions are all described in this
/// class. It is also possible to construct the different limiters here.
/// Furthermore, you can choose here between the different types of friction
/// that are given in SavageHutter1DBase. The other function in here is
/// computePhysicalFlux, which is the function F in (h,hu)_t + F(h,hu)_x =
/// S(h,hu). This is described here because it is different for this system than
/// for example for the basic system or the width-averaged system This system is
/// not stable, so this application is not usable yet.
class SavageHutter1DBidispersed : public SavageHutter1DBase {
   public:
    ///\brief Constructor: initialise parent classes and set parameters.
    SavageHutter1DBidispersed(std::size_t polyOrder, std::string meshName);

    ///\brief Put the initial solution in here.
    LinearAlgebra::MiddleSizeVector getInitialSolution(
        const PointPhysicalT &pPhys, const double &startTime,
        const std::size_t orderTimeDerivative = 0) override final;

    ///\brief Put the analytical solution of your system in here. If there is no
    /// analytical solution, put in anything and set the flag in main::solve to
    /// false.
    LinearAlgebra::MiddleSizeVector getExactSolution(
        const PointPhysicalT &pPhys, const double &time,
        const std::size_t orderTimeDerivative = 0) override final;

    ///\brief Set which functions should be written in the VTK output file
    void registerVTKWriteFunctions();

    ///\brief Construct the slope limiter that will be used in this application.
    SlopeLimiter *createSlopeLimiter() override final;

    ///\brief Construct the non-negativity limiter that will be used in this
    /// application.
    HeightLimiter *createHeightLimiter() override final;

    ///\brief Compute S in (h,hu)_t + F(h,hu)_x = S(h,hu)
    LinearAlgebra::MiddleSizeVector computeSourceTerm(
        const LinearAlgebra::MiddleSizeVector &numericalSolution,
        const PointPhysicalT &pPhys, const double time) override final;

    ///\brief Compute F in (h,hu)_t + F(h,hu)_x = S(h,hu)
    LinearAlgebra::MiddleSizeVector computePhysicalFlux(
        const LinearAlgebra::MiddleSizeVector &numericalSolution,
        const PointPhysicalT &pPhys);

    void setInflowBC(double time) override final;

   private:
    double computeFrictionBidispersed(
        const LinearAlgebra::MiddleSizeVector &numericalSolution);
    double computeFrictionExponentialBidispersed(
        const LinearAlgebra::MiddleSizeVector &numericalSolution);

    void tasksAfterTimeStep() override final;

    void tasksAfterSolving() override final;

    /// shape factor of the velocity of the flow, 0<=alpha_<=1 .
    /// This is a different alpha_ than in the basic application!
    double alpha_;
    std::vector<double> maximumHeights_;
};
#endif // HPGEM_APP_SAVAGEHUTTER1DBIDISPERSED_H
