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

#ifndef SAVAGEHUTTER2DINTERIORTEST_H
#define SAVAGEHUTTER2DINTERIORTEST_H

#include "SavageHutter2DBase.h"

///Class meant to test the convergence of the code. It is mostly copied from
///SavageHutter2DBasic, but with a different source term and boundary conditions
///and without unnecessary features like width-averaging.
class SavageHutter2DInteriorTest : public SavageHutter2DBase
{
public:

    ///\brief Constructor: initialise parent classes and set parameters.
    SavageHutter2DInteriorTest(std::size_t polyOrder, std::size_t numberOfElements);

    ///\brief Create the description of the domain and the mesh.
    Base::RectangularMeshDescriptor<2> createMeshDescription(const std::size_t numOfElementsPerDirection);

    ///\brief Put the initial solution in here.
    LinearAlgebra::MiddleSizeVector getInitialSolution(const PointPhysicalT &pPhys, const double &startTime, const std::size_t orderTimeDerivative = 0) override final;

    ///\brief Put the analytical solution of your system in here. If there is no analytical solution, put in anything and set the flag in main::solve to false.
    LinearAlgebra::MiddleSizeVector getExactSolution(const PointPhysicalT &pPhys, const double &time, const std::size_t orderTimeDerivative = 0) override final;

    ///\brief Compute S in (h,hu,hv)_t + F(h,hu,hv)_x = S(h,hu,hv)
    LinearAlgebra::MiddleSizeVector computeSourceTerm(const LinearAlgebra::MiddleSizeVector &numericalSolution, const PointPhysicalT& pPhys, const double time) override final;

    ///\brief Compute F in (h,hu,hv)_t + F(h,hu,hv)_x = S(h,hu,hv)
    LinearAlgebra::MiddleSizeVector computePhysicalFlux(const LinearAlgebra::MiddleSizeVector &numericalSolution);

private:

};

#endif	/* SAVAGEHUTTER2DBASIC_H */

