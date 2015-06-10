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

#ifndef INVISCID_H_
#define INVISCID_H_

#include "CompressibleDimension.h"

class Inviscid
{
public:
	Inviscid(const CompressibleNavierStokes& instance);

    /// *****************************************
    /// ***   Element integration functions   ***
    /// *****************************************

    /// Compute integrand of righthandside on an element
    LinearAlgebra::MiddleSizeVector integrandAtElement(Base::PhysicalElement<DIM>& element, const double &time, const double pressureTerm, const LinearAlgebra::MiddleSizeVector &qSolution);

    /// *****************************************
    /// ***    face integration functions     ***
    /// *****************************************

    /// \brief Compute the Roe Riemann Flux.
    LinearAlgebra::MiddleSizeVector RoeRiemannFluxFunction(const LinearAlgebra::MiddleSizeVector &qReconstructionLeft, const LinearAlgebra::MiddleSizeVector &qReconstructionRight, const LinearAlgebra::SmallVector<DIM>& normal);

    /// \brief Compute the integrand for the right hand side for the reference face corresponding to a boundary face.
	LinearAlgebra::MiddleSizeVector integrandAtFace(Base::PhysicalFace<DIM>& face, const double &time, const LinearAlgebra::MiddleSizeVector &solutionCoefficients);

	/// \brief Compute the integrand for the right hand side for the reference face corresponding to an internal face.
	LinearAlgebra::MiddleSizeVector integrandAtFace(Base::PhysicalFace<DIM>& face, const double &time, const Base::Side &iSide, const LinearAlgebra::MiddleSizeVector &qSolutionInternal, const LinearAlgebra::MiddleSizeVector &qSolutionExternal);

private:
	const CompressibleNavierStokes& instance_;

};

#endif /* INVISCID_H_ */
