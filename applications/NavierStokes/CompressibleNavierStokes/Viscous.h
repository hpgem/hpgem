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

#ifndef VISCOUS_H_
#define VISCOUS_H_

class Viscous
{
public:
	Viscous(const CompressibleNavierStokes& instance);

    /// *****************************************
    /// ***      flux function functions      ***
    /// *****************************************

	/// Computes the temperature based on the pressure
	double computeTemperature(const LinearAlgebra::NumericalVector qSolution, const double pressure);

	/// Computes the viscosity as function of temperature, based on Sutherlands law.
	double computeViscosity(double temperature);

	/// Computes ATensor_ for a given partialState (containing velocities and total Energy), viscosity, kappa and c_v
	void computeATensor(const LinearAlgebra::NumericalVector partialState, const double viscosity);

	LinearAlgebra::Matrix computeATensorMatrixContraction(LinearAlgebra::Matrix matrix);

    /// *****************************************
    /// ***   Element integration functions   ***
    /// *****************************************

	/// Computes the integrand used in the element integration routine.
	LinearAlgebra::NumericalVector integrandAtElement(const Base::Element *ptrElement, const LinearAlgebra::NumericalVector qSolution, const LinearAlgebra::Matrix qSolutionJacobian, const double pressure, const LinearAlgebra::NumericalVector partialState, const Geometry::PointReference &pRef);

    /// *****************************************
    /// ***    face integration functions     ***
    /// *****************************************

	LinearAlgebra::NumericalVector integrandViscousAtFace(const Base::Face *ptrFace, const Base::Side &iSide, LinearAlgebra::NumericalVector qSolutionInternal, LinearAlgebra::NumericalVector qSolutionExternal, double pressure, LinearAlgebra::NumericalVector partialState, const LinearAlgebra::NumericalVector normal, const Geometry::PointReference &pRef);

private:
	const CompressibleNavierStokes& instance_;

	const double PrInv_; //Inverse of Pr
	const double cp_;
	std::vector<LinearAlgebra::Matrix> ATensor_;
};

#endif /* VISCOUS_H_ */
