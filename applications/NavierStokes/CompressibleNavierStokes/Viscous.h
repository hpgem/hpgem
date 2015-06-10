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

#include "CompressibleDimension.h"

class Viscous
{
public:

	Viscous(CompressibleNavierStokes& instance);

	void setInverseStabilityMassMatrix(LinearAlgebra::MiddleSizeMatrix &inverseStabilityMassMatrix);

    /// *****************************************
    /// ***      flux function functions      ***
    /// *****************************************

	/// Computes the temperature based on the pressure
	double computeTemperature(const LinearAlgebra::MiddleSizeVector qSolution, const double pressure);

	/// Computes the viscosity as function of temperature, based on Sutherlands law.
	double computeViscosity(double temperature);

	/// Computes ATensor_ for a given partialState (containing velocities and total Energy), viscosity, kappa and c_v
	std::vector<LinearAlgebra::MiddleSizeMatrix> computeATensor( const LinearAlgebra::MiddleSizeVector partialState, const double viscosity);

	LinearAlgebra::MiddleSizeMatrix computeATensorMatrixContraction(const std::vector<LinearAlgebra::MiddleSizeMatrix> ATensor, const LinearAlgebra::MiddleSizeMatrix matrix);

    /// *****************************************
    /// ***   Element integration functions   ***
    /// *****************************************

	/// Computes the integrand used in the element integration routine.
	LinearAlgebra::MiddleSizeVector integrandAtElement(Base::PhysicalElement<DIM>& element, const LinearAlgebra::MiddleSizeVector qSolution, const LinearAlgebra::MiddleSizeMatrix qSolutionJacobian, const double pressure, const LinearAlgebra::MiddleSizeVector partialState);

    /// *****************************************
    /// ***    face integration functions     ***
    /// *****************************************

	//Computes the fluxFunction in the stability parameter calculation used for the integrand
	void computeStabilityFluxFunction(const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorInternal, const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorExternal, const LinearAlgebra::MiddleSizeVector &stateInternal, const LinearAlgebra::MiddleSizeVector &stateExternal, const LinearAlgebra::MiddleSizeVector &normalInternal);

	//Computes the integrand required for the stability parameter calculations
	LinearAlgebra::MiddleSizeVector integrandStabilityRightHandSideOnRefFace(Base::PhysicalFace<DIM>& face, const Base::Side side, const LinearAlgebra::MiddleSizeMatrix stabilityFluxFunction, const std::size_t iD);

	//Computes the rhs by integrating the rhs integrand
	LinearAlgebra::MiddleSizeVector computeRhs(const Base::Face* face, const Base::Side side, const LinearAlgebra::MiddleSizeMatrix stabilityFluxFunction, const std::size_t iD);

	//Computes the stability parameters
	LinearAlgebra::MiddleSizeMatrix computeStabilityParameters(const Base::Face *ptrFace, Base::Side iSide, const std::vector<LinearAlgebra::MiddleSizeMatrix> ATensorInternal,	const std::vector<LinearAlgebra::MiddleSizeMatrix> ATensorExternal, const LinearAlgebra::MiddleSizeVector stateInternal, const LinearAlgebra::MiddleSizeVector stateExternal, const LinearAlgebra::SmallVector<DIM> normalInternal, const Geometry::PointReference<DIM - 1> &pRef);

	/// computes the fluxFunction for the auxilliary values at the face
	LinearAlgebra::MiddleSizeMatrix computeAuxilliaryFlux(const Base::Face *ptrFace, Base::Side iSide, const LinearAlgebra::MiddleSizeVector stateInternal, const LinearAlgebra::MiddleSizeVector stateExternal, const double pressureInternal, const double pressureExternal, const LinearAlgebra::MiddleSizeVector partialStateInternal, const LinearAlgebra::MiddleSizeVector partialStateExternal, const LinearAlgebra::MiddleSizeMatrix stateJacobianInternal, const LinearAlgebra::MiddleSizeMatrix stateJacobianExternal,  const LinearAlgebra::MiddleSizeVector normalInternal, const Geometry::PointReference<DIM - 1> &pRef);

	/// computes the viscous integral at the face, based on the viscous fluxes
	LinearAlgebra::MiddleSizeVector integrandViscousAtFace(Base::PhysicalFace<DIM>& face, const Base::Side &iSide, LinearAlgebra::MiddleSizeVector qSolutionInternal, LinearAlgebra::MiddleSizeVector qSolutionExternal, double pressure, LinearAlgebra::MiddleSizeVector partialState);

	/// computes the auxilliary integral at the face, based on the auxilliary values
	LinearAlgebra::MiddleSizeVector integrandAuxilliaryAtFace(Base::PhysicalFace<DIM>& face, const Base::Side &iSide, LinearAlgebra::MiddleSizeVector stateInternal, const LinearAlgebra::MiddleSizeVector stateExternal, const double pressureInternal, const double pressureExternal, const LinearAlgebra::MiddleSizeVector partialStateInternal, const LinearAlgebra::MiddleSizeVector partialStateExternal, const LinearAlgebra::MiddleSizeMatrix stateJacobianInternal, const LinearAlgebra::MiddleSizeMatrix stateJacobianExternal);
private:
	CompressibleNavierStokes& instance_;

	const double PrInv_; //Inverse of Pr
	const double cp_;
	std::vector<LinearAlgebra::MiddleSizeMatrix> ATensorInternal_;
	std::vector<LinearAlgebra::MiddleSizeMatrix> ATensorExternal_;
	LinearAlgebra::MiddleSizeMatrix stabilityMassMatrix_; //Note: this breaks down if p is not the same in all elements.
	LinearAlgebra::MiddleSizeMatrix stabilityFluxFunctionInternal_;
	LinearAlgebra::MiddleSizeMatrix stabilityFluxFunctionExternal_;
};

#endif /* VISCOUS_H_ */
