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

    /// ****************************************
    /// ***      Constitutive functions      ***
    /// ****************************************

	/// \brief Computes the temperature based on the pressure
	double computeTemperature(const LinearAlgebra::MiddleSizeVector &state, const double pressure);

	/// \brief Computes the viscosity as function of temperature, based on Sutherlands law.
	double computeViscosity(const double temperature);


    /// *******************************************
    /// ***      Elliptic Tensor functions      ***
    /// *******************************************

	/// \brief Computes ATensor_ for a given partialState, viscosity, kappa and c_v
	std::vector<LinearAlgebra::MiddleSizeMatrix> computeATensor( const LinearAlgebra::MiddleSizeVector &partialState, const double viscosity);

	/// \brief Computes a second order contraction between an ATensor and a matrix resulting in a matrix
	LinearAlgebra::MiddleSizeMatrix computeATensorMatrixContraction(const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensor, const LinearAlgebra::MiddleSizeMatrix &matrix);

	LinearAlgebra::MiddleSizeMatrix computeATensorMatrixContractionFast(const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensor, const LinearAlgebra::MiddleSizeMatrix &matrix);


    /// *****************************************
    /// ***   Element integration functions   ***
    /// *****************************************

	/// \brief Computes the integrand used in the element integration routine.
	LinearAlgebra::MiddleSizeVector integrandAtElement(Base::PhysicalElement<DIM> &element, const LinearAlgebra::MiddleSizeVector &state, const LinearAlgebra::MiddleSizeMatrix &stateJacobian, const double pressure, const LinearAlgebra::MiddleSizeVector &partialState);


    /// **************************************************
    /// ***    External face integration functions     ***
    /// **************************************************

	//todo: ATensor matrix might be a small matrix.
	//todo: check if the order of boundary and internal state inputs are correct
	/// \brief Computes the fluxFunction in the stability parameter calculation used for the external stability parameter
	void computeStabilityFluxFunction(const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorExternal, const LinearAlgebra::MiddleSizeVector &stateInternal, const LinearAlgebra::MiddleSizeVector &stateExternal, const LinearAlgebra::SmallVector<DIM> &normalInternal);

	/// \brief Computes the stability parameters used in the external auxilliary integrand
	LinearAlgebra::MiddleSizeMatrix computeStabilityParameters(Base::PhysicalFace<DIM> &face, const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorExternal, const LinearAlgebra::MiddleSizeVector &stateInternal, const LinearAlgebra::MiddleSizeVector &stateExternal, const LinearAlgebra::SmallVector<DIM> &normalInternal);

	/// \brief Computes the flux function for the auxilliary value at an external face
	LinearAlgebra::MiddleSizeMatrix computeAuxilliaryFlux(Base::PhysicalFace<DIM> &face, const double temperatureExternal, const LinearAlgebra::MiddleSizeVector &stateInternal, const LinearAlgebra::MiddleSizeVector &stateExternal, const LinearAlgebra::MiddleSizeVector &partialStateExternal, const LinearAlgebra::MiddleSizeMatrix &stateJacobianInternal, const LinearAlgebra::SmallVector<DIM> &normalInternal);

	/// \brief Computes the auxilliary integrand at the face, based on the auxilliary flux
	LinearAlgebra::MiddleSizeVector integrandAuxilliaryAtFace(Base::PhysicalFace<DIM> &face, const double temperatureExternal, const LinearAlgebra::MiddleSizeVector &stateInternal, const LinearAlgebra::MiddleSizeVector &stateBoundary, const LinearAlgebra::MiddleSizeVector &partialStateExternal, const LinearAlgebra::SmallVector<DIM> &normalInternal, const LinearAlgebra::MiddleSizeMatrix &stateJacobianInternal);

	/// \brief
	///todo: integrand viscous at boundary write it here

    /// **************************************************
    /// ***    Internal face integration functions     ***
    /// **************************************************

	//Computes the fluxFunction in the stability parameter calculation used for the integrand
	void computeStabilityFluxFunction(const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorInternal, const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorExternal, const LinearAlgebra::MiddleSizeVector &stateInternal, const LinearAlgebra::MiddleSizeVector &stateExternal, const LinearAlgebra::SmallVector<DIM> &normalInternal);

	/// \brief Computes the stability parameters used in the auxilliary integrand for an internal face
	LinearAlgebra::MiddleSizeMatrix computeStabilityParameters(Base::PhysicalFace<DIM> &face, const Base::Side &iSide, const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorInternal,	const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorExternal, const LinearAlgebra::MiddleSizeVector &stateInternal, const LinearAlgebra::MiddleSizeVector &stateExternal, const LinearAlgebra::SmallVector<DIM> &normalInternal);

	/// \brief Computes the fluxFunction for the auxilliary values at an internal face
	LinearAlgebra::MiddleSizeMatrix computeAuxilliaryFlux(Base::PhysicalFace<DIM> &face, const Base::Side &iSide, const LinearAlgebra::MiddleSizeVector &stateInternal, const LinearAlgebra::MiddleSizeVector &stateExternal, const double pressureInternal, const double pressureExternal, const LinearAlgebra::MiddleSizeVector &partialStateInternal, const LinearAlgebra::MiddleSizeVector &partialStateExternal, const LinearAlgebra::MiddleSizeMatrix &stateJacobianInternal, const LinearAlgebra::MiddleSizeMatrix &stateJacobianExternal,  const LinearAlgebra::SmallVector<DIM> &unitNormalInternal);

	/// \brief Computes the auxilliary integrand at an internal face
	LinearAlgebra::MiddleSizeVector integrandAuxilliaryAtFace(Base::PhysicalFace<DIM> &face, const Base::Side &iSide, const LinearAlgebra::MiddleSizeVector &stateInternal, const LinearAlgebra::MiddleSizeVector &stateExternal, const double pressureInternal, const double pressureExternal, const LinearAlgebra::MiddleSizeVector &partialStateInternal, const LinearAlgebra::MiddleSizeVector &partialStateExternal, const LinearAlgebra::MiddleSizeMatrix &stateJacobianInternal, const LinearAlgebra::MiddleSizeMatrix &stateJacobianExternal, const LinearAlgebra::SmallVector<DIM> &unitNormalInternal);

	/// \brief Computes the viscous integral at the face, based on the viscous fluxes
	LinearAlgebra::MiddleSizeVector integrandViscousAtFace(Base::PhysicalFace<DIM> &face, const Base::Side &iSide, const LinearAlgebra::MiddleSizeVector &stateInternal, const LinearAlgebra::MiddleSizeVector &stateExternal, const double pressure, const LinearAlgebra::MiddleSizeVector &partialStateInternal, const LinearAlgebra::SmallVector<DIM> &unitNormalInternal);

	/// \brief Compute both the auxilliary face integral for the left element as the right element at the same time. (reducing flux calculations)
	std::pair<LinearAlgebra::MiddleSizeVector,LinearAlgebra::MiddleSizeVector> integrandsAuxilliaryAtFace(Base::PhysicalFace<DIM> &face, const LinearAlgebra::MiddleSizeVector &stateLeft, const LinearAlgebra::MiddleSizeVector &stateRight, const double pressureLeft, const double pressureRight, const LinearAlgebra::MiddleSizeVector &partialStateLeft, const LinearAlgebra::MiddleSizeVector &partialStateRight, const LinearAlgebra::MiddleSizeMatrix &stateJacobianLeft, const LinearAlgebra::MiddleSizeMatrix &stateJacobianRight);

	/// \brief Compute both the viscous face integral for the left element as the right element at the same time. (reducing flux calculations)
	std::pair<LinearAlgebra::MiddleSizeVector,LinearAlgebra::MiddleSizeVector> integrandsViscousAtFace(Base::PhysicalFace<DIM> &face, const LinearAlgebra::MiddleSizeVector &stateLeft, const LinearAlgebra::MiddleSizeVector &stateRight, double pressureLeft, double pressureRight, const LinearAlgebra::MiddleSizeVector &partialStateLeft, const LinearAlgebra::MiddleSizeVector &partialStateRight);

    /// **************************************************
    /// ***    general face integration functions      ***
    /// **************************************************

	/// \brief Computes the integrand required for the stability parameter calculations
	LinearAlgebra::MiddleSizeVector integrandStabilityRightHandSideOnFace(Base::PhysicalFace<DIM> &face, const Base::Side &side, const LinearAlgebra::MiddleSizeMatrix &stabilityFluxFunction, const std::size_t iD, const LinearAlgebra::SmallVector<DIM> &normalInternal);

	/// \brief Computes the rhs, for the system of equations solving the stability parameters, by integrating the rhs  stability parameter integrand
	LinearAlgebra::MiddleSizeVector computeRhs(const Base::Face *ptrFace, const Base::Side &side, const LinearAlgebra::MiddleSizeMatrix &stabilityFluxFunction, const LinearAlgebra::SmallVector<DIM> &normalInternal, const std::size_t iD);

	/// \brief Sets the mass matrix used in the computation of the stability parameters
	void setStabilityMassMatrix(LinearAlgebra::MiddleSizeMatrix &StabilityMassMatrix);


private:
	CompressibleNavierStokes& instance_;

	/// \var Integrator voor the stability parameters
	Integration::FaceIntegral<DIM> stabilityFaceIntegrator_;

	const double PrInv_; //Inverse of Pr
	const double cp_;
	std::vector<LinearAlgebra::MiddleSizeMatrix> ATensorInternal_;
	std::vector<LinearAlgebra::MiddleSizeMatrix> ATensorExternal_;
	std::vector<LinearAlgebra::MiddleSizeMatrix> ATensorLeft_;
	std::vector<LinearAlgebra::MiddleSizeMatrix> ATensorRight_;
	LinearAlgebra::MiddleSizeMatrix stabilityMassMatrix_; //Note: this breaks down if p is not the same in all elements.
	LinearAlgebra::MiddleSizeMatrix stabilityFluxFunctionInternal_;
	LinearAlgebra::MiddleSizeMatrix stabilityFluxFunctionExternal_;
};

#endif /* VISCOUS_H_ */
