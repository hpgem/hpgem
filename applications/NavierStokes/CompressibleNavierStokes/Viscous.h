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

#ifndef VISCOUS_H_
#define VISCOUS_H_

#include "CompressibleDimension.h"

class Viscous {
   public:
    Viscous(CompressibleNavierStokes &instance);

    /// ****************************************
    /// ***      Constitutive functions      ***
    /// ****************************************

    /// \brief Computes the temperature based on the pressure
    double computeTemperature(const LinearAlgebra::MiddleSizeVector &state,
                              const double pressure);

    /// \brief Computes the viscosity as function of temperature, based on
    /// Sutherlands law.
    double computeViscosity(const double temperature);

    /// *******************************************
    /// ***      Elliptic Tensor functions      ***
    /// *******************************************

    /// \brief Computes ATensor_ for a given partialState, viscosity, kappa and
    /// c_v
    std::vector<LinearAlgebra::MiddleSizeMatrix> computeATensor(
        const LinearAlgebra::MiddleSizeVector &partialState,
        const double viscosity);

    /// \brief Computes a second order contraction between an ATensor and a
    /// matrix resulting in a matrix
    LinearAlgebra::MiddleSizeMatrix computeATensorMatrixContraction(
        const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensor,
        const LinearAlgebra::MiddleSizeMatrix &matrix);

    LinearAlgebra::MiddleSizeMatrix computeATensorMatrixContractionFast(
        const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensor,
        const LinearAlgebra::MiddleSizeMatrix &matrix);

    /// *****************************************
    /// ***   Element integration functions   ***
    /// *****************************************

    /// \brief Computes the integrand used in the element integration routine.
    LinearAlgebra::MiddleSizeVector integrandAtElement(
        Base::PhysicalElement<DIM> &element,
        const LinearAlgebra::MiddleSizeVector &state,
        const LinearAlgebra::MiddleSizeMatrix &stateJacobian,
        const double pressure,
        const LinearAlgebra::MiddleSizeVector &partialState);

    /// **************************************************
    /// ***    External face integration functions     ***
    /// **************************************************

    // Computes the flux function required for the stability parameter
    // calculation
    LinearAlgebra::MiddleSizeMatrix computeStabilityFluxFunctionBoundary(
        const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorBoundary,
        const LinearAlgebra::MiddleSizeVector &stateLeft,
        const LinearAlgebra::MiddleSizeVector &stateBoundary,
        const LinearAlgebra::SmallVector<DIM> &unitNormalLeft);

    /// \brief Computes the stability parameters used in the auxilliary
    /// integrand for an external face
    LinearAlgebra::MiddleSizeMatrix computeStabilityParameters(
        Base::PhysicalFace<DIM> &face,
        const LinearAlgebra::MiddleSizeVector &stateLeft,
        const LinearAlgebra::MiddleSizeVector &stateBoundary,
        const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorBoundary,
        const LinearAlgebra::SmallVector<DIM> &unitNormalLeft);

    /// \brief Computes the flux function for the auxilliary value at an
    /// external face
    LinearAlgebra::MiddleSizeMatrix computeAuxilliaryFlux(
        Base::PhysicalFace<DIM> &face,
        const LinearAlgebra::MiddleSizeVector &stateLeft,
        const LinearAlgebra::MiddleSizeVector &stateBoundary,
        const LinearAlgebra::MiddleSizeMatrix &stateJacobianLeft,
        const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorBoundary,
        const LinearAlgebra::SmallVector<DIM> &unitNormalLeft);

    /// \brief Computes the viscous integral at the face, based on the viscous
    /// fluxes for an external face
    LinearAlgebra::MiddleSizeVector integrandViscousAtFace(
        Base::PhysicalFace<DIM> &face,
        const LinearAlgebra::MiddleSizeVector &stateLeft,
        const LinearAlgebra::MiddleSizeVector &stateBoundary,
        const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorLeft,
        const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorBoundary);

    /// \brief Computes the auxilliary integrand at an internal face
    LinearAlgebra::MiddleSizeVector integrandAuxilliaryAtFace(
        Base::PhysicalFace<DIM> &face,
        const LinearAlgebra::MiddleSizeVector &stateLeft,
        const LinearAlgebra::MiddleSizeVector &stateBoundary,
        const LinearAlgebra::MiddleSizeMatrix &stateJacobianLeft,
        const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorBoundary);

    /// **************************************************
    /// ***    Internal face integration functions     ***
    /// **************************************************

    // Computes the fluxFunction in the stability parameter calculation used for
    // the integrand
    LinearAlgebra::MiddleSizeMatrix computeStabilityFluxFunction(
        const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensor,
        const LinearAlgebra::MiddleSizeVector &stateInternal,
        const LinearAlgebra::MiddleSizeVector &stateExternal,
        const LinearAlgebra::SmallVector<DIM> &normalInternal);

    /// \brief Computes the stability parameters used in the auxilliary
    /// integrand for an internal face
    LinearAlgebra::MiddleSizeMatrix computeStabilityParameters(
        Base::PhysicalFace<DIM> &face,
        const LinearAlgebra::MiddleSizeVector &stateLeft,
        const LinearAlgebra::MiddleSizeVector &stateRight,
        const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorLeft,
        const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorRight,
        const LinearAlgebra::SmallVector<DIM> &normalInternal);

    /// \brief Computes the fluxFunction for the auxilliary values at an
    /// internal face
    LinearAlgebra::MiddleSizeMatrix computeAuxilliaryFlux(
        Base::PhysicalFace<DIM> &face,
        const LinearAlgebra::MiddleSizeVector &stateLeft,
        const LinearAlgebra::MiddleSizeVector &stateRight,
        const LinearAlgebra::MiddleSizeMatrix &stateJacobianLeft,
        const LinearAlgebra::MiddleSizeMatrix &stateJacobianRight,
        const LinearAlgebra::SmallVector<DIM> &unitNormalInternal,
        const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorLeft,
        const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorRight);

    /// \brief Compute both the auxilliary face integral for the left element as
    /// the right element at the same time. (reducing flux calculations)
    std::pair<LinearAlgebra::MiddleSizeVector, LinearAlgebra::MiddleSizeVector>
        integrandsAuxilliaryAtFace(
            Base::PhysicalFace<DIM> &face,
            const LinearAlgebra::MiddleSizeVector &stateLeft,
            const LinearAlgebra::MiddleSizeVector &stateRight,
            const LinearAlgebra::MiddleSizeMatrix &stateJacobianLeft,
            const LinearAlgebra::MiddleSizeMatrix &stateJacobianRight,
            const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorLeft,
            const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorRight);

    /// \brief Compute both the viscous face integral for the left element as
    /// the right element at the same time. (reducing flux calculations)
    std::pair<LinearAlgebra::MiddleSizeVector, LinearAlgebra::MiddleSizeVector>
        integrandsViscousAtFace(
            Base::PhysicalFace<DIM> &face,
            const LinearAlgebra::MiddleSizeVector &stateLeft,
            const LinearAlgebra::MiddleSizeVector &stateRight,
            const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorLeft,
            const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorRight);

    /// **************************************************
    /// ***    general face integration functions      ***
    /// **************************************************

    /// \brief Computes the integrand required for the stability parameter
    /// calculations
    LinearAlgebra::MiddleSizeVector integrandStabilityRightHandSideOnFace(
        Base::PhysicalFace<DIM> &face, const Base::Side &side,
        const LinearAlgebra::MiddleSizeMatrix &stabilityFluxFunction,
        const std::size_t iD);

    /// \brief Computes the rhs, for the system of equations solving the
    /// stability parameters, by integrating the rhs  stability parameter
    /// integrand
    LinearAlgebra::MiddleSizeVector computeRhs(
        const Base::Face *ptrFace, const Base::Side &side,
        const LinearAlgebra::MiddleSizeMatrix &stabilityFluxFunction,
        const std::size_t iD);

    /// \brief Sets the mass matrix used in the computation of the stability
    /// parameters
    void setStabilityMassMatrix(
        LinearAlgebra::MiddleSizeMatrix &StabilityMassMatrix);

   private:
    CompressibleNavierStokes &instance_;

    /// \var Integrator voor the stability parameters
    Integration::FaceIntegral<DIM> stabilityFaceIntegrator_;

    LinearAlgebra::MiddleSizeMatrix stabilityMassMatrix_;  // Note: this breaks
                                                           // down if p is not
                                                           // the same in all
                                                           // elements.
};

#endif /* VISCOUS_H_ */
