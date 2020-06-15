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

#ifndef VISCOUSTERMS_H_
#define VISCOUSTERMS_H_

#include "NavierStokesConstants.h"
#include "UnsteadyNavierStokesAPI.h"

template <std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
class ViscousTerms {
   public:
    ViscousTerms(UnsteadyNavierStokesAPI<DIM, NUMBER_OF_VARIABLES> *instance)
        : instance_(instance) {}

    virtual ~ViscousTerms() {}

    /// *****************************************
    /// ***   Element integration functions   ***
    /// *****************************************

    /// \brief Computes the integrand used in the element integration routine.
    LinearAlgebra::MiddleSizeVector integrandAtElement(
        Base::PhysicalElement<DIM> &element,
        const StateCoefficientsStruct<DIM, NUMBER_OF_VARIABLES>
            &elementStateStruct,
        const double time);

    /// **************************************************
    /// ***    External face integration functions     ***
    /// **************************************************

    /// \brief Computes the flux function for the auxilliary value at an
    /// external face
    LinearAlgebra::MiddleSizeMatrix fluxAuxilliaryBoundary(
        Base::PhysicalFace<DIM> &face,
        const StateCoefficientsStruct<DIM, NUMBER_OF_VARIABLES>
            &faceStateStructBoundary,
        const StateCoefficientsStruct<DIM, NUMBER_OF_VARIABLES>
            &faceStateStructLeft,
        const double &time);

    /// \brief Computes the auxilliary integrand at an internal face
    LinearAlgebra::MiddleSizeVector integrandAuxilliaryAtBoundaryFace(
        Base::PhysicalFace<DIM> &face, const BoundaryType boundaryType,
        const StateCoefficientsStruct<DIM, NUMBER_OF_VARIABLES>
            &faceStateStructBoundary,
        const StateCoefficientsStruct<DIM, NUMBER_OF_VARIABLES>
            &faceStateStructLeft,
        const double &time);

    /// \brief Computes the viscous integral at the face, based on the viscous
    /// fluxes for an external face
    LinearAlgebra::MiddleSizeVector integrandViscousAtBoundaryFace(
        Base::PhysicalFace<DIM> &face, const BoundaryType boundaryType,
        const StateCoefficientsStruct<DIM, NUMBER_OF_VARIABLES>
            &faceStateStructBoundary,
        const StateCoefficientsStruct<DIM, NUMBER_OF_VARIABLES>
            &faceStateStructLeft,
        const double &time);

    /// ***************************************************
    /// ***   External Stability Parameter functions    ***
    /// ***************************************************

    /// \brief Computes the fluxFunction in the stability parameter calculation
    /// used for the integrand
    LinearAlgebra::MiddleSizeMatrix fluxStabilityParameters(
        const StateCoefficientsStruct<DIM, NUMBER_OF_VARIABLES>
            &faceStateStructBoundary,
        const LinearAlgebra::MiddleSizeVector &stateLeft,
        const LinearAlgebra::SmallVector<DIM> &normalInternal);

    /// \brief Computes the integrand required for the stability parameter
    /// calculations
    LinearAlgebra::MiddleSizeMatrix
        integrandStabilityRightHandSideOnBoundaryFace(
            Base::PhysicalFace<DIM> &face,
            const LinearAlgebra::MiddleSizeVector stateCoefficientsLeft,
            const double time);

    /// \brief Computes the rhs, for the system of equations solving the
    /// stability parameters, by integrating the rhs  stability parameter
    /// integrand
    LinearAlgebra::MiddleSizeMatrix computeRhsStabilityParametersBoundary(
        const Base::Face *ptrFace,
        const LinearAlgebra::MiddleSizeVector stateCoefficientsLeft,
        const double time);

    /// \brief Computes the stability parameters used in the auxilliary
    /// integrand for an internal face for given stabilityParameterflux
    LinearAlgebra::MiddleSizeMatrix computeStabilityParametersBoundary(
        Base::PhysicalFace<DIM> &face,
        const LinearAlgebra::MiddleSizeVector stateCoefficientsLeft,
        const double time);

    /// **************************************************
    /// ***    Internal face integration functions     ***
    /// **************************************************

    /// \brief Computes the fluxFunction for the auxilliary values at an
    /// internal face
    LinearAlgebra::MiddleSizeMatrix fluxAuxilliary(
        Base::PhysicalFace<DIM> &face,
        const StateCoefficientsStruct<DIM, NUMBER_OF_VARIABLES>
            &faceStateStructLeft,
        const StateCoefficientsStruct<DIM, NUMBER_OF_VARIABLES>
            &faceStateStructRight);

    /// \brief Compute both the auxilliary face integral for the left element as
    /// the right element at the same time. (reducing flux calculations)
    std::pair<LinearAlgebra::MiddleSizeVector, LinearAlgebra::MiddleSizeVector>
        integrandsAuxilliaryAtFace(
            Base::PhysicalFace<DIM> &face,
            const StateCoefficientsStruct<DIM, NUMBER_OF_VARIABLES>
                &faceStateStructLeft,
            const StateCoefficientsStruct<DIM, NUMBER_OF_VARIABLES>
                &faceStateStructRight,
            const double time);

    /// \brief Compute both the viscous face integral for the left element as
    /// the right element at the same time.
    std::pair<LinearAlgebra::MiddleSizeVector, LinearAlgebra::MiddleSizeVector>
        integrandsViscousAtFace(
            Base::PhysicalFace<DIM> &face,
            const StateCoefficientsStruct<DIM, NUMBER_OF_VARIABLES>
                &faceStateStructLeft,
            const StateCoefficientsStruct<DIM, NUMBER_OF_VARIABLES>
                &faceStateStructRight,
            const double time);

    /// ***************************************************
    /// ***   Internal Stability Parameter functions    ***
    /// ***************************************************

    /// \brief Computes the fluxFunction in the stability parameter calculation
    /// used for the integrand
    LinearAlgebra::MiddleSizeMatrix fluxStabilityParameters(
        const StateCoefficientsStruct<DIM, NUMBER_OF_VARIABLES> faceStateStruct,
        const LinearAlgebra::MiddleSizeVector &stateInternal,
        const LinearAlgebra::MiddleSizeVector &stateExternal,
        const LinearAlgebra::SmallVector<DIM> &normalInternal);

    /// \brief Computes the integrand required for the stability parameter
    /// calculations
    LinearAlgebra::MiddleSizeMatrix integrandStabilityRightHandSideOnFace(
        Base::PhysicalFace<DIM> &face,
        const LinearAlgebra::MiddleSizeVector stateCoefficientsLeft,
        const LinearAlgebra::MiddleSizeVector stateCoefficientsRight,
        const Base::Side &side);

    /// \brief Computes the rhs, for the system of equations solving the
    /// stability parameters, by integrating the rhs  stability parameter
    /// integrand
    LinearAlgebra::MiddleSizeMatrix computeRhsStabilityParameters(
        const Base::Face *ptrFace,
        const LinearAlgebra::MiddleSizeVector stateCoefficientsLeft,
        const LinearAlgebra::MiddleSizeVector stateCoefficientsRight,
        const Base::Side &side);

    /// \brief Computes the stability parameters used in the auxilliary
    /// integrand for an internal face for given stabilityParameterflux
    LinearAlgebra::MiddleSizeMatrix computeStabilityParameters(
        Base::PhysicalFace<DIM> &face,
        const LinearAlgebra::MiddleSizeVector stateCoefficientsLeft,
        const LinearAlgebra::MiddleSizeVector stateCoefficientsRight);

    /// \brief Sets the mass matrix used in the computation of the stability
    /// parameters
    void setStabilityMassMatrix(
        LinearAlgebra::MiddleSizeMatrix &stabilityMassMatrix);

    /// ************************************************
    /// ***    Jacobian Exact Derivative Functions   ***
    /// ************************************************
    /*  //todo: Complete this section at some point. Many functions are not yet
       finished. The question is, do I require the exact Jacobian or is a
       difference approach good enough?

            /// \brief computes the derivative of the viscous element integral
       flux with respect to the solutionCoefficients
            LinearAlgebra::MiddleSizeMatrix
       computeViscousFluxDiv(Base::PhysicalElement<DIM> &element, const
       std::size_t iV, const std::size_t iB, const
       LinearAlgebra::MiddleSizeVector &state, const
       LinearAlgebra::MiddleSizeVector &partialState, const
       LinearAlgebra::MiddleSizeMatrix &stateJacobian);

            /// \brief computes the stabilityfluxFunction derivative, used for
       calculating the auxilliary jacobian face integrand
            LinearAlgebra::MiddleSizeMatrix
       computeStabilityFluxFunctionDiv(const
       std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensor, const
       std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorDiv, const
       LinearAlgebra::MiddleSizeVector &stateLeft, const
       LinearAlgebra::MiddleSizeVector &stateLeftDiv, const
       LinearAlgebra::MiddleSizeVector &stateRight, const
       LinearAlgebra::MiddleSizeVector &stateRightDiv, const
       LinearAlgebra::SmallVector<DIM> &unitNormalLeft);

            /// This function computes the Auxilliary flux derivative for the
       auxilliary jacobian face integral LinearAlgebra::MiddleSizeMatrix
       computeAuxilliaryFluxDiv(const LinearAlgebra::MiddleSizeVector
       &stateLeft, const LinearAlgebra::MiddleSizeVector &stateLeftDiv, const
       LinearAlgebra::MiddleSizeVector &stateRight, const
       LinearAlgebra::MiddleSizeVector &stateRightDiv, const
       std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorLeft, const
       std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorLeftDiv, const
       std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorRight,	const
       std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorRightDiv, const
       LinearAlgebra::MiddleSizeMatrix &stateJacobianLeft, const
       LinearAlgebra::MiddleSizeMatrix &stateJacobianLeftDiv, const
       LinearAlgebra::MiddleSizeMatrix &stateJacobianRight, const
       LinearAlgebra::MiddleSizeMatrix &stateJacobianRightDiv, const
       LinearAlgebra::SmallVector<DIM> &unitNormalLeft) ;
    */

    /// ************************************************
    /// ***    Jacobian Element Matrix Functions     ***
    /// ************************************************

    /// \brief Compute the ViscousElementIntegrand contribution to the Jacobian
    /// matrix
    LinearAlgebra::MiddleSizeMatrix integrandJacobianViscousElement(
        Base::PhysicalElement<DIM> &element,
        const StateCoefficientsStruct<DIM, NUMBER_OF_VARIABLES>
            &elementStateStruct);

    /// *********************************************
    /// ***    Jacobian Face Matrix Functions     ***
    /// *********************************************

    /// \brief This function computes the integrand of the Jacobian integral for
    /// the viscous part of the equations. The Jacobian is approximated by a
    /// difference approach
    LinearAlgebra::MiddleSizeMatrix integrandJacobianViscousFace(
        Base::PhysicalFace<DIM> &face,
        const StateCoefficientsStruct<DIM, NUMBER_OF_VARIABLES>
            &faceStateStructLeft,
        const StateCoefficientsStruct<DIM, NUMBER_OF_VARIABLES>
            &faceStateStructRight,
        const Base::Side &elementSide, const Base::Side &derivativeSide);

    /// \brief This function computes the integrand of the Jacobian integral for
    /// the Auxilliary part of the equations. The Jacobian is approximated by a
    /// difference approach
    LinearAlgebra::MiddleSizeMatrix integrandJacobianAuxilliaryFace(
        Base::PhysicalFace<DIM> &face,
        const StateCoefficientsStruct<DIM, NUMBER_OF_VARIABLES>
            &faceStateStructLeft,
        const StateCoefficientsStruct<DIM, NUMBER_OF_VARIABLES>
            &faceStateStructRight,
        const Base::Side &elementSide, const Base::Side &derivativeSide);

   private:
    /// \var Integrator voor the stability parameters
    Integration::FaceIntegral<DIM> stabilityFaceIntegrator_;

    LinearAlgebra::MiddleSizeMatrix stabilityMassMatrix_;  // Note: this breaks
                                                           // down if p is not
                                                           // the same in all
                                                           // elements.

    UnsteadyNavierStokesAPI<DIM, NUMBER_OF_VARIABLES> *instance_;
};

#include "ViscousTerms_Impl.h"

#endif /* VISCOUSTERMS_H_ */
