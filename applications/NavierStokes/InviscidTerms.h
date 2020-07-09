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

#ifndef HPGEM_APP_INVISCIDTERMS_H
#define HPGEM_APP_INVISCIDTERMS_H

#include "NavierStokesConstants.h"
#include "UnsteadyNavierStokesAPI.h"
#include "StateCoefficientsStruct.h"

template <std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
class InviscidTerms {
   public:
    InviscidTerms(UnsteadyNavierStokesAPI<DIM, NUMBER_OF_VARIABLES> *instance)
        : instance_(instance) {}

    virtual ~InviscidTerms() {}

    /// *****************************************
    /// ***   Element integration functions   ***
    /// *****************************************

    /// \brief Compute integrand of righthandside on an element
    LinearAlgebra::MiddleSizeVector integrandAtElement(
        Base::PhysicalElement<DIM> &element,
        const StateCoefficientsStruct<DIM, NUMBER_OF_VARIABLES>
            &elementStateStruct,
        const double time);

    /// **************************************************
    /// ***    General face integration functions      ***
    /// **************************************************

    /// \brief Compute the local Lax-Friedrichs flux
    LinearAlgebra::MiddleSizeVector computeLLFFluxFunction(
        const LinearAlgebra::MiddleSizeVector stateLeft,
        const LinearAlgebra::MiddleSizeVector stateRight,
        const double pressureLeft, const double pressureRight,
        const LinearAlgebra::SmallVector<DIM> &unitNormalLeft);

    /// \brief Compute the local Lax-Friedrichs flux
    LinearAlgebra::MiddleSizeVector computeLLFFluxFunction(
        const StateCoefficientsStruct<DIM, NUMBER_OF_VARIABLES>
            faceStateStructLeft,
        const StateCoefficientsStruct<DIM, NUMBER_OF_VARIABLES>
            faceStateStructRight,
        const LinearAlgebra::SmallVector<DIM> &unitNormalLeft);

    /// \brief Compute the Roe Rieman flux
    LinearAlgebra::MiddleSizeVector computeRoeFluxFunction(
        const StateCoefficientsStruct<DIM, NUMBER_OF_VARIABLES>
            faceStateStructLeft,
        const StateCoefficientsStruct<DIM, NUMBER_OF_VARIABLES>
            faceStateStructRight,
        const LinearAlgebra::SmallVector<DIM> &unitNormalLeft);

    /// \brief Compute the HLLC Flux function, based on Klaij et al. 2006
    LinearAlgebra::MiddleSizeVector computeHLLCFluxFunction(
        const StateCoefficientsStruct<DIM, NUMBER_OF_VARIABLES>
            &faceStateStructLeft,
        const StateCoefficientsStruct<DIM, NUMBER_OF_VARIABLES>
            &faceStateStructRight,
        const LinearAlgebra::SmallVector<DIM> &unitNormalLeft);

    /// **************************************************
    /// ***    External face integration functions     ***
    /// **************************************************

    /// \brief Compute the integrand for the right hand side for the reference
    /// face corresponding to a external face.
    LinearAlgebra::MiddleSizeVector integrandAtBoundaryFace(
        Base::PhysicalFace<DIM> &face, const BoundaryType boundaryType,
        const StateCoefficientsStruct<DIM, NUMBER_OF_VARIABLES>
            &faceStateStructBoundary,
        const StateCoefficientsStruct<DIM, NUMBER_OF_VARIABLES>
            &faceStateStructLeft,
        const double &time);

    /// **************************************************
    /// ***    Internal face integration functions     ***
    /// **************************************************

    /// \brief Compute both the face integral for the left element as the right
    /// element at the same time. (reducing flux calculations)
    std::pair<LinearAlgebra::MiddleSizeVector, LinearAlgebra::MiddleSizeVector>
        integrandsAtFace(Base::PhysicalFace<DIM> &face,
                         const StateCoefficientsStruct<DIM, NUMBER_OF_VARIABLES>
                             &faceStateStructLeft,
                         const StateCoefficientsStruct<DIM, NUMBER_OF_VARIABLES>
                             &faceStateStructRight,
                         const double &time);

    /// ************************************************
    /// ***    Jacobian Exact Derivative Functions   ***
    /// ************************************************

    /*
            /// \brief Compute the Jacobian with respect to the solution
       coefficients of the element std::vector<LinearAlgebra::MiddleSizeMatrix>
       computeFluxJacobian2D(const LinearAlgebra::MiddleSizeVector &state);

            /// \brief Computes the derivative of the InviscidFace flux
            LinearAlgebra::MiddleSizeVector
       computeInviscidFaceFluxDiv(Base::PhysicalFace<DIM> &face, const
       LinearAlgebra::MiddleSizeVector &stateLeft, const
       LinearAlgebra::MiddleSizeVector &stateRight, const
       LinearAlgebra::SmallVector<DIM> &unitNormalLeft, const Base::Side
       derivativeSide, const std::size_t iV2, const std::size_t iB2);

            //LinearAlgebra::MiddleSizeVector
       computeInviscidFaceLLFFluxDiv(Base::PhysicalFace<DIM> &face, const
       LinearAlgebra::MiddleSizeVector &stateLeft, const
       LinearAlgebra::MiddleSizeVector &stateRight, const
       LinearAlgebra::SmallVector<DIM> &unitNormalLeft, const Base::Side
       derivativeSide, const std::size_t iV2, const std::size_t iB2);

            //LinearAlgebra::MiddleSizeMatrix
       computeFluxTimesNormalJacobian(const LinearAlgebra::MiddleSizeVector
       &state, const LinearAlgebra::SmallVector<DIM> &unitNormal);


    */

    /// ************************************************
    /// ***    Jacobian Element Matrix Functions     ***
    /// ************************************************

    /// \brief Compute the inviscid element integral integrand contribution to
    /// the Jacobian matrix
    LinearAlgebra::MiddleSizeMatrix integrandJacobianInviscidElement(
        Base::PhysicalElement<DIM> &element,
        const StateCoefficientsStruct<DIM, NUMBER_OF_VARIABLES>
            &elementStateStruct);

    /// *********************************************
    /// ***    Jacobian Face Matrix Functions     ***
    /// *********************************************

    /// \brief Computes the inviscid face integral integrand contribution to the
    /// Jacobian matrix
    LinearAlgebra::MiddleSizeMatrix integrandJacobianInviscidFace(
        Base::PhysicalFace<DIM> &face,
        const StateCoefficientsStruct<DIM, NUMBER_OF_VARIABLES>
            &faceStateStructLeft,
        const StateCoefficientsStruct<DIM, NUMBER_OF_VARIABLES>
            &faceStateStructRight,
        const Base::Side elementSide, const Base::Side derivativeSide);

   private:
    UnsteadyNavierStokesAPI<DIM, NUMBER_OF_VARIABLES> *instance_;
};

#include "InviscidTerms_Impl.h"

#endif // HPGEM_APP_INVISCIDTERMS_H
