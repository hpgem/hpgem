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

#ifndef COMPRESSIBLENAVIERSTOKES_H_
#define COMPRESSIBLENAVIERSTOKES_H_

#include "Base/HpgemAPISimplified.h"

class CompressibleNavierStokes;

#include "Inviscid.h"
#include "Viscous.h"

class CompressibleNavierStokes : public Base::HpgemAPISimplified<DIM>
{
public:
	CompressibleNavierStokes
	(
		const std::size_t numOfVariables,
		const double endTime,
		const std::size_t polynomialOrder,
		const TimeIntegration::ButcherTableau * const ptrButcherTableau,
		const bool computeBothFaces
	);

    /// \brief Create a domain
    Base::RectangularMeshDescriptor<DIM> createMeshDescription(const std::size_t numOfElementPerDirection) override;

    void setStabilityMassMatrix();

    /// \brief Computes pressure for a given state
    double computePressure(const LinearAlgebra::MiddleSizeVector &state);

    /// *************************************************
    /// ***   Element integration support functions   ***
    /// *************************************************

    /// \brief Compute state at an element
    LinearAlgebra::MiddleSizeVector computeStateOnElement(Base::PhysicalElement<DIM> &element, const LinearAlgebra::MiddleSizeVector &stateCoefficients);

    /// \brief Compute state derivatives at an element.
    LinearAlgebra::MiddleSizeMatrix computeStateJacobianAtElement(Base::PhysicalElement<DIM> &element, const LinearAlgebra::MiddleSizeVector &stateCoefficients);

/*
	//todo: Remove this function, it is outdated
    /// Compute the Jacobian of the velocities
    LinearAlgebra::MiddleSizeMatrix computePartialStateJacobian(const LinearAlgebra::MiddleSizeMatrix qSolutionGradient, const LinearAlgebra::MiddleSizeVector qSolution);
*/

    /// \brief Computes the partial states: all states except the density are divided by the density
    LinearAlgebra::MiddleSizeVector computePartialState(const LinearAlgebra::MiddleSizeVector &state);

    /// *****************************************
    /// ***   Element integration functions   ***
    /// *****************************************

    /// \brief Compute source function at an element
    LinearAlgebra::MiddleSizeVector integrandSourceAtElement(Base::PhysicalElement<DIM> &element, const LinearAlgebra::MiddleSizeVector &state, const double &pressureTerm, const double &time);

    /// \brief Compute integrand of righthandside on an element
    LinearAlgebra::MiddleSizeVector integrandRightHandSideOnElement(Base::PhysicalElement<DIM>& element, const double &time, const LinearAlgebra::MiddleSizeVector &stateCoefficients);

    /// \brief Compute the right hand side on an element
    LinearAlgebra::MiddleSizeVector computeRightHandSideAtElement(Base::Element *ptrElement, LinearAlgebra::MiddleSizeVector &solutionCoefficients, const double time) override final;

    /// *************************************************
    /// ***    face integration support functions     ***
    /// *************************************************

    /// \brief Compute state at a face
    LinearAlgebra::MiddleSizeVector computeStateOnFace(Base::PhysicalFace<DIM> &face, const Base::Side &iSide, const LinearAlgebra::MiddleSizeVector &stateCoefficients) const;

    /// \brief Compute state at a face
    //LinearAlgebra::MiddleSizeVector computeStateOnFace(Base::PhysicalFace<DIM> &face, const Base::Side &iSide, LinearAlgebra::MiddleSizeVector &stateCoefficients);

    /// \brief Compute state Jacobian on face
    LinearAlgebra::MiddleSizeMatrix computeStateJacobianAtFace(Base::PhysicalFace<DIM> &face, const Base::Side &iSide, const LinearAlgebra::MiddleSizeVector &stateCoefficients);

    /// **************************************************
    /// ***    external face integration functions     ***
    /// **************************************************

/*    /// \brief Compute the integrand for the right hand side for the face corresponding to an external face.
    LinearAlgebra::MiddleSizeVector integrandRightHandSideOnFace(Base::PhysicalFace<DIM> &face, const double &time, const LinearAlgebra::MiddleSizeVector &stateCoefficients);

    /// \brief Compute the right-hand side corresponding to an external face
    LinearAlgebra::MiddleSizeVector computeRightHandSideAtFace(Base::Face *ptrFace, LinearAlgebra::MiddleSizeVector &solutionCoefficients, const double time) override final;*/

    /// **************************************************
    /// ***    internal face integration functions     ***
    /// **************************************************


/*
    /// \brief Compute the integrand for the right hand side for the face corresponding to an internal face.
    LinearAlgebra::MiddleSizeVector integrandRightHandSideOnFace(Base::PhysicalFace<DIM>& face, const double &time, const Base::Side &iSide, const LinearAlgebra::MiddleSizeVector &stateCoefficientsLeft, const LinearAlgebra::MiddleSizeVector &stateCoefficientsRight);
*/

    /// \brief Compute the right-hand side corresponding to an internal face
    /// Note: this line is only required to avoid a designed safety error in hpGEM (hack)
    LinearAlgebra::MiddleSizeVector computeRightHandSideAtFace(Base::Face *ptrFace, const Base::Side side, LinearAlgebra::MiddleSizeVector &solutionCoefficientsLeft, LinearAlgebra::MiddleSizeVector &solutionCoefficientsRight, const double time) override final
    		{
    		LinearAlgebra::MiddleSizeVector empty;
    		return empty;
    		}


    /// \brief Compute integrands on both sides of an internal face
    std::pair<LinearAlgebra::MiddleSizeVector,LinearAlgebra::MiddleSizeVector> integrandsRightHandSideOnFace(Base::PhysicalFace<DIM>& face, const double &time, const LinearAlgebra::MiddleSizeVector &stateCoefficientsLeft, const LinearAlgebra::MiddleSizeVector &stateCoefficientsRight);

    /// \brief Compute the right-hand side corresponding to an internal face. Computing both left and right face integrals
    std::pair<LinearAlgebra::MiddleSizeVector,LinearAlgebra::MiddleSizeVector> computeBothRightHandSidesAtFace(Base::Face *ptrFace, LinearAlgebra::MiddleSizeVector &solutionCoefficientsLeft, LinearAlgebra::MiddleSizeVector &solutionCoefficientsRight, const double time) override final;

    /// *****************************************
    /// ***    		Various Functions         ***
    /// *****************************************

    /// \brief Computes the exact solution of the given problem (if it is avaiable)
    LinearAlgebra::MiddleSizeVector getExactSolution(const PointPhysicalT &pPhys, const double &time, const std::size_t orderTimeDerivative) override final;

    /// \brief Compute the initial solution at a given point in space and time.
	LinearAlgebra::MiddleSizeVector getInitialSolution(const PointPhysicalT &pPhys, const double &startTime, const std::size_t orderTimeDerivative = 0) override final;

	/// \brief Computes the Error at the end of the simulation, compared to the exact solution
	LinearAlgebra::MiddleSizeVector Error(const double time);

	/// \brief Shows the progress in the terminal as output
	void showProgress(const double time, const std::size_t timeStepID) override final;

/*	deprecated
	/// \brief Ensure that referenceIntegrate works correctly with the newly templated API. This is a hack.
	void beforeTimeIntegration()
	{
	    faceIntegrator_.setTransformation(std::shared_ptr<Base::CoordinateTransformation<DIM> >(new Base::DoNotScaleIntegrands<DIM>(new Base::H1ConformingTransformation<DIM>())));
	}*/

private:
    /// \var Dimension of the domain
    const std::size_t DIM_;

    /// \var Specific heat ratio
    const double gamma_ = 1.4;

    /// \var specific gas constant
    const double Rs_ = (1 - 1/gamma_)*1000; // cp-cv;

    /// \var specific heat at constant pressure
    //todo: remove this from viscous.cpp
    const double cp_ = 1000;

	/// \var Number of variables
	const std::size_t numOfVariables_;

	/// \var Inviscid class, treating the inviscid part of the NS equations
	Inviscid inviscidTerms_;

	/// \var Viscous class, treating the viscosity part of the NS equations
	Viscous viscousTerms_;

	/// Plate characteristics
	const double uPlateTop_ = 0.0; 			//Velocity of the top plate
	const double uPlateBottom_ = 0.0;		//Velocity of the bottom plate
	const double tPlateTop_ = 288;			//temperature of the top plate
	const double tPlateBottom_ = 288;		//temperature of the bottom plate

	/// Simulation parameters
	const double Tc_ = 3.0/5.0;

	/// non-Dimensionless parameters
	const double rho0_ = 1.225;
	const double U0_ = 1.0;
	const double E0_ = 1.225*1000/1.4*288;
	const double H_ = 1.0;

	const double nonDIM1_ = 1.0;//E0_/(U0_*U0_);// //E0_/(uWall_ *uWall_); //Used for scaling the pressure correctly

	friend class Inviscid;
	friend class Viscous;




};

#endif /* COMPRESSIBLENAVIERSTOKES_H_ */
