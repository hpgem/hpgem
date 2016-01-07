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

#ifndef UNSTEADYNAVIERSTOKESAPI_H_
#define UNSTEADYNAVIERSTOKESAPI_H_

//#include "Base/HpgemAPINonLinearSteadyState.h"
#include "Base/HpgemAPISimplified.h"

template<std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
class UnsteadyNavierStokesAPI;

#include "InviscidTerms.h"
#include "ViscousTerms.h"


//class Thomas : public Base::HpgemAPINonLinearSteadyState<DIM>
template<std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
class UnsteadyNavierStokesAPI : public Base::HpgemAPISimplified<DIM>
{
public:

	UnsteadyNavierStokesAPI(const std::size_t numOfVariables, const double endTime, const std::size_t polynomialOrder, const TimeIntegration::ButcherTableau * const ptrButcherTableau, const bool computeBothFaces);

	virtual ~UnsteadyNavierStokesAPI()
	{
	}

	/// *****************************************
	/// ***   Element integration functions   ***
	/// *****************************************

	/// \brief Compute source function at an element
	virtual LinearAlgebra::MiddleSizeVector integrandSourceAtElement(Base::PhysicalElement<DIM> &element, const StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> &elementStateStruct, const double &time)
	{
		std::size_t numOfBasisFunctions = element.getNumOfBasisFunctions();

		LinearAlgebra::MiddleSizeVector integrandSource(NUMBER_OF_VARIABLES * numOfBasisFunctions);

		//Convert pRef to pPhys
		Geometry::PointPhysical<DIM> pPhys = element.getPointPhysical();

		std::cout << "There is no source function implemented." << std::endl;
		return integrandSource;
	}


	/// \brief Compute integrand of righthandside on an element
	LinearAlgebra::MiddleSizeVector integrandRightHandSideOnElement(Base::PhysicalElement<DIM>& element, const double &time, const LinearAlgebra::MiddleSizeVector &stateCoefficients);

	/// \brief Compute the right hand side on an element
	LinearAlgebra::MiddleSizeVector computeRightHandSideAtElement(Base::Element *ptrElement, LinearAlgebra::MiddleSizeVector &solutionCoefficients, const double time) override final;

	/// **************************************************
	/// ***    external face integration functions     ***
	/// **************************************************
	/*

	    /// \brief Compute the integrand for the right hand side for the face corresponding to an external face.
	    	LinearAlgebra::MiddleSizeVector integrandRightHandSideOnFace(Base::PhysicalFace<DIM> &face, const double &time, const LinearAlgebra::MiddleSizeVector &stateCoefficients);

	    /// \brief Compute the right-hand side corresponding to an external face
	    LinearAlgebra::MiddleSizeVector computeRightHandSideAtFace(Base::Face *ptrFace, LinearAlgebra::MiddleSizeVector &solutionCoefficients, const double time) override final;
	*/

	/// **************************************************
	/// ***    internal face integration functions     ***
	/// **************************************************

	/// \brief Compute integrands on both sides of an internal face
	std::pair<LinearAlgebra::MiddleSizeVector,LinearAlgebra::MiddleSizeVector> integrandsRightHandSideOnFace(Base::PhysicalFace<DIM>& face, const double &time, const LinearAlgebra::MiddleSizeVector &stateCoefficientsLeft, const LinearAlgebra::MiddleSizeVector &stateCoefficientsRight);

	/// \brief Compute the right-hand side corresponding to an internal face. Computing both left and right face integrals
	std::pair<LinearAlgebra::MiddleSizeVector,LinearAlgebra::MiddleSizeVector> computeBothRightHandSidesAtFace(Base::Face *ptrFace, LinearAlgebra::MiddleSizeVector &solutionCoefficientsLeft, LinearAlgebra::MiddleSizeVector &solutionCoefficientsRight, const double time) override final;

	/// ************************************************
	/// ***     Jacobian Element Matrix Functions    ***
	/// ************************************************

	/// \brief This function computes the integrand of local element integral Jacobian contributions
	LinearAlgebra::MiddleSizeMatrix integrandJacobianAtElement(Base::PhysicalElement<DIM>& element, const LinearAlgebra::MiddleSizeVector &solutionCoefficients, const double time);

	/// \brief This function computes the local Jacobian contributions from the element solutionCoefficients
	LinearAlgebra::MiddleSizeMatrix computeJacobianAtElement(Base::Element *ptrElement, const LinearAlgebra::MiddleSizeVector &solutionCoefficients, const double time);

	/// ************************************************
	/// ***      Jacobian  Face Matrix Functions     ***
	/// ************************************************

	/// \brief This function computes the integrand of local and non local face integral Jacobian contributions
	LinearAlgebra::MiddleSizeMatrix integrandJacobianAtFace(Base::PhysicalFace<DIM>& face, const LinearAlgebra::MiddleSizeVector &solutionCoefficientsLeft, const LinearAlgebra::MiddleSizeVector &solutionCoefficientsRight, const  Base::Side elementSide, const Base::Side derivativeSide, const double time);

	/// \brief This functions computes the Jacobian face matrix of element elementSide with respect to variables of derivativeSide
	LinearAlgebra::MiddleSizeMatrix computeJacobianAtFace(Base::Face *ptrFace, LinearAlgebra::MiddleSizeVector &solutionCoefficientsLeft, LinearAlgebra::MiddleSizeVector &solutionCoefficientsRight, const double time, Base::Side elementSide, Base::Side derivativeSide);

	/// *****************************************
	/// ***    		Constitutive Relations    ***
	/// *****************************************

	/// \brief This equation obtains the functions to compute constitutive relations
	virtual StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> constutiveRelationsStruct() = 0;

	/// \brief This equation couples the constitutive relations to the Navier-Stokes API for an Element
	virtual StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> computeElementStateStruct(Base::PhysicalElement<DIM> &element, const LinearAlgebra::MiddleSizeVector &stateCoefficients, const double time) = 0;

	/// \brief This equation couples the constitutive relations to the Navier-Stokes API for a Face
	virtual StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> computeFaceStateStruct(Base::PhysicalFace<DIM> &face, const LinearAlgebra::MiddleSizeVector &stateCoefficients, const Base::Side side, const double time) = 0;

	/// *****************************************
	/// ***    		Various Functions         ***
	/// *****************************************

	/// \brief Create a domain
	Base::RectangularMeshDescriptor<DIM> createMeshDescription(const std::size_t numOfElementPerDirection) override final;

	void tasksBeforeSolving() override final;

	//todo: this can be removed
	/// \brief Computes the Error at the end of the simulation, compared to the exact solution
	LinearAlgebra::MiddleSizeVector Error(const double time);

	/// \brief Shows the progress in the terminal as output
	void showProgress(const double time, const std::size_t timeStepID) override final;

private:
		/// \var Inviscid class, treating the inviscid part of the NS equations
		InviscidTerms<DIM,NUMBER_OF_VARIABLES> inviscidTerms_;

		/// \var Viscous class, treating the viscosity part of the NS equations
		ViscousTerms<DIM,NUMBER_OF_VARIABLES> viscousTerms_;
};

#include "UnsteadyNavierStokesAPI_Impl.h"

#endif /* UNSTEADYNAVIERSTOKESAPI */
