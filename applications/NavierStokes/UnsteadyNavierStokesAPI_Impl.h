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

//#include "Base/HpgemAPINonLinearSteadyState.h"
/*
#include "../StateCoefficientsStruct.h"
#include "Base/HpgemAPISimplified.h"
#include "ThomasConstants.h"
#include "TInviscid.h"
#include "TViscous.h"
*/

#include "Base/HpgemAPISimplified.h"
#include "UnsteadyNavierStokesAPI.h"

//class Thomas : public Base::HpgemAPINonLinearSteadyState<DIM>


/*	/// \brief Steady State constructor
	Thomas
	(
		const double endTime,
		const std::size_t polynomialOrder,
		const bool computeBothFaces
	);*/

	/// \brief Unsteady Constructor
template<std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
UnsteadyNavierStokesAPI<DIM,NUMBER_OF_VARIABLES>::UnsteadyNavierStokesAPI
(
 const std::size_t numOfVariables,
 const double endTime,
 const std::size_t polynomialOrder,
 const TimeIntegration::ButcherTableau * const ptrButcherTableau,
 const bool computeBothFaces
) :
Base::HpgemAPISimplified<DIM>(numOfVariables, polynomialOrder, ptrButcherTableau, 1, computeBothFaces),
inviscidTerms_(this),
viscousTerms_(this),
time_(0)
{
}

/// *****************************************
/// ***   Element integration functions   ***
/// *****************************************

/// \brief Compute integrand of righthandside on an element
template<std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
LinearAlgebra::MiddleSizeVector UnsteadyNavierStokesAPI<DIM,NUMBER_OF_VARIABLES>::integrandRightHandSideOnElement
(
 Base::PhysicalElement<DIM>& element,
 const double &time,
 const LinearAlgebra::MiddleSizeVector &stateCoefficients)
{
	//Compute all functions that depend on the state coefficients, required to calculate the integrands
	StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> elementStateStruct = computeElementStateStruct(element, stateCoefficients, time);

	//Compute inviscid terms
	LinearAlgebra::MiddleSizeVector integrandInviscid = inviscidTerms_.integrandAtElement(element, elementStateStruct, time);

	//Compute viscous terms
	//LinearAlgebra::MiddleSizeVector integrandViscous = viscousTerms_.integrandAtElement(element, elementStateStruct, time);

	//Compute source terms
	LinearAlgebra::MiddleSizeVector integrandSource = integrandSourceAtElement(element, elementStateStruct, time);

/*	std::cout << "===Element integral===" << std::endl;
	std::cout << "Inviscid: " << integrandInviscid << std::endl;
	std::cout << "Source: " << integrandSource << std::endl;
	std::cout << "Viscous: " << integrandViscous << std::endl;*/

	return integrandInviscid + integrandSource;// + integrandViscous;
}

/// \brief Compute the right hand side on an element
template<std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
LinearAlgebra::MiddleSizeVector UnsteadyNavierStokesAPI<DIM,NUMBER_OF_VARIABLES>::computeRightHandSideAtElement
(
 Base::Element *ptrElement,
 LinearAlgebra::MiddleSizeVector &solutionCoefficients,
 const double time)
{
	std::function<LinearAlgebra::MiddleSizeVector(Base::PhysicalElement<DIM> &)> integrandFunction = [=](Base::PhysicalElement<DIM>& element) -> LinearAlgebra::MiddleSizeVector
			{   return this->integrandRightHandSideOnElement(element, time, solutionCoefficients);};

	return this->elementIntegrator_.integrate(ptrElement, integrandFunction, ptrElement->getGaussQuadratureRule());
}

/// **************************************************
/// ***    external face integration functions     ***
/// **************************************************

/// \brief Compute the integrand for the right hand side for the face corresponding to an external face.
template<std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
LinearAlgebra::MiddleSizeVector UnsteadyNavierStokesAPI<DIM,NUMBER_OF_VARIABLES>::integrandRightHandSideOnFace(Base::PhysicalFace<DIM> &face, const double &time, const LinearAlgebra::MiddleSizeVector &stateCoefficients)
{
	//Compute the datastructures

	StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> faceStateStructLeft = computeFaceStateStruct(face, stateCoefficients, Base::Side::LEFT, time);

	//todo: in future: return the boundary condition + if it is dirichlet or neumann
	//NOTE: Jacobian information of a dirichlet boundary is not available (!!!!)
	LinearAlgebra::MiddleSizeVector stateBoundary = computeBoundaryState(face, faceStateStructLeft, time);
	//std::cout << "boundaryState: " << stateBoundary << std::endl;
	StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> faceStateStructBoundary = computeBoundaryFaceStateStruct(stateBoundary, time);

	//Compute inviscid terms
	LinearAlgebra::MiddleSizeVector integrandInviscid = inviscidTerms_.integrandAtBoundaryFace(face, faceStateStructBoundary, faceStateStructLeft, time);

	//Compute viscous terms
	//LinearAlgebra::MiddleSizeVector integrandViscous = viscousTerms_.integrandViscousAtBoundaryFace(face, faceStateStructBoundary, faceStateStructLeft, time);

	//Compute support variable terms
	//LinearAlgebra::MiddleSizeVector integrandAuxilliary = viscousTerms_.integrandAuxilliaryAtBoundaryFace(face, faceStateStructBoundary, faceStateStructLeft, time);

/*		std::cout << "===External Face Integral===" << std::endl;
        std::cout << "inv: " << integrandInviscid << std::endl;
        std::cout << "vis: " << integrandViscous << std::endl;
        std::cout << "aux: " << integrandAuxilliary << std::endl;*/
	
	return  integrandInviscid;// + integrandViscous + integrandAuxilliary;
}

/// \brief Compute the right-hand side corresponding to an external face
template<std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
LinearAlgebra::MiddleSizeVector UnsteadyNavierStokesAPI<DIM,NUMBER_OF_VARIABLES>::computeRightHandSideAtFace(Base::Face *ptrFace, LinearAlgebra::MiddleSizeVector &solutionCoefficients, const double time)
{
	    // Define the integrand function for the right hand side for the reference face.
	    std::function<LinearAlgebra::MiddleSizeVector(Base::PhysicalFace<DIM>&)> integrandFunction = [=](Base::PhysicalFace<DIM>& face) -> LinearAlgebra::MiddleSizeVector
	    {   return this->integrandRightHandSideOnFace(face, time, solutionCoefficients);};

	    return this->faceIntegrator_.integrate(ptrFace, integrandFunction);
}

/// **************************************************
/// ***    internal face integration functions     ***
/// **************************************************

/// \brief Compute integrands on both sides of an internal face
template<std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
std::pair<LinearAlgebra::MiddleSizeVector,LinearAlgebra::MiddleSizeVector> UnsteadyNavierStokesAPI<DIM,NUMBER_OF_VARIABLES>::integrandsRightHandSideOnFace
(
 Base::PhysicalFace<DIM>& face,
 const double &time,
 const LinearAlgebra::MiddleSizeVector &stateCoefficientsLeft,
 const LinearAlgebra::MiddleSizeVector &stateCoefficientsRight)
{
	//Compute all functions that depend on the state coefficients, required to calculate the integrands
	StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> faceStateStructLeft = computeFaceStateStruct(face, stateCoefficientsLeft, Base::Side::LEFT, time);
	StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> faceStateStructRight = computeFaceStateStruct(face, stateCoefficientsRight, Base::Side::RIGHT, time);

/*	if (face.getID() == 10000000000)
	{
		std::cout << "============" << std::endl;
		std::cout << "State Left: " << faceStateStructLeft.getState() << std::endl;
		std::cout << "partialState Left: " << faceStateStructLeft.getPartialState() << std::endl;
		//std::cout << "JacobianState: " << faceStateStructLeft.getStateJacobian() << std::endl;
		std::cout << "hyperbolic matrix Left: " << faceStateStructLeft.getHyperbolicMatrix() << std::endl;
		std::cout << "pressure Left: " << faceStateStructLeft.getPressure() << std::endl;
		std::cout << "Speed of Sound Left: " << faceStateStructLeft.getSpeedOfSound() << std::endl;
		std::cout << "+++++++++++++++++++" << std::endl;
		std::cout << "State Right: " << faceStateStructRight.getState() << std::endl;
		std::cout << "partialState Right: " << faceStateStructRight.getPartialState() << std::endl;
		//std::cout << "JacobianState: " << faceStateStructLeft.getStateJacobian() << std::endl;
		std::cout << "hyperbolic matrix Right: " << faceStateStructRight.getHyperbolicMatrix() << std::endl;
		std::cout << "pressure Right: " << faceStateStructRight.getPressure() << std::endl;
		std::cout << "Speed of Sound Right: " << faceStateStructRight.getSpeedOfSound() << std::endl;
	}*/

	//Compute inviscid terms
	std::pair<LinearAlgebra::MiddleSizeVector,LinearAlgebra::MiddleSizeVector> integrandsInviscid = inviscidTerms_.integrandsAtFace(face, faceStateStructLeft, faceStateStructRight, time);

	//Compute viscous terms
	//std::pair<LinearAlgebra::MiddleSizeVector,LinearAlgebra::MiddleSizeVector> integrandsViscous = viscousTerms_.integrandsViscousAtFace(face, faceStateStructLeft, faceStateStructRight, time);

	//Compute support variable terms
	//std::pair<LinearAlgebra::MiddleSizeVector,LinearAlgebra::MiddleSizeVector> integrandsAuxilliary = viscousTerms_.integrandsAuxilliaryAtFace(face, faceStateStructLeft, faceStateStructRight, time);

/*	std::cout << "===Internal Face Integral===" << std::endl;
	std::cout << "integrandsInviscid Left: " << integrandsInviscid.first << std::endl;
	std::cout << "integrandsInvisecid Right: " << integrandsInviscid.second << std::endl;
	std::cout << "integrandVsicous Left: " << integrandsViscous.first << std::endl;
	std::cout << "integrandVsicous Right: " << integrandsViscous.second << std::endl;
	std::cout << "integrandAuxilliary Left: " << integrandsAuxilliary.first << std::endl;
	std::cout << "integrandAuxilliary Right: " << integrandsAuxilliary.second << std::endl;*/

	//combine all integrands into one pair
	std::pair<LinearAlgebra::MiddleSizeVector,LinearAlgebra::MiddleSizeVector> integrands;
	integrands.first = integrandsInviscid.first;// + integrandsViscous.first + integrandsAuxilliary.first;
	integrands.second = integrandsInviscid.second;// + integrandsViscous.second + integrandsAuxilliary.second;

	return  integrands;
}

/// \brief Compute the right-hand side corresponding to an internal face. Computing both left and right face integrals
template<std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
std::pair<LinearAlgebra::MiddleSizeVector,LinearAlgebra::MiddleSizeVector> UnsteadyNavierStokesAPI<DIM,NUMBER_OF_VARIABLES>::computeBothRightHandSidesAtFace
(
 Base::Face *ptrFace,
 LinearAlgebra::MiddleSizeVector &solutionCoefficientsLeft,
 LinearAlgebra::MiddleSizeVector &solutionCoefficientsRight,
 const double time
)
{
	// Define the integrand function for the right hand side for the face.
	std::function<std::pair<LinearAlgebra::MiddleSizeVector,LinearAlgebra::MiddleSizeVector>(Base::PhysicalFace<DIM>&)> integrandFunction = [=](Base::PhysicalFace<DIM>& face) -> std::pair<LinearAlgebra::MiddleSizeVector,LinearAlgebra::MiddleSizeVector>
	{   return this->integrandsRightHandSideOnFace(face, time, solutionCoefficientsLeft, solutionCoefficientsRight);};

	return this->faceIntegrator_.integratePair(ptrFace, integrandFunction);
}

/// ************************************************
/// ***     Jacobian Element Matrix Functions    ***
/// ************************************************

/// \brief This function computes the integrand of local element integral Jacobian contributions
template<std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
LinearAlgebra::MiddleSizeMatrix UnsteadyNavierStokesAPI<DIM,NUMBER_OF_VARIABLES>::integrandJacobianAtElement
(
 Base::PhysicalElement<DIM>& element,
 const LinearAlgebra::MiddleSizeVector &solutionCoefficients,
 const double time
)
{
	//Create integrand matrix
	std::size_t numberOfBasisFunctions =  element.getNumberOfBasisFunctions();
	LinearAlgebra::MiddleSizeMatrix integrand(NUMBER_OF_VARIABLES*numberOfBasisFunctions,NUMBER_OF_VARIABLES*numberOfBasisFunctions);

	//reconstruct the solution and all other solution structures
	StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> elementStateStruct = computeElementStateStruct(element, solutionCoefficients, time);

	// Inviscid element integral contribution
	integrand += inviscidTerms_.integrandJacobianInviscidElement(element, elementStateStruct);

	// Viscous element integral contribution
	integrand += viscousTerms_.integrandJacobianViscousElement(element, elementStateStruct);

	return integrand;
}

/// \brief This function computes the local Jacobian contributions from the element solutionCoefficients
template<std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
LinearAlgebra::MiddleSizeMatrix UnsteadyNavierStokesAPI<DIM,NUMBER_OF_VARIABLES>::computeJacobianAtElement
(
 Base::Element *ptrElement,
 const LinearAlgebra::MiddleSizeVector &solutionCoefficients,
 const double time
)
{
	std::function<LinearAlgebra::MiddleSizeMatrix(Base::PhysicalElement<DIM> &)> integrandFunction = [=](Base::PhysicalElement<DIM>& element) -> LinearAlgebra::MiddleSizeMatrix
			{   return integrandJacobianAtElement(element, solutionCoefficients, time);};

	return this->elementIntegrator_.integrate(ptrElement, integrandFunction, ptrElement->getGaussQuadratureRule());
}

/// ************************************************
/// ***      Jacobian  Face Matrix Functions     ***
/// ************************************************

/// \brief This function computes the integrand of local and non local face integral Jacobian contributions
template<std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
LinearAlgebra::MiddleSizeMatrix UnsteadyNavierStokesAPI<DIM,NUMBER_OF_VARIABLES>::integrandJacobianAtFace
(
 Base::PhysicalFace<DIM>& face,
 const LinearAlgebra::MiddleSizeVector &solutionCoefficientsLeft,
 const LinearAlgebra::MiddleSizeVector &solutionCoefficientsRight,
 const  Base::Side elementSide,
 const Base::Side derivativeSide,
 const double time
 )
{
	//Compute all functions that depend on the state coefficients, required to calculate the integrands
	StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> faceStateStructLeft = computeFaceStateStruct(face, solutionCoefficientsLeft, Base::Side::LEFT,time);
	StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> faceStateStructRight = computeFaceStateStruct(face, solutionCoefficientsRight, Base::Side::RIGHT,time);

	//Compute inviscid contribution
	LinearAlgebra::MiddleSizeMatrix integrandInviscid =  inviscidTerms_.integrandJacobianInviscidFace(face, faceStateStructLeft, faceStateStructRight, elementSide, derivativeSide);

	//Compute viscous contribution
	LinearAlgebra::MiddleSizeMatrix integrandViscous =  viscousTerms_.integrandJacobianViscousFace(face, faceStateStructLeft, faceStateStructRight, elementSide, derivativeSide);

	//Compute auxilliary contribution
	LinearAlgebra::MiddleSizeMatrix integrandAuxilliary = viscousTerms_.integrandJacobianAuxilliaryFace(face, faceStateStructLeft, faceStateStructRight, elementSide, derivativeSide);

	return integrandInviscid + integrandViscous + 0.0*integrandAuxilliary;
}

/// \brief This functions computes the Jacobian face matrix of element elementSide with respect to variables of derivativeSide
template<std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
LinearAlgebra::MiddleSizeMatrix UnsteadyNavierStokesAPI<DIM,NUMBER_OF_VARIABLES>::computeJacobianAtFace
(
 Base::Face *ptrFace,
 LinearAlgebra::MiddleSizeVector &solutionCoefficientsLeft,
 LinearAlgebra::MiddleSizeVector &solutionCoefficientsRight,
 const double time,
 Base::Side elementSide,
 Base::Side derivativeSide
)
{
	// Define the integrand function for the right hand side for the face.
	std::function<LinearAlgebra::MiddleSizeMatrix(Base::PhysicalFace<DIM>&)> integrandFunction = [=](Base::PhysicalFace<DIM>& face) -> LinearAlgebra::MiddleSizeMatrix
			{   return integrandJacobianAtFace(face, solutionCoefficientsLeft, solutionCoefficientsRight, elementSide, derivativeSide, time);};

	return this->faceIntegrator_.integrate(ptrFace, integrandFunction);
}

/// *****************************************
/// ***    		Various Functions         ***
/// *****************************************

template<std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
void UnsteadyNavierStokesAPI<DIM,NUMBER_OF_VARIABLES>::tasksBeforeSolving
(
 )
{
	//Stability parameters in the auxiliary flux needs the mass matrix
	//For a single element create the mass matrix: note this breaks down with p-refinement or limiters
    //the result of getElementsList() no longer provides a subscript-operator, so I traced the iterator to the beginning instead -FB
	LinearAlgebra::MiddleSizeMatrix stabilityMassMatrix  = this->computeMassMatrixAtElement(*this->meshes_[0]->elementColBegin());
	viscousTerms_.setStabilityMassMatrix(stabilityMassMatrix);
}

/// \brief Shows the progress in the terminal as output
template<std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
void UnsteadyNavierStokesAPI<DIM,NUMBER_OF_VARIABLES>::showProgress
(
 const double time,
 const std::size_t timeStepID
)
{
	time_ = time;
 	std::cout << "Time is " << time << std::endl;
 	std::cout << "timeStepID is " << timeStepID << std::endl;
}
