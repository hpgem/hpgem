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

#include "CompressibleNavierStokes.h"
#include "Inviscid.h"
#include "Viscous.h"

#include <cmath>

CompressibleNavierStokes::CompressibleNavierStokes
(
const std::size_t numOfVariables,
const double endTime,
const std::size_t polynomialOrder,
const Base::ButcherTableau * const ptrButcherTableau
) :
HpgemAPISimplified<DIM>(numOfVariables, polynomialOrder, ptrButcherTableau),
DIM_(DIM),
numOfVariables_(numOfVariables),
inviscidTerms_(*this),
viscousTerms_(*this)
{
}

// todo: rewrite all the normal vector queries etc.

/// \brief General mesh description
Base::RectangularMeshDescriptor<DIM> CompressibleNavierStokes::createMeshDescription(const std::size_t numOfElementPerDirection)
{
    // Create the domain. In this case the domain is the square [0,1]^DIM and periodic.
    Base::RectangularMeshDescriptor<DIM> description;
    for (std::size_t i = 0; i < DIM_; ++i)
    {
        description.bottomLeft_[i] = 0;
        description.topRight_[i] = 1;
        description.numElementsInDIM_[i] = numOfElementPerDirection;
    }
    description.boundaryConditions_[0] = Base::BoundaryType::PERIODIC;
    description.boundaryConditions_[1] = Base::BoundaryType::PERIODIC;

    return description;
}

void CompressibleNavierStokes::setStabilityMassMatrix()
{
	//For a single element create the mass matrix: note this breaks down with p-refinement or limiters
	LinearAlgebra::MiddleSizeMatrix stabilityMassMatrix  = computeMassMatrixAtElement(meshes_[0]->getElementsList()[0]);
	viscousTerms_.setStabilityMassMatrix(stabilityMassMatrix);
}

double CompressibleNavierStokes::computePressure(const LinearAlgebra::MiddleSizeVector &state)
{
	//Compute pressure term
	double pressure = 0.0;
	double q1Inverse = 1.0/state(0);
	for(std::size_t iD = 0; iD < DIM_; iD++)
	{
		pressure += state(iD+1)*state(iD+1); // (u^2 + v^2 + w^2)*rho^2
	}

	pressure = (gamma_ -1)*(state(DIM_+1)*nonDIM1_ - 0.5*q1Inverse*(pressure)); // (gamma-1)*rho*(e*nondim1- (u^2 + v^2 + w^2)/2), where nondim1 is scaling such that it is dimensionless

	logger.assert(pressure > 0, "Negative pressure.");

	return pressure;
}


/// *************************************************
/// ***   Element integration support functions   ***
/// *************************************************

///  \brief Constructs the solution based on the solutionCoefficients.
LinearAlgebra::MiddleSizeVector CompressibleNavierStokes::computeStateOnElement(Base::PhysicalElement<DIM> &element, const LinearAlgebra::MiddleSizeVector &solutionCoefficients)
{
		std::size_t numberOfBasisFunctions =  element.getNumOfBasisFunctions();
		LinearAlgebra::MiddleSizeVector elementState(numOfVariables_);
		std::size_t iVB; // Index in solution coefficients for variable i and basisfunction j

		for(std::size_t iV = 0; iV < numOfVariables_; iV++)
		{
			elementState(iV) = 0.0;
			for (std::size_t iB = 0; iB < numberOfBasisFunctions; iB++)
			{
				iVB = element.convertToSingleIndex(iB,iV);
				elementState(iV) += solutionCoefficients(iVB)*element.basisFunction(iB);
			}
		}

		return elementState;
}

LinearAlgebra::MiddleSizeMatrix CompressibleNavierStokes::computeStateJacobianAtElement(Base::PhysicalElement<DIM> &element, const LinearAlgebra::MiddleSizeVector &solutionCoefficients)
{
		std::size_t numberOfBasisFunctions =  element.getNumOfBasisFunctions();
		std::size_t iVB; // Index in solution coefficients for variable i and basisfunction j

		LinearAlgebra::MiddleSizeMatrix solutionGradient(numOfVariables_,DIM_);
		LinearAlgebra::SmallVector<DIM> gradientBasisFunction;

		for (std::size_t iB = 0; iB < numberOfBasisFunctions; iB++) //Note, For the current purpose the derivative of rhoE is not required.
		{
			gradientBasisFunction = element.basisFunctionDeriv(iB);
			for (std::size_t iV = 0; iV < numOfVariables_; iV++)
			{
				iVB = element.convertToSingleIndex(iB,iV);
				for (std::size_t iD = 0; iD < DIM_; iD++)
				{
					solutionGradient(iV,iD) += solutionCoefficients(iVB)*gradientBasisFunction(iD);
				}
			}
		}

		return solutionGradient;
}

LinearAlgebra::MiddleSizeVector CompressibleNavierStokes::computePartialState(const LinearAlgebra::MiddleSizeVector &state)
{
	LinearAlgebra::MiddleSizeVector partialState(DIM_ + 2);
	double q1Inverse = 1.0/state(0);

	partialState(0) = state(0);

	for (std::size_t iD = 0; iD < DIM_ + 1; iD++)
	{
		partialState(iD+1) = state(iD+1)*q1Inverse;
	}
	return partialState;
}


/// *****************************************
/// ***   Element integration functions   ***
/// *****************************************

/// \brief computes the source at an element
LinearAlgebra::MiddleSizeVector CompressibleNavierStokes::integrandSourceAtElement(Base::PhysicalElement<DIM>& element, const LinearAlgebra::MiddleSizeVector &state, const double &pressureTerm, const double &time)
{
	std::size_t numOfBasisFunctions = element.getNumOfBasisFunctions();

	LinearAlgebra::MiddleSizeVector integrandSource(numOfVariables_ * numOfBasisFunctions);

	//Empty source function: Put your own source function here.

	return integrandSource;
}

LinearAlgebra::MiddleSizeVector CompressibleNavierStokes::integrandRightHandSideOnRefElement(Base::PhysicalElement<DIM>& element, const double &time, const LinearAlgebra::MiddleSizeVector &stateCoefficients)
{

	//reconstruct the solution, partial state jacobian and pressure at pRef
	//todo: check if partialState and its Jacobian can somehow be computed efficiently together
	const LinearAlgebra::MiddleSizeVector state = computeStateOnElement(element, stateCoefficients);
	const LinearAlgebra::MiddleSizeMatrix stateJacobian = computeStateJacobianAtElement(element, stateCoefficients);
	const LinearAlgebra::MiddleSizeVector partialState = computePartialState(state);
	const double pressure = computePressure(state);

	//Compute inviscid terms
	LinearAlgebra::MiddleSizeVector integrandInviscid = inviscidTerms_.integrandAtElement(element, time, pressure, state);

	//Compute viscous terms
	LinearAlgebra::MiddleSizeVector integrandViscous = viscousTerms_.integrandAtElement(element,  state, stateJacobian, pressure, partialState);

	//Compute source terms
	LinearAlgebra::MiddleSizeVector integrandSource = integrandSourceAtElement(element, state, pressure, time);

    return (integrandInviscid + integrandViscous + integrandSource);

}

LinearAlgebra::MiddleSizeVector CompressibleNavierStokes::computeRightHandSideAtElement(Base::Element *ptrElement, LinearAlgebra::MiddleSizeVector &stateCoefficients, const double time)
{
	std::function<LinearAlgebra::MiddleSizeVector(Base::PhysicalElement<DIM> &)> integrandFunction = [=](Base::PhysicalElement<DIM>& element) -> LinearAlgebra::MiddleSizeVector
	    {   return this->integrandRightHandSideOnRefElement(element, time, stateCoefficients);};

    return elementIntegrator_.integrate(ptrElement, integrandFunction, ptrElement->getGaussQuadratureRule());

}


/// *************************************************
/// ***    face integration support functions     ***
/// *************************************************

LinearAlgebra::MiddleSizeVector CompressibleNavierStokes::computeStateOnFace(Base::PhysicalFace<DIM> &face, const Base::Side &iSide, const LinearAlgebra::MiddleSizeVector &stateCoefficients) const
{
	std::size_t numOfBasisFunctions =  face.getPhysicalElement(iSide).getNumOfBasisFunctions();
	LinearAlgebra::MiddleSizeVector elementState(numOfVariables_);
	std::size_t iVB; // Index in solution coefficients for variable i and basisfunction j

	for(std::size_t iV = 0; iV < numOfVariables_; iV++)
	{
		elementState(iV) = 0.0;
		for (std::size_t iB = 0; iB < numOfBasisFunctions; iB++)
		{
			iVB = face.getPhysicalElement(iSide).convertToSingleIndex(iB,iV);
			elementState(iV) += stateCoefficients(iVB)*face.basisFunction(iSide, iB); //basisFunction returns physical value
		}
	}

	return elementState;
}


LinearAlgebra::MiddleSizeMatrix CompressibleNavierStokes::computeStateJacobianAtFace(Base::PhysicalFace<DIM> &face, const Base::Side &iSide, const LinearAlgebra::MiddleSizeVector &stateCoefficients)
{
	std::size_t numberOfBasisFunctions =  face.getPhysicalElement(iSide).getNumOfBasisFunctions();
	std::size_t iVB; // Index in solution coefficients for variable i and basisfunction j

	LinearAlgebra::MiddleSizeMatrix stateJacobian(numOfVariables_,DIM_);
	LinearAlgebra::SmallVector<DIM> gradientBasisFunction;

	for (std::size_t iB = 0; iB < numberOfBasisFunctions; iB++)
	{
		gradientBasisFunction = face.basisFunctionDeriv(iSide,iB);
		for (std::size_t iV = 0; iV < numOfVariables_; iV++)
		{
			iVB = face.getPhysicalElement(iSide).convertToSingleIndex(iB,iV);
			for (std::size_t iD = 0; iD < DIM_; iD++)
			{
				stateJacobian(iV,iD) += stateCoefficients(iVB)*gradientBasisFunction(iD);
			}
		}
	}

	return stateJacobian;
}


/// **************************************************
/// ***    external face integration functions     ***
/// **************************************************

/// \brief Compute the integrand for the right hand side for the reference face corresponding to an external face.
//This function is written for the Couette type flow, a plate at top and a plate at bottom
//The state is reconstructed based on the plate's temperature and movement speed, and therefore this is a dirichlet type of BC
//It is handled like an internal element.
LinearAlgebra::MiddleSizeVector CompressibleNavierStokes::integrandRightHandSideOnRefFace(Base::PhysicalFace<DIM> &face, const double &time, const LinearAlgebra::MiddleSizeVector &stateCoefficients)
{
	//Compute the internal state
	const LinearAlgebra::MiddleSizeVector stateInternal = computeStateOnFace(face, Base::Side::LEFT, stateCoefficients);
	const LinearAlgebra::SmallVector<DIM> normalInternal = face.getNormalVector();
	//const double pressureInternal = computePressure(stateInternal);
	const LinearAlgebra::MiddleSizeVector partialStateInternal = computePartialState(stateInternal);
	const LinearAlgebra::MiddleSizeMatrix stateJacobianInternal = computeStateJacobianAtFace(face, Base::Side::LEFT, stateCoefficients);

	//Determine if this is the top or bottom boundary and set boundary state
	LinearAlgebra::MiddleSizeVector stateBoundary(DIM_+2);
	double temperatureBoundary;
	double velocityBoundary;
	double exponent;
	const Geometry::PointPhysical<DIM> pPhys = face.getPointPhysical();
	if(pPhys[1] > 0.8)
	{
		exponent = time*time/(Tc_*Tc_);
		velocityBoundary = uPlateTop_ - std::exp(-exponent);							// Spatial blending factor for the transient phase
		temperatureBoundary = tPlateTop_;
		stateBoundary(0) = (gamma_ - 1.0)*stateInternal(DIM_ + 1)/tPlateTop_/Rs_;	 	// density
		stateBoundary(1) = stateBoundary(0)*velocityBoundary; 								// velocity u on plate is uPlateTop_
		stateBoundary(2) = 0; 															//velocity v on plate is zero
		stateBoundary(DIM_+1) = stateInternal(DIM_ + 1);								// energy is copied from internal
	}
	else
	{
		exponent = time*time/(Tc_*Tc_);
		velocityBoundary = uPlateBottom_ - std::exp(-exponent);
		temperatureBoundary = tPlateBottom_;
		stateBoundary(0) = (gamma_ - 1.0)*stateInternal(DIM_ + 1)/tPlateBottom_/Rs_;	// density
		stateBoundary(1) = stateBoundary(0)*uPlateBottom_; 								// velocity u on plate is uPlateBottom_
		stateBoundary(2) = 0; 															// velocity v on plate is zero
		stateBoundary(DIM_+1) = stateInternal(DIM_ + 1);								// energy is copied from internal
	}
	const LinearAlgebra::MiddleSizeVector partialStateBoundary = computePartialState(stateBoundary);

	//Compute inviscid terms
	//todo:FIX ERROR: unitNormal goes in here
	LinearAlgebra::MiddleSizeVector integrandInviscid = inviscidTerms_.integrandAtFace(face, time, Base::Side::LEFT, stateInternal, stateBoundary, normalInternal);

	//Compute viscous terms
	LinearAlgebra::MiddleSizeVector integrandViscous = integrandInviscid; //viscousTerms_.integrandViscousAtFace(face, Base::Side::LEFT, stateInternal, stateBoundary, pressureInternal, partialStateInternal, normalInternal);


	//Compute support variable terms
	LinearAlgebra::MiddleSizeVector integrandAuxilliary = integrandInviscid; //viscousTerms_.integrandAuxilliaryAtFace(face, temperatureBoundary, stateInternal, stateBoundary, partialStateBoundary, normalInternal, stateJacobianInternal);

	return  integrandInviscid + integrandViscous + integrandAuxilliary;
}

/// \brief Compute the right-hand side corresponding to a boundary face
LinearAlgebra::MiddleSizeVector CompressibleNavierStokes::computeRightHandSideAtFace(Base::Face *ptrFace, LinearAlgebra::MiddleSizeVector &solutionCoefficients, const double time)
{
	    // Define the integrand function for the right hand side for the reference face.
	    std::function<LinearAlgebra::MiddleSizeVector(Base::PhysicalFace<DIM>&)> integrandFunction = [=](Base::PhysicalFace<DIM>& face) -> LinearAlgebra::MiddleSizeVector
	    {   return this->integrandRightHandSideOnRefFace(face, time, solutionCoefficients);};

	    return faceIntegrator_.integrate(ptrFace, integrandFunction);
}


/// **************************************************
/// ***    internal face integration functions     ***
/// **************************************************

/// \brief Compute the integrand for the right hand side for the reference face corresponding to an internal face.
   LinearAlgebra::MiddleSizeVector CompressibleNavierStokes::integrandRightHandSideOnRefFace(Base::PhysicalFace<DIM>& face, const double &time, const Base::Side &iSide, const LinearAlgebra::MiddleSizeVector &stateCoefficientsLeft, const LinearAlgebra::MiddleSizeVector &stateCoefficientsRight)
   {
       //reconstruct the solution, partial state Jacobian and pressure at pRef, left and right of the interface
       const LinearAlgebra::MiddleSizeVector stateLeft = computeStateOnFace(face, Base::Side::LEFT, stateCoefficientsLeft);
       const LinearAlgebra::MiddleSizeVector stateRight = computeStateOnFace(face, Base::Side::RIGHT, stateCoefficientsRight);
       const LinearAlgebra::MiddleSizeMatrix stateJacobianLeft = computeStateJacobianAtFace(face, Base::Side::LEFT, stateCoefficientsLeft);
       const LinearAlgebra::MiddleSizeMatrix stateJacobianRight = computeStateJacobianAtFace(face, Base::Side::RIGHT, stateCoefficientsRight);

		//Determine internal and external solutions
		LinearAlgebra::MiddleSizeVector stateInternal;
		LinearAlgebra::MiddleSizeVector stateExternal;
		LinearAlgebra::MiddleSizeMatrix stateJacobianInternal;
		LinearAlgebra::MiddleSizeMatrix stateJacobianExternal;
		LinearAlgebra::SmallVector<DIM> normalInternal;
		LinearAlgebra::SmallVector<DIM> normalExternal;
		LinearAlgebra::SmallVector<DIM> unitNormalInternal;
		LinearAlgebra::SmallVector<DIM> unitNormalExternal;
		LinearAlgebra::SmallVector<DIM> normal = face.getNormalVector();
		double area = Base::L2Norm(normal);

		if (iSide == Base::Side::RIGHT)
		{
			stateInternal = stateRight;
			stateExternal = stateLeft;
			stateJacobianInternal = stateJacobianRight;
			stateJacobianExternal = stateJacobianLeft;
			normalInternal = -normal;
			unitNormalInternal = normalInternal/area;
		}
		else
		{
			stateInternal = stateLeft;
			stateExternal = stateRight;
			stateJacobianInternal = stateJacobianLeft;
			stateJacobianExternal = stateJacobianRight;
			normalInternal = normal;
			unitNormalInternal = normalInternal/area;
		}

		const double pressureInternal = computePressure(stateInternal);
		const double pressureExternal = computePressure(stateExternal);
		const LinearAlgebra::MiddleSizeVector partialStateInternal = computePartialState(stateInternal);
		const LinearAlgebra::MiddleSizeVector partialStateExternal = computePartialState(stateExternal);

		//Compute inviscid terms

		LinearAlgebra::MiddleSizeVector integrandInviscid = inviscidTerms_.integrandAtFace(face, time, iSide, stateInternal, stateExternal, unitNormalInternal);

		//Compute viscous terms
		LinearAlgebra::MiddleSizeVector integrandViscous = viscousTerms_.integrandViscousAtFace(face, iSide, stateInternal, stateExternal, pressureInternal, partialStateInternal);

		//Compute support variable terms
		//todo: write out the integral as summ, see if things cancel
		//todo: Integrate over whole face in one go? More efficient?
		LinearAlgebra::MiddleSizeVector integrandAuxilliary = viscousTerms_.integrandAuxilliaryAtFace(face, iSide, stateInternal, stateExternal, pressureInternal, pressureExternal, partialStateInternal, partialStateExternal, stateJacobianInternal, stateJacobianExternal);

		return  area*integrandInviscid + integrandViscous + integrandAuxilliary;

   }

   /// \brief Compute the right-hand side corresponding to an internal face
   LinearAlgebra::MiddleSizeVector CompressibleNavierStokes::computeRightHandSideAtFace(Base::Face *ptrFace, const Base::Side side, LinearAlgebra::MiddleSizeVector &solutionCoefficientsLeft, LinearAlgebra::MiddleSizeVector &solutionCoefficientsRight, const double time)
   {
	    // Define the integrand function for the right hand side for the reference face.
	    std::function<LinearAlgebra::MiddleSizeVector(Base::PhysicalFace<DIM>&)> integrandFunction = [=](Base::PhysicalFace<DIM>& face) -> LinearAlgebra::MiddleSizeVector
	    {   return this->integrandRightHandSideOnRefFace(face, time, side, solutionCoefficientsLeft, solutionCoefficientsRight);};
	    return faceIntegrator_.integrate(ptrFace, integrandFunction);
   }



   /// *****************************************
   /// ***    		Various Functions        ***
   /// *****************************************

   LinearAlgebra::MiddleSizeVector CompressibleNavierStokes::getExactSolution(const PointPhysicalT &pPhys, const double &time, const std::size_t orderTimeDerivative)
   {
		LinearAlgebra::MiddleSizeVector exactSolution(numOfVariables_);

		double amplitude = 0.1;
		double frequency = 2.0*M_PI;
		double function = amplitude*std::cos(frequency*time);

		for (std::size_t iD = 0; iD < DIM_; iD++)
		{
			function *= std::cos(frequency*pPhys[iD]);
		}

		exactSolution(0) = 1.225;// + function;

		exactSolution(1) = exactSolution(0)*function;

		exactSolution(2) = exactSolution(0)*function;

		exactSolution(DIM_ + 1) = exactSolution(0)*718.0*288.0; // + 0.5*exactSolution(1)*exactSolution(1)/exactSolution(0); // + function;

	    return exactSolution;
   }

   /// \brief Compute the initial solution at a given point in space and time.
   LinearAlgebra::MiddleSizeVector CompressibleNavierStokes::getInitialSolution(const PointPhysicalT &pPhys, const double &startTime, const std::size_t orderTimeDerivative)
   {
       return getExactSolution(pPhys, startTime, orderTimeDerivative);
   }

   /// \brief Computes the error for output purposes
   LinearAlgebra::MiddleSizeVector CompressibleNavierStokes::Error(const double time)
   {
	   return computeMaxError(solutionTimeLevel_, time);
   }

   /// \brief Show the progress of the time integration.
   void CompressibleNavierStokes::showProgress(const double time, const std::size_t timeStepID)
   {
	   std::cout << "Time is " << time << std::endl;
	   std::cout << "timeStepID is " << timeStepID << std::endl;
   }




































