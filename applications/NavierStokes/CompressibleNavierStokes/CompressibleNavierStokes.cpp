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

CompressibleNavierStokes::CompressibleNavierStokes
(
const std::size_t dimension,
const std::size_t numOfVariables,
const double endTime,
const std::size_t polynomialOrder,
const Base::ButcherTableau * const ptrButcherTableau
) :
HpgemAPISimplified(dimension, numOfVariables, polynomialOrder, ptrButcherTableau),
DIM_(dimension),
numOfVariables_(numOfVariables),
inviscidTerms_(*this),
viscousTerms_(*this)
{
}

// todo: Opmaak van deze file
// todo: check for the use of const in function names
// todo: check if thigns are passed by reference
// todo: rewrite qSolution to state
// todo: think beter about function names, reorganise code

/// \brief General mesh description
Base::RectangularMeshDescriptor CompressibleNavierStokes::createMeshDescription(const std::size_t numOfElementPerDirection)
{
    // Create the domain. In this case the domain is the square [0,1]^DIM and periodic.
    Base::RectangularMeshDescriptor description(DIM_);
    for (std::size_t i = 0; i < DIM_; ++i)
    {
        description.bottomLeft_[i] = 0;
        description.topRight_[i] = 1;
        description.numElementsInDIM_[i] = numOfElementPerDirection;
        description.boundaryConditions_[i] = Base::BoundaryType::PERIODIC;
    }

    return description;
}

void CompressibleNavierStokes::setStabilityMassMatrix()
{
	//For a single element create the mass matrix: note this breaks down with p-refinement or limiters
	LinearAlgebra::Matrix stabilityMassMatrix  = computeMassMatrixAtElement(meshes_[0]->getElementsList()[0]);
	viscousTerms_.setInverseStabilityMassMatrix(stabilityMassMatrix);
}

double CompressibleNavierStokes::computePressure(const LinearAlgebra::NumericalVector &qSolution)
{
	//Compute pressure term
	double pressure = 0.0;
	double q1Inverse = 1.0/qSolution(0);
	for(std::size_t iD = 0; iD < DIM_; iD++)
	{
		pressure += qSolution(iD+1)*qSolution(iD+1); // (u^2 + v^2 + w^2)*rho^2
	}
	pressure = (gamma_ -1)*(qSolution(DIM_+1) - 0.5*q1Inverse*(pressure)); // (gamma-1)*rho*(e- (u^2 + v^2 + w^2)/2)

	logger.assert(pressure > 0, "Negative pressure.");

	return pressure;
}

/// *****************************************
/// ***   Element integration functions   ***
/// *****************************************

///  \brief Constructs the solution based on the solutionCoefficients.
LinearAlgebra::NumericalVector CompressibleNavierStokes::computeSolutionOnElement(const Base::Element *ptrElement, const LinearAlgebra::NumericalVector &solutionCoefficients, const Geometry::PointReference &pRef)
{
		std::size_t numOfBasisFunctions =  ptrElement->getNrOfBasisFunctions();
		LinearAlgebra::NumericalVector elementSolution(numOfVariables_);
		std::size_t iVB; // Index in solution coefficients for variable i and basisfunction j

		for(std::size_t iV = 0; iV < numOfVariables_; iV++)
		{
			elementSolution(iV) = 0.0;
			for (std::size_t iB = 0; iB < numOfBasisFunctions; iB++)
			{
				iVB = ptrElement->convertToSingleIndex(iB,iV);
				elementSolution(iV) += solutionCoefficients(iVB)*ptrElement->basisFunction(iB, pRef); //basisFunction returns physical value
			}
		}

		return elementSolution;
}

LinearAlgebra::Matrix CompressibleNavierStokes::computeSolutionJacobianAtElement(const Base::Element *ptrElement, const LinearAlgebra::NumericalVector &solutionCoefficients, const Geometry::PointReference &pRef)
{
		std::size_t numOfBasisFunctions =  ptrElement->getNrOfBasisFunctions();
		std::size_t iVB; // Index in solution coefficients for variable i and basisfunction j

		LinearAlgebra::Matrix solutionGradient(numOfVariables_,DIM_);
		LinearAlgebra::NumericalVector gradientBasisFunction(DIM_);

		for (std::size_t iB = 0; iB < numOfBasisFunctions; iB++) //Note, For the current purpose the derivative of rhoE is not required.
		{
			gradientBasisFunction = ptrElement->basisFunctionDeriv(iB, pRef);
			for (std::size_t iV = 0; iV < numOfVariables_; iV++)
			{
				iVB = ptrElement->convertToSingleIndex(iB,iV);
				for (std::size_t iD = 0; iD < DIM_; iD++)
				{
					solutionGradient(iV,iD) += solutionCoefficients(iVB)*gradientBasisFunction(iD);
				}
			}
		}

		return solutionGradient;
}

LinearAlgebra::NumericalVector CompressibleNavierStokes::computePartialState(const LinearAlgebra::NumericalVector qSolution)
{
	LinearAlgebra::NumericalVector partialState(DIM_ + 2);
	double q1Inverse = 1.0/qSolution(0);

	partialState(0) = qSolution(0);

	for (std::size_t iD = 0; iD < DIM_ + 1; iD++)
	{
		partialState(iD+1) = qSolution(iD+1)*q1Inverse;
	}

	return partialState;
}

/// \brief computes the source at an element
LinearAlgebra::NumericalVector CompressibleNavierStokes::integrandSourceAtElement(const Base::Element *ptrElement, const LinearAlgebra::NumericalVector qSolution, const double pressureTerm, const double &time, const Geometry::PointReference &pRef)
{
	std::size_t numOfBasisFunctions = ptrElement->getNrOfBasisFunctions();
	std::size_t iVB;

	LinearAlgebra::NumericalVector integrandSource(numOfVariables_ * numOfBasisFunctions);

	//Empty source function: Put your own source function here.

	return integrandSource;
}

LinearAlgebra::NumericalVector CompressibleNavierStokes::integrandRightHandSideOnRefElement(const Base::Element *ptrElement, const double &time, const Geometry::PointReference &pRef, const LinearAlgebra::NumericalVector &solutionCoefficients)
{
	//reconstruct the solution, partial state jacobian and pressure at pRef
	//todo: check if partialState and its Jacobian can somehow be computed efficiently together
	const LinearAlgebra::NumericalVector qSolution = computeSolutionOnElement(ptrElement, solutionCoefficients, pRef);
	const LinearAlgebra::Matrix qSolutionJacobian = computeSolutionJacobianAtElement(ptrElement, solutionCoefficients, pRef);
	const LinearAlgebra::NumericalVector partialState = computePartialState(qSolution);
	const double pressure = computePressure(qSolution);

	//Compute inviscid terms
	LinearAlgebra::NumericalVector integrandInviscid = inviscidTerms_.integrandAtElement(ptrElement, time, pRef, pressure, qSolution);

	//Compute viscous terms
	LinearAlgebra::NumericalVector integrandViscous = viscousTerms_.integrandAtElement(ptrElement,  qSolution, qSolutionJacobian, pressure, partialState, pRef);

	//Compute source terms
	LinearAlgebra::NumericalVector integrandSource = integrandSourceAtElement(ptrElement, qSolution, pressure, time, pRef);

    return integrandInviscid + integrandViscous + integrandSource;

}

LinearAlgebra::NumericalVector CompressibleNavierStokes::computeRightHandSideAtElement(Base::Element *ptrElement, LinearAlgebra::NumericalVector &solutionCoefficients, const double time)
{
	std::function<LinearAlgebra::NumericalVector(const Base::Element*, const Geometry::PointReference &)> integrandFunction = [=](const Base::Element *El, const Geometry::PointReference & pRef) -> LinearAlgebra::NumericalVector
	    {   return this->integrandRightHandSideOnRefElement(El, time, pRef, solutionCoefficients);};

    return elementIntegrator_.integrate(ptrElement, integrandFunction, ptrElement->getGaussQuadratureRule());

}

/// *****************************************
/// ***    face integration functions     ***
/// *****************************************

LinearAlgebra::Matrix CompressibleNavierStokes::computeSolutionJacobianAtFace(const Base::Face *ptrFace, const Base::Side &iSide, const LinearAlgebra::NumericalVector &solutionCoefficients, const Geometry::PointReference &pRef)
{
	std::size_t numOfBasisFunctions =  ptrFace->getPtrElement(iSide)->getNrOfBasisFunctions();
	std::size_t iVB; // Index in solution coefficients for variable i and basisfunction j

	LinearAlgebra::Matrix solutionGradient(numOfVariables_,DIM_);
	LinearAlgebra::NumericalVector gradientBasisFunction(DIM_);

	for (std::size_t iB = 0; iB < numOfBasisFunctions; iB++) //Note, For the current purpose the derivative of rhoE is not required.
	{
		gradientBasisFunction = ptrFace->basisFunctionDeriv(iSide, iB, pRef);
		for (std::size_t iV = 0; iV < numOfVariables_; iV++)
		{
			iVB = ptrFace->getPtrElement(iSide)->convertToSingleIndex(iB,iV);
			for (std::size_t iD = 0; iD < DIM_; iD++)
			{
				solutionGradient(iV,iD) += solutionCoefficients(iVB)*gradientBasisFunction(iD);
			}
		}
	}

	return solutionGradient;
}

LinearAlgebra::NumericalVector CompressibleNavierStokes::computeSolutionOnFace(const Base::Face *ptrFace, const Base::Side &iSide, const LinearAlgebra::NumericalVector &solutionCoefficients, const Geometry::PointReference &pRef) const
{
	std::size_t numOfBasisFunctions =  ptrFace->getPtrElement(iSide)->getNrOfBasisFunctions();
	LinearAlgebra::NumericalVector elementSolution(numOfVariables_);
	std::size_t iVB; // Index in solution coefficients for variable i and basisfunction j

	for(std::size_t iV = 0; iV < numOfVariables_; iV++)
	{
		elementSolution(iV) = 0.0;
		for (std::size_t iB = 0; iB < numOfBasisFunctions; iB++)
		{
			iVB = ptrFace->getPtrElement(iSide)->convertToSingleIndex(iB,iV);
			elementSolution(iV) += solutionCoefficients(iVB)*ptrFace->basisFunction(iSide, iB, pRef); //basisFunction returns physical value
		}
	}

	return elementSolution;
}

/// \brief Compute the integrand for the right hand side for the reference face corresponding to an external face.
LinearAlgebra::NumericalVector CompressibleNavierStokes::integrandRightHandSideOnRefFace(const Base::Face *ptrFace, const double &time, const Geometry::PointReference &pRef, const LinearAlgebra::NumericalVector &solutionCoefficients)
{

	//Compute inviscid terms
	LinearAlgebra::NumericalVector integrandInviscid = inviscidTerms_.integrandAtFace(ptrFace, time, pRef, solutionCoefficients);

	//Compute viscous terms
	//todo: write viscousTerms_.integrandAtFace()
	LinearAlgebra::NumericalVector integrandViscous = integrandInviscid; //integrandViscousAtFace();

	//Compute support variable terms
	//todo: write viscousTerms_.SupportAtFace()
	LinearAlgebra::NumericalVector integrandSupport = integrandInviscid; // integrandSupportAtFace();

	return  integrandInviscid + integrandViscous + integrandSupport;
}

/// \brief Compute the integrand for the right hand side for the reference face corresponding to an internal face.
   LinearAlgebra::NumericalVector CompressibleNavierStokes::integrandRightHandSideOnRefFace(const Base::Face *ptrFace, const double &time, const Geometry::PointReference &pRef, const Base::Side &iSide, const LinearAlgebra::NumericalVector &solutionCoefficientsLeft, const LinearAlgebra::NumericalVector &solutionCoefficientsRight)
   {
		//reconstruct the solution, partial state Jacobian and pressure at pRef, left and right of the interface
		const LinearAlgebra::NumericalVector stateLeft = computeSolutionOnFace(ptrFace, Base::Side::LEFT, solutionCoefficientsLeft, pRef);
		const LinearAlgebra::NumericalVector stateRight = computeSolutionOnFace(ptrFace, Base::Side::RIGHT, solutionCoefficientsRight, pRef);
		const LinearAlgebra::Matrix stateJacobianLeft = computeSolutionJacobianAtFace(ptrFace, Base::Side::LEFT, solutionCoefficientsLeft, pRef);
		const LinearAlgebra::Matrix stateJacobianRight = computeSolutionJacobianAtFace(ptrFace, Base::Side::RIGHT, solutionCoefficientsRight, pRef);

		//Determine internal and external solutions
		LinearAlgebra::NumericalVector stateInternal;
		LinearAlgebra::NumericalVector stateExternal;
		LinearAlgebra::Matrix stateJacobianInternal;
		LinearAlgebra::Matrix stateJacobianExternal;
		LinearAlgebra::NumericalVector normalInternal;
		LinearAlgebra::NumericalVector normalExternal;
		LinearAlgebra::NumericalVector normal = ptrFace->getNormalVector(pRef);
		double area = Base::L2Norm(normal);
		if (iSide == Base::Side::RIGHT)
		{
			stateInternal = stateRight;
			stateExternal = stateLeft;
			stateJacobianInternal = stateJacobianRight;
			stateJacobianExternal = stateJacobianLeft;
			normalInternal = -normal/area;
		}
		else
		{
			stateInternal = stateLeft;
			stateExternal = stateRight;
			stateJacobianInternal = stateJacobianLeft;
			stateJacobianExternal = stateJacobianRight;
			normalInternal = normal/area;
		}

		const double pressureInternal = computePressure(stateInternal);
		const double pressureExternal = computePressure(stateExternal);
		const LinearAlgebra::NumericalVector partialStateInternal = computePartialState(stateInternal);
		const LinearAlgebra::NumericalVector partialStateExternal = computePartialState(stateExternal);

		//Compute inviscid terms
		LinearAlgebra::NumericalVector integrandInviscid = inviscidTerms_.integrandAtFace(ptrFace, time, pRef, iSide, stateInternal, stateExternal, normalInternal);

		//Compute viscous terms
		LinearAlgebra::NumericalVector integrandViscous = viscousTerms_.integrandViscousAtFace(ptrFace, iSide, stateInternal, stateExternal, pressureInternal, partialStateInternal, normalInternal, pRef);

		//Compute support variable terms
		//todo: write out the integral as summ, see if things cancel
		//todo: Integrate over whole face in one go? More efficient?
		LinearAlgebra::NumericalVector integrandAuxilliary = viscousTerms_.integrandAuxilliaryAtFace(ptrFace, iSide, stateInternal, stateExternal, pressureInternal, pressureExternal, partialStateInternal, partialStateExternal,  normalInternal, stateJacobianInternal, stateJacobianExternal, pRef);

		// Note: correct with area, because the integration uses area, this might be avoided in future
		return  (integrandInviscid + integrandViscous + integrandAuxilliary)*area;

   }

   /// \brief Compute the right-hand side corresponding to a boundary face
   LinearAlgebra::NumericalVector CompressibleNavierStokes::computeRightHandSideAtFace(Base::Face *ptrFace, LinearAlgebra::NumericalVector &solutionCoefficients, const double time)
   {
	    // Define the integrand function for the right hand side for the reference face.
	    std::function<LinearAlgebra::NumericalVector(const Geometry::PointReference &)> integrandFunction = [=](const Geometry::PointReference &pRef) -> LinearAlgebra::NumericalVector
	    {   return this->integrandRightHandSideOnRefFace(ptrFace, time, pRef, solutionCoefficients);};

	    return faceIntegrator_.referenceFaceIntegral(ptrFace->getGaussQuadratureRule(), integrandFunction);
   }

   /// \brief Compute the right-hand side corresponding to an internal face
   LinearAlgebra::NumericalVector CompressibleNavierStokes::computeRightHandSideAtFace(Base::Face *ptrFace, const Base::Side side, LinearAlgebra::NumericalVector &solutionCoefficientsLeft, LinearAlgebra::NumericalVector &solutionCoefficientsRight, const double time)
   {
	    // Define the integrand function for the right hand side for the reference face.
	    std::function<LinearAlgebra::NumericalVector(const Geometry::PointReference &)> integrandFunction = [=](const Geometry::PointReference &pRef) -> LinearAlgebra::NumericalVector
	    {   return this->integrandRightHandSideOnRefFace(ptrFace, time, pRef, side, solutionCoefficientsLeft, solutionCoefficientsRight);};
	    return faceIntegrator_.referenceFaceIntegral(ptrFace->getGaussQuadratureRule(), integrandFunction);
   }



   /// *****************************************
   /// ***    		Various Functions        ***
   /// *****************************************

   LinearAlgebra::NumericalVector CompressibleNavierStokes::getExactSolution(const PointPhysicalT &pPhys, const double &time, const std::size_t orderTimeDerivative)
   {
		LinearAlgebra::NumericalVector exactSolution(numOfVariables_);

		double amplitude = 0.2;
		double frequency = 2.0*M_PI;
		double function = amplitude*std::cos(frequency*time);

		for (std::size_t iD = 0; iD < DIM_; iD++)
		{
			function *= std::cos(frequency*pPhys[iD]);
		}

		exactSolution(0) = 1.5;// + function;

		for (std::size_t iD = 0; iD < DIM_; iD++)
		{
			exactSolution(iD+1) = 0.2; //function;
		}

		exactSolution(DIM_ + 1) = 30.0;// + function;

	    return exactSolution;
   }

   /// \brief Compute the initial solution at a given point in space and time.
   LinearAlgebra::NumericalVector CompressibleNavierStokes::getInitialSolution(const PointPhysicalT &pPhys, const double &startTime, const std::size_t orderTimeDerivative)
   {
       return getExactSolution(pPhys, startTime, orderTimeDerivative);
   }

   /// \brief Computes the error for output purposes
   LinearAlgebra::NumericalVector CompressibleNavierStokes::Error(const double time)
   {
	   return computeMaxError(solutionTimeLevel_, time);
   }






































