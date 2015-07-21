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

// todo: Opmaak van deze file
// todo: DIM_
// todo: check for the use of const in function names
// todo: check if thigns are passed by reference
// todo: rewrite qSolution to state
// todo: think beter about function names, reorganise code

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
        description.boundaryConditions_[i] = Base::BoundaryType::PERIODIC;
    }

    return description;
}

void CompressibleNavierStokes::setStabilityMassMatrix()
{
	//For a single element create the mass matrix: note this breaks down with p-refinement or limiters
	LinearAlgebra::MiddleSizeMatrix stabilityMassMatrix  = computeMassMatrixAtElement(meshes_[0]->getElementsList()[0]);
	viscousTerms_.setInverseStabilityMassMatrix(stabilityMassMatrix);
}

double CompressibleNavierStokes::computePressure(const LinearAlgebra::MiddleSizeVector &qSolution)
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
LinearAlgebra::MiddleSizeVector CompressibleNavierStokes::computeSolutionOnElement(const Base::Element *ptrElement, const LinearAlgebra::MiddleSizeVector &solutionCoefficients, const Geometry::PointReference<DIM> &pRef)
{
		std::size_t numOfBasisFunctions =  ptrElement->getNumberOfBasisFunctions();
		LinearAlgebra::MiddleSizeVector elementSolution(numOfVariables_);
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

LinearAlgebra::MiddleSizeMatrix CompressibleNavierStokes::computeSolutionJacobianAtElement(const Base::Element *ptrElement, const LinearAlgebra::MiddleSizeVector &solutionCoefficients, const Geometry::PointReference<DIM> &pRef)
{
		std::size_t numOfBasisFunctions =  ptrElement->getNumberOfBasisFunctions();
		std::size_t iVB; // Index in solution coefficients for variable i and basisfunction j

		LinearAlgebra::MiddleSizeMatrix solutionGradient(numOfVariables_,DIM_);
		LinearAlgebra::MiddleSizeVector gradientBasisFunction(DIM_);

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

LinearAlgebra::MiddleSizeVector CompressibleNavierStokes::computePartialState(const LinearAlgebra::MiddleSizeVector qSolution)
{
	LinearAlgebra::MiddleSizeVector partialState(DIM_ + 2);
	double q1Inverse = 1.0/qSolution(0);

	partialState(0) = qSolution(0);

	for (std::size_t iD = 0; iD < DIM_ + 1; iD++)
	{
		partialState(iD+1) = qSolution(iD+1)*q1Inverse;
	}

	return partialState;
}

/// \brief computes the source at an element
LinearAlgebra::MiddleSizeVector CompressibleNavierStokes::integrandSourceAtElement(Base::PhysicalElement<DIM>& element, const LinearAlgebra::MiddleSizeVector qSolution, const double pressureTerm, const double &time)
{
	std::size_t numOfBasisFunctions = element.getElement()->getNumberOfBasisFunctions();
	std::size_t iVB;

	LinearAlgebra::MiddleSizeVector integrandSource(numOfVariables_ * numOfBasisFunctions);

	//Empty source function: Put your own source function here.

	return integrandSource;
}

LinearAlgebra::MiddleSizeVector CompressibleNavierStokes::integrandRightHandSideOnRefElement(Base::PhysicalElement<DIM>& element, const double &time, const LinearAlgebra::MiddleSizeVector &solutionCoefficients)
{
    const Base::Element* ptrElement = element.getElement();
    const Geometry::PointReference<DIM>& pRef = element.getPointReference();
	//reconstruct the solution, partial state jacobian and pressure at pRef
	//todo: check if partialState and its Jacobian can somehow be computed efficiently together
	const LinearAlgebra::MiddleSizeVector qSolution = computeSolutionOnElement(ptrElement, solutionCoefficients, pRef);
	const LinearAlgebra::MiddleSizeMatrix qSolutionJacobian = computeSolutionJacobianAtElement(ptrElement, solutionCoefficients, pRef);
	const LinearAlgebra::MiddleSizeVector partialState = computePartialState(qSolution);
	const double pressure = computePressure(qSolution);

	//Compute inviscid terms
	LinearAlgebra::MiddleSizeVector integrandInviscid = inviscidTerms_.integrandAtElement(element, time, pressure, qSolution);

	//Compute viscous terms
	LinearAlgebra::MiddleSizeVector integrandViscous = viscousTerms_.integrandAtElement(element,  qSolution, qSolutionJacobian, pressure, partialState);

	//Compute source terms
	LinearAlgebra::MiddleSizeVector integrandSource = integrandSourceAtElement(element, qSolution, pressure, time);

    return integrandInviscid + integrandViscous + integrandSource;

}

LinearAlgebra::MiddleSizeVector CompressibleNavierStokes::computeRightHandSideAtElement(Base::Element *ptrElement, LinearAlgebra::MiddleSizeVector &solutionCoefficients, const double time)
{
	std::function<LinearAlgebra::MiddleSizeVector(Base::PhysicalElement<DIM> &)> integrandFunction = [=](Base::PhysicalElement<DIM>& element) -> LinearAlgebra::MiddleSizeVector
	    {   return this->integrandRightHandSideOnRefElement(element, time, solutionCoefficients);};

    return elementIntegrator_.integrate(ptrElement, integrandFunction, ptrElement->getGaussQuadratureRule());

}

/// *****************************************
/// ***    face integration functions     ***
/// *****************************************

LinearAlgebra::MiddleSizeMatrix CompressibleNavierStokes::computeSolutionJacobianAtFace(const Base::Face *ptrFace, const Base::Side &iSide, const LinearAlgebra::MiddleSizeVector &solutionCoefficients, const Geometry::PointReference<DIM - 1> &pRef)
{
	std::size_t numOfBasisFunctions =  ptrFace->getPtrElement(iSide)->getNumberOfBasisFunctions();
	std::size_t iVB; // Index in solution coefficients for variable i and basisfunction j

	LinearAlgebra::MiddleSizeMatrix solutionGradient(numOfVariables_,DIM_);
	LinearAlgebra::MiddleSizeVector gradientBasisFunction(DIM_);

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

LinearAlgebra::MiddleSizeVector CompressibleNavierStokes::computeSolutionOnFace(const Base::Face *ptrFace, const Base::Side &iSide, const LinearAlgebra::MiddleSizeVector &solutionCoefficients, const Geometry::PointReference<DIM - 1> &pRef) const
{
	std::size_t numOfBasisFunctions =  ptrFace->getPtrElement(iSide)->getNumberOfBasisFunctions();
	LinearAlgebra::MiddleSizeVector elementSolution(numOfVariables_);
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
LinearAlgebra::MiddleSizeVector CompressibleNavierStokes::integrandRightHandSideOnRefFace(Base::PhysicalFace<DIM>& face, const double &time, const LinearAlgebra::MiddleSizeVector &solutionCoefficients)
{

	//Compute inviscid terms
	LinearAlgebra::MiddleSizeVector integrandInviscid = inviscidTerms_.integrandAtFace(face, time, solutionCoefficients);

	//Compute viscous terms
	//todo: write viscousTerms_.integrandAtFace()
	LinearAlgebra::MiddleSizeVector integrandViscous = integrandInviscid; //integrandViscousAtFace();

	//Compute support variable terms
	//todo: write viscousTerms_.SupportAtFace()
	LinearAlgebra::MiddleSizeVector integrandSupport = integrandInviscid; // integrandSupportAtFace();

	return  integrandInviscid + integrandViscous + integrandSupport;
}

/// \brief Compute the integrand for the right hand side for the reference face corresponding to an internal face.
   LinearAlgebra::MiddleSizeVector CompressibleNavierStokes::integrandRightHandSideOnRefFace(Base::PhysicalFace<DIM>& face, const double &time, const Base::Side &iSide, const LinearAlgebra::MiddleSizeVector &solutionCoefficientsLeft, const LinearAlgebra::MiddleSizeVector &solutionCoefficientsRight)
   {
       const Base::Face* ptrFace = face.getFace();
       const Geometry::PointReference<DIM - 1>& pRef = face.getPointReference();
		//reconstruct the solution, partial state Jacobian and pressure at pRef, left and right of the interface
		const LinearAlgebra::MiddleSizeVector stateLeft = computeSolutionOnFace(ptrFace, Base::Side::LEFT, solutionCoefficientsLeft, pRef);
		const LinearAlgebra::MiddleSizeVector stateRight = computeSolutionOnFace(ptrFace, Base::Side::RIGHT, solutionCoefficientsRight, pRef);
		const LinearAlgebra::MiddleSizeMatrix stateJacobianLeft = computeSolutionJacobianAtFace(ptrFace, Base::Side::LEFT, solutionCoefficientsLeft, pRef);
		const LinearAlgebra::MiddleSizeMatrix stateJacobianRight = computeSolutionJacobianAtFace(ptrFace, Base::Side::RIGHT, solutionCoefficientsRight, pRef);

		//Determine internal and external solutions
		LinearAlgebra::MiddleSizeVector stateInternal;
		LinearAlgebra::MiddleSizeVector stateExternal;
		LinearAlgebra::MiddleSizeMatrix stateJacobianInternal;
		LinearAlgebra::MiddleSizeMatrix stateJacobianExternal;
		LinearAlgebra::MiddleSizeVector normalInternal;
		LinearAlgebra::MiddleSizeVector normalExternal;
		LinearAlgebra::MiddleSizeVector normal = ptrFace->getNormalVector(pRef);
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
		const LinearAlgebra::MiddleSizeVector partialStateInternal = computePartialState(stateInternal);
		const LinearAlgebra::MiddleSizeVector partialStateExternal = computePartialState(stateExternal);

		//Compute inviscid terms
		LinearAlgebra::MiddleSizeVector integrandInviscid = inviscidTerms_.integrandAtFace(face, time, iSide, stateInternal, stateExternal);

		//Compute viscous terms
		LinearAlgebra::MiddleSizeVector integrandViscous = viscousTerms_.integrandViscousAtFace(face, iSide, stateInternal, stateExternal, pressureInternal, partialStateInternal);

		//Compute support variable terms
		//todo: write out the integral as summ, see if things cancel
		//todo: Integrate over whole face in one go? More efficient?
		LinearAlgebra::MiddleSizeVector integrandAuxilliary = viscousTerms_.integrandAuxilliaryAtFace(face, iSide, stateInternal, stateExternal, pressureInternal, pressureExternal, partialStateInternal, partialStateExternal, stateJacobianInternal, stateJacobianExternal);

		// Note: correct with area, because the integration uses area, this might be avoided in future
		return  (integrandInviscid + integrandViscous + integrandAuxilliary)*area;

   }

   /// \brief Compute the right-hand side corresponding to a boundary face
   LinearAlgebra::MiddleSizeVector CompressibleNavierStokes::computeRightHandSideAtFace(Base::Face *ptrFace, LinearAlgebra::MiddleSizeVector &solutionCoefficients, const double time)
   {
	    // Define the integrand function for the right hand side for the reference face.
	    std::function<LinearAlgebra::MiddleSizeVector(Base::PhysicalFace<DIM>&)> integrandFunction = [=](Base::PhysicalFace<DIM>& face) -> LinearAlgebra::MiddleSizeVector
	    {   return this->integrandRightHandSideOnRefFace(face, time, solutionCoefficients);};

	    return faceIntegrator_.integrate(ptrFace, integrandFunction);
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
   LinearAlgebra::MiddleSizeVector CompressibleNavierStokes::getInitialSolution(const PointPhysicalT &pPhys, const double &startTime, const std::size_t orderTimeDerivative)
   {
       return getExactSolution(pPhys, startTime, orderTimeDerivative);
   }

   /// \brief Computes the error for output purposes
   LinearAlgebra::MiddleSizeVector CompressibleNavierStokes::Error(const double time)
   {
	   return computeMaxError(solutionTimeLevel_, time);
   }






































