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

#include "ViscousTerms.h"

/// *****************************************
/// ***   Element integration functions   ***
/// *****************************************

/// \brief Computes the integrand used in the element integration routine.
template<std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
LinearAlgebra::MiddleSizeVector ViscousTerms<DIM,NUMBER_OF_VARIABLES>::integrandAtElement
(
 Base::PhysicalElement<DIM> &element,
 const StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> &elementStateStruct,
 const double time
)
{
	// Get the number of basis functions in an element.
	std::size_t numberOfBasisFunctions = element.getNumberOfBasisFunctions();

	// Create data structures for calculating the integrand
	LinearAlgebra::MiddleSizeVector integrand(NUMBER_OF_VARIABLES * numberOfBasisFunctions);
	LinearAlgebra::SmallVector<DIM> gradientBasisFunction;
	std::size_t iVB;

	//Calculate flux function
	LinearAlgebra::MiddleSizeMatrix fluxFunction = elementStateStruct.computeEllipticTensorMatrixContractionFast(elementStateStruct.getStateJacobian());

/*	std::cout << "stateJacobian: " << elementStateStruct.getStateJacobian() << std::endl;
	std::cout << "fluxFunction: " << fluxFunction << std::endl;*/

	for (std::size_t iB = 0; iB < numberOfBasisFunctions; iB++) // for all basis functions
	{
		gradientBasisFunction = element.basisFunctionDeriv(iB);
		for (std::size_t iD = 0; iD < DIM; iD++) // for all DIMs
		{
			for (std::size_t iV = 0; iV < NUMBER_OF_VARIABLES; iV++) // for all variables
			{
				iVB = element.convertToSingleIndex(iB,iV);
				integrand(iVB) += fluxFunction(iV,iD)*gradientBasisFunction(iD);
			}
		}
	}

	return -integrand; //Integral is on right hand side of the equation, hence the minus sign
}


/// **************************************************
/// ***    External face integration functions     ***
/// **************************************************

/// \brief Computes the flux function for the auxilliary value at an external face
template<std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
LinearAlgebra::MiddleSizeMatrix ViscousTerms<DIM,NUMBER_OF_VARIABLES>::fluxAuxilliaryBoundary
(
 Base::PhysicalFace<DIM> &face,
 const StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> &faceStateStructBoundary,
 const StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> &faceStateStructLeft,
 const double &time
 )
{
	LinearAlgebra::MiddleSizeMatrix fluxBoundary  = faceStateStructBoundary.computeEllipticTensorMatrixContractionFast(faceStateStructLeft.getStateJacobian());
	LinearAlgebra::MiddleSizeMatrix stabilityParametersBoundary = computeStabilityParametersBoundary(face, faceStateStructLeft.getStateCoefficients(), time);

	//std::cout << "fluxBoundary: " << fluxBoundary << std::endl;
	//std::cout << "stabilityParametersBoundary: " << stabilityParametersBoundary << std::endl;

	return (fluxBoundary - ETA_BOUNDARY*stabilityParametersBoundary);
}

/// \brief Computes the auxilliary integrand at an internal face
template<std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
LinearAlgebra::MiddleSizeVector ViscousTerms<DIM,NUMBER_OF_VARIABLES>::integrandAuxilliaryAtBoundaryFace
(
 Base::PhysicalFace<DIM> &face,
 const StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> &faceStateStructBoundary,
 const StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> &faceStateStructLeft,
 const double &time
 )
{
	//Data structures for left and right integrand
	std::size_t numberOfTestBasisFunctionsLeft = face.getPhysicalElement(Base::Side::LEFT).getNumberOfBasisFunctions();
	LinearAlgebra::MiddleSizeVector integrand(NUMBER_OF_VARIABLES*numberOfTestBasisFunctionsLeft);
	LinearAlgebra::SmallVector<DIM> unitNormalLeft = face.getUnitNormalVector();
	std::size_t iVB;

	//Compute Left flux
	//note: left flux is the same as the right flux
	LinearAlgebra::MiddleSizeMatrix fluxFunctionBoundary = fluxAuxilliaryBoundary(face, faceStateStructBoundary, faceStateStructLeft, time);

	//Compute Left integrand
	LinearAlgebra::MiddleSizeVector partialIntegrand(NUMBER_OF_VARIABLES);
	for (std::size_t iV = 0; iV < NUMBER_OF_VARIABLES; iV++)
	{
		for (std::size_t iD = 0; iD < DIM; iD++)
		{
			partialIntegrand(iV) += fluxFunctionBoundary(iV,iD)*unitNormalLeft(iD);
		}
	}

//	std::cout << "partialIntegrand: " << partialIntegrand << std::endl;

	for (std::size_t iV = 0; iV < NUMBER_OF_VARIABLES; iV++)
	{
		for (std::size_t iB = 0; iB < numberOfTestBasisFunctionsLeft; iB++)
		{
			iVB = face.getPhysicalElement(Base::Side::LEFT).convertToSingleIndex(iB,iV);
			integrand(iVB) += face.basisFunction(Base::Side::LEFT, iB)*partialIntegrand(iV);
		}

	}

	return integrand;
}

/// \brief Computes the viscous integral at the face, based on the viscous fluxes for an external face
template<std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
LinearAlgebra::MiddleSizeVector ViscousTerms<DIM,NUMBER_OF_VARIABLES>::integrandViscousAtBoundaryFace
(
 Base::PhysicalFace<DIM> &face,
 const StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> &faceStateStructBoundary,
 const StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> &faceStateStructLeft,
 const double &time
 )
{
 	std::size_t numberOfTestBasisFunctionsLeft = face.getPhysicalElement(Base::Side::LEFT).getNumberOfBasisFunctions();
 	LinearAlgebra::MiddleSizeVector integrand(NUMBER_OF_VARIABLES*numberOfTestBasisFunctionsLeft);
 	LinearAlgebra::SmallVector<DIM> gradientBasisFunction;
 	std::size_t iVB;

 	//todo: this function can be written slightly more efficient. Now A^{+}*(u^{+}*n) - A^{b}*(u^{b}*n) is computed while (A^{+}*u^{+} - A^{b}*u^{b})*n is faster.
 	// It however takes a more sophisticated function of computing a 4th rank tensor with a vector and is not very important at the moment.
 	//Compute state times normal matrix for the Left face, this is the same as for the Right side
 	//todo: move this to the StateCoefficientsFunctions class
 	LinearAlgebra::MiddleSizeMatrix stateNormalBoundary(NUMBER_OF_VARIABLES,DIM);
 	LinearAlgebra::MiddleSizeMatrix stateNormalLeft(NUMBER_OF_VARIABLES,DIM);
 	LinearAlgebra::SmallVector<DIM> unitNormalLeft = face.getUnitNormalVector();
 	for (std::size_t iV = 0; iV < NUMBER_OF_VARIABLES; iV++)
 	{
 		for (std::size_t iD = 0; iD < DIM; iD++)
 		{
 			stateNormalBoundary(iV,iD) = faceStateStructBoundary.getState()(iV)*unitNormalLeft(iD);
 			stateNormalLeft(iV,iD) = faceStateStructLeft.getState()(iV)*unitNormalLeft(iD);
 		}
 	}

 	//Compute flux
 	LinearAlgebra::MiddleSizeMatrix fluxBoundary = faceStateStructLeft.computeEllipticTensorMatrixContractionFast(stateNormalLeft) - faceStateStructBoundary.computeEllipticTensorMatrixContractionFast(stateNormalBoundary);


/*	std::cout << "viscosity Boundary: " << faceStateStructBoundary.getViscosity() << std::endl;
	std::cout << "viscosity left: " << faceStateStructLeft.getViscosity() << std::endl;
	std::cout << "elliptic11 boundary: " << faceStateStructBoundary.getEllipticTensor()[0] << std::endl;
	std::cout << "elliptic22 left: " << faceStateStructLeft.getEllipticTensor()[0] << std::endl; 
	std::cout << "statenormalboundary: " << stateNormalBoundary << std::endl;
	std::cout << "stateNormalLeft: " << stateNormalLeft << std::endl;
	std::cout << "Boundary contraction: " << faceStateStructBoundary.computeEllipticTensorMatrixContractionFast(stateNormalBoundary) << std::endl;
	std::cout << "left contraction: " << faceStateStructLeft.computeEllipticTensorMatrixContractionFast(stateNormalLeft) << std::endl;
*/
 	//Compute integrand Left
 	for (std::size_t iB = 0; iB < numberOfTestBasisFunctionsLeft; iB++) // for all basis functions
 	{
 		gradientBasisFunction = face.basisFunctionDeriv(Base::Side::LEFT,iB);
 		for (std::size_t iD = 0; iD < DIM; iD++) // for all DIMs
 		{
 			for (std::size_t iV = 0; iV < NUMBER_OF_VARIABLES; iV++) // for all equations
 			{
 				iVB = face.getPhysicalElement(Base::Side::LEFT).convertToSingleIndex(iB,iV);
 				integrand(iVB) += fluxBoundary(iV,iD)*gradientBasisFunction(iD);
 			}
 		}
 	}

 	return integrand;
}

/// ***************************************************
/// ***   External Stability Parameter functions    ***
/// ***************************************************
// todo: it is weird that left is always first in the arguments, but incase of aboundary, boundary is first
/// \brief Computes the fluxFunction in the stability parameter calculation used for the integrand
template<std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
LinearAlgebra::MiddleSizeMatrix ViscousTerms<DIM,NUMBER_OF_VARIABLES>::fluxStabilityParameters
(
 const StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> &faceStateStructBoundary,
 const LinearAlgebra::MiddleSizeVector &stateLeft,
 const LinearAlgebra::SmallVector<DIM> &normalInternal
 )
{
	//Compute velocity normal matrix
	LinearAlgebra::MiddleSizeMatrix stateNormal(NUMBER_OF_VARIABLES,DIM);
	LinearAlgebra::MiddleSizeVector stateDifference = stateLeft - faceStateStructBoundary.getState();
	for (std::size_t iV = 0; iV < NUMBER_OF_VARIABLES; iV++)
	{
		for (std::size_t iD = 0; iD < DIM; iD++)
		{
			stateNormal(iV,iD) = stateDifference(iV)*normalInternal(iD);
		}
	}

	return faceStateStructBoundary.computeEllipticTensorMatrixContractionFast(stateNormal);
}

/// \brief Computes the integrand required for the stability parameter calculations
template<std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
LinearAlgebra::MiddleSizeMatrix ViscousTerms<DIM,NUMBER_OF_VARIABLES>::integrandStabilityRightHandSideOnBoundaryFace
(
 Base::PhysicalFace<DIM> &face,
 const LinearAlgebra::MiddleSizeVector stateCoefficientsLeft,
 const double time
)
{
	std::size_t numberOfBasisFunctions = face.getPhysicalElement(Base::Side::LEFT).getNumberOfBasisFunctions();
	LinearAlgebra::MiddleSizeMatrix integrand(numberOfBasisFunctions*NUMBER_OF_VARIABLES,DIM);
	std::size_t iVB;

	//compute stabilityFluxFunction
	//Note: This can be done more efficiently
	const StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> faceStateStructLeft = instance_->computeFaceStateStruct(face, stateCoefficientsLeft, Base::Side::LEFT, time);
	const LinearAlgebra::MiddleSizeVector stateBoundary = instance_->computeBoundaryState(face, faceStateStructLeft, time);
	const StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> faceStateStructBoundary = instance_->computeBoundaryFaceStateStruct(stateBoundary, time);
	const LinearAlgebra::MiddleSizeVector stateLeft = computeStateOnFace<DIM,NUMBER_OF_VARIABLES>(face, Base::Side::LEFT, stateCoefficientsLeft);
	LinearAlgebra::MiddleSizeMatrix stabilityFluxFunction = fluxStabilityParameters(faceStateStructBoundary, stateLeft, face.getUnitNormalVector());

	//Compute integrand
	for (std::size_t iV = 0; iV < NUMBER_OF_VARIABLES; iV++)
	{
		for (std::size_t iB = 0; iB < numberOfBasisFunctions; iB++)
		{
			iVB = face.getPhysicalElement(Base::Side::LEFT).convertToSingleIndex(iB,iV);
			for (std::size_t iD = 0; iD < DIM; iD++)
			{
				integrand(iVB,iD) += stabilityFluxFunction(iV,iD)*face.basisFunction(Base::Side::LEFT, iB);
			}
		}
	}
	//std::cout << "-----------------" << std::endl;
	//std::cout << "stabilityflux: " << stabilityFluxFunction << std::endl;
	//std::cout << "integrand: " << integrand << std::endl;

	return integrand;
}

/// \brief Computes the rhs, for the system of equations solving the stability parameters, by integrating the rhs  stability parameter integrand
template<std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
LinearAlgebra::MiddleSizeMatrix ViscousTerms<DIM,NUMBER_OF_VARIABLES>::computeRhsStabilityParametersBoundary
(
 const Base::Face *ptrFace,
 const LinearAlgebra::MiddleSizeVector stateCoefficientsLeft,
 const double time
 )
{
	std::function<LinearAlgebra::MiddleSizeMatrix(Base::PhysicalFace<DIM>&)> integrandFunction = [=](Base::PhysicalFace<DIM> &face) -> LinearAlgebra::MiddleSizeMatrix
			{   return this->integrandStabilityRightHandSideOnBoundaryFace(face, stateCoefficientsLeft, time);};
	return stabilityFaceIntegrator_.integrate(ptrFace, integrandFunction);
}


/// \brief Computes the stability parameters used in the auxilliary integrand for an internal face for given stabilityParameterflux
template<std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
LinearAlgebra::MiddleSizeMatrix ViscousTerms<DIM,NUMBER_OF_VARIABLES>::computeStabilityParametersBoundary
(
 Base::PhysicalFace<DIM> &face,
 const LinearAlgebra::MiddleSizeVector stateCoefficientsLeft,
 const double time
 )
{
	//Datastructures
	std::size_t numberOfTestBasisFunctionsLeft = face.getPhysicalElement(Base::Side::LEFT).getNumberOfBasisFunctions();

	//TODO: NOTE that it is quicker to determin the inverse once and then multiply it forever with the inverse.
	//Compute stability parameters on the Left side
	LinearAlgebra::MiddleSizeMatrix rhsBoundary = computeRhsStabilityParametersBoundary(face.getFace(), stateCoefficientsLeft, time);
	this->stabilityMassMatrix_.solve(rhsBoundary);


	//Compute the average of the stability parameters
	return computeStateOnFace<DIM,NUMBER_OF_VARIABLES>(face, Base::Side::LEFT, rhsBoundary);
}


/// **************************************************
/// ***    Internal face integration functions     ***
/// **************************************************

/// \brief Computes the fluxFunction for the auxilliary values at an internal face
template<std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
LinearAlgebra::MiddleSizeMatrix ViscousTerms<DIM,NUMBER_OF_VARIABLES>::fluxAuxilliary
(
 Base::PhysicalFace<DIM> &face,
 const StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> &faceStateStructLeft,
 const StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> &faceStateStructRight
)
{
	LinearAlgebra::MiddleSizeMatrix fluxLeft  = faceStateStructLeft.computeEllipticTensorMatrixContractionFast(faceStateStructLeft.getStateJacobian());
	LinearAlgebra::MiddleSizeMatrix fluxRight = faceStateStructRight.computeEllipticTensorMatrixContractionFast(faceStateStructRight.getStateJacobian());

	//std::cout << "stateJacobian" << faceStateStructLeft.getStateJacobian() << std::endl;

	LinearAlgebra::MiddleSizeMatrix stabilityParameters = computeStabilityParameters(face, faceStateStructLeft.getStateCoefficients(), faceStateStructRight.getStateCoefficients());

	//std::cout << "==============" << std::endl;
	//std::cout << "State Left: " << faceStateStructLeft.getState() << std::endl;
	//std::cout << "State Right: " << faceStateStructRight.getState() << std::endl;
	//std::cout << "stabilityParameters: " << stabilityParameters << std::endl;

	return (0.5*(fluxLeft + fluxRight) - ETA*stabilityParameters);
}

/// \brief Compute both the auxilliary face integral for the left element as the right element at the same time. (reducing flux calculations)
template<std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
std::pair<LinearAlgebra::MiddleSizeVector,LinearAlgebra::MiddleSizeVector> ViscousTerms<DIM,NUMBER_OF_VARIABLES>::integrandsAuxilliaryAtFace
(
 Base::PhysicalFace<DIM> &face,
 const StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> &faceStateStructLeft,
 const StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> &faceStateStructRight,
 const double time
 )
{

	//Data structures for left and right integrand
	std::size_t numberOfTestBasisFunctionsLeft = face.getPhysicalElement(Base::Side::LEFT).getNumberOfBasisFunctions();
	std::size_t numberOfTestBasisFunctionsRight = face.getPhysicalElement(Base::Side::RIGHT).getNumberOfBasisFunctions();
	std::pair<LinearAlgebra::MiddleSizeVector,LinearAlgebra::MiddleSizeVector> integrands(
			std::piecewise_construct,
			std::forward_as_tuple(NUMBER_OF_VARIABLES*numberOfTestBasisFunctionsLeft),
			std::forward_as_tuple(NUMBER_OF_VARIABLES*numberOfTestBasisFunctionsRight));
	LinearAlgebra::SmallVector<DIM> unitNormalLeft = face.getUnitNormalVector();
	std::size_t iVB;

	//Compute Left flux
	//note: left flux is the same as the right flux
	LinearAlgebra::MiddleSizeMatrix fluxFunctionLeft = fluxAuxilliary(face, faceStateStructLeft, faceStateStructRight);

	//Compute Left integrand
	LinearAlgebra::MiddleSizeVector partialIntegrand(NUMBER_OF_VARIABLES);
	for (std::size_t iV = 0; iV < NUMBER_OF_VARIABLES; iV++)
	{
		for (std::size_t iD = 0; iD < DIM; iD++)
		{
			partialIntegrand(iV) += fluxFunctionLeft(iV,iD)*unitNormalLeft(iD);
		}
	}

	for (std::size_t iV = 0; iV < NUMBER_OF_VARIABLES; iV++)
	{
		for (std::size_t iB = 0; iB < numberOfTestBasisFunctionsLeft; iB++)
		{
			iVB = face.getPhysicalElement(Base::Side::LEFT).convertToSingleIndex(iB,iV);
			integrands.first(iVB) += face.basisFunction(Base::Side::LEFT, iB)*partialIntegrand(iV);
		}

		for (std::size_t iB = 0; iB < numberOfTestBasisFunctionsRight; iB++)
		{
			iVB = face.getPhysicalElement(Base::Side::RIGHT).convertToSingleIndex(iB,iV);
			integrands.second(iVB) += -face.basisFunction(Base::Side::RIGHT, iB)*partialIntegrand(iV);
		}
	}

	return integrands;
}

/// \brief Compute both the viscous face integral for the left element as the right element at the same time.
template<std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
std::pair<LinearAlgebra::MiddleSizeVector,LinearAlgebra::MiddleSizeVector> ViscousTerms<DIM,NUMBER_OF_VARIABLES>::integrandsViscousAtFace
(
 Base::PhysicalFace<DIM> &face,
 const StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> &faceStateStructLeft,
 const StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> &faceStateStructRight,
 const double time
 )
{
 	std::size_t numOfTestBasisFunctionsLeft = face.getPhysicalElement(Base::Side::LEFT).getNumberOfBasisFunctions();
 	std::size_t numOfTestBasisFunctionsRight = face.getPhysicalElement(Base::Side::RIGHT).getNumberOfBasisFunctions();
 	std::pair<LinearAlgebra::MiddleSizeVector,LinearAlgebra::MiddleSizeVector> integrands(
 			std::piecewise_construct,
 			std::forward_as_tuple(NUMBER_OF_VARIABLES*numOfTestBasisFunctionsLeft),
 			std::forward_as_tuple(NUMBER_OF_VARIABLES*numOfTestBasisFunctionsRight));
 	LinearAlgebra::SmallVector<DIM> unitNormalLeft = face.getUnitNormalVector();
 	LinearAlgebra::SmallVector<DIM> gradientBasisFunctionLeft;
 	LinearAlgebra::SmallVector<DIM> gradientBasisFunctionRight;
 	std::size_t iVB;

 	//Compute state times normal matrix for the Left face, this is the same as for the Right side
 	//todo: move this to the StateCoefficientsFunctions class
 	LinearAlgebra::MiddleSizeMatrix stateNormalLeft(NUMBER_OF_VARIABLES,DIM);
 	LinearAlgebra::MiddleSizeVector stateDifferenceLeft;
 	stateDifferenceLeft = faceStateStructLeft.getState() - faceStateStructRight.getState();
 	for (std::size_t iV = 0; iV < NUMBER_OF_VARIABLES; iV++)
 	{
 		for (std::size_t iD = 0; iD < DIM; iD++)
 		{
 			stateNormalLeft(iV,iD) = 0.5*stateDifferenceLeft(iV)*unitNormalLeft(iD);
 		}
 	}

 	//Compute fluxLeft
 	LinearAlgebra::MiddleSizeMatrix fluxLeft = faceStateStructLeft.computeEllipticTensorMatrixContractionFast(stateNormalLeft);

 	//compute fluxRight
 	LinearAlgebra::MiddleSizeMatrix fluxRight = faceStateStructLeft.computeEllipticTensorMatrixContractionFast(stateNormalLeft);

 	//Compute integrand Left
 	//note: when p differ in both elements, this cannot be done in the same loop
 	//todo: generalise this such that it works for p-refinement

 	for (std::size_t iV = 0; iV < NUMBER_OF_VARIABLES; iV++) // for all equations
 	{
 		for (std::size_t iB = 0; iB < numOfTestBasisFunctionsLeft; iB++) // for all basis functions
 		{
 			gradientBasisFunctionLeft = face.basisFunctionDeriv(Base::Side::LEFT,iB);
 			for (std::size_t iD = 0; iD < DIM; iD++) // for all DIMs
 			{
 				iVB = face.getPhysicalElement(Base::Side::LEFT).convertToSingleIndex(iB,iV);
 				integrands.first(iVB) += fluxLeft(iV,iD)*gradientBasisFunctionLeft(iD);
 			}
 		}

 		for (std::size_t iB = 0; iB < numOfTestBasisFunctionsRight; iB++) // for all basis functions
 		{
 			gradientBasisFunctionRight = face.basisFunctionDeriv(Base::Side::RIGHT,iB);
/* 			std::cout << "---" << std::endl;
 			std::cout << "fluxLeft: " << fluxLeft << std::endl;
 			std::cout << "griadient: " << gradientBasisFunctionRight << std::endl;*/
 			for (std::size_t iD = 0; iD < DIM; iD++) // for all DIMs
 			{
 				iVB = face.getPhysicalElement(Base::Side::RIGHT).convertToSingleIndex(iB,iV);
 				integrands.second(iVB) += fluxRight(iV,iD)*gradientBasisFunctionRight(iD);
 			}
 		}
 	}
 	return integrands;
}

/// ***************************************************
/// ***   Internal Stability Parameter functions    ***
/// ***************************************************

/// \brief Computes the fluxFunction in the stability parameter calculation used for the integrand
template<std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
LinearAlgebra::MiddleSizeMatrix ViscousTerms<DIM,NUMBER_OF_VARIABLES>::fluxStabilityParameters
(
 const StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> faceStateStruct,
 const LinearAlgebra::MiddleSizeVector &stateInternal,
 const LinearAlgebra::MiddleSizeVector &stateExternal,
 const LinearAlgebra::SmallVector<DIM> &normalInternal
 )
{
	//Compute velocity normal matrix
	LinearAlgebra::MiddleSizeMatrix stateNormal(NUMBER_OF_VARIABLES,DIM);
	LinearAlgebra::MiddleSizeVector stateDifference = stateInternal - stateExternal;
	for (std::size_t iV = 0; iV < NUMBER_OF_VARIABLES; iV++)
	{
		for (std::size_t iD = 0; iD < DIM; iD++)
		{
			stateNormal(iV,iD) = stateDifference(iV)*normalInternal(iD);
		}
	}
	/*
	std::cout << "========Computing stablity flux=========" << std::endl;
	std::cout << "stateInternal: " << stateInternal << std::endl;
	std::cout << "stateExternal: " << stateExternal << std::endl;
	std::cout << "stateNormal: " << stateNormal << std::endl;

	std::cout << "EllipticTensors: " << std::endl;
	std::cout << "11: " << faceStateStruct.getEllipticTensor()[0] << std::endl;
	std::cout << "12: " << faceStateStruct.getEllipticTensor()[1] << std::endl;
	std::cout << "21: " << faceStateStruct.getEllipticTensor()[2] << std::endl;
	std::cout << "22: " << faceStateStruct.getEllipticTensor()[3] << std::endl;

	std::cout << "result: " << 0.5*faceStateStruct.computeEllipticTensorMatrixContractionFast(stateNormal) << std::endl;
*/


	return 0.5*faceStateStruct.computeEllipticTensorMatrixContractionFast(stateNormal);
}

/// \brief Computes the integrand required for the stability parameter calculations
template<std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
LinearAlgebra::MiddleSizeMatrix ViscousTerms<DIM,NUMBER_OF_VARIABLES>::integrandStabilityRightHandSideOnFace
(
 Base::PhysicalFace<DIM> &face,
 const LinearAlgebra::MiddleSizeVector stateCoefficientsLeft,
 const LinearAlgebra::MiddleSizeVector stateCoefficientsRight,
 const Base::Side &side
 )
{

	std::size_t numberOfBasisFunctions = face.getPhysicalElement(side).getNumberOfBasisFunctions();
	std::size_t iVB;
	LinearAlgebra::MiddleSizeMatrix integrand(numberOfBasisFunctions*NUMBER_OF_VARIABLES,DIM);


	//compute stabilityFluxFunction
	LinearAlgebra::MiddleSizeMatrix stabilityFluxFunction;
		if (side == Base::Side::LEFT)
		{
			//From the given state coefficients, and coordinate in face, compute the new StateStruct for the left side (Elliptic tensor is required)
			const StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> faceStateStructLeft = instance_->computeFaceStateStruct(face, stateCoefficientsLeft, Base::Side::LEFT, 0);
			//From the given state coefficients reconstruct the state on the right side
			const LinearAlgebra::MiddleSizeVector stateRight = computeStateOnFace<DIM,NUMBER_OF_VARIABLES>(face, Base::Side::RIGHT, stateCoefficientsRight);
			stabilityFluxFunction = fluxStabilityParameters(faceStateStructLeft, faceStateStructLeft.getState(), stateRight, face.getUnitNormalVector());
		}
		else
		{
			const StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> faceStateStructRight = instance_->computeFaceStateStruct(face, stateCoefficientsRight, Base::Side::RIGHT, 0);
			const LinearAlgebra::MiddleSizeVector stateLeft = computeStateOnFace<DIM,NUMBER_OF_VARIABLES>(face, Base::Side::LEFT, stateCoefficientsLeft);
			stabilityFluxFunction = fluxStabilityParameters(faceStateStructRight, stateLeft, faceStateStructRight.getState(), face.getUnitNormalVector());
		}

	//Compute integrand
	for (std::size_t iV = 0; iV < NUMBER_OF_VARIABLES; iV++)
	{
		for (std::size_t iB = 0; iB < numberOfBasisFunctions; iB++)
		{
			iVB = face.getPhysicalElement(side).convertToSingleIndex(iB,iV);
			for (std::size_t iD = 0; iD < DIM; iD++)
			{
				integrand(iVB,iD) += stabilityFluxFunction(iV,iD)*face.basisFunction(side, iB);
			}
		}
	}
	//std::cout << "-----------------" << std::endl;
	//std::cout << "stabilityflux: " << stabilityFluxFunction << std::endl;
	//std::cout << "integrand: " << integrand << std::endl;
	
	return integrand;
}

/// \brief Computes the rhs, for the system of equations solving the stability parameters, by integrating the rhs  stability parameter integrand
template<std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
LinearAlgebra::MiddleSizeMatrix ViscousTerms<DIM,NUMBER_OF_VARIABLES>::computeRhsStabilityParameters
(
 const Base::Face *ptrFace,
 const LinearAlgebra::MiddleSizeVector stateCoefficientsLeft,
 const LinearAlgebra::MiddleSizeVector stateCoefficientsRight,
 const Base::Side &side
 )
{
	std::function<LinearAlgebra::MiddleSizeMatrix(Base::PhysicalFace<DIM>&)> integrandFunction = [=](Base::PhysicalFace<DIM> &face) -> LinearAlgebra::MiddleSizeMatrix
			{   return this->integrandStabilityRightHandSideOnFace(face, stateCoefficientsLeft, stateCoefficientsRight, side);};
	return stabilityFaceIntegrator_.integrate(ptrFace, integrandFunction);
}


/// \brief Computes the stability parameters used in the auxilliary integrand for an internal face for given stabilityParameterflux
template<std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
LinearAlgebra::MiddleSizeMatrix ViscousTerms<DIM,NUMBER_OF_VARIABLES>::computeStabilityParameters
(
 Base::PhysicalFace<DIM> &face,
 const LinearAlgebra::MiddleSizeVector stateCoefficientsLeft,
 const LinearAlgebra::MiddleSizeVector stateCoefficientsRight
 )
{
	//Datastructures
	std::size_t numberOfTestBasisFunctionsLeft = face.getPhysicalElement(Base::Side::LEFT).getNumberOfBasisFunctions();
	std::size_t numberOfTestBasisFunctionsRight = face.getPhysicalElement(Base::Side::RIGHT).getNumberOfBasisFunctions();

        //LinearAlgebra::MiddleSizeMatrix MassMatrixLeft = this->computeMassMatrixAtFace(face.getFace(),Base::Side::LEFT);
        //LinearAlgebra::MiddleSizeMatrix MassMatrixRight = this->computeMassMatrixAtFace(face.getFace(),Base::Side::RIGHT);

	//std::cout << "========COMPUTING RHSLEFT STABILITY========" << std::endl;
	//Compute stability parameters on the Left side
	LinearAlgebra::MiddleSizeMatrix rhsLeft = computeRhsStabilityParameters(face.getFace(), stateCoefficientsLeft, stateCoefficientsRight, Base::Side::LEFT);
	//std::cout << "rhsLeft: " << rhsLeft << std::endl;
	this->stabilityMassMatrix_.solve(rhsLeft);
	//LinearAlgebra::MiddleSizeMatrix result = stabilityMassMatrix_*test;
	//std::cout << "identity?" << result << std::end;
	
	//MassMatrixLeft.solve(rhsLeft);
	//std::cout << "CoefficientsLeft: " << rhsLeft << std::endl;
	//std::cout << "MassMatrix: " << stabilityMassMatrix_ << std::endl;

//	std::cout << "size matrix: " << stabilityMassMatrix_.getNRows() << std::endl;
//	std::cout << "size vector: " << rhsLeft.getNRows() << std::endl;

        //std::cout << "========COMPUTING RHSRIGHT STABILITY========" << std::endl;
	//Compute stability parameters on the Right side
	LinearAlgebra::MiddleSizeMatrix rhsRight = computeRhsStabilityParameters(face.getFace(), stateCoefficientsLeft, stateCoefficientsRight, Base::Side::RIGHT);
	//std::cout << "rhsRight: " << rhsRight << std::endl;
	//TODO: NOTE that it is quicker to determin the inverse once and then multiply it forever with the inverse.
	//std::cout << "rhsRight: " << rhsRight << std::endl;
	stabilityMassMatrix_.solve(rhsRight);
	//MassMatrixRight.solve(rhsRight);
        //std::cout << "CoefficientsRight: " << rhsRight << std::endl;
        //std::cout << "MassMatrix: " << stabilityMassMatrix_ << std::endl;

        //std::cout << "face mass matrix left: " << MassMatrixLeft << std::endl;
        //std::cout << "face mass matrix right: " << MassMatrixRight << std::endl;


	

	//Compute the average of the stability parameters
	return 0.5*(computeStateOnFace<DIM,NUMBER_OF_VARIABLES>(face, Base::Side::LEFT, rhsLeft) + computeStateOnFace<DIM,NUMBER_OF_VARIABLES>(face, Base::Side::RIGHT, rhsRight));
}

/// \brief Sets the mass matrix used in the computation of the stability parameters
template<std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
void ViscousTerms<DIM,NUMBER_OF_VARIABLES>::setStabilityMassMatrix
(
 LinearAlgebra::MiddleSizeMatrix &stabilityMassMatrix
 )
{
	stabilityMassMatrix_ = stabilityMassMatrix;
}


	/// ************************************************
	/// ***    Jacobian Exact Derivative Functions   ***
	/// ************************************************
/*  //todo: Complete this section at some point. Many functions are not yet finished. The question is, do I require the exact Jacobian or is a difference approach good enough?

	/// \brief computes the derivative of the viscous element integral flux with respect to the solutionCoefficients
	LinearAlgebra::MiddleSizeMatrix computeViscousFluxDiv(Base::PhysicalElement<DIM> &element, const std::size_t iV, const std::size_t iB, const LinearAlgebra::MiddleSizeVector &state, const LinearAlgebra::MiddleSizeVector &partialState, const LinearAlgebra::MiddleSizeMatrix &stateJacobian);

	/// \brief computes the stabilityfluxFunction derivative, used for calculating the auxilliary jacobian face integrand
	LinearAlgebra::MiddleSizeMatrix computeStabilityFluxFunctionDiv(const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensor, const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorDiv, const LinearAlgebra::MiddleSizeVector &stateLeft, const LinearAlgebra::MiddleSizeVector &stateLeftDiv, const LinearAlgebra::MiddleSizeVector &stateRight, const LinearAlgebra::MiddleSizeVector &stateRightDiv, const LinearAlgebra::SmallVector<DIM> &unitNormalLeft);

	/// This function computes the Auxilliary flux derivative for the auxilliary jacobian face integral
	LinearAlgebra::MiddleSizeMatrix computeAuxilliaryFluxDiv(const LinearAlgebra::MiddleSizeVector &stateLeft, const LinearAlgebra::MiddleSizeVector &stateLeftDiv, const LinearAlgebra::MiddleSizeVector &stateRight, const LinearAlgebra::MiddleSizeVector &stateRightDiv, const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorLeft, const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorLeftDiv, const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorRight,	const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorRightDiv, const LinearAlgebra::MiddleSizeMatrix &stateJacobianLeft, const LinearAlgebra::MiddleSizeMatrix &stateJacobianLeftDiv, const LinearAlgebra::MiddleSizeMatrix &stateJacobianRight, const LinearAlgebra::MiddleSizeMatrix &stateJacobianRightDiv, const LinearAlgebra::SmallVector<DIM> &unitNormalLeft) ;
*/

/// ************************************************
/// ***    Jacobian Element Matrix Functions     ***
/// ************************************************

/// \brief Compute the ViscousElementIntegrand contribution to the Jacobian matrix
template<std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
LinearAlgebra::MiddleSizeMatrix ViscousTerms<DIM,NUMBER_OF_VARIABLES>::integrandJacobianViscousElement
(
 Base::PhysicalElement<DIM> &element,
 const StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> &elementStateStruct
 )
{
	//Create integrand vector
	std::size_t numberOfBasisFunctions =  element.getNumOfBasisFunctions();
	LinearAlgebra::MiddleSizeMatrix integrand(NUMBER_OF_VARIABLES*numberOfBasisFunctions,NUMBER_OF_VARIABLES*numberOfBasisFunctions);
	LinearAlgebra::SmallVector<DIM> gradientBasisFunction; //Gradient function based on the number of DIMs
	LinearAlgebra::SmallVector<DIM> flux;

	//First flux function
	LinearAlgebra::MiddleSizeMatrix fluxFunction = elementStateStruct.computeEllipticTensorMatrixContractionFast(elementStateStruct.getStateJacobian());
	LinearAlgebra::MiddleSizeMatrix fluxFirst = fluxFunction;

	//compute the integrand matrix
	std::size_t iVB1, iVB2;
	for(std::size_t iV2 = 0; iV2 < NUMBER_OF_VARIABLES; iV2++) // This is the variable in the integrand vector
	{
		for (std::size_t iB2 = 0; iB2 < numberOfBasisFunctions; iB2++) // This is the basis function in the integrand vector
		{
			iVB2 = element.convertToSingleIndex(iB2,iV2); // the integrand vector index

			//Compute second state
			LinearAlgebra::MiddleSizeVector stateCoefficients = elementStateStruct.getStateCoefficients();
			stateCoefficients(iVB2) += EPSILON;
			StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> elementStateStructSecond(element,stateCoefficients,0);

			//Compute second flux
			LinearAlgebra::MiddleSizeMatrix fluxSecond = elementStateStructSecond.computeEllipticTensorMatrixContractionFast(elementStateStructSecond.getStateJacobian());

			//Compute the difference
			LinearAlgebra::MiddleSizeMatrix fluxFunctionDiv = (1.0/EPSILON)*(fluxSecond - fluxFirst);

			for (std::size_t iV1 = 0; iV1 < NUMBER_OF_VARIABLES; iV1++) // This is the variable in the coefficient vector
			{
				for (std::size_t iB1 = 0; iB1 < numberOfBasisFunctions; iB1++) // This is the basis function in the coefficient vector
				{
					gradientBasisFunction = element.basisFunctionDeriv(iB1);
					iVB1 = element.convertToSingleIndex(iB1,iV1);
					for (std::size_t iD = 0; iD < DIM; iD++)
					{
						integrand(iVB1,iVB2) += fluxFunctionDiv(iV1,iD)*gradientBasisFunction(iD);
					}
				}
			}
		}
	}
	return -integrand; //Minus sign because it is on the right hand side
}

/// *********************************************
/// ***    Jacobian Face Matrix Functions     ***
/// *********************************************

/// \brief This function computes the integrand of the Jacobian integral for the viscous part of the equations. The Jacobian is approximated by a difference approach
template<std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
LinearAlgebra::MiddleSizeMatrix ViscousTerms<DIM,NUMBER_OF_VARIABLES>::integrandJacobianViscousFace
(
 Base::PhysicalFace<DIM> &face,
 const StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> &faceStateStructLeft,
 const StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> &faceStateStructRight,
 const Base::Side &elementSide,
 const Base::Side &derivativeSide)
{
	//Data structures for left and right integrand
	std::size_t numberOfBasisFunctionsElementSide = face.getPhysicalElement(elementSide).getNumberOfBasisFunctions();
	std::size_t numberOfBasisFunctionsDerivativeSide = face.getPhysicalElement(derivativeSide).getNumOfBasisFunctions();
	LinearAlgebra::MiddleSizeMatrix integrand(NUMBER_OF_VARIABLES*numberOfBasisFunctionsElementSide,NUMBER_OF_VARIABLES*numberOfBasisFunctionsDerivativeSide);
	LinearAlgebra::SmallVector<DIM> unitNormalLeft = face.getUnitNormalVector();
	LinearAlgebra::SmallVector<DIM> gradientBasisFunction;

	//Compute First Flux
	//todo: move this stateNormal function to StatecoefficientsFunctions or let BLAS do it
	LinearAlgebra::MiddleSizeMatrix fluxFirst;
	LinearAlgebra::MiddleSizeMatrix stateNormalFirst(NUMBER_OF_VARIABLES,DIM);
	LinearAlgebra::MiddleSizeVector stateDifferenceFirst = faceStateStructLeft.getState() - faceStateStructRight.getState();
	for (std::size_t iV = 0; iV < NUMBER_OF_VARIABLES; iV++)
	{
		for (std::size_t iD = 0; iD < DIM; iD++)
		{
			stateNormalFirst(iV,iD) = 0.5*stateDifferenceFirst(iV)*unitNormalLeft(iD);
		}
	}

	if (elementSide == Base::Side::LEFT)
	{
		fluxFirst = faceStateStructLeft.computeEllipticTensorMatrixContractionFast(stateNormalFirst);
	}
	else
	{
		fluxFirst = faceStateStructRight.computeEllipticTensorMatrixContractionFast(stateNormalFirst);
	}

	//Start computing the Jacobian
	std::size_t iVB1, iVB2;
	for (std::size_t iV2 = 0; iV2 < NUMBER_OF_VARIABLES; iV2++)
	{
		for (std::size_t iB2 = 0; iB2 < numberOfBasisFunctionsDerivativeSide; iB2++)
		{
			iVB2 = face.getPhysicalElement(derivativeSide).convertToSingleIndex(iB2,iV2);

			//Compute new state by perturbing one of the coefficients and compute the second flux
			LinearAlgebra::MiddleSizeVector solutionCoefficientsPerturb;
			LinearAlgebra::MiddleSizeMatrix fluxSecond;
			LinearAlgebra::MiddleSizeMatrix stateNormalSecond(NUMBER_OF_VARIABLES,DIM);
			LinearAlgebra::MiddleSizeVector stateDifferenceSecond;
			if (derivativeSide == Base::Side::LEFT)
			{
				//Compute new state
				solutionCoefficientsPerturb = faceStateStructLeft.getStateCoefficients();
				solutionCoefficientsPerturb(iVB2) += EPSILON;
				StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> faceStateStructLeftNew(face, solutionCoefficientsPerturb, Base::Side::LEFT,0);

				//New state normal
				//todo: move this to a different file or let BLAS do it
				stateDifferenceSecond = faceStateStructLeftNew.getState() - faceStateStructRight.getState();
				for (std::size_t iV = 0; iV < NUMBER_OF_VARIABLES; iV++)
				{
					for (std::size_t iD = 0; iD < DIM; iD++)
					{
						stateNormalSecond(iV,iD) = 0.5*stateDifferenceSecond(iV)*unitNormalLeft(iD);
					}
				}

				//New fluxFunction
				if (elementSide == Base::Side::LEFT)
				{
					fluxSecond = faceStateStructLeftNew.computeEllipticTensorMatrixContractionFast(stateNormalSecond);
				}
				else
				{
					fluxSecond = faceStateStructRight.computeEllipticTensorMatrixContractionFast(stateNormalSecond);
				}
			}
			else
			{
				//Compute new state
				solutionCoefficientsPerturb = faceStateStructRight.getStateCoefficients();
				solutionCoefficientsPerturb(iVB2) += EPSILON;
				StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> faceStateStructRightNew(face, solutionCoefficientsPerturb, Base::Side::RIGHT,0);

				//New state normal
				//todo: move this to a different file
				stateDifferenceSecond = faceStateStructLeft.getState() - faceStateStructRightNew.getState();
				for (std::size_t iV = 0; iV < NUMBER_OF_VARIABLES; iV++)
				{
					for (std::size_t iD = 0; iD < DIM; iD++)
					{
						stateNormalSecond(iV,iD) = 0.5*stateDifferenceSecond(iV)*unitNormalLeft(iD);
					}
				}

				//New fluxFunction
				if (elementSide == Base::Side::LEFT)
				{
					fluxSecond = faceStateStructLeft.computeEllipticTensorMatrixContractionFast(stateNormalSecond);
				}
				else
				{
					fluxSecond = faceStateStructRightNew.computeEllipticTensorMatrixContractionFast(stateNormalSecond);
				}
			}

			//Compute fluxdifference
			LinearAlgebra::MiddleSizeMatrix fluxDifference = (1.0/EPSILON)*(fluxSecond - fluxFirst);

			for (std::size_t iB1 = 0; iB1 < numberOfBasisFunctionsDerivativeSide; iB1++) // for all basis functions
			{
				gradientBasisFunction = face.basisFunctionDeriv(elementSide,iB1);
				for (std::size_t iD = 0; iD < DIM; iD++) // for all DIMs
				{
					for (std::size_t iV1 = 0; iV1 < NUMBER_OF_VARIABLES; iV1++) // for all equations
					{
						iVB1 = face.getPhysicalElement(elementSide).convertToSingleIndex(iB1,iV1);
						integrand(iVB1,iVB2) += fluxDifference(iV1,iD)*gradientBasisFunction(iD);
					}
				}
			}
		}
	}

	return integrand;
}

/// \brief This function computes the integrand of the Jacobian integral for the Auxilliary part of the equations. The Jacobian is approximated by a difference approach
template<std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
LinearAlgebra::MiddleSizeMatrix ViscousTerms<DIM,NUMBER_OF_VARIABLES>::integrandJacobianAuxilliaryFace
(
 Base::PhysicalFace<DIM> &face,
 const StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> &faceStateStructLeft,
 const StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> &faceStateStructRight,
 const Base::Side &elementSide,
 const Base::Side &derivativeSide
 )
{
	//Data structures for left and right integrand
	std::size_t numberOfBasisFunctionsElementSide = face.getPhysicalElement(elementSide).getNumberOfBasisFunctions();
	std::size_t numberOfBasisFunctionsDerivativeSide = face.getPhysicalElement(derivativeSide).getNumOfBasisFunctions();
	LinearAlgebra::MiddleSizeMatrix integrand(NUMBER_OF_VARIABLES*numberOfBasisFunctionsElementSide,NUMBER_OF_VARIABLES*numberOfBasisFunctionsDerivativeSide);
	LinearAlgebra::SmallVector<DIM> unitNormalLeft = face.getUnitNormalVector();

	//Compute First Flux
	LinearAlgebra::MiddleSizeMatrix fluxFunctionFirst = fluxAuxilliary(face, faceStateStructLeft, faceStateStructRight);

	//Flux sign, depending on the elementSide
	double sign;
	if (elementSide == Base::Side::LEFT)
	{
		sign = 1.0;
	}
	else
	{
		sign = -1.0;
	}

	std::size_t iVB1, iVB2; // iVB1 stands for the index in the RHS vector, and iVB2 stands for the perturbed index in the coefficient vector
	for (std::size_t iV2 = 0; iV2 < NUMBER_OF_VARIABLES; iV2++)
	{
		for (std::size_t iB2 = 0; iB2 < numberOfBasisFunctionsDerivativeSide; iB2++)
		{
			iVB2 = face.getPhysicalElement(derivativeSide).convertToSingleIndex(iB2,iV2);

			//compute new state
			LinearAlgebra::MiddleSizeVector solutionCoefficientsPerturb;
			LinearAlgebra::MiddleSizeMatrix fluxFunctionSecond;
			if (derivativeSide == Base::Side::LEFT)
			{
				solutionCoefficientsPerturb = faceStateStructLeft.getStateCoefficients();
				solutionCoefficientsPerturb(iVB2) += EPSILON;
				StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> faceStateStructLeftNew(face, solutionCoefficientsPerturb, Base::Side::LEFT,0);
				//Compute second Flux
				fluxFunctionSecond = fluxAuxilliary(face, faceStateStructLeftNew, faceStateStructRight);
			}
			else
			{
				solutionCoefficientsPerturb = faceStateStructRight.getStateCoefficients();
				solutionCoefficientsPerturb(iVB2) += EPSILON;
				StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> faceStateStructRightNew(face, solutionCoefficientsPerturb, Base::Side::RIGHT,0);
				//Compute second Flux
				fluxFunctionSecond = fluxAuxilliary(face, faceStateStructLeft, faceStateStructRightNew);
			}

			//Compute fluxdifference
			LinearAlgebra::MiddleSizeMatrix fluxFunctionDifference = (1.0/EPSILON)*(fluxFunctionSecond - fluxFunctionFirst);

			//Compute partial integrand second
			//todo: this is a simple matrix vector multiplciation: let BLAS do it
			LinearAlgebra::MiddleSizeVector partialIntegrandDifference(NUMBER_OF_VARIABLES);
			for (std::size_t iV = 0; iV < NUMBER_OF_VARIABLES; iV++)
			{
				for (std::size_t iD = 0; iD < DIM; iD++)
				{
					partialIntegrandDifference(iV) += fluxFunctionDifference(iV,iD)*unitNormalLeft(iD);
				}
			}

			//Compute the integrand
			for (std::size_t iB1 = 0; iB1 < numberOfBasisFunctionsElementSide; iB1++)
			{
				for (std::size_t iV1 = 0; iV1 < NUMBER_OF_VARIABLES; iV1++)
				{
					iVB1 = face.getPhysicalElement(elementSide).convertToSingleIndex(iB1,iV1);
					integrand(iVB1,iVB2) += sign*face.basisFunction(elementSide, iB1)*partialIntegrandDifference(iV1);
				}
			}
		}
	}

	return integrand;
}

