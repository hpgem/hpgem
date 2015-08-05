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
#include "Viscous.h"
#include <tuple>
#include <chrono> //todo: remove this, not required

Viscous::Viscous(CompressibleNavierStokes& instance) : instance_(instance), PrInv_(1/0.71), cp_(1000)
{
}


/// ****************************************
/// ***      Constitutive functions      ***
/// ****************************************

double Viscous::computeTemperature(const LinearAlgebra::MiddleSizeVector &state, const double pressure)
{
	return pressure*instance_.gamma_/((instance_.gamma_ - 1)*cp_*state(0)); //T = p/(R*rho);
}


double Viscous::computeViscosity(double temperature)
{
	double temperatureS =  110;
	double temperatureRef = 288.16;
	double muRef = 0.000017894;

	double temp = temperature/temperatureRef;
	return muRef*(temp)*sqrt(temp)*(temperatureRef + temperatureS)/(temperature + temperatureS);
}


/// *******************************************
/// ***      Elliptic Tensor functions      ***
/// *******************************************

std::vector<LinearAlgebra::MiddleSizeMatrix> Viscous::computeATensor(const LinearAlgebra::MiddleSizeVector &partialState, const double viscosity)
{
	//todo: note that the kinetic velocity must be computed for a wide range of problems. Fix this.
	//todo: remove the kinetic velocity in this computation.
	//todo: Check if this also works correctly in 3D

	std::vector<LinearAlgebra::MiddleSizeMatrix> ATensor(instance_.DIM_*instance_.DIM_);
	double velocityNormSquared = 0.0;
	double thermalFactor = instance_.gamma_*viscosity/PrInv_;

	double pos1;
	double pos2;

	double factor43 = 4.0/3.0;
	double vis23 = 2.0/3.0*viscosity;
	double inverseRho = 1.0/partialState(0);

	for (std::size_t iD = 0; iD < instance_.DIM_; iD++)
	{
		velocityNormSquared += partialState(iD+1)*partialState(iD+1);
	}

	LinearAlgebra::MiddleSizeMatrix APartial1(instance_.DIM_+2,instance_.DIM_+2);
	LinearAlgebra::MiddleSizeMatrix APartial2(instance_.DIM_+2,instance_.DIM_+2);

	//A11 A22 en A33: For documentation see the full matrix in Klaij et al. 2006
	for (std::size_t iD = 0; iD < instance_.DIM_; iD++)
	{
		//Reset matrix
		APartial1 *= 0.0;

		//Tensor index
		pos1 = (instance_.DIM_)*iD + iD;

		// (1) viscosity contributions
		for (std::size_t iD2 = 0; iD2 < instance_.DIM_; iD2++)
		{
			APartial1(iD2+1,iD2+1) = viscosity;
			APartial1(iD2+1,0) = -viscosity*partialState(iD2+1);
			APartial1(instance_.DIM_+1,iD2+1) = viscosity;
		}

		// Multiply by correct value of the dominant component
		APartial1(iD+1,0) *= factor43;
		APartial1(iD+1,iD+1) *= factor43;
		APartial1(instance_.DIM_+1,iD+1) *= factor43;

		// (2) temperature contribution
		for (std::size_t iD2 = 0; iD2 < instance_.DIM_; iD2++)
		{
			APartial1(instance_.DIM_+1,iD2+1) += -thermalFactor;
		}
		APartial1(instance_.DIM_+1,instance_.DIM_+1) = thermalFactor;

		// Complete energy part by multiplying by velocity.
		for (std::size_t iD2 = 0; iD2 < instance_.DIM_; iD2++)
		{
			APartial1(instance_.DIM_+1,iD2+1) *= partialState(iD2+1);
		}

		APartial1(instance_.DIM_+1,0) = -(1.0/3.0)*partialState(iD+1)*partialState(iD+1) - viscosity*velocityNormSquared - thermalFactor*(partialState(instance_.DIM_+1) - velocityNormSquared);

		//Divide by rho
		APartial1 *= inverseRho;
		ATensor[pos1] = APartial1;

	}

	//A12 A13 A23 ? Note: A12 --> A(iD1)(iD2)
	for (std::size_t iD1 = 0; iD1 < instance_.DIM_ - 1; iD1++)
	{
		for (std::size_t iD2 = iD1 + 1; iD2 < instance_.DIM_; iD2++)
		{
			//Reset matrix
			APartial1 *= 0.0;
			APartial2 *= 0.0;

			//Tensor index
			pos1 = (instance_.DIM_)*iD1 + iD2;
			pos2 = (instance_.DIM_)*iD2 + iD1;

			// viscosity contributions for A(iD1)(iD2)
			APartial1(iD1+1,0) = vis23*partialState(iD2+1);
			APartial1(iD2+1,0) = -viscosity*partialState(iD1+1);

			APartial1(iD1+1,iD2+1) = -vis23;
			APartial1(iD2+1,iD1+1) =  viscosity;

			APartial1(instance_.DIM_ + 1, 0) = -1.0/3.0*viscosity*partialState(iD1+1)*partialState(iD2+1);
			APartial1(instance_.DIM_ + 1, iD2+1) = -vis23*partialState(iD1+1);
			APartial1(instance_.DIM_ + 1, iD1+1) = viscosity*partialState(iD2+1);

			// viscosity contributions for A(iD2)(iD1)
			APartial2(iD1+1,0) = -viscosity*partialState(iD2+1);
			APartial2(iD2+1,0) = vis23*partialState(iD1+1);

			APartial2(iD1+1,iD2+1) = viscosity;
			APartial2(iD2+1,iD1+1) =  -vis23;

			APartial2(instance_.DIM_ + 1, 0) = -1.0/3.0*viscosity*partialState(iD1+1)*partialState(iD2+1);
			APartial2(instance_.DIM_ + 1, iD2+1) = viscosity*partialState(iD1+1);
			APartial2(instance_.DIM_ + 1, iD1+1) = -vis23*partialState(iD2+1);

			//Divide by rho
			APartial1 *= inverseRho;
			APartial2 *= inverseRho;

			//Assign matrices to the tensor vector
			ATensor[pos1] = APartial1;
			ATensor[pos2] = APartial2;

		}
	}
	return ATensor;
}


//todo: fix this part, faster integration
//A_ikrs = A_(iV)(iD)(iVm)(iDm)
//This is a naive implementation
LinearAlgebra::MiddleSizeMatrix Viscous::computeATensorMatrixContraction(
		const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensor,
		const LinearAlgebra::MiddleSizeMatrix &matrix)
{
	LinearAlgebra::MiddleSizeMatrix result(instance_.numOfVariables_,instance_.DIM_);
	double pos;

	for (std::size_t iD = 0; iD < instance_.DIM_; iD++)
	{
		for (std::size_t iV = 0; iV < instance_.numOfVariables_; iV++)
		{
			for (std::size_t iDm = 0; iDm < instance_.DIM_; iDm++)
			{
				for (std::size_t iVm = 0; iVm < instance_.numOfVariables_; iVm++)
				{
					pos = (instance_.DIM_)*iD + iDm;
					result(iV,iD) += ATensor[pos](iV,iVm)*matrix(iVm,iDm);
				}
			}
		}
	}

	return result;
}

//A_ikrs = A_(iV)(iD)(iVm)(iDm)
//matrix_rs = matrix(iVm,iDm)
LinearAlgebra::MiddleSizeMatrix Viscous::computeATensorMatrixContractionFast(
		const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensor,
		const LinearAlgebra::MiddleSizeMatrix &matrix)
{
	LinearAlgebra::MiddleSizeMatrix result(instance_.numOfVariables_,instance_.DIM_);
	double pos;

	//A11 A22 and A33
	for (std::size_t iD = 0; iD < instance_.DIM_; iD++)
	{
		std::size_t iDm = iD;
		pos = (instance_.DIM_)*iD + iDm;
		//Velocity part
		for (std::size_t it1 = 0; it1 < instance_.DIM_; it1++)
		{
			result(it1+1,iD) += ATensor[pos](it1+1,0)*matrix(0,iDm);
			result(it1+1,iD) += ATensor[pos](it1+1,it1+1)*matrix(it1+1,iDm);
		}
		//Energy part
		for (std::size_t iVm = 0; iVm < instance_.numOfVariables_; iVm++)
		{
			result(instance_.DIM_+1,iD) += ATensor[pos](instance_.DIM_ + 1,iVm)*matrix(iVm,iDm);
		}
	}

	//Other
	for (std::size_t iD = 0; iD < instance_.DIM_; iD++)
	{
		for ( std::size_t iDm = iD + 1; iDm < instance_.DIM_; iDm++)
		{

			//A12 A13 and A23
			pos = (instance_.DIM_)*iD + iDm;
			//Velocity part
			result(iD+1,iD) += ATensor[pos](iD+1,0)*matrix(0,iDm);
			result(iD+1,iD) += ATensor[pos](iD+1,iDm+1)*matrix(iDm+1,iDm);
			result(iDm+1,iD) += ATensor[pos](iDm+1,0)*matrix(0,iDm);
			result(iDm+1,iD) += ATensor[pos](iDm+1,iD+1)*matrix(iD+1,iDm);
			//Energy part
			for (std::size_t iVm = 0; iVm < instance_.numOfVariables_; iVm++)
			{
				result(instance_.DIM_+1,iD) += ATensor[pos](instance_.DIM_ + 1,iVm)*matrix(iVm,iDm);
			}

			//A21 A31 and A32
			//Note: Same structure as above, however iD and iDm are switched
			pos = (instance_.DIM_)*iDm + iD;
			//Velocity part
			result(iDm+1,iDm) += ATensor[pos](iDm+1,0)*matrix(0,iD);
			result(iDm+1,iDm) += ATensor[pos](iDm+1,iD+1)*matrix(iD+1,iD);
			result(iD+1,iDm) += ATensor[pos](iD+1,0)*matrix(0,iD);
			result(iD+1,iDm) += ATensor[pos](iD+1,iDm+1)*matrix(iDm+1,iD);
			//Energy part
			for (std::size_t iVm = 0; iVm < instance_.numOfVariables_; iVm++)
			{
				result(instance_.DIM_+1,iDm) += ATensor[pos](instance_.DIM_ + 1,iVm)*matrix(iVm,iD);
			}

		}
	}

	return result;
}


/// *****************************************
/// ***   Element integration functions   ***
/// *****************************************

LinearAlgebra::MiddleSizeVector Viscous::integrandAtElement(
		Base::PhysicalElement<DIM> &element,
		const LinearAlgebra::MiddleSizeVector &state,
		const LinearAlgebra::MiddleSizeMatrix &stateJacobian,
		const double pressure,
		const LinearAlgebra::MiddleSizeVector &partialState)
{
	// Get the number of basis functions in an element.
	std::size_t numberOfBasisFunctions = element.getNumOfBasisFunctions();

	// Create data structures for calculating the integrand
	LinearAlgebra::MiddleSizeVector integrand(instance_.numOfVariables_ * numberOfBasisFunctions);
	LinearAlgebra::SmallVector<DIM> gradientBasisFunction;
	std::size_t iVB;


	// Compute A tensor
	double temperature = computeTemperature(state, pressure);
	double viscosity = computeViscosity(temperature);
	std::vector<LinearAlgebra::MiddleSizeMatrix> ATensor = computeATensor(partialState, viscosity);

	// Compute flux (A_ikrs*Ur,s)
	LinearAlgebra::MiddleSizeMatrix fluxFunction = computeATensorMatrixContractionFast(ATensor, stateJacobian);

	for (std::size_t iB = 0; iB < numberOfBasisFunctions; iB++) // for all basis functions
	{
		gradientBasisFunction = element.basisFunctionDeriv(iB);
		for (std::size_t iD = 0; iD < instance_.DIM_; iD++) // for all dimensions
		{
			for (std::size_t iV = 0; iV < instance_.numOfVariables_; iV++) // for all variables
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

/*
void Viscous::computeStabilityFluxFunction(
		const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorExternal,
		const LinearAlgebra::MiddleSizeVector &stateInternal,
		const LinearAlgebra::MiddleSizeVector &stateBoundary,
		const LinearAlgebra::SmallVector<DIM> &normalInternal)
{
	//Compute velocity normal matrix
	LinearAlgebra::MiddleSizeMatrix velocityNormal(instance_.DIM_+2,instance_.DIM_);
	LinearAlgebra::MiddleSizeVector stateDifference = stateInternal - stateBoundary;
	for (std::size_t iV = 0; iV < instance_.numOfVariables_; iV++)
	{
		for (std::size_t iD = 0; iD < instance_.DIM_; iD++)
		{
			velocityNormal(iV,iD) = 0.5*stateDifference(iV)*normalInternal(iD);
		}
	}

	//Compute fluxFunction
	stabilityFluxFunctionExternal_ = computeATensorMatrixContraction(ATensorExternal, velocityNormal);
}

//todo: check if the integrand is correctly calculated, see the table in fidkowski
LinearAlgebra::MiddleSizeMatrix Viscous::computeStabilityParameters(
		Base::PhysicalFace<DIM> &face,
		const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorExternal,
		const LinearAlgebra::MiddleSizeVector &stateInternal,
		const LinearAlgebra::MiddleSizeVector &stateExternal,
		const LinearAlgebra::SmallVector<DIM> &normalInternal)
{
	std::size_t numberOfTestBasisFunctionsInternal = face.getPhysicalElement(Base::Side::LEFT).getNumOfBasisFunctions();

	//Create datastructures
	LinearAlgebra::MiddleSizeMatrix stabilityParameters(instance_.numOfVariables_,instance_.DIM_);
	LinearAlgebra::MiddleSizeVector stabilityParametersExternal(instance_.numOfVariables_);
	LinearAlgebra::MiddleSizeVector rhsInternal(numberOfTestBasisFunctionsInternal);
	LinearAlgebra::MiddleSizeVector stabilityParameterCoefficientsInternal;

	//Compute the flux functions as function of iV and iD
	computeStabilityFluxFunction(ATensorExternal, stateInternal, stateExternal, normalInternal);

	for (std::size_t iD = 0; iD < instance_.DIM_; iD++)
	{
			//Integrate the correct fluxfunction
			rhsInternal = computeRhs(face.getFace(), Base::Side::LEFT, stabilityFluxFunctionExternal_, normalInternal, iD);

			//Solve the stabilityParameter coefficients
			stabilityMassMatrix_.solve(rhsInternal);


			//Reconstruct the stabilityparameter coefficients for a given iD
			stabilityParametersExternal = instance_.computeStateOnFace(face, Base::Side::LEFT, rhsInternal);

			//Put result into matrix
			for (std::size_t iV = 0; iV < instance_.numOfVariables_; iV++)
			{
				stabilityParameters(iV,iD) = stabilityParametersExternal(iV);
			}
	}

	return stabilityParameters;
}

//todo: check if the integrand is correctly calculated, see the table in fidkowski
LinearAlgebra::MiddleSizeMatrix Viscous::computeAuxilliaryFlux(
		Base::PhysicalFace<DIM> &face,
		const double temperatureExternal,
		const LinearAlgebra::MiddleSizeVector &stateInternal,
		const LinearAlgebra::MiddleSizeVector &stateExternal,
		const LinearAlgebra::MiddleSizeVector &partialStateExternal,
		const LinearAlgebra::MiddleSizeMatrix &stateJacobianInternal,
		const LinearAlgebra::SmallVector<DIM>   &normalInternal)
{
	double eta = 3.0; //stability value hard coded
	LinearAlgebra::MiddleSizeMatrix fluxBoundary;
	LinearAlgebra::MiddleSizeMatrix stabilityFluxBoundary;
	LinearAlgebra::MiddleSizeMatrix AuxilliaryFlux;

	//Compute A and A contraction with velocitynormal matrix
	double viscosityBoundary = computeViscosity(temperatureExternal);
	ATensorExternal_ = computeATensor(partialStateExternal, viscosityBoundary);
	fluxBoundary  = computeATensorMatrixContraction(ATensorExternal_, stateJacobianInternal);

	stabilityFluxBoundary = computeStabilityParameters(face, ATensorExternal_, stateInternal, stateExternal, normalInternal);

	AuxilliaryFlux = fluxBoundary - eta*stabilityFluxBoundary;

	return AuxilliaryFlux;
}

LinearAlgebra::MiddleSizeVector Viscous::integrandAuxilliaryAtFace(
		Base::PhysicalFace<DIM> &face,
		const double temperatureExternal,
		const LinearAlgebra::MiddleSizeVector &stateInternal,
		const LinearAlgebra::MiddleSizeVector &stateExternal,
		const LinearAlgebra::MiddleSizeVector &partialStateExternal,
		const LinearAlgebra::SmallVector<DIM> &normalInternal,
		const LinearAlgebra::MiddleSizeMatrix &stateJacobianInternal)
{
	LinearAlgebra::MiddleSizeMatrix fluxFunction = computeAuxilliaryFlux(face, temperatureExternal, stateInternal, stateExternal, partialStateExternal, stateJacobianInternal, normalInternal);

	//Compute integrand
	std::size_t numOfTestBasisFunctions = face.getPhysicalElement(Base::Side::LEFT).getNumOfBasisFunctions();
	LinearAlgebra::MiddleSizeVector integrand(instance_.numOfVariables_*numOfTestBasisFunctions);
	LinearAlgebra::MiddleSizeVector partialIntegrand(instance_.numOfVariables_);
	std::size_t iVB;

	for (std::size_t iV = 0; iV < instance_.numOfVariables_; iV++)
	{
		for (std::size_t iD = 0; iD < instance_.DIM_; iD++)
		{
			partialIntegrand(iV) += fluxFunction(iV,iD)*normalInternal(iD);
		}
	}

	for (std::size_t iB = 0; iB < numOfTestBasisFunctions; iB++)
	{
		for (std::size_t iV = 0; iV < instance_.numOfVariables_; iV++)
		{
			iVB = face.getPhysicalElement(Base::Side::LEFT).convertToSingleIndex(iB,iV);
			integrand(iVB) += face.basisFunction(Base::Side::LEFT, iB)*partialIntegrand(iV);
		}
	}

	return integrand; //no minus sign because integral is on rhs
}
*/


/// **************************************************
/// ***    Internal face integration functions     ***
/// **************************************************

//todo: stateNormal needs to be computed in the main integrand function - required by two fluxes, also the contraction
LinearAlgebra::MiddleSizeMatrix Viscous::computeStabilityFluxFunction(
		const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensor,
		const LinearAlgebra::MiddleSizeVector &stateInternal,
		const LinearAlgebra::MiddleSizeVector &stateExternal,
		const LinearAlgebra::SmallVector<DIM> &normalInternal)
{
	//Compute velocity normal matrix
	LinearAlgebra::MiddleSizeMatrix stateNormal(instance_.DIM_+2,instance_.DIM_);
	LinearAlgebra::MiddleSizeVector stateDifference = stateInternal - stateExternal;
	for (std::size_t iV = 0; iV < instance_.numOfVariables_; iV++)
	{
		for (std::size_t iD = 0; iD < instance_.DIM_; iD++)
		{
			stateNormal(iV,iD) = stateDifference(iV)*normalInternal(iD);
		}
	}

	//note: this is already computed in the viscous part
	return computeATensorMatrixContractionFast(ATensor, stateNormal);
}

//todo: rewrite such that it reduces the vector allocations by two
LinearAlgebra::MiddleSizeMatrix Viscous::computeStabilityParameters(
		Base::PhysicalFace<DIM> &face,
		const LinearAlgebra::MiddleSizeVector &stateLeft,
		const LinearAlgebra::MiddleSizeVector &stateRight,
		const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorLeft,
		const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorRight,
		const LinearAlgebra::SmallVector<DIM> &unitNormalLeft)
{

	//Datastructures
	LinearAlgebra::MiddleSizeMatrix stabilityParametersAverage(instance_.numOfVariables_,instance_.DIM_);
	std::size_t numberOfTestBasisFunctionsLeft = face.getPhysicalElement(Base::Side::LEFT).getNumOfBasisFunctions();
	std::size_t numberOfTestBasisFunctionsRight = face.getPhysicalElement(Base::Side::RIGHT).getNumOfBasisFunctions();
	LinearAlgebra::MiddleSizeVector rhsLeft(numberOfTestBasisFunctionsLeft);
	LinearAlgebra::MiddleSizeVector rhsRight(numberOfTestBasisFunctionsRight);
	LinearAlgebra::MiddleSizeVector resultLeft(instance_.numOfVariables_);
	LinearAlgebra::MiddleSizeVector resultRight(instance_.numOfVariables_);

	//Compute the flux functions as function of iV and iD
	//note: this is also computed in the viscous integrand: should not have to be computed twice.
	LinearAlgebra::MiddleSizeMatrix fluxLeft  = computeStabilityFluxFunction(ATensorLeft, stateLeft, stateRight, unitNormalLeft);
	LinearAlgebra::MiddleSizeMatrix fluxRight = computeStabilityFluxFunction(ATensorRight, stateLeft, stateRight, unitNormalLeft);

	//Compute the average of the stability parameters
	for (std::size_t iD = 0; iD < instance_.DIM_; iD++)
	{
			//Left
			rhsLeft = 0.5*computeRhs(face.getFace(), Base::Side::LEFT, fluxLeft, iD);
			stabilityMassMatrix_.solve(rhsLeft);
			resultLeft = instance_.computeStateOnFace(face, Base::Side::LEFT, rhsLeft);

			//Right
			rhsRight = 0.5*computeRhs(face.getFace(), Base::Side::RIGHT, fluxRight, iD);
			stabilityMassMatrix_.solve(rhsRight);
			resultRight = instance_.computeStateOnFace(face, Base::Side::RIGHT, rhsRight);

			//Combine
			//Put result into matrix
			for (std::size_t iV = 0; iV < instance_.numOfVariables_; iV++)
			{
				stabilityParametersAverage(iV,iD) = 0.5*(resultLeft(iV) + resultRight(iV));
			}
	}

	return stabilityParametersAverage;
}

LinearAlgebra::MiddleSizeMatrix Viscous::computeAuxilliaryFlux(
		Base::PhysicalFace<DIM> &face,
		const LinearAlgebra::MiddleSizeVector &stateInternal,
		const LinearAlgebra::MiddleSizeVector &stateExternal,
		const LinearAlgebra::MiddleSizeMatrix &stateJacobianLeft,
		const LinearAlgebra::MiddleSizeMatrix &stateJacobianRight,
		const LinearAlgebra::SmallVector<DIM> &unitNormalInternal,
		const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorLeft,
		const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorRight)
{

	double eta = 3.0; //stability value hard coded

	LinearAlgebra::MiddleSizeMatrix fluxLeft  = computeATensorMatrixContractionFast(ATensorLeft, stateJacobianLeft);
	LinearAlgebra::MiddleSizeMatrix fluxRight = computeATensorMatrixContractionFast(ATensorRight, stateJacobianRight);
	LinearAlgebra::MiddleSizeMatrix stabilityParameters = computeStabilityParameters(face, stateInternal, stateExternal, ATensorLeft, ATensorRight, unitNormalInternal);

	return (0.5*(fluxLeft + fluxRight) - eta*stabilityParameters);
}

std::pair<LinearAlgebra::MiddleSizeVector,LinearAlgebra::MiddleSizeVector> Viscous::integrandsAuxilliaryAtFace(
		Base::PhysicalFace<DIM> &face,
		const LinearAlgebra::MiddleSizeVector &stateLeft,
		const LinearAlgebra::MiddleSizeVector &stateRight,
		const LinearAlgebra::MiddleSizeMatrix &stateJacobianLeft,
		const LinearAlgebra::MiddleSizeMatrix &stateJacobianRight,
		const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorLeft,
		const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorRight)
{

	//Data structures for left and right integrand
	std::size_t numOfTestBasisFunctionsLeft = face.getPhysicalElement(Base::Side::LEFT).getNumOfBasisFunctions();
	std::size_t numOfTestBasisFunctionsRight = face.getPhysicalElement(Base::Side::RIGHT).getNumOfBasisFunctions();
	std::pair<LinearAlgebra::MiddleSizeVector,LinearAlgebra::MiddleSizeVector> integrands(
			std::piecewise_construct,
			std::forward_as_tuple(instance_.numOfVariables_*numOfTestBasisFunctionsLeft),
			std::forward_as_tuple(instance_.numOfVariables_*numOfTestBasisFunctionsRight));
	LinearAlgebra::SmallVector<DIM> unitNormalLeft = face.getUnitNormalVector();
	std::size_t iVB;

	//Compute Left flux
	//note: left flux is the same as the right flux
	LinearAlgebra::MiddleSizeMatrix fluxFunctionLeft = computeAuxilliaryFlux(face, stateLeft, stateRight, stateJacobianLeft, stateJacobianRight, unitNormalLeft, ATensorLeft, ATensorRight);

	//Compute Left integrand
	LinearAlgebra::MiddleSizeVector partialIntegrand(instance_.numOfVariables_);
	for (std::size_t iV = 0; iV < instance_.numOfVariables_; iV++)
	{
		for (std::size_t iD = 0; iD < instance_.DIM_; iD++)
		{
			partialIntegrand(iV) += fluxFunctionLeft(iV,iD)*unitNormalLeft(iD);
		}
	}

	for (std::size_t iB = 0; iB < numOfTestBasisFunctionsLeft; iB++)
	{
		for (std::size_t iV = 0; iV < instance_.numOfVariables_; iV++)
		{
			iVB = face.getPhysicalElement(Base::Side::LEFT).convertToSingleIndex(iB,iV);
			integrands.first(iVB) += face.basisFunction(Base::Side::LEFT, iB)*partialIntegrand(iV);
			integrands.second(iVB) += -face.basisFunction(Base::Side::RIGHT, iB)*partialIntegrand(iV);
		}
	}

	return integrands;
}

std::pair<LinearAlgebra::MiddleSizeVector,LinearAlgebra::MiddleSizeVector> Viscous::integrandsViscousAtFace(
		Base::PhysicalFace<DIM> &face,
		const LinearAlgebra::MiddleSizeVector &stateLeft,
		const LinearAlgebra::MiddleSizeVector &stateRight,
		const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorLeft,
		const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorRight
		)
{

	std::size_t numOfTestBasisFunctionsLeft = face.getPhysicalElement(Base::Side::LEFT).getNumOfBasisFunctions();
	std::size_t numOfTestBasisFunctionsRight = face.getPhysicalElement(Base::Side::RIGHT).getNumOfBasisFunctions();
	std::pair<LinearAlgebra::MiddleSizeVector,LinearAlgebra::MiddleSizeVector> integrands(
			std::piecewise_construct,
			std::forward_as_tuple(instance_.numOfVariables_*numOfTestBasisFunctionsLeft),
			std::forward_as_tuple(instance_.numOfVariables_*numOfTestBasisFunctionsRight));
	LinearAlgebra::SmallVector<DIM> unitNormalLeft = face.getUnitNormalVector();
	LinearAlgebra::SmallVector<DIM> gradientBasisFunctionLeft;
	LinearAlgebra::SmallVector<DIM> gradientBasisFunctionRight;
	std::size_t iVB;

	//Compute state times normal matrix for the Left face, this is the same as for the Right sight times a minus
	LinearAlgebra::MiddleSizeMatrix stateNormalLeft(instance_.numOfVariables_,instance_.DIM_);
	LinearAlgebra::MiddleSizeVector stateDifferenceLeft;
	stateDifferenceLeft = stateLeft - stateRight;
	for (std::size_t iV = 0; iV < instance_.numOfVariables_; iV++)
	{
		for (std::size_t iD = 0; iD < instance_.DIM_; iD++)
		{
			stateNormalLeft(iV,iD) = 0.5*stateDifferenceLeft(iV)*unitNormalLeft(iD);
		}
	}

	//Compute fluxLeft
	LinearAlgebra::MiddleSizeMatrix fluxLeft = computeATensorMatrixContractionFast(ATensorLeft, stateNormalLeft);

	//compute fluxRight
	LinearAlgebra::MiddleSizeMatrix fluxRight = computeATensorMatrixContractionFast(ATensorRight, stateNormalLeft); // negative sign on matrix, see comment above

 	//Compute integrand Left
	for (std::size_t iB = 0; iB < numOfTestBasisFunctionsLeft; iB++) // for all basis functions
	{
		gradientBasisFunctionLeft = face.basisFunctionDeriv(Base::Side::LEFT,iB);
		gradientBasisFunctionRight = face.basisFunctionDeriv(Base::Side::RIGHT,iB);
		for (std::size_t iD = 0; iD < instance_.DIM_; iD++) // for all dimensions
		{
			for (std::size_t iV = 0; iV < instance_.numOfVariables_; iV++) // for all equations
			{
				iVB = face.getPhysicalElement(Base::Side::LEFT).convertToSingleIndex(iB,iV);
				integrands.first(iVB) += fluxLeft(iV,iD)*gradientBasisFunctionLeft(iD);
				integrands.second(iVB) += fluxRight(iV,iD)*gradientBasisFunctionRight(iD);
			}
		}
	}

	return integrands;
}


/// **************************************************
/// ***    general face integration functions      ***
/// **************************************************

LinearAlgebra::MiddleSizeVector Viscous::integrandStabilityRightHandSideOnFace(
		Base::PhysicalFace<DIM> &face,
		const Base::Side &side,
		const LinearAlgebra::MiddleSizeMatrix &stabilityFluxFunction,
		const std::size_t iD)

{
	std::size_t numBasisFunctions = face.getPhysicalElement(side).getNumOfBasisFunctions();
	std::size_t iVB;
	LinearAlgebra::MiddleSizeVector integrand(numBasisFunctions*instance_.numOfVariables_);

	for (std::size_t iV = 0; iV < instance_.numOfVariables_; iV++)
	{
		for (std::size_t iB = 0; iB < numBasisFunctions; iB++)
		{
			iVB = face.getPhysicalElement(side).convertToSingleIndex(iB,iV);
			integrand(iVB) = stabilityFluxFunction(iV,iD)*face.basisFunction(side, iB);
		}
	}
	return integrand;
}

LinearAlgebra::MiddleSizeVector Viscous::computeRhs(
		const Base::Face *ptrFace,
		const Base::Side &side,
		const LinearAlgebra::MiddleSizeMatrix &stabilityFluxFunction,
		const std::size_t iD)
{
    std::function<LinearAlgebra::MiddleSizeVector(Base::PhysicalFace<DIM>&)> integrandFunction = [=](Base::PhysicalFace<DIM> &face) -> LinearAlgebra::MiddleSizeVector
    {   return this->integrandStabilityRightHandSideOnFace(face, side, stabilityFluxFunction, iD);};
	return stabilityFaceIntegrator_.integrate(ptrFace, integrandFunction);
}

void Viscous::setStabilityMassMatrix(LinearAlgebra::MiddleSizeMatrix &stabilityMassMatrix)
{
	stabilityMassMatrix_ = stabilityMassMatrix;
}






















