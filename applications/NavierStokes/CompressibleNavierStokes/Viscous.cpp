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
//todo: build another face integrator
Viscous::Viscous(CompressibleNavierStokes& instance) : instance_(instance), PrInv_(1/0.71), cp_(1000), ATensorInternal_(instance_.DIM_*instance_.DIM_), ATensorExternal_(instance_.DIM_*instance_.DIM_)
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
	//todo: division by rho has to be computed an aweful lot of times
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


//A_ikrs = A_(iV)(iD)(iVm)(iDm)
//This is a naive implementation, Many terms could be zero in this computation
//todo: A quicker specialised function can be written to improve speed
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
	ATensorInternal_ = computeATensor(partialState, viscosity);

	// Compute flux (A_ikrs*Ur,s)
	LinearAlgebra::MiddleSizeMatrix fluxFunction = computeATensorMatrixContraction(ATensorInternal_, stateJacobian);

	for (std::size_t iB = 0; iB < numberOfBasisFunctions; iB++) // for all basis functions
	{
		gradientBasisFunction = element.basisFunctionDeriv(iB);
		for (std::size_t iD = 0; iD < instance_.DIM_; iD++) // for all dimensions
		{
			for (std::size_t iV = 0; iV < instance_.DIM_ + 2; iV++) // for all variables
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
		for (std::size_t iV = 0; iV < instance_.DIM_; iV++)
		{
			iVB = face.getPhysicalElement(Base::Side::LEFT).convertToSingleIndex(iB,iV);
			integrand(iVB) += face.basisFunction(Base::Side::LEFT, iB)*partialIntegrand(iV);
		}
	}

	return integrand; //no minus sign because integral is on rhs
}


/// **************************************************
/// ***    Internal face integration functions     ***
/// **************************************************

void Viscous::computeStabilityFluxFunction(
		const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorInternal,
		const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorExternal,
		const LinearAlgebra::MiddleSizeVector &stateInternal,
		const LinearAlgebra::MiddleSizeVector &stateExternal,
		const LinearAlgebra::SmallVector<DIM> &normalInternal)
{
	//Compute velocity normal matrix
	LinearAlgebra::MiddleSizeMatrix velocityNormal(instance_.DIM_+2,instance_.DIM_);
	LinearAlgebra::MiddleSizeVector stateDifference = stateInternal - stateExternal;
	for (std::size_t iV = 0; iV < instance_.numOfVariables_; iV++)
	{
		for (std::size_t iD = 0; iD < instance_.DIM_; iD++)
		{
			velocityNormal(iV,iD) = stateDifference(iV)*normalInternal(iD);
		}
	}

	//Compute fluxFunction
	stabilityFluxFunctionInternal_ = computeATensorMatrixContraction(ATensorInternal, velocityNormal);
	stabilityFluxFunctionExternal_ = computeATensorMatrixContraction(ATensorExternal, velocityNormal);
}


LinearAlgebra::MiddleSizeMatrix Viscous::computeStabilityParameters(
		Base::PhysicalFace<DIM> &face,
		const Base::Side &iSide,
		const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorInternal,
		const std::vector<LinearAlgebra::MiddleSizeMatrix> &ATensorExternal,
		const LinearAlgebra::MiddleSizeVector &stateInternal,
		const LinearAlgebra::MiddleSizeVector &stateExternal,
		const LinearAlgebra::SmallVector<DIM> &normalInternal)
{
	Base::Side eSide;
	std::size_t numberOfTestBasisFunctionsInternal;
	std::size_t numberOfTestBasisFunctionsExternal;

	//todo: Investigate if this can be written down easier
	if (iSide == Base::Side::RIGHT)
	{
		eSide = Base::Side::LEFT;
		numberOfTestBasisFunctionsInternal = face.getPhysicalElement(iSide).getNumOfBasisFunctions();
		numberOfTestBasisFunctionsExternal= face.getPhysicalElement(eSide).getNumOfBasisFunctions();
	}
	else
	{
		eSide = Base::Side::RIGHT;
		numberOfTestBasisFunctionsInternal = face.getPhysicalElement(iSide).getNumOfBasisFunctions();
		numberOfTestBasisFunctionsExternal = face.getPhysicalElement(eSide).getNumOfBasisFunctions();
	}

	//Create datastructures
	LinearAlgebra::MiddleSizeMatrix stabilityParametersAverage(instance_.numOfVariables_,instance_.DIM_);
	LinearAlgebra::MiddleSizeVector stabilityParametersInternal(instance_.numOfVariables_);
	LinearAlgebra::MiddleSizeVector stabilityParametersExternal(instance_.numOfVariables_);
	LinearAlgebra::MiddleSizeVector rhsInternal(numberOfTestBasisFunctionsInternal);
	LinearAlgebra::MiddleSizeVector rhsExternal(numberOfTestBasisFunctionsExternal);
	LinearAlgebra::MiddleSizeVector stabilityParameterCoefficientsInternal;
	LinearAlgebra::MiddleSizeVector stabilityParameterCoefficientsExternal;



	//Compute the flux functions as function of iE and iD
	computeStabilityFluxFunction(ATensorInternal, ATensorExternal, stateInternal, stateExternal, normalInternal);

	for (std::size_t iD = 0; iD < instance_.DIM_; iD++)
	{
			//Integrate the correct fluxfunction
			rhsInternal = 0.5*computeRhs(face.getFace(), iSide, stabilityFluxFunctionInternal_, normalInternal, iD); //Add 0.5 factor from the average of the Left and Right A tensor
			rhsExternal = 0.5*computeRhs(face.getFace(), eSide, stabilityFluxFunctionExternal_, normalInternal, iD);

			//Solve the stabilityParameter coefficients
			stabilityMassMatrix_.solve(rhsInternal);
			stabilityMassMatrix_.solve(rhsExternal);

			//Reconstruct the stabilityparameter coefficients for a given iD
			stabilityParametersInternal = instance_.computeStateOnFace(face, iSide, rhsInternal); // rhsInternal are now the reconstruction parameters
			stabilityParametersExternal = instance_.computeStateOnFace(face, eSide, rhsExternal);

			//Put result into matrix
			for (std::size_t iE = 0; iE < instance_.numOfVariables_; iE++)
			{
				stabilityParametersAverage(iE,iD) = 0.5*(stabilityParametersInternal(iE) + stabilityParametersExternal(iE));
			}
	}

	return stabilityParametersAverage;
}

LinearAlgebra::MiddleSizeMatrix Viscous::computeAuxilliaryFlux(
		Base::PhysicalFace<DIM> &face,
		const Base::Side &iSide,
		const LinearAlgebra::MiddleSizeVector &stateInternal,
		const LinearAlgebra::MiddleSizeVector &stateExternal,
		const double pressureInternal,
		const double pressureExternal,
		const LinearAlgebra::MiddleSizeVector &partialStateInternal,
		const LinearAlgebra::MiddleSizeVector &partialStateExternal,
		const LinearAlgebra::MiddleSizeMatrix &stateJacobianInternal,
		const LinearAlgebra::MiddleSizeMatrix &stateJacobianExternal,
		const LinearAlgebra::SmallVector<DIM> &normalInternal)
{
	double eta = 3.0; //stability value hard coded
	LinearAlgebra::MiddleSizeMatrix fluxInternal;
	LinearAlgebra::MiddleSizeMatrix fluxExternal;
	LinearAlgebra::MiddleSizeMatrix stabilityFlux;
	LinearAlgebra::MiddleSizeMatrix AuxilliaryFlux;

	//Compute A and A contraction with Jacobian
	double temperatureInternal = computeTemperature(stateInternal, pressureInternal);
	double viscosityInternal = computeViscosity(temperatureInternal);
	double temperatureExternal = computeTemperature(stateExternal, pressureExternal);
	double viscosityExternal = computeViscosity(temperatureExternal);

	ATensorInternal_ = computeATensor(partialStateInternal, viscosityInternal);
	ATensorExternal_ = computeATensor(partialStateExternal, viscosityExternal);

	fluxInternal  = computeATensorMatrixContraction(ATensorInternal_, stateJacobianInternal);
	fluxExternal = computeATensorMatrixContraction(ATensorExternal_, stateJacobianExternal);
	stabilityFlux = computeStabilityParameters(face, iSide, ATensorInternal_, ATensorExternal_, stateInternal, stateExternal, normalInternal);

	AuxilliaryFlux = 0.5*(fluxInternal + fluxExternal) - eta*stabilityFlux;

	return AuxilliaryFlux;
}

//todo:: Check if the normal vector here is used correctly
LinearAlgebra::MiddleSizeVector Viscous::integrandAuxilliaryAtFace(
		Base::PhysicalFace<DIM> &face,
		const Base::Side &iSide,
		const LinearAlgebra::MiddleSizeVector &stateInternal,
		const LinearAlgebra::MiddleSizeVector &stateExternal,
		const double pressureInternal,
		const double pressureExternal,
		const LinearAlgebra::MiddleSizeVector &partialStateInternal,
		const LinearAlgebra::MiddleSizeVector &partialStateExternal,
		const LinearAlgebra::MiddleSizeMatrix &stateJacobianInternal,
		const LinearAlgebra::MiddleSizeMatrix &stateJacobianExternal)
{
    const LinearAlgebra::SmallVector<DIM> normalInternal = face.getUnitNormalVector();
	LinearAlgebra::MiddleSizeMatrix fluxFunction = computeAuxilliaryFlux(face, iSide, stateInternal, stateExternal, pressureInternal, pressureExternal, partialStateInternal, partialStateExternal, stateJacobianInternal, stateJacobianExternal, normalInternal);

	//Compute integrand
	std::size_t numOfTestBasisFunctions = face.getPhysicalElement(iSide).getNumOfBasisFunctions();
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
		for (std::size_t iV = 0; iV < instance_.DIM_; iV++)
		{
			iVB = face.getPhysicalElement(iSide).convertToSingleIndex(iB,iV);
			integrand(iVB) += face.basisFunction(iSide, iB)*partialIntegrand(iV);
		}
	}

	return -integrand;
}

LinearAlgebra::MiddleSizeVector Viscous::integrandViscousAtFace(
		Base::PhysicalFace<DIM> &face,
		const Base::Side &iSide,
		const LinearAlgebra::MiddleSizeVector &stateInternal,
		const LinearAlgebra::MiddleSizeVector &stateExternal,
		double pressure,
		const LinearAlgebra::MiddleSizeVector &partialStateInternal)
{
	LinearAlgebra::SmallVector<DIM> normal = face.getUnitNormalVector();

	//Compute velocity normal matrix
	LinearAlgebra::MiddleSizeMatrix velocityNormal(instance_.DIM_+2,instance_.DIM_);
	LinearAlgebra::MiddleSizeVector stateDifference;
	stateDifference = stateInternal - stateExternal;
	for (std::size_t iV = 0; iV < instance_.numOfVariables_; iV++)
	{
		for (std::size_t iD = 0; iD < instance_.DIM_; iD++)
		{
			velocityNormal(iV,iD) = 0.5*stateDifference(iV)*normal(iD);
		}
	}

	//Compute A and A contraction with velocitynormal matrix
	double temperature = computeTemperature(stateInternal, pressure);
	double viscosity = computeViscosity(temperature);

	ATensorInternal_ = computeATensor(partialStateInternal, viscosity);
	LinearAlgebra::MiddleSizeMatrix fluxFunction = computeATensorMatrixContraction(ATensorInternal_, velocityNormal);

	//Compute integrand
	std::size_t numOfTestBasisFunctions = face.getPhysicalElement(iSide).getNumOfBasisFunctions();
	LinearAlgebra::MiddleSizeVector integrand(instance_.numOfVariables_*numOfTestBasisFunctions);
	LinearAlgebra::SmallVector<DIM> gradientBasisFunction;
	std::size_t iVB; // Index for both variable and basis function.

	for (std::size_t iB = 0; iB < numOfTestBasisFunctions; iB++) // for all basis functions
	{
		gradientBasisFunction = face.basisFunctionDeriv(iSide,iB);
		for (std::size_t iD = 0; iD < instance_.DIM_; iD++) // for all dimensions
		{
			for (std::size_t iV = 0; iV < instance_.DIM_ + 2; iV++) // for all equations
			{
				iVB = face.getPhysicalElement(iSide).convertToSingleIndex(iB,iV);
				integrand(iVB) += fluxFunction(iV,iD)*gradientBasisFunction(iD);
			}
		}
	}
	return integrand; //no negative sign: integral is on right hand side of equation
}


/// **************************************************
/// ***    general face integration functions      ***
/// **************************************************

LinearAlgebra::MiddleSizeVector Viscous::integrandStabilityRightHandSideOnRefFace(
		Base::PhysicalFace<DIM> &face,
		const Base::Side &side,
		const LinearAlgebra::MiddleSizeMatrix &stabilityFluxFunction,
		const std::size_t iD,
		const LinearAlgebra::SmallVector<DIM> &normalInternal)

{
	std::size_t numBasisFunctions = face.getPhysicalElement(side).getNumOfBasisFunctions();
	std::size_t pos;
	LinearAlgebra::MiddleSizeVector integrand(numBasisFunctions*instance_.numOfVariables_);

	for (std::size_t iV = 0; iV < instance_.numOfVariables_; iV++)
	{
		for (std::size_t iB = 0; iB < numBasisFunctions; iB++)
		{
			pos = numBasisFunctions*iV +iB;
			integrand(pos) = stabilityFluxFunction(iV,iD)*face.basisFunction(side, iB);
		}
	}
	return integrand;
}

//todo: I use integrate here, not referenceIntegrate: does it make a problem (?)
LinearAlgebra::MiddleSizeVector Viscous::computeRhs(
		const Base::Face *ptrFace,
		const Base::Side &side,
		const LinearAlgebra::MiddleSizeMatrix &stabilityFluxFunction,
		const LinearAlgebra::SmallVector<DIM> &normalInternal,
		const std::size_t iD)
{
    std::function<LinearAlgebra::MiddleSizeVector(Base::PhysicalFace<DIM>&)> integrandFunction = [=](Base::PhysicalFace<DIM> &face) -> LinearAlgebra::MiddleSizeVector
    {   return this->integrandStabilityRightHandSideOnRefFace(face, side, stabilityFluxFunction, iD, normalInternal);};
	return instance_.faceIntegrator_.integrate(ptrFace, integrandFunction);
}

void Viscous::setStabilityMassMatrix(LinearAlgebra::MiddleSizeMatrix &stabilityMassMatrix)
{
	stabilityMassMatrix_ = stabilityMassMatrix;
}






















