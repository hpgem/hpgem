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

#include "InviscidTerms.h"

/// *****************************************
/// ***   Element integration functions   ***
/// *****************************************

/// \brief Compute integrand of righthandside on an element
template<std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
LinearAlgebra::MiddleSizeVector InviscidTerms<DIM,NUMBER_OF_VARIABLES>::integrandAtElement
(
 Base::PhysicalElement<DIM> &element,
 const StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> &elementStateStruct,
 const double time
 )
{
	// Get the number of basis functions in an element.
	std::size_t numberOfBasisFunctions =  element.getNumberOfBasisFunctions();

	//Create data structures for calculating the integrand
	LinearAlgebra::MiddleSizeVector integrand(NUMBER_OF_VARIABLES * numberOfBasisFunctions); //The final integrand value will be stored in this vector
	LinearAlgebra::SmallVector<DIM> gradientBasisFunction; //Gradient function based on the number of DIMs

	//Create iteration values
	std::size_t iVB; // Index in solution coefficients for variable i and basisfunction j

	//Compute the flux
	LinearAlgebra::MiddleSizeMatrix fluxMatrix = elementStateStruct.getHyperbolicMatrix();

	// Compute the integrand for all equations
	for (std::size_t iB = 0; iB < numberOfBasisFunctions; iB++) // For every basis function
	{
		gradientBasisFunction = element.basisFunctionDeriv(iB); // Compute the gradient of the ith test function

		for (std::size_t iV = 0; iV < NUMBER_OF_VARIABLES; iV++) // For all variables
		{
			for (std::size_t iD = 0; iD < DIM; iD++) // For all DIMs
			{

				iVB = element.convertToSingleIndex(iB,iV);
				integrand(iVB) += fluxMatrix(iV,iD)*gradientBasisFunction(iD);
			}
		}
	}

	return integrand;
}

/// **************************************************
/// ***    General face integration functions      ***
/// **************************************************

/// \brief Compute the local Lax-Friedrichs flux
template<std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
LinearAlgebra::MiddleSizeVector InviscidTerms<DIM,NUMBER_OF_VARIABLES>::computeLLFFluxFunction
(
 const LinearAlgebra::MiddleSizeVector stateLeft,
 const LinearAlgebra::MiddleSizeVector stateRight,
 const double pressureLeft,
 const double pressureRight,
 const LinearAlgebra::SmallVector<DIM> &unitNormalLeft
 )
{
	const StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> constitutiveRelations = instance_->constitutiveRelationsStruct;

	//compute the speed normal to the face
	double inverseDensityLeft = 1.0/stateLeft(0);
	double inverseDensityRight = 1.0/stateRight(0);
	double normalSpeedLeft = 0.0;
	double normalSpeedRight = 0.0;
	for (std::size_t iD = 0; iD < DIM; iD++)
	{
		normalSpeedLeft += stateLeft(iD+1)*inverseDensityLeft*unitNormalLeft(iD);
		normalSpeedRight += stateRight(iD+1)*inverseDensityRight*unitNormalLeft(iD);
	}

	///compute speed of sound, left and right of the face
	double aLeft = constitutiveRelations.computeSpeedOfSound(stateLeft,pressureLeft);
	double aRight = constitutiveRelations.computeSpeedOfSound(stateRight, pressureRight);

	///compute wave speeds
	double waveSpeedLeft = std::min(normalSpeedLeft - aLeft, normalSpeedRight - aRight);
	double waveSpeedRight = std::max(normalSpeedLeft + aLeft, normalSpeedRight + aRight);

	double lMax = std::max(std::abs(waveSpeedLeft), std::abs(waveSpeedRight));

	//todo: remove this, temporary hack
	LinearAlgebra::MiddleSizeMatrix fluxMatrixCombined = constitutiveRelations.computeHyperbolicMatrix(stateLeft, pressureLeft)+ constitutiveRelations.computeHyperbolicMatrix(stateLeft, pressureLeft);
	LinearAlgebra::MiddleSizeVector flux(NUMBER_OF_VARIABLES);
	for (std::size_t iV = 0; iV < NUMBER_OF_VARIABLES; iV++)
	{
		for (std::size_t iD = 0; iD < DIM; iD++)
		{
			flux(iV) = fluxMatrixCombined(iV,iD)*unitNormalLeft(iD);
		}
	}

	return 0.5*flux	- 0.5*lMax*(stateRight- stateLeft);
}

/// \brief Compute the local Lax-Friedrichs flux
template<std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
LinearAlgebra::MiddleSizeVector InviscidTerms<DIM,NUMBER_OF_VARIABLES>::computeLLFFluxFunction
(
 const StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> faceStateStructLeft,
 const StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> faceStateStructRight,
 const LinearAlgebra::SmallVector<DIM> &unitNormalLeft
 )
{
	//For convenience:
	LinearAlgebra::MiddleSizeVector stateLeft = faceStateStructLeft.getState();
	LinearAlgebra::MiddleSizeVector stateRight = faceStateStructRight.getState();
	double pressureLeft = faceStateStructLeft.getPressure();
	double pressureRight = faceStateStructRight.getPressure();

	//compute the speed normal to the face
	double inverseDensityLeft = 1.0/stateLeft(0);
	double inverseDensityRight = 1.0/stateRight(0);
	double normalSpeedLeft = 0.0;
	double normalSpeedRight = 0.0;
	for (std::size_t iD = 0; iD < DIM; iD++)
	{
		normalSpeedLeft += stateLeft(iD+1)*inverseDensityLeft*unitNormalLeft(iD);
		normalSpeedRight += stateRight(iD+1)*inverseDensityRight*unitNormalLeft(iD);
	}

	///compute speed of sound, left and right of the face
	double aLeft = faceStateStructLeft.getSpeedOfSound();
	double aRight = faceStateStructRight.getSpeedOfSound();

	///compute wave speeds
	double waveSpeedLeft = std::min(normalSpeedLeft - aLeft, normalSpeedRight - aRight);
	double waveSpeedRight = std::max(normalSpeedLeft + aLeft, normalSpeedRight + aRight);

/*	std::cout << "aLeft: " << aLeft << std::endl;
	std::cout << "aRight: " << aRight << std::endl;
	std::cout << "waveSpeedLeft: " << waveSpeedLeft << std::endl;
	std::cout << "waveSpeedright: "  << waveSpeedRight << std::endl;
	std::exit(-1);*/

	double lMax = std::max(std::abs(waveSpeedLeft), std::abs(waveSpeedRight));

	//todo: remove this, temporary hack
	LinearAlgebra::MiddleSizeMatrix fluxMatrixCombined = faceStateStructLeft.getHyperbolicMatrix() + faceStateStructRight.getHyperbolicMatrix();
	LinearAlgebra::MiddleSizeVector flux(NUMBER_OF_VARIABLES);
	for (std::size_t iV = 0; iV < NUMBER_OF_VARIABLES; iV++)
	{
		for (std::size_t iD = 0; iD < DIM; iD++)
		{
			flux(iV) += fluxMatrixCombined(iV,iD)*unitNormalLeft(iD);
		}
	}

/*	std::cout << "====" << std::endl;
	std::cout << "flux: " << 0.5*flux << std::endl;
	std::cout << "stability: " << - 0.5*lMax*(stateRight - stateLeft) << std::endl;
	std::cout << "----" << std::endl;*/
	return 0.5*flux	- 0.5*lMax*(stateRight - stateLeft);
}

//Compute the Roe Riemann Flux function
template<std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
LinearAlgebra::MiddleSizeVector InviscidTerms<DIM,NUMBER_OF_VARIABLES>::computeRoeFluxFunction
(
		 const StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> faceStateStructLeft,
		 const StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> faceStateStructRight,
		 const LinearAlgebra::SmallVector<DIM> &unitNormalLeft
		 )
{
	//Compute correct normal direction and difference vector
	LinearAlgebra::MiddleSizeVector stateLeft = faceStateStructLeft.getState();
	LinearAlgebra::MiddleSizeVector stateRight = faceStateStructRight.getState();
	LinearAlgebra::MiddleSizeVector stateDifference = stateRight- stateLeft;

	//Compute the Roe average state
	LinearAlgebra::MiddleSizeVector stateAverage(DIM+1);
	double zL = std::sqrt(stateLeft(0));
	double zR = std::sqrt(stateRight(0));
	double tmp1 = 1.0/(stateLeft(0) + zL*zR);
	double tmp2 = 1.0/(stateRight(0) + zL*zR);
	double ruSquaredLeft = 0.0;
	double ruSquaredRight = 0.0;
	double rhoInverseLeft = 1.0/stateLeft(0);
	double rhoInverseRight = 1.0/stateRight(0);
	double pressureLeft;
	double pressureRight;


	for (std::size_t iD = 0; iD < DIM; iD++)
	{
		stateAverage(iD) = (stateLeft(iD+1)*tmp1 + stateRight(iD+1)*tmp2); 	// u_average
		ruSquaredLeft += stateLeft(iD+1)*stateLeft(iD+1); 											// Kinetic part of the left pressure term
		ruSquaredRight += stateRight(iD+1)*stateRight(iD+1);											// Kinetic part of the right pressure term
	}

	//todo: function layer above contains this information already
	pressureLeft = (GAMMA - 1)*(stateLeft(DIM + 1) - 0.5*ruSquaredLeft*rhoInverseLeft);
	pressureRight = (GAMMA - 1)*(stateRight(DIM + 1) - 0.5*ruSquaredRight*rhoInverseRight);

	stateAverage(DIM) = (stateLeft(DIM+1) + pressureLeft)*tmp1 + (stateRight(DIM+1) + pressureRight)*tmp2; //H_average

	//Compute useful variables for constructing the flux function
	double alphaAvg = 0.0;
	double unAvg = 0.0;
	for (std::size_t iD = 0; iD < DIM; iD++)
	{
		alphaAvg += (stateAverage(iD)*stateAverage(iD));
		unAvg += stateAverage(iD)*unitNormalLeft(iD);
	}
	alphaAvg *= 0.5;

	const double a2Avg = std::abs((GAMMA -1)*(stateAverage(DIM) - alphaAvg));
	const double aAvg = std::sqrt(a2Avg);
	const double ovaAvg = 1.0/aAvg;
	const double ova2Avg = 1.0/a2Avg;

	//Compute eigenvalues
	double lam1 = std::abs(unAvg + aAvg);
	double lam2 = std::abs(unAvg - aAvg);
	double lam3 = std::abs(unAvg);

	//Add entropy correction
	double epsilon = 0.01;
	if (lam1 < epsilon)
	{
		lam1 = (lam1*lam1 + epsilon*epsilon)/(2.0*epsilon);
	}

	if (lam2 < epsilon)
	{
		lam2 = (lam2*lam2 + epsilon*epsilon)/(2.0*epsilon);
	}

	if (lam3 < epsilon)
	{
		lam3 = (lam3*lam3 + epsilon*epsilon)/(2.0*epsilon);
	}

	//Compute useful abbreviations for constructing the flux function
	const double abv1 = 0.5*(lam1 + lam2);
	const double abv2 = 0.5*(lam1 - lam2);
	const double abv3 = abv1 - lam3;

	double abv4 = alphaAvg*stateDifference(0);
	double abv5 = -unAvg*stateDifference(0);
	for (std::size_t iD = 0; iD < DIM; iD++)
	{
		abv4 += -stateAverage(iD)*stateDifference(iD+1);
		abv5 += unitNormalLeft(iD)*stateDifference(iD+1);
	}
	abv4 += stateDifference(DIM+1);
	abv4 *= (GAMMA -1);

	const double abv6 = abv3*abv4*ova2Avg + abv2*abv5*ovaAvg;
	const double abv7 = abv2*abv4*ovaAvg + abv3*abv5;

	//Compute the Roe Riemann Flux function: h(u_L,u_R) = 0.5*(F(u_L) + F(u_R) - |A|(u_R - u_L))
	LinearAlgebra::MiddleSizeVector flux(DIM + 2);

	double pLR = pressureLeft + pressureRight;

	double runR = 0.0;
	double runL = 0.0;
	for (std::size_t iD = 0; iD < DIM; iD++)
	{
		runL += stateLeft(iD+1)*unitNormalLeft(iD);
		runR += stateRight(iD+1)*unitNormalLeft(iD);
	}

	double unL = runL*rhoInverseLeft;
	double unR = runR*rhoInverseRight;

	//continuity equation
	flux(0) = (runL + runR - (lam3*stateDifference(0) + abv6));

	//momentum equations
	for (std::size_t iD = 0; iD < DIM; iD++)
	{
		flux(iD+1) = runL*stateLeft(iD+1)*rhoInverseLeft + runR*stateRight(iD+1)*rhoInverseRight + pLR*unitNormalLeft(iD) - (lam3*stateDifference(iD+1) + stateAverage(iD)*abv6 + unitNormalLeft(iD)*abv7);
	}

	//energy equation
	flux(DIM+1) = (unL*(stateLeft(DIM+1) + pressureLeft) + unR*(stateRight(DIM+1) + pressureRight) - (lam3*stateDifference(DIM+1) + stateAverage(DIM)*abv6 + unAvg*abv7));

	//Note: Twice the flux is computed above, hence the factor 0.5 in front of the equation
	return 0.5*flux;
}

/// \brief Compute the HLLC Flux function, based on Klaij et al. 2006
template<std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
LinearAlgebra::MiddleSizeVector InviscidTerms<DIM,NUMBER_OF_VARIABLES>::computeHLLCFluxFunction
(
 const StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> &faceStateStructLeft,
 const StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> &faceStateStructRight,
 const LinearAlgebra::SmallVector<DIM> &unitNormalLeft
)
{
	//For convenience:
	LinearAlgebra::MiddleSizeVector stateLeft = faceStateStructLeft.getState();
	LinearAlgebra::MiddleSizeVector stateRight = faceStateStructRight.getState();
	double pressureLeft = faceStateStructLeft.getPressure();
	double pressureRight = faceStateStructRight.getPressure();

	//compute the speed normal to the face
	double inverseDensityLeft = 1.0/stateLeft(0);
	double inverseDensityRight = 1.0/stateRight(0);
	double normalSpeedLeft = 0.0;
	double normalSpeedRight = 0.0;
	for (std::size_t iD = 0; iD < DIM; iD++)
	{
		normalSpeedLeft += stateLeft(iD+1)*inverseDensityLeft*unitNormalLeft(iD);
		normalSpeedRight += stateRight(iD+1)*inverseDensityRight*unitNormalLeft(iD);
	}

	///compute speed of sound, left and right of the face
	double aLeft = faceStateStructLeft.getSpeedOfSound();
	double aRight = faceStateStructRight.getSpeedOfSound();

	///compute wave speeds
	double waveSpeedLeft = std::min(normalSpeedLeft - aLeft, normalSpeedRight - aRight);
	double waveSpeedRight = std::max(normalSpeedLeft + aLeft, normalSpeedRight + aRight);



	///compute middle wave speed
	double tempLeft = stateLeft(0)*normalSpeedLeft*(waveSpeedLeft - normalSpeedLeft) - pressureLeft;
	double tempRight = stateRight(0)*normalSpeedRight*(waveSpeedRight - normalSpeedRight) - pressureRight;
	double tempBottom = stateRight(0)*(waveSpeedRight - normalSpeedRight) - stateLeft(0)*(waveSpeedLeft - normalSpeedLeft);
	double waveSpeedMiddle = (tempRight - tempLeft)/tempBottom;

/*	std::cout << "=========" << std::endl;
	std::cout << "aLeft: " << aLeft << std::endl;
	std::cout << "aRight: " << aRight << std::endl;
	std::cout << "normalSpeedLeft: " << normalSpeedLeft << std::endl;
	std::cout << "normalSpeedright: " << normalSpeedRight << std::endl;
	std::cout << "inverseDensityLeft: " << inverseDensityLeft << std::endl;
	std::cout << "inverseDensityRight: " << inverseDensityRight << std::endl;
	std::cout << "--------------" << std::endl;
	std::cout << "left: " << waveSpeedLeft << std::endl;
	std::cout << "right: " << waveSpeedRight << std::endl;
	std::cout << "middle: " << waveSpeedMiddle << std::endl;*/

	///compute intermediate pressure
	double pIntermediate = stateLeft(0)*(waveSpeedLeft - normalSpeedLeft)*(waveSpeedMiddle - normalSpeedLeft) + pressureLeft;

	///compute intermediate states
	double factorLeft = 1.0/(waveSpeedLeft - waveSpeedMiddle);
	double factorRight = 1.0/(waveSpeedRight - waveSpeedMiddle);
	LinearAlgebra::MiddleSizeVector stateIntermediateLeft(NUMBER_OF_VARIABLES);
	LinearAlgebra::MiddleSizeVector stateIntermediateRight(NUMBER_OF_VARIABLES);

	//compute intermediate states: part 1
	for (std::size_t iV = 0; iV < NUMBER_OF_VARIABLES; iV++)
	{
		stateIntermediateLeft(iV) = factorLeft*(waveSpeedLeft - normalSpeedLeft)*stateLeft(iV);
		stateIntermediateRight(iV) = factorRight*(waveSpeedRight - normalSpeedRight)*stateRight(iV);
	}

	//compute intermediate states: part 2
	for (std::size_t iD = 0; iD < DIM; iD++)
	{
		stateIntermediateLeft(iD+1) +=  factorLeft*(pIntermediate - pressureLeft)*unitNormalLeft(iD);
		stateIntermediateRight(iD+1) += factorRight*(pIntermediate - pressureRight)*unitNormalLeft(iD);
	}
	//TEMP energy addition: FIX THIS SHIT IN FUTURE
	//todo: see line above
	stateIntermediateLeft(DIM+1) += factorLeft*(pIntermediate*waveSpeedMiddle - pressureLeft*normalSpeedLeft);
	stateIntermediateRight(DIM+1) += factorRight*(pIntermediate*waveSpeedMiddle - pressureRight*normalSpeedRight);

	///compute normalFluxLeft and normalFluxRight
	LinearAlgebra::MiddleSizeMatrix fluxMatrixLeft = faceStateStructLeft.getHyperbolicMatrix();
	LinearAlgebra::MiddleSizeMatrix fluxMatrixRight = faceStateStructRight.getHyperbolicMatrix();
	LinearAlgebra::MiddleSizeVector normalFluxLeft(NUMBER_OF_VARIABLES);
	LinearAlgebra::MiddleSizeVector normalFluxRight(NUMBER_OF_VARIABLES);
	for (std::size_t iV = 0; iV < NUMBER_OF_VARIABLES; iV++)
	{
		for (std::size_t iD = 0; iD < DIM; iD++)
		{
			normalFluxLeft(iV) += fluxMatrixLeft(iV,iD)*unitNormalLeft(iD);
			normalFluxRight(iV) += fluxMatrixRight(iV,iD)*unitNormalLeft(iD);
		}
	}

	//compute flux function
	LinearAlgebra::MiddleSizeVector fluxFunction(NUMBER_OF_VARIABLES);

	fluxFunction = 0.5*(normalFluxLeft + normalFluxRight
				 + (std::abs(waveSpeedMiddle) - std::abs(waveSpeedLeft))*stateIntermediateLeft + std::abs(waveSpeedLeft)*stateLeft
				 + (std::abs(waveSpeedRight) - std::abs(waveSpeedMiddle))*stateIntermediateRight - std::abs(waveSpeedRight)*stateRight);

	return fluxFunction;
}

/// **************************************************
/// ***    External face integration functions     ***
/// **************************************************

/// \brief Compute the integrand for the right hand side for the reference face corresponding to a external face.
template<std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
LinearAlgebra::MiddleSizeVector InviscidTerms<DIM,NUMBER_OF_VARIABLES>::integrandAtBoundaryFace
(
 Base::PhysicalFace<DIM> &face,
 const BoundaryType boundaryType,
 const StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> &faceStateStructBoundary,
 const StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> &faceStateStructLeft,
 const double &time
 )
{
	//Get the number of basis functions
	std::size_t numOfTestBasisFunctions = face.getPhysicalElement(Base::Side::LEFT).getNumberOfBasisFunctions();

	LinearAlgebra::MiddleSizeVector integrand(NUMBER_OF_VARIABLES*numOfTestBasisFunctions);

	//Compute the flux. Based on the type of boundary condition this is either done exact, or using a flux function
	LinearAlgebra::MiddleSizeVector flux(NUMBER_OF_VARIABLES);
	if (boundaryType == BoundaryType::FULL_STATE)
	{
		 flux = computeRoeFluxFunction(faceStateStructLeft, faceStateStructBoundary, face.getUnitNormalVector());
	}
	else
	{
		//todo: move this part to a separate function
		LinearAlgebra::MiddleSizeMatrix fluxMatrix = faceStateStructBoundary.getHyperbolicMatrix();
		for (std::size_t iV = 0; iV < NUMBER_OF_VARIABLES; iV++)
		{
			for (std::size_t iD = 0; iD < DIM; iD++)
			{
				flux(iV) += fluxMatrix(iV,iD)*face.getUnitNormalVector()(iD);
			}
		}
	}

	// Compute integrand on the reference element.
	std::size_t iVB; // Index for both variable and basis function.
	for (std::size_t iB = 0; iB < numOfTestBasisFunctions; iB++)
	{
		for (std::size_t iV = 0; iV < NUMBER_OF_VARIABLES; iV++) // Index for direction
		{
			iVB = face.getPhysicalElement(Base::Side::LEFT).convertToSingleIndex(iB, iV);
			integrand(iVB) = -flux(iV)*face.basisFunction(Base::Side::LEFT, iB); // Minus sign because the integral is on the right hand side
		}
	}

	return  integrand;
}


/// **************************************************
/// ***    Internal face integration functions     ***
/// **************************************************

/// \brief Compute both the face integral for the left element as the right element at the same time. (reducing flux calculations)
template<std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
std::pair<LinearAlgebra::MiddleSizeVector,LinearAlgebra::MiddleSizeVector> InviscidTerms<DIM,NUMBER_OF_VARIABLES>::integrandsAtFace
(
 Base::PhysicalFace<DIM> &face,
 const StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> &faceStateStructLeft,
 const StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> &faceStateStructRight,
 const double &time
 )
{
	//Data structures for left and right integrand
	std::size_t numOfTestBasisFunctionsLeft = face.getPhysicalElement(Base::Side::LEFT).getNumberOfBasisFunctions();
	std::size_t numOfTestBasisFunctionsRight = face.getPhysicalElement(Base::Side::RIGHT).getNumberOfBasisFunctions();
	std::pair<LinearAlgebra::MiddleSizeVector,LinearAlgebra::MiddleSizeVector> integrands(
	 			std::piecewise_construct,
	 			std::forward_as_tuple(NUMBER_OF_VARIABLES*numOfTestBasisFunctionsLeft),
	 			std::forward_as_tuple(NUMBER_OF_VARIABLES*numOfTestBasisFunctionsRight));
	std::size_t iVB;

	//Compute left flux
	LinearAlgebra::SmallVector<DIM> unitNormalLeft = face.getUnitNormalVector(); //todo: remove this line
	LinearAlgebra::MiddleSizeVector flux = computeRoeFluxFunction(faceStateStructLeft, faceStateStructRight, unitNormalLeft);

	//Compute integrands
	for (std::size_t iV = 0; iV < NUMBER_OF_VARIABLES; iV++) // Index for direction
	{
		//Compute left integrand
		for (std::size_t iB = 0; iB < numOfTestBasisFunctionsLeft; iB++)
		{
			iVB = face.getPhysicalElement(Base::Side::LEFT).convertToSingleIndex(iB, iV);
			integrands.first(iVB) = -flux(iV)*face.basisFunction(Base::Side::LEFT, iB); // Minus sign because the integral is on the right hand side
		}
		//compute right integrand
		for (std::size_t iB = 0; iB < numOfTestBasisFunctionsRight; iB++)
		{
			iVB = face.getPhysicalElement(Base::Side::RIGHT).convertToSingleIndex(iB, iV);
			integrands.second(iVB) = flux(iV)*face.basisFunction(Base::Side::RIGHT,iB);
		}

	}

	return integrands;
}


/// ************************************************
/// ***    Jacobian Exact Derivative Functions   ***
/// ************************************************

/*
	/// \brief Compute the Jacobian with respect to the solution coefficients of the element
	std::vector<LinearAlgebra::MiddleSizeMatrix> computeFluxJacobian2D(const LinearAlgebra::MiddleSizeVector &state);

	/// \brief Computes the derivative of the InviscidFace flux
	LinearAlgebra::MiddleSizeVector computeInviscidFaceFluxDiv(Base::PhysicalFace<DIM> &face, const LinearAlgebra::MiddleSizeVector &stateLeft, const LinearAlgebra::MiddleSizeVector &stateRight, const LinearAlgebra::SmallVector<DIM> &unitNormalLeft, const Base::Side derivativeSide, const std::size_t iV2, const std::size_t iB2);

	//LinearAlgebra::MiddleSizeVector computeInviscidFaceLLFFluxDiv(Base::PhysicalFace<DIM> &face, const LinearAlgebra::MiddleSizeVector &stateLeft, const LinearAlgebra::MiddleSizeVector &stateRight, const LinearAlgebra::SmallVector<DIM> &unitNormalLeft, const Base::Side derivativeSide, const std::size_t iV2, const std::size_t iB2);

  	//LinearAlgebra::MiddleSizeMatrix computeFluxTimesNormalJacobian(const LinearAlgebra::MiddleSizeVector &state, const LinearAlgebra::SmallVector<DIM> &unitNormal);


*/

/// ************************************************
/// ***    Jacobian Element Matrix Functions     ***
/// ************************************************

/// \brief Compute the inviscid element integral integrand contribution to the Jacobian matrix
template<std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
LinearAlgebra::MiddleSizeMatrix InviscidTerms<DIM,NUMBER_OF_VARIABLES>::integrandJacobianInviscidElement
(
 Base::PhysicalElement<DIM> &element,
 const StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> &elementStateStruct
 )
{
	//Create integrand vector
	std::size_t numberOfBasisFunctions =  element.getNumOfBasisFunctions();
	LinearAlgebra::MiddleSizeMatrix integrand(NUMBER_OF_VARIABLES*numberOfBasisFunctions,NUMBER_OF_VARIABLES*numberOfBasisFunctions);
	LinearAlgebra::SmallVector<DIM> gradientBasisFunction; //Gradient function based on the number of DIMs
	StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> newState;

	//Compute the jacobian of all the flux components. In 2D f, g and in 3D f,g,h
	//std::vector<LinearAlgebra::MiddleSizeMatrix> fluxJacobian = computeFluxJacobian2D(elementStateFunctions.getState());

	//Compute the first stateflux
	LinearAlgebra::MiddleSizeMatrix fluxFirst = computeFluxMatrix(elementStateStruct.getState(), elementStateStruct.getPressure());

	//compute the integrand matrix
	std::size_t iVB1, iVB2;
	for(std::size_t iV2 = 0; iV2 < NUMBER_OF_VARIABLES; iV2++) // This is the variable in the integrand vector
	{
		for (std::size_t iB2 = 0; iB2 < numberOfBasisFunctions; iB2++) // This is the basis function in the integrand vector
		{
			iVB2 = element.convertToSingleIndex(iB2,iV2); // the integrand vector index

			//Compute the peturbed state
			LinearAlgebra::MiddleSizeVector stateCoefficientsPerturb = elementStateStruct.getStateCoefficients();
			//add pertubation to the correct solutionCoefficient
			stateCoefficientsPerturb(iVB2) += EPSILON;
			LinearAlgebra::MiddleSizeVector stateNew = computeStateOnElement<DIM,NUMBER_OF_VARIABLES>(element, stateCoefficientsPerturb);
			double pressureNew = newState.computePressure(stateNew);

			//Compute new state flux
			LinearAlgebra::MiddleSizeMatrix fluxSecond = newState.computeHyperbolicMatrix(stateNew,pressureNew);

			//Compute the difference
			LinearAlgebra::MiddleSizeMatrix fluxDifference = (1.0/EPSILON)*(fluxSecond - fluxFirst);

			for (std::size_t iV1 = 0; iV1 < NUMBER_OF_VARIABLES; iV1++) // This is the variable in the coefficient vector
			{
				//Here I already know that all derivatives with respect to other variables than iV2 are zero
				//so we only need to focus on the derivative of the flux function with respect to iV2
				// dF/dU = dF/dq * dq/dU
				for (std::size_t iB1 = 0; iB1 < numberOfBasisFunctions; iB1++) // This is the basis function in the coefficient vector
				{
					//The derivative of variable iV2 with respect to the state coefficient is the basis function
					// i.e.: u = u0*phi0 + u1*phi1 + u2*phi2 ==> du/du1 = phi1
					iVB1 = element.convertToSingleIndex(iB1,iV1);
					gradientBasisFunction = element.basisFunctionDeriv(iB1);
					//iVB1 gives the value in the residual integrand vector corresponding to equation V1 and basisfunction B1
					//iVB2 gives the value in the coefficient vector corresponding to state value V2 and basisfunction B2
					for (std::size_t iD = 0; iD < DIM; iD++)
					{
						integrand(iVB1,iVB2) += fluxDifference(iV1,iD)*gradientBasisFunction(iD);
					}
				}
			}
		}
	}

	return integrand;
}

/// *********************************************
/// ***    Jacobian Face Matrix Functions     ***
/// *********************************************

/// \brief Computes the inviscid face integral integrand contribution to the Jacobian matrix
template<std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
LinearAlgebra::MiddleSizeMatrix InviscidTerms<DIM,NUMBER_OF_VARIABLES>::integrandJacobianInviscidFace
(
 Base::PhysicalFace<DIM> &face,
 const StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> &faceStateStructLeft,
 const StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> &faceStateStructRight,
 const Base::Side elementSide,
 const Base::Side derivativeSide
 )
{
	//Datastructures
	std::size_t numberOfBasisFunctionsElementSide = face.getPhysicalElement(elementSide).getNumberOfBasisFunctions();
	std::size_t numberOfBasisFunctionsDerivativeSide = face.getPhysicalElement(derivativeSide).getNumOfBasisFunctions();
	LinearAlgebra::MiddleSizeMatrix integrand(NUMBER_OF_VARIABLES*numberOfBasisFunctionsElementSide,NUMBER_OF_VARIABLES*numberOfBasisFunctionsDerivativeSide);
	LinearAlgebra::MiddleSizeVector jacobianFluxFunction(NUMBER_OF_VARIABLES);
	StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> newState;

	//compute the integrand matrix
	//iVB1 gives the value in the residual vector corresponding to state variable V1 and basisfunction B1
	//iVB2 gives the value in the coefficient vector corresponding to state variable V2 and basisfunction B2
	std::size_t iVB1, iVB2;

	//Set the correct sign for the flux
	double sign;
	if (elementSide == Base::Side::LEFT)
	{
		sign = -1.0;
	}
	else
	{
		sign = 1.0;
	}

	//Compute the first flux function
	LinearAlgebra::MiddleSizeVector fluxFirst = computeLLFFluxFunction(faceStateStructLeft, faceStateStructRight.getState(),
																		faceStateStructLeft.getPressure(), faceStateStructRight.getPressure(), face.getUnitNormalVector());
	//First two for loops compute the derivative of the flux function with respect to iVB2
	//Second two loops computes the integrand matrix
	for (std::size_t iV2 = 0; iV2 < NUMBER_OF_VARIABLES; iV2++) // This is the variable in the coefficient vector
	{
		for (std::size_t iB2 = 0; iB2 < numberOfBasisFunctionsDerivativeSide; iB2++) // This is the basis function in the coefficient vector
		{
			//Location of the solutionCoefficient that the derivative is taken to
			iVB2 = face.getPhysicalElement(derivativeSide).convertToSingleIndex(iB2,iV2);
				//Compute new perturbed state and flux
			LinearAlgebra::MiddleSizeVector solutionCoefficientsPerturb;
			LinearAlgebra::MiddleSizeVector fluxSecond;
			LinearAlgebra::MiddleSizeVector stateLeftNew, stateRightNew;
			if (derivativeSide == Base::Side::LEFT)
			{
				solutionCoefficientsPerturb = faceStateStructLeft.getStateCoefficients();
				solutionCoefficientsPerturb(iVB2) += EPSILON;
				stateLeftNew = computeStateOnFace<DIM,NUMBER_OF_VARIABLES>(face, Base::Side::LEFT, solutionCoefficientsPerturb);
				StateCoefficientsStruct<DIM,NUMBER_OF_VARIABLES> newState;
				double pressureLeftNew = newState.computePressure(stateLeftNew);
				fluxSecond  = computeLLFFluxFunction(stateLeftNew, faceStateStructRight.getState(), pressureLeftNew, faceStateStructRight.getPressure(), face.getUnitNormalVector());
			}
			else
			{
				solutionCoefficientsPerturb = faceStateStructRight.getStateCoefficients();
				solutionCoefficientsPerturb(iVB2) += EPSILON;
				stateRightNew = computeStateOnFace<DIM,NUMBER_OF_VARIABLES>(face, Base::Side::RIGHT, solutionCoefficientsPerturb);
				double pressureRightNew = newState.computePressure(stateRightNew);
				fluxSecond  = computeLLFFluxFunction(faceStateStructRight.getState(), stateRightNew, faceStateStructLeft.getPressure(), pressureRightNew, face.getUnitNormalVector());
			}

			//Compute the difference
			LinearAlgebra::MiddleSizeVector differenceFlux = (fluxSecond-fluxFirst)/EPSILON;
			//Now compute the normal integrand, with the correct derivative flux function
			for(std::size_t iV1 = 0; iV1 < NUMBER_OF_VARIABLES; iV1++) // This is the variable in the residual vector
			{
				for (std::size_t iB1 = 0; iB1 < numberOfBasisFunctionsElementSide; iB1++) // This is the basis function in the residual vector
				{
					iVB1 = face.getPhysicalElement(elementSide).convertToSingleIndex(iB1,iV1);
					integrand(iVB1,iVB2) = sign*differenceFlux(iV1)*face.basisFunction(elementSide, iB1);
				}
			}
		}
	}
	return integrand;
}

