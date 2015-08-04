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
#include <cmath>

Inviscid::Inviscid(const CompressibleNavierStokes& instance) : instance_(instance)
{
}

LinearAlgebra::MiddleSizeVector Inviscid::integrandAtElement(Base::PhysicalElement<DIM> &element, const double &time, const double &pressureTerm, const LinearAlgebra::MiddleSizeVector &state)
{
	// Get the number of basis functions in an element.
	std::size_t numberOfBasisFunctions =  element.getNumOfBasisFunctions();

	//Create data structures for calculating the integrand
	LinearAlgebra::MiddleSizeVector integrand(instance_.numOfVariables_ * numberOfBasisFunctions); //The final integrand value will be stored in this vector
	LinearAlgebra::SmallVector<DIM> gradientBasisFunction; //Gradient function based on the number of dimensions

	//Create temporary result values
	double integrandTerm = 0.0;
	double q1Inverse = 1.0/state(0);

	//Create iteration values
	std::size_t iVB; // Index in solution coefficients for variable i and basisfunction j

	// Compute the integrand for all equations
	for (std::size_t iB = 0; iB < numberOfBasisFunctions; iB++) // For every basis function
	{
		gradientBasisFunction = element.basisFunctionDeriv(iB); // Compute the gradient of the ith test function

		for (std::size_t iD = 0; iD < instance_.DIM_; iD++) // For every dimension
		{

			//Calculate convection terms
			integrandTerm = state(iD+1)*gradientBasisFunction(iD)*q1Inverse; // (u)*d(phi_iB)/dx or (v)*d(phi_iB)/dy

			//Density integrand
			iVB = element.convertToSingleIndex(iB,0);
			integrand(iVB) += state(iD+1)*gradientBasisFunction(iD);//(iD); // rho*u*d(phi)/dx + rho*v*d(phi)/dy


			for (std::size_t jD = 0; jD < instance_.DIM_; jD++)
			{
				//Momentum integrand
				iVB = element.convertToSingleIndex(iB,jD+1);
				integrand(iVB) += state(jD+1)*integrandTerm; //rho*u*u*d(phi)/dx + rho*u*v*d(phi)/dy
			}

			//Energy integrand
			iVB = element.convertToSingleIndex(iB,instance_.DIM_+1);
			integrand(iVB) += (state(instance_.DIM_+1) + pressureTerm)*integrandTerm; // rho*u*h*d(phi_iB)/dx or rho*v*h*d(phi_iB)/dy

			//Calculate pressure terms in momentum equations
			iVB = element.convertToSingleIndex(iB,iD+1);
			integrand(iVB) += pressureTerm*gradientBasisFunction(iD);
		}

	}

    return integrand;
}

//Compute the Roe Riemann Flux function
LinearAlgebra::MiddleSizeVector Inviscid::RoeRiemannFluxFunction(const LinearAlgebra::MiddleSizeVector &stateLeft, const LinearAlgebra::MiddleSizeVector &stateRight, const LinearAlgebra::SmallVector<DIM> &unitNormalInternal)
{

	   //Compute correct normal direction and difference vector
	   LinearAlgebra::MiddleSizeVector stateDifference = stateRight - stateLeft;

	   //Compute the Roe average state
	   LinearAlgebra::MiddleSizeVector stateAverage(instance_.DIM_+1);
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


	   for (std::size_t iD = 0; iD < instance_.DIM_; iD++)
	   {
		   stateAverage(iD) = stateLeft(iD+1)*tmp1 + stateRight(iD+1)*tmp2; // u_average
		   ruSquaredLeft += stateLeft(iD+1)*stateLeft(iD+1); 			// Kinetic part of the left pressure term
		   ruSquaredRight += stateRight(iD+1)*stateRight(iD+1);			// Kinetic part of the right pressure term
	   }

	   //todo: function layer above contains this information already
	   pressureLeft = (instance_.gamma_ - 1)*(stateLeft(instance_.DIM_ + 1) - 0.5*ruSquaredLeft*rhoInverseLeft);
	   pressureRight = (instance_.gamma_ - 1)*(stateRight(instance_.DIM_ + 1) - 0.5*ruSquaredRight*rhoInverseRight);

 	   stateAverage(instance_.DIM_) = (stateLeft(instance_.DIM_+1) + pressureLeft)*tmp1 + (stateRight(instance_.DIM_+1) + pressureRight)*tmp2; //H_average

	   //Compute useful variables for constructing the flux function
	   double alphaAvg = 0.0;
	   double unAvg = 0.0;
	   for (std::size_t iD = 0; iD < instance_.DIM_; iD++)
	   {
		   alphaAvg += (stateAverage(iD)*stateAverage(iD));
		   unAvg += stateAverage(iD)*unitNormalInternal(iD);
	   }
	   alphaAvg *= 0.5;

	   const double a2Avg = std::abs((instance_.gamma_ -1)*(stateAverage(instance_.DIM_) - alphaAvg));
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
	   for (std::size_t iD = 0; iD < instance_.DIM_; iD++)
	   {
		   abv4 += -stateAverage(iD)*stateDifference(iD+1);
		   abv5 += unitNormalInternal(iD)*stateDifference(iD+1);
	   }
	   abv4 += stateDifference(instance_.DIM_+1);
	   abv4 *= (instance_.gamma_ -1);

	   const double abv6 = abv3*abv4*ova2Avg + abv2*abv5*ovaAvg;
	   const double abv7 = abv2*abv4*ovaAvg + abv3*abv5;

	   //Compute the Roe Riemann Flux function: h(u_L,u_R) = 0.5*(F(u_L) + F(u_R) - |A|(u_R - u_L))
	   LinearAlgebra::MiddleSizeVector flux(instance_.DIM_ + 2);

	    double pLR = pressureLeft + pressureRight;

	    double runR = 0.0;
	    double runL = 0.0;
	    for (std::size_t iD = 0; iD < instance_.DIM_; iD++)
	    {
	    	runL += stateLeft(iD+1)*unitNormalInternal(iD);
	    	runR += stateRight(iD+1)*unitNormalInternal(iD);
	    }

	    double unL = runL*rhoInverseLeft;
	    double unR = runR*rhoInverseRight;

	    //continuity equation
	    flux(0) = (runL + runR - (lam3*stateDifference(0) + abv6));

	    //momentum equations
	    for (std::size_t iD = 0; iD < instance_.DIM_; iD++)
	    {
	    	flux(iD+1) = runL*stateLeft(iD+1)*rhoInverseLeft + runR*stateRight(iD+1)*rhoInverseRight + pLR*unitNormalInternal(iD) - (lam3*stateDifference(iD+1) + stateAverage(iD)*abv6 + unitNormalInternal(iD)*abv7);
	    }

	    //energy equation
	    flux(instance_.DIM_+1) = (unL*(stateLeft(instance_.DIM_+1)+pressureLeft) + unR*(stateRight(instance_.DIM_+1)+pressureRight) - (lam3*stateDifference(instance_.DIM_+1) + stateAverage(instance_.DIM_)*abv6 + unAvg*abv7));

	    //Note: Twice the flux is computed above, hence the factor 0.5 in front of the equation
     return 0.5*flux;
}

/// \brief Compute the integrand for the right hand side for the reference face corresponding to an external face.
LinearAlgebra::MiddleSizeVector Inviscid::integrandAtFace(Base::PhysicalFace<DIM> &face, const double &time, const LinearAlgebra::MiddleSizeVector &stateInternal)
{
	   //Get the number of basis functions
	   std::size_t numOfBasisFunctionsLeft= face.getPhysicalElement(Base::Side::LEFT).getNumOfBasisFunctions();

	   LinearAlgebra::MiddleSizeVector integrand(instance_.numOfVariables_*numOfBasisFunctionsLeft);
	   LinearAlgebra::MiddleSizeVector stateExternal(instance_.numOfVariables_);

/*	   //Deprecated code:
	   //Compute left and right states
	    std::size_t jVB; // Index for both variable and basis function.
	    for (std::size_t jV = 0; jV < instance_.numOfVariables_; jV++)
	    {
	        qReconstructionLeft(jV) = 0;
	        for (std::size_t jB = 0; jB < numOfBasisFunctionsLeft; jB++)
	        {
	            jVB = ptrFace->getPtrElementLeft()->convertToSingleIndex(jB, jV);
	            qReconstructionLeft(jV) += ptrFace->basisFunction(Base::Side::LEFT, jB, pRef) * solutionCoefficients(jVB);
	        }
	    }*/

	   // Boundary face is assumed to be a solid wall: set reflective solution on the other side
	   stateExternal(0) = stateInternal(0);
	   for (std::size_t iD = 0; iD < instance_.DIM_; iD++)
	   {
		   stateExternal(iD+1) = -stateInternal(iD+1);
	   }
	   stateExternal(instance_.DIM_+1) = stateInternal(instance_.DIM_+1);

	   //WARNING: not a unit normal vector here
	   // Compute normal vector, with size of the ref-to-phys face scale, pointing outward of the left element.
	   LinearAlgebra::MiddleSizeVector normalInternal = face.getNormalVector();
	   LinearAlgebra::MiddleSizeVector unitNormalInternal = normalInternal/Base::L2Norm(normalInternal);

	   //Compute flux
	   //todo: check if this function is called correctly
	   LinearAlgebra::MiddleSizeVector flux = RoeRiemannFluxFunction(stateExternal, stateInternal, unitNormalInternal);

	   // Compute integrand on the reference element.
	   std::size_t iVB; // Index for both variable and basis function.
	   for (std::size_t iB = 0; iB < numOfBasisFunctionsLeft; iB++)
	   {
	       for (std::size_t iV = 0; iV < instance_.numOfVariables_; iV++) // Index for direction
	       {
	           iVB = face.getPhysicalElement(Base::Side::LEFT).convertToSingleIndex(iB, iV);
	           integrand(iVB) = -flux(iV)*face.basisFunction(Base::Side::LEFT, iB);
	       }
	   }
	   std::cout << "I should not be here" << std::endl;

	   return  integrand;
}

/// \brief Compute the integrand for the right hand side for the reference face corresponding to an internal face.
/*
LinearAlgebra::MiddleSizeVector Inviscid::integrandAtFace(
		Base::PhysicalFace<DIM> &face,
		const double &time,
		const Base::Side &iSide,
		const LinearAlgebra::MiddleSizeVector &stateInternal,
		const LinearAlgebra::MiddleSizeVector &stateExternal,
		const LinearAlgebra::SmallVector<DIM> &unitNormalInternal)
{

	   //Get the number of basis functions
	   std::size_t numOfTestBasisFunctions = face.getPhysicalElement(iSide).getNumOfBasisFunctions();

	   LinearAlgebra::MiddleSizeVector integrand(instance_.numOfVariables_*numOfTestBasisFunctions);

	   //Compute flux
	   LinearAlgebra::MiddleSizeVector flux = RoeRiemannFluxFunction(stateInternal, stateExternal, unitNormalInternal);

 	   // Compute integrand on the reference element.
	   std::size_t iVB; // Index for both variable and basis function.
	   for (std::size_t iB = 0; iB < numOfTestBasisFunctions; iB++)
	   {
	       for (std::size_t iV = 0; iV < instance_.numOfVariables_; iV++) // Index for direction
	       {
	           iVB = face.getPhysicalElement(iSide).convertToSingleIndex(iB, iV);
	           integrand(iVB) = flux(iV)*face.basisFunction(iSide, iB);
	       }
	   }

	   return  -integrand;
}*/

std::pair<LinearAlgebra::MiddleSizeVector,LinearAlgebra::MiddleSizeVector> Inviscid::integrandsAtFace(
		Base::PhysicalFace<DIM> &face,
		const double &time,
		const LinearAlgebra::MiddleSizeVector &stateLeft,
		const LinearAlgebra::MiddleSizeVector &stateRight)
{
	std::pair<LinearAlgebra::MiddleSizeVector,LinearAlgebra::MiddleSizeVector> Integrands;

	//Data structures for left and right integrand
	std::size_t numOfTestBasisFunctionsLeft = face.getPhysicalElement(Base::Side::LEFT).getNumOfBasisFunctions();
	LinearAlgebra::MiddleSizeVector integrandLeft(instance_.numOfVariables_*numOfTestBasisFunctionsLeft);

	std::size_t numOfTestBasisFunctionsRight = face.getPhysicalElement(Base::Side::RIGHT).getNumOfBasisFunctions();
	LinearAlgebra::MiddleSizeVector integrandRight(instance_.numOfVariables_*numOfTestBasisFunctionsRight);

	std::size_t iVB;

	//Compute left flux
	LinearAlgebra::SmallVector<DIM> unitNormalLeft = face.getUnitNormalVector();
	LinearAlgebra::MiddleSizeVector flux = RoeRiemannFluxFunction(stateLeft, stateRight, unitNormalLeft);

	//Compute left integrand
	for (std::size_t iB = 0; iB < numOfTestBasisFunctionsLeft; iB++)
	{
		for (std::size_t iV = 0; iV < instance_.numOfVariables_; iV++) // Index for direction
		{
			iVB = face.getPhysicalElement(Base::Side::LEFT).convertToSingleIndex(iB, iV);
			integrandLeft(iVB) = -flux(iV)*face.basisFunction(Base::Side::LEFT, iB); // Minus sign because the integral is on the right hand side
			integrandRight(iVB) = flux(iV)*face.basisFunction(Base::Side::RIGHT,iB);
		}
	}

	//Assign integrand values to pair
	Integrands.first = integrandLeft;
	Integrands.second = integrandRight;

	return Integrands;
}

