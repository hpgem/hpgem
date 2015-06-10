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

LinearAlgebra::MiddleSizeVector Inviscid::integrandAtElement(Base::PhysicalElement<DIM>& element, const double &time, const double pressureTerm, const LinearAlgebra::MiddleSizeVector &qSolution)
{
    const Base::Element* ptrElement = element.getElement();
    const Geometry::PointReference<DIM>& pRef = element.getPointReference();
	// Get the number of basis functions in an element.
	std::size_t numOfBasisFunctions =  ptrElement->getNrOfBasisFunctions();

	//Create data structures for calculating the integrand
	LinearAlgebra::MiddleSizeVector integrand(instance_.numOfVariables_ * numOfBasisFunctions); //The final integrand value will be stored in this vector
	LinearAlgebra::MiddleSizeVector gradientBasisFunction(instance_.DIM_); //Gradient function based on the number of dimensions

	//Create temporary result values
	double integrandTerm = 0.0;
	double q1Inverse = 1.0/qSolution(0);

	//Create iteration values
	std::size_t iVB; // Index in solution coefficients for variable i and basisfunction j

	// Compute the integrand for all equations
	for (std::size_t iB = 0; iB < numOfBasisFunctions; iB++) // For every basis function
	{
		gradientBasisFunction = ptrElement->basisFunctionDeriv(iB, pRef); // Compute the gradient of the ith test function

		for (std::size_t iD = 0; iD < instance_.DIM_; iD++) // For every dimension
		{

			//Calculate convection terms
			integrandTerm = qSolution(iD+1)*gradientBasisFunction(iD)*q1Inverse; // (u)*d(phi_iB)/dx or (v)*d(phi_iB)/dy

			//Density integrand
			iVB = ptrElement->convertToSingleIndex(iB,0);
			integrand(iVB) += qSolution(iD+1)*gradientBasisFunction(iD);//(iD); // rho*u*d(phi)/dx + rho*v*d(phi)/dy


			for (std::size_t jD = 0; jD < instance_.DIM_; jD++)
			{
				//Momentum integrand
				iVB = ptrElement->convertToSingleIndex(iB,jD+1);
				integrand(iVB) += qSolution(jD+1)*integrandTerm; //rho*u*u*d(phi)/dx + rho*u*v*d(phi)/dy
			}

			//Energy integrand
			iVB = ptrElement->convertToSingleIndex(iB,instance_.DIM_+1);
			integrand(iVB) += (qSolution(instance_.DIM_+1) + pressureTerm)*integrandTerm; // rho*u*h*d(phi_iB)/dx or rho*v*h*d(phi_iB)/dy

			//Calculate pressure terms in momentum equations
			iVB = ptrElement->convertToSingleIndex(iB,iD+1);
			integrand(iVB) += pressureTerm*gradientBasisFunction(iD);
		}

	}

    return integrand;
}

//Compute the Roe Riemann Flux function
LinearAlgebra::MiddleSizeVector Inviscid::RoeRiemannFluxFunction(const LinearAlgebra::MiddleSizeVector &qSolutionLeft, const LinearAlgebra::MiddleSizeVector &qSolutionRight, const LinearAlgebra::SmallVector<DIM> &normal)
{

	   LinearAlgebra::MiddleSizeVector qDifference = qSolutionRight - qSolutionLeft;

	   //Compute the Roe average state
	   LinearAlgebra::MiddleSizeVector qAverage(instance_.DIM_+1);
	   double zL = std::sqrt(qSolutionLeft(0));
	   double zR = std::sqrt(qSolutionRight(0));
	   double tmp1 = 1.0/(qSolutionLeft(0) + zL*zR);
	   double tmp2 = 1.0/(qSolutionRight(0) + zL*zR);
	   double ruSquaredLeft = 0.0;
	   double ruSquaredRight = 0.0;
	   double rhoInverseLeft = 1.0/qSolutionLeft(0);
	   double rhoInverseRight = 1.0/qSolutionRight(0);
	   double pressureLeft;
	   double pressureRight;

	   for (std::size_t iD = 0; iD < instance_.DIM_; iD++)
	   {
		   qAverage(iD) = qSolutionLeft(iD+1)*tmp1 + qSolutionRight(iD+1)*tmp2; // u_average
		   ruSquaredLeft += qSolutionLeft(iD+1)*qSolutionLeft(iD+1); 			// Kinetic part of the left pressure term
		   ruSquaredRight += qSolutionRight(iD+1)*qSolutionRight(iD+1);			// Kinetic part of the right pressure term
	   }

	   pressureLeft = (instance_.gamma_ - 1)*(qSolutionLeft(instance_.DIM_ + 1) - 0.5*ruSquaredLeft*rhoInverseLeft);
	   pressureRight = (instance_.gamma_ - 1)*(qSolutionRight(instance_.DIM_ + 1) - 0.5*ruSquaredRight*rhoInverseRight);

	   qAverage(instance_.DIM_) = (qSolutionLeft(instance_.DIM_+1) + pressureLeft)*tmp1 + (qSolutionRight(instance_.DIM_+1) + pressureRight)*tmp2; //H_average

	   //Compute useful variables for constructing the flux function
	   double alphaAvg = 0.0;
	   double unAvg = 0.0;
	   for (std::size_t iD = 0; iD < instance_.DIM_; iD++)
	   {
		   alphaAvg += (qAverage(iD)*qAverage(iD));
		   unAvg += qAverage(iD)*normal(iD);
	   }
	   alphaAvg *= 0.5;

	   const double a2Avg = std::abs((instance_.gamma_ -1)*(qAverage(instance_.DIM_) - alphaAvg));
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

	   double abv4 = alphaAvg*qDifference(0);
	   double abv5 = -unAvg*qDifference(0);
	   for (std::size_t iD = 0; iD < instance_.DIM_; iD++)
	   {
		   abv4 += -qAverage(iD)*qDifference(iD+1);
		   abv5 += normal(iD)*qDifference(iD+1);
	   }
	   abv4 += qDifference(instance_.DIM_+1);
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
	    	runL += qSolutionLeft(iD+1)*normal(iD);
	    	runR += qSolutionRight(iD+1)*normal(iD);
	    }

	    double unL = runL*rhoInverseLeft;
	    double unR = runR*rhoInverseRight;

	    //continuity equation
	    flux(0) = (runL + runR - (lam3*qDifference(0) + abv6));

	    //momentum equations
	    for (std::size_t iD = 0; iD < instance_.DIM_; iD++)
	    {
	    	flux(iD+1) = runL*qSolutionLeft(iD+1)*rhoInverseLeft + runR*qSolutionRight(iD+1)*rhoInverseRight + pLR*normal(iD) - (lam3*qDifference(iD+1) + qAverage(iD)*abv6 + normal(iD)*abv7);
	    }

	    //energy equation
	    flux(instance_.DIM_+1) = (unL*(qSolutionLeft(instance_.DIM_+1)+pressureLeft) + unR*(qSolutionRight(instance_.DIM_+1)+pressureRight) - (lam3*qDifference(instance_.DIM_+1) + qAverage(instance_.DIM_)*abv6 + unAvg*abv7));

	    //Note: Twice the flux is computed above, hence the factor 0.5 in front of the equation
     return 0.5*flux;
}

/// \brief Compute the integrand for the right hand side for the reference face corresponding to an external face.
LinearAlgebra::MiddleSizeVector Inviscid::integrandAtFace(Base::PhysicalFace<DIM>& face, const double &time, const LinearAlgebra::MiddleSizeVector &solutionCoefficients)
{
    const Base::Face* ptrFace = face.getFace();
    const Geometry::PointReference<DIM - 1>& pRef = face.getPointReference();
	   //Get the number of basis functions
	   std::size_t numOfBasisFunctionsLeft= ptrFace->getPtrElementLeft()->getNrOfBasisFunctions(); //Get the number of basis functions on the left

	   LinearAlgebra::MiddleSizeVector integrand(instance_.numOfVariables_*numOfBasisFunctionsLeft);
	   LinearAlgebra::MiddleSizeVector qReconstructionLeft(instance_.numOfVariables_);
	   LinearAlgebra::MiddleSizeVector qReconstructionRight(instance_.numOfVariables_);

	   //todo: this has been computed already in previous step: remove it to reduce code duplication
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
	    }

	   // Boundary face is assumed to be a solid wall: set reflective solution on the other side
	   qReconstructionRight(0) = qReconstructionLeft(0);
	   for (std::size_t iD = 0; iD < instance_.DIM_; iD++)
	   {
		   qReconstructionRight(iD+1) = -qReconstructionLeft(iD+1);
	   }
	   qReconstructionRight(instance_.DIM_+1) = qReconstructionLeft(instance_.DIM_+1);

	   // Compute normal vector, with size of the ref-to-phys face scale, pointing outward of the left element.
	   LinearAlgebra::MiddleSizeVector normal = ptrFace->getNormalVector(pRef);


	   //Compute flux
	   LinearAlgebra::MiddleSizeVector flux = RoeRiemannFluxFunction(qReconstructionRight, qReconstructionLeft, normal);

	   //todo: Other implementation is simply to put the flux to zero. No mass, momentum and energy should get out. This saves computational time

	   // Compute integrand on the reference element.
	   std::size_t iVB; // Index for both variable and basis function.
	   for (std::size_t iB = 0; iB < numOfBasisFunctionsLeft; iB++)
	   {
	       for (std::size_t iV = 0; iV < instance_.numOfVariables_; iV++) // Index for direction
	       {
	           iVB = ptrFace->getPtrElementLeft()->convertToSingleIndex(iB, iV);
	           integrand(iVB) = -flux(iV)*ptrFace->basisFunction(Base::Side::LEFT, iB, pRef);
	       }
	   }

	   return  integrand;
}

/// \brief Compute the integrand for the right hand side for the reference face corresponding to an internal face.
LinearAlgebra::MiddleSizeVector Inviscid::integrandAtFace(Base::PhysicalFace<DIM>& face, const double &time, const Base::Side &iSide, const LinearAlgebra::MiddleSizeVector &qSolutionInternal, const LinearAlgebra::MiddleSizeVector &qSolutionExternal)
{
    const Base::Face* ptrFace = face.getFace();
    const Geometry::PointReference<DIM - 1>& pRef = face.getPointReference();
    LinearAlgebra::SmallVector<DIM> normal = face.getUnitNormalVector();
	   //Get the number of basis functions
	   std::size_t numOfTestBasisFunctions = ptrFace->getPtrElement(iSide)->getNrOfBasisFunctions(); // Get the number of test basis functions on a given side, iSide

	   LinearAlgebra::MiddleSizeVector integrand(instance_.numOfVariables_*numOfTestBasisFunctions);

	   //Compute flux
	   LinearAlgebra::MiddleSizeVector flux = RoeRiemannFluxFunction(qSolutionInternal, qSolutionExternal, normal);

	   // Compute integrand on the reference element.
	   std::size_t iVB; // Index for both variable and basis function.
	   for (std::size_t iB = 0; iB < numOfTestBasisFunctions; iB++)
	   {
	       for (std::size_t iV = 0; iV < instance_.numOfVariables_; iV++) // Index for direction
	       {
	           iVB = ptrFace->getPtrElement(iSide)->convertToSingleIndex(iB, iV);
	           integrand(iVB) = -flux(iV)*ptrFace->basisFunction(iSide, iB, pRef);
	       }
	   }

	    if (integrand[0] != integrand[0])
	    {
	    	logger(ERROR,"integrand on ref face is nan");
	    }
	   return  integrand;
}

