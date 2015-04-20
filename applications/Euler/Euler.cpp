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

#include "Euler.h"
#include <iomanip>
#include <cmath>

Euler::Euler
(
 const std::size_t dimension,
 const std::size_t numOfVariables,
 const double endTime,
 const std::size_t polynomialOrder,
 const Base::ButcherTableau * const ptrButcherTableau
) :
HpgemAPISimplified(dimension, numOfVariables, polynomialOrder, ptrButcherTableau),
DIM_(dimension),
numOfVariables_(numOfVariables)
{
}

/// \brief General mesh description
Base::RectangularMeshDescriptor Euler::createMeshDescription(const std::size_t numOfElementPerDirection)
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

/// *****************************************
/// ***   Element integration functions   ***
/// *****************************************

///  \brief computes the initial solution at an element
LinearAlgebra::NumericalVector Euler::computeSolutionAtElement(const Base::Element *ptrElement, const LinearAlgebra::NumericalVector &solutionCoefficients, const Geometry::PointReference &pRef)
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

/// \brief computes the source at an element
LinearAlgebra::NumericalVector Euler::integrandSourceAtElement(const Base::Element *ptrElement, const LinearAlgebra::NumericalVector qSolution, const double pressureTerm, const double &time, const Geometry::PointReference &pRef)
{
	std::size_t numOfBasisFunctions = ptrElement->getNrOfBasisFunctions();
	std::size_t iVB;

	LinearAlgebra::NumericalVector integrandSource(numOfVariables_ * numOfBasisFunctions);

	//Convert pRef to pPhys
	Geometry::PointPhysical pPhys = ptrElement->referenceToPhysical(pRef);


	//*********************************************************
	//***	Calculate derivative terms for source function	***
	//*********************************************************

	//Calculate base source functions: S_t, S_x, S_y, S_z, see manual
	double amplitude = 0.2;
	double frequency = 2.0*M_PI;
	LinearAlgebra::NumericalVector sourceValue(DIM_ + 1);

	//Add the time-dependent part of the function
	sourceValue(0) = -amplitude*frequency*std::sin(frequency*time);
	for (std::size_t iD = 0; iD < DIM_; iD++)
	{
		sourceValue(iD + 1) = -amplitude*frequency*std::cos(frequency*time);
	}

	//Add the space-dependent part of the function
	for (std::size_t iD = 0; iD < DIM_; iD++)
	{
		sourceValue(0) *= std::cos(frequency*pPhys[iD]);
		sourceValue(iD+1) *= std::sin(frequency*pPhys[iD]);
	}
	for (std::size_t iD = 0; iD < DIM_; iD++)
	{
		for (std::size_t iD2 = 0; iD2 < DIM_; iD2++)
		{
			if (iD != iD2)
			{
				sourceValue(iD+1) *= cos(frequency*pPhys[iD2]);
			}
		}
	}

	//*****************************************************
	//***	Calculate values for various source terms	***
	//*****************************************************

	double q1Inverse = 1.0/qSolution(0);

	//Calculate Source term values: momentum convection
	LinearAlgebra::Matrix sourceConvection(DIM_,DIM_); // d(rho*u^2)/dx or d(rho*v^2)/dy or d(rho*w^2)/dz on the diagonal and terms like d(rho*u*v)/dx and d(rho*w*u)/dz	on the off-diagonal

	for (std::size_t iD = 0; iD < DIM_; iD++)
	{
		sourceConvection(iD,iD) = (2.0*qSolution(iD+1) - qSolution(iD+1)*qSolution(iD+1)*q1Inverse)*sourceValue(iD+1)*q1Inverse; // d(rho*u^2)/dx or d(rho*v^2)/dy or d(rho*w^2)/dz
	}
	//off diagonal convection
	for (std::size_t iD = 0; iD < DIM_; iD++)
	{
		for (std::size_t iD2 = 0; iD2 < DIM_; iD2++)
		{
			if (iD!=iD2)
			{
				sourceConvection(iD,iD2) = (qSolution(iD+1) + qSolution(iD2+1) - qSolution(iD+1)*qSolution(iD2+1)*q1Inverse)*sourceValue(iD2+1)*q1Inverse; // terms like d(rho*u*v)/dx and d(rho*w*u)/dz
			}
		}
	}

	//Calculate Source term values: pressure
	LinearAlgebra::NumericalVector sourcePressure(DIM_); // dp/dx or dp/dy or dp/dz;
	double kineticPressure = 0.0; // This is the kinetic part of the pressure term;

	for (std::size_t iD = 0; iD < DIM_; iD++)
	{
		kineticPressure += (-qSolution(iD+1) + 0.5*qSolution(iD+1)*qSolution(iD+1)*q1Inverse)*q1Inverse;// part 1 of calculation
	}
	for (std::size_t iD = 0; iD < DIM_; iD++)
	{
		sourcePressure(iD) = (gamma_ -1)*(1 + kineticPressure)*sourceValue(iD+1); // dp/dx or dp/dy or dp/dz
	}

	//Calculate Source term values: Enthalpy convection
	LinearAlgebra::NumericalVector sourceEnthalpy(DIM_);
	for (std::size_t iD = 0; iD < DIM_; iD++)
	{
		sourceEnthalpy(iD) = (qSolution(iD+1) + qSolution(DIM_+1) - qSolution(iD+1)*qSolution(DIM_+1)*q1Inverse + pressureTerm - qSolution(iD+1)*pressureTerm*q1Inverse)*sourceValue(iD+1)*q1Inverse + qSolution(iD+1)*q1Inverse*sourcePressure(iD); // d(rho*u*h)/dx or d(rho*v*h)/dy or d(rho*w*h)/dz
	}

	//*************************************************************************
	//***	Calculate the complete source functions used for integration	***
	//*************************************************************************
	double sDensity;
	LinearAlgebra::NumericalVector sMomentum(DIM_);
	double sEnergy;

	// Add time derivative term
	sDensity = sourceValue(0);
	for (std::size_t iD = 0; iD < DIM_; iD++)
	{
		sMomentum(iD) = sourceValue(0);
	}
	sEnergy = sourceValue(0);

	// Add space terms
	for (std::size_t iD = 0; iD < DIM_ ; iD++)
	{
		sDensity += sourceValue(iD+1);
		for (std::size_t iD2 = 0; iD2 < DIM_; iD2++)
		{
			sMomentum(iD) += sourceConvection(iD,iD2);
		}
		sMomentum(iD) += sourcePressure(iD);
		sEnergy += sourceEnthalpy(iD);
	}



	//*********************************************************
	//*** Calculate the integrand of the Source integral	***
	//*********************************************************

	for(std::size_t iB = 0; iB < numOfBasisFunctions; iB++)
	{
		// Density
		iVB = ptrElement->convertToSingleIndex(iB,0);
		integrandSource(iVB) = sDensity*ptrElement->basisFunction(iB, pRef);

		// Momentum
		for (std::size_t iD = 0; iD < DIM_; iD++)
		{
			iVB = ptrElement->convertToSingleIndex(iB,iD+1);
			integrandSource(iVB) = sMomentum(iD)*ptrElement->basisFunction(iB, pRef);
		}

		// Energy
		iVB = ptrElement->convertToSingleIndex(iB,DIM_+1);
		integrandSource(iVB) = sEnergy*ptrElement->basisFunction(iB, pRef);
	}

	return integrandSource;
}

//NOTE: this function says RefElement, but it is on a physical element
LinearAlgebra::NumericalVector Euler::integrandRightHandSideOnRefElement(const Base::Element *ptrElement, const double &time, const Geometry::PointReference &pRef, const LinearAlgebra::NumericalVector &solutionCoefficients)
{
	// Get the number of basis functions in an element.
	std::size_t numOfBasisFunctions =  ptrElement->getNrOfBasisFunctions();

	//Create data structures for calculating the integrand
	LinearAlgebra::NumericalVector integrand(numOfVariables_ * numOfBasisFunctions); //The final integrand value will be stored in this vector
	LinearAlgebra::NumericalVector qSolution = computeSolutionAtElement(ptrElement, solutionCoefficients, pRef);
	LinearAlgebra::NumericalVector gradientBasisFunction(DIM_); //Gradient function based on the number of dimensions

	//Create temporary result values
	double integrandTerm = 0.0;
	double pressureTerm = 0.0;

	//Create iteration values
	std::size_t iVB; // Index in solution coefficients for variable i and basisfunction j

	//Compute pressure term
	double q1Inverse = 1.0/qSolution(0);
	for(std::size_t iD = 0; iD < DIM_; iD++)
	{
		pressureTerm += qSolution(iD+1)*qSolution(iD+1); // (u^2 + v^2 + w^2)*rho^2
	}
	pressureTerm = (gamma_ -1)*(qSolution(DIM_+1) - 0.5*q1Inverse*(pressureTerm)); // (gamma-1)*rho*(e- (u^2 + v^2 + w^2)/2)

	logger.assert(pressureTerm > 0, "Negative pressure.");

	// Compute the integrand for all equations
	for (std::size_t iB = 0; iB < numOfBasisFunctions; iB++) // For every basis function
	{
		gradientBasisFunction = ptrElement->basisFunctionDeriv(iB, pRef); // Compute the gradient of the ith test function

		for (std::size_t iD = 0; iD < DIM_; iD++) // For every dimension
		{

			//Calculate convection terms
			integrandTerm = qSolution(iD+1)*gradientBasisFunction(iD)*q1Inverse; // (u)*d(phi_iB)/dx or (v)*d(phi_iB)/dy

			//Density integrand
			iVB = ptrElement->convertToSingleIndex(iB,0);
			integrand(iVB) += qSolution(iD+1)*gradientBasisFunction(iD);//(iD); // rho*u*d(phi)/dx + rho*v*d(phi)/dy


			for (std::size_t jD = 0; jD < DIM_; jD++)
			{
				//Momentum integrand
				iVB = ptrElement->convertToSingleIndex(iB,jD+1);
				integrand(iVB) += qSolution(jD+1)*integrandTerm; //rho*u*u*d(phi)/dx + rho*u*v*d(phi)/dy
			}

			//Energy integrand
			iVB = ptrElement->convertToSingleIndex(iB,DIM_+1);
			integrand(iVB) += (qSolution(DIM_+1) + pressureTerm)*integrandTerm; // rho*u*h*d(phi_iB)/dx or rho*v*h*d(phi_iB)/dy

			//Calculate pressure terms in momentum equations
			iVB = ptrElement->convertToSingleIndex(iB,iD+1);
			integrand(iVB) += pressureTerm*gradientBasisFunction(iD);
		}

	}

	//Compute source terms
	LinearAlgebra::NumericalVector integrandSource = integrandSourceAtElement(ptrElement, qSolution, pressureTerm, time, pRef);

    return integrand + integrandSource;

}

LinearAlgebra::NumericalVector Euler::computeRightHandSideAtElement(Base::Element *ptrElement, LinearAlgebra::NumericalVector &solutionCoefficients, const double time)
{
	std::function<LinearAlgebra::NumericalVector(const Base::Element*, const Geometry::PointReference &)> integrandFunction = [&](const Base::Element *El, const Geometry::PointReference & pRef) -> LinearAlgebra::NumericalVector
	    {   return this->integrandRightHandSideOnRefElement(El, time, pRef, solutionCoefficients);};

    return elementIntegrator_.integrate(ptrElement, integrandFunction, ptrElement->getGaussQuadratureRule());

}

/// *****************************************
/// ***    face integration functions     ***
/// *****************************************

	// todo: Write function that computes the qSolution at the interface

   //Compute the Roe Riemann Flux function
   LinearAlgebra::NumericalVector Euler::RoeRiemannFluxFunction(const LinearAlgebra::NumericalVector &qSolutionLeft, const LinearAlgebra::NumericalVector &qSolutionRight, LinearAlgebra::NumericalVector &normal)
   {

	   //Compute correct normal direction and difference vector
	   double area = Base::L2Norm(normal);
	   normal = normal/area;

	   LinearAlgebra::NumericalVector qDifference = qSolutionRight - qSolutionLeft;

	   //Compute the Roe average state
	   LinearAlgebra::NumericalVector qAverage(DIM_+1);
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

	   for (std::size_t iD = 0; iD < DIM_; iD++)
	   {
		   qAverage(iD) = qSolutionLeft(iD+1)*tmp1 + qSolutionRight(iD+1)*tmp2; // u_average
		   ruSquaredLeft += qSolutionLeft(iD+1)*qSolutionLeft(iD+1); 			// Kinetic part of the left pressure term
		   ruSquaredRight += qSolutionRight(iD+1)*qSolutionRight(iD+1);			// Kinetic part of the right pressure term
	   }

	   pressureLeft = (gamma_ - 1)*(qSolutionLeft(DIM_ + 1) - 0.5*ruSquaredLeft*rhoInverseLeft);
	   pressureRight = (gamma_ - 1)*(qSolutionRight(DIM_ + 1) - 0.5*ruSquaredRight*rhoInverseRight);

	   qAverage(DIM_) = (qSolutionLeft(DIM_+1) + pressureLeft)*tmp1 + (qSolutionRight(DIM_+1) + pressureRight)*tmp2; //H_average

	   //Compute useful variables for constructing the flux function
	   double alphaAvg = 0.0;
	   double unAvg = 0.0;
	   for (std::size_t iD = 0; iD < DIM_; iD++)
	   {
		   alphaAvg += (qAverage(iD)*qAverage(iD));
		   unAvg += qAverage(iD)*normal(iD);
	   }
	   alphaAvg *= 0.5;

	   const double a2Avg = std::abs((gamma_ -1)*(qAverage(DIM_) - alphaAvg));
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
	   for (std::size_t iD = 0; iD < DIM_; iD++)
	   {
		   abv4 += -qAverage(iD)*qDifference(iD+1);
		   abv5 += normal(iD)*qDifference(iD+1);
	   }
	   abv4 += qDifference(DIM_+1);
	   abv4 *= (gamma_ -1);

	   const double abv6 = abv3*abv4*ova2Avg + abv2*abv5*ovaAvg;
	   const double abv7 = abv2*abv4*ovaAvg + abv3*abv5;

	   //Compute the Roe Riemann Flux function: h(u_L,u_R) = 0.5*(F(u)_L + F(u)_R - |A|(u_R - u_L))
	   LinearAlgebra::NumericalVector flux(DIM_ + 2);

	    double pLR = pressureLeft + pressureRight;

	    double runR = 0.0;
	    double runL = 0.0;
	    for (std::size_t iD = 0; iD < DIM_; iD++)
	    {
	    	runL += qSolutionLeft(iD+1)*normal(iD);
	    	runR += qSolutionRight(iD+1)*normal(iD);
	    }

	    double unL = runL*rhoInverseLeft;
	    double unR = runR*rhoInverseRight;

	    //continuity equation
	    flux(0) = (runL + runR - (lam3*qDifference(0) + abv6));

	    //momentum equations
	    for (std::size_t iD = 0; iD < DIM_; iD++)
	    {
	    	flux(iD+1) = runL*qSolutionLeft(iD+1)*rhoInverseLeft + runR*qSolutionRight(iD+1)*rhoInverseRight + pLR*normal(iD) - (lam3*qDifference(iD+1) + qAverage(iD)*abv6 + normal(iD)*abv7);
	    }

	    //energy equation
	    flux(DIM_+1) = (unL*(qSolutionLeft(DIM_+1)+pressureLeft) + unR*(qSolutionRight(DIM_+1)+pressureRight) - (lam3*qDifference(DIM_+1) + qAverage(DIM_)*abv6 + unAvg*abv7));

	    //Note: Twice the flux is computed above, hence the factor 0.5 in front of the equation
	    //Note: Correction is made to the flux since F*n was computed above, where n is the normal unit vector
	    //However the face integral function does not use a normalised vector.
        return 0.5*flux*area;
   }

   /// \brief Compute the integrand for the right hand side for the reference face corresponding to an external face.
   LinearAlgebra::NumericalVector Euler::integrandRightHandSideOnRefFace(const Base::Face *ptrFace, const double &time, const Geometry::PointReference &pRef, const LinearAlgebra::NumericalVector &solutionCoefficients)
   {
	   //Get the number of basis functions
	   std::size_t numOfBasisFunctionsLeft= ptrFace->getPtrElementLeft()->getNrOfBasisFunctions(); //Get the number of basis functions on the left

	   LinearAlgebra::NumericalVector integrand(numOfVariables_*numOfBasisFunctionsLeft);
	   LinearAlgebra::NumericalVector qReconstructionLeft(numOfVariables_);
	   LinearAlgebra::NumericalVector qReconstructionRight(numOfVariables_);

	   //Compute left and right states
	    std::size_t jVB; // Index for both variable and basis function.
	    for (std::size_t jV = 0; jV < numOfVariables_; jV++)
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
	   for (std::size_t iD = 0; iD < DIM_; iD++)
	   {
		   qReconstructionRight(iD+1) = -qReconstructionLeft(iD+1);
	   }
	   qReconstructionRight(DIM_+1) = qReconstructionLeft(DIM_+1);

	   // Compute normal vector, with size of the ref-to-phys face scale, pointing outward of the left element.
	   LinearAlgebra::NumericalVector normal = ptrFace->getNormalVector(pRef);


	   //Compute flux
	   LinearAlgebra::NumericalVector flux = RoeRiemannFluxFunction(qReconstructionRight, qReconstructionLeft, normal);

	   // Compute integrand on the reference element.
	   std::size_t iVB; // Index for both variable and basis function.
	   for (std::size_t iB = 0; iB < numOfBasisFunctionsLeft; iB++)
	   {
	       for (std::size_t iV = 0; iV < numOfVariables_; iV++) // Index for direction
	       {
	           iVB = ptrFace->getPtrElementLeft()->convertToSingleIndex(iB, iV);
	           integrand(iVB) = -flux(iV)*ptrFace->basisFunction(Base::Side::LEFT, iB, pRef);
	       }
	   }

	   return  integrand;
   }

   /// \brief Compute the integrand for the right hand side for the reference face corresponding to an internal face.
   LinearAlgebra::NumericalVector Euler::integrandRightHandSideOnRefFace(const Base::Face *ptrFace, const double &time, const Geometry::PointReference &pRef, const Base::Side &iSide, const LinearAlgebra::NumericalVector &solutionCoefficientsLeft, const LinearAlgebra::NumericalVector &solutionCoefficientsRight)
   {
	   //Get the number of basis functions
	   std::size_t numOfTestBasisFunctions = ptrFace->getPtrElement(iSide)->getNrOfBasisFunctions(); // Get the number of test basis functions on a given side, iSide
	   std::size_t numOfSolutionBasisFunctionsLeft = ptrFace->getPtrElementLeft()->getNrOfBasisFunctions(); //Get the number of basis functions on the left
	   std::size_t numOfSolutionBasisFunctionsRight = ptrFace->getPtrElementRight()->getNrOfBasisFunctions(); //Get the number of basis functions on the right side


	   LinearAlgebra::NumericalVector integrand(numOfVariables_*numOfTestBasisFunctions);
	   LinearAlgebra::NumericalVector qReconstructionLeft(numOfVariables_);
	   LinearAlgebra::NumericalVector qReconstructionRight(numOfVariables_);

	   // /todo: Remove this to a seperate function to reduce code duplication
	   //Compute left and right states
	    std::size_t jVB; // Index for both variable and basis function.
	    for (std::size_t jV = 0; jV < numOfVariables_; jV++)
	    {
	        qReconstructionLeft(jV) = 0;
	        qReconstructionRight(jV) = 0;
	        for (std::size_t jB = 0; jB < numOfSolutionBasisFunctionsLeft; jB++)
	        {
	            jVB = ptrFace->getPtrElementLeft()->convertToSingleIndex(jB, jV);
	            qReconstructionLeft(jV) += ptrFace->basisFunction(Base::Side::LEFT, jB, pRef) * solutionCoefficientsLeft(jVB);
	        }
	        for (std::size_t jB = 0; jB < numOfSolutionBasisFunctionsRight; jB++)
	        {
	            jVB = ptrFace->getPtrElementRight()->convertToSingleIndex(jB, jV);
	            qReconstructionRight(jV) += ptrFace->basisFunction(Base::Side::RIGHT, jB, pRef) * solutionCoefficientsRight(jVB);
	        }
	    }

	   // Compute normal vector, with size of the ref-to-phys face scale, pointing outward of the left element.
	   LinearAlgebra::NumericalVector normal = ptrFace->getNormalVector(pRef);

	   //Compute flux
	   LinearAlgebra::NumericalVector flux;
	   if (iSide == Base::Side::RIGHT)
	   {
		   flux = -RoeRiemannFluxFunction(qReconstructionLeft, qReconstructionRight, normal);
	   }
	   else
	   {
		   flux = RoeRiemannFluxFunction(qReconstructionLeft, qReconstructionRight, normal);
	   }

	   // Compute integrand on the reference element.
	   std::size_t iVB; // Index for both variable and basis function.
	   for (std::size_t iB = 0; iB < numOfTestBasisFunctions; iB++)
	   {
	       for (std::size_t iV = 0; iV < numOfVariables_; iV++) // Index for direction
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

   /// \brief Compute the right-hand side corresponding to a boundary face
   LinearAlgebra::NumericalVector Euler::computeRightHandSideAtFace(Base::Face *ptrFace, LinearAlgebra::NumericalVector &solutionCoefficients, const double time)
   {
	    // Define the integrand function for the right hand side for the reference face.
	    std::function<LinearAlgebra::NumericalVector(const Geometry::PointReference &)> integrandFunction = [&](const Geometry::PointReference &pRef) -> LinearAlgebra::NumericalVector
	    {   return this->integrandRightHandSideOnRefFace(ptrFace, time, pRef, solutionCoefficients);};

	    return faceIntegrator_.referenceFaceIntegral(ptrFace->getGaussQuadratureRule(), integrandFunction);
   }

   /// \brief Compute the right-hand side corresponding to an internal face
   LinearAlgebra::NumericalVector Euler::computeRightHandSideAtFace(Base::Face *ptrFace, const Base::Side side, LinearAlgebra::NumericalVector &solutionCoefficientsLeft, LinearAlgebra::NumericalVector &solutionCoefficientsRight, const double time)
   {
	    // Define the integrand function for the right hand side for the reference face.
	    std::function<LinearAlgebra::NumericalVector(const Geometry::PointReference &)> integrandFunction = [&](const Geometry::PointReference &pRef) -> LinearAlgebra::NumericalVector
	    {   return this->integrandRightHandSideOnRefFace(ptrFace, time, pRef, side, solutionCoefficientsLeft, solutionCoefficientsRight);};
	    return faceIntegrator_.referenceFaceIntegral(ptrFace->getGaussQuadratureRule(), integrandFunction);
   }


   /// *****************************************
   /// ***    		Various Functions        ***
   /// *****************************************

   LinearAlgebra::NumericalVector Euler::getExactSolution(const PointPhysicalT &pPhys, const double &time, const std::size_t orderTimeDerivative)
   {
		LinearAlgebra::NumericalVector exactSolution(numOfVariables_);

		double amplitude = 0.2;
		double frequency = 2.0*M_PI;
		double function = amplitude*std::cos(frequency*time);

		for (std::size_t iD = 0; iD < DIM_; iD++)
		{
			function *= std::cos(frequency*pPhys[iD]);
		}

		exactSolution(0) = 1.5 + function;

		for (std::size_t iD = 0; iD < DIM_; iD++)
		{
			exactSolution(iD+1) = function;
		}

		exactSolution(DIM_ + 1) = 30.0 + function;

	    return exactSolution;
   }

   /// \brief Compute the initial solution at a given point in space and time.
   LinearAlgebra::NumericalVector Euler::getInitialSolution(const PointPhysicalT &pPhys, const double &startTime, const std::size_t orderTimeDerivative)
   {
       return getExactSolution(pPhys, startTime, orderTimeDerivative);
   }

   /// \brief Computes the error for output purposes
   LinearAlgebra::NumericalVector Euler::Error(const double time)
   {
	   return computeMaxError(solutionTimeLevel_, time);
   }
















