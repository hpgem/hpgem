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

template<std::size_t DIM>
Euler<DIM>::Euler
(
 const std::size_t numOfVariables,
 const double endTime,
 const std::size_t polynomialOrder,
 const Base::ButcherTableau * const ptrButcherTableau
) :
Base::HpgemAPISimplified<DIM>(numOfVariables, polynomialOrder, ptrButcherTableau),
numOfVariables_(numOfVariables)
{
}

/// \brief General mesh description
template<std::size_t DIM>
Base::RectangularMeshDescriptor<DIM> Euler<DIM>::createMeshDescription(const std::size_t numOfElementPerDirection)
{
    // Create the domain. In this case the domain is the square [0,1]^DIM and periodic.
    Base::RectangularMeshDescriptor<DIM> description;
    for (std::size_t i = 0; i < DIM; ++i)
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
template<std::size_t DIM>
LinearAlgebra::MiddleSizeVector Euler<DIM>::computeSolutionAtElement(const Base::Element *ptrElement, const LinearAlgebra::MiddleSizeVector &solutionCoefficients, const PointReferenceT &pRef)
{
		std::size_t numOfBasisFunctions =  ptrElement->getNrOfBasisFunctions();
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

/// \brief computes the source at an element
template<std::size_t DIM>
LinearAlgebra::MiddleSizeVector Euler<DIM>::integrandSourceAtElement(Base::PhysicalElement<DIM> &ptrElement, const LinearAlgebra::MiddleSizeVector qSolution, const double pressureTerm, const double &time)
{
	std::size_t numOfBasisFunctions = ptrElement.getElement()->getNrOfBasisFunctions();
	std::size_t iVB;

	//getResultVector already contains partial computation from integrandRightHandSideOnRefElement
	LinearAlgebra::MiddleSizeVector integrandSource(numOfVariables_ * numOfBasisFunctions);

	return integrandSource;
}

//NOTE: this function says RefElement, but it is on a physical element
template<std::size_t DIM>
LinearAlgebra::MiddleSizeVector Euler<DIM>::integrandRightHandSideOnRefElement(Base::PhysicalElement<DIM> &ptrElement, const double &time, const LinearAlgebra::MiddleSizeVector &solutionCoefficients)
{
	// Get the number of basis functions in an element.
	std::size_t numOfBasisFunctions =  ptrElement.getElement()->getNrOfBasisFunctions();

	//Create data structures for calculating the integrand
	LinearAlgebra::MiddleSizeVector& integrand = ptrElement.getResultVector();
	const LinearAlgebra::MiddleSizeVector& qSolution = ptrElement.getSolution();
	LinearAlgebra::SmallVector<DIM> gradientBasisFunction; //Gradient function based on the number of dimensions

	//Create temporary result values
	double integrandTerm = 0.0;
	double pressureTerm = 0.0;

	//Create iteration values
	std::size_t iVB; // Index in solution coefficients for variable i and basisfunction j

	//Compute pressure term
	double q1Inverse = 1.0/qSolution(0);
	for(std::size_t iD = 0; iD < DIM; iD++)
	{
		pressureTerm += qSolution(iD+1)*qSolution(iD+1); // (u^2 + v^2 + w^2)*rho^2
	}
	pressureTerm = (gamma_ -1)*(qSolution(DIM+1) - 0.5*q1Inverse*(pressureTerm)); // (gamma-1)*rho*(e- (u^2 + v^2 + w^2)/2)

	logger.assert(pressureTerm > 0, "Negative pressure.");

	// Compute the integrand for all equations
	for (std::size_t iB = 0; iB < numOfBasisFunctions; iB++) // For every basis function
	{
		gradientBasisFunction = ptrElement.basisFunctionDeriv(iB); // Compute the gradient of the ith test function

		for (std::size_t iD = 0; iD < DIM; iD++) // For every dimension
		{

			//Calculate convection terms
			integrandTerm = qSolution(iD+1)*gradientBasisFunction(iD)*q1Inverse; // (u)*d(phi_iB)/dx or (v)*d(phi_iB)/dy

			//Density integrand
			iVB = ptrElement.getElement()->convertToSingleIndex(iB,0);
			integrand(iVB) += qSolution(iD+1)*gradientBasisFunction(iD);//(iD); // rho*u*d(phi)/dx + rho*v*d(phi)/dy


			for (std::size_t jD = 0; jD < DIM; jD++)
			{
				//Momentum integrand
				iVB = ptrElement.getElement()->convertToSingleIndex(iB,jD+1);
				integrand(iVB) += qSolution(jD+1)*integrandTerm; //rho*u*u*d(phi)/dx + rho*u*v*d(phi)/dy
			}

			//Energy integrand
			iVB = ptrElement.getElement()->convertToSingleIndex(iB,DIM+1);
			integrand(iVB) += (qSolution(DIM+1) + pressureTerm)*integrandTerm; // rho*u*h*d(phi_iB)/dx or rho*v*h*d(phi_iB)/dy

			//Calculate pressure terms in momentum equations
			iVB = ptrElement.getElement()->convertToSingleIndex(iB,iD+1);
			integrand(iVB) += pressureTerm*gradientBasisFunction(iD);
		}

	}

	//Compute source terms
	LinearAlgebra::MiddleSizeVector integrandSource = integrandSourceAtElement(ptrElement, qSolution, pressureTerm, time);

    return integrand + integrandSource;

}

template<std::size_t DIM>
LinearAlgebra::MiddleSizeVector Euler<DIM>::computeRightHandSideAtElement(Base::Element *ptrElement, LinearAlgebra::MiddleSizeVector &solutionCoefficients, const double time)
{
	std::function<LinearAlgebra::MiddleSizeVector(Base::PhysicalElement<DIM>&)> integrandFunction = [&](Base::PhysicalElement<DIM>& El) -> LinearAlgebra::MiddleSizeVector
	    {   return this->integrandRightHandSideOnRefElement(El, time, solutionCoefficients);};

    return this->elementIntegrator_.integrate(ptrElement, integrandFunction, ptrElement->getGaussQuadratureRule());

}

/// *****************************************
/// ***    face integration functions     ***
/// *****************************************

	// todo: Write function that computes the qSolution at the interface

   //Compute the Roe Riemann Flux function
    template<std::size_t DIM>
   LinearAlgebra::MiddleSizeVector Euler<DIM>::RoeRiemannFluxFunction(const LinearAlgebra::MiddleSizeVector &qSolutionLeft, const LinearAlgebra::MiddleSizeVector &qSolutionRight, LinearAlgebra::SmallVector<DIM> &normal)
   {

	   //Compute correct normal direction and difference vector
	   double area = Base::L2Norm(normal);
	   normal = normal/area;

	   LinearAlgebra::MiddleSizeVector qDifference = qSolutionRight - qSolutionLeft;

	   //Compute the Roe average state
	   LinearAlgebra::MiddleSizeVector qAverage(DIM+1);
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

	   for (std::size_t iD = 0; iD < DIM; iD++)
	   {
		   qAverage(iD) = qSolutionLeft(iD+1)*tmp1 + qSolutionRight(iD+1)*tmp2; // u_average
		   ruSquaredLeft += qSolutionLeft(iD+1)*qSolutionLeft(iD+1); 			// Kinetic part of the left pressure term
		   ruSquaredRight += qSolutionRight(iD+1)*qSolutionRight(iD+1);			// Kinetic part of the right pressure term
	   }

	   pressureLeft = (gamma_ - 1)*(qSolutionLeft(DIM + 1) - 0.5*ruSquaredLeft*rhoInverseLeft);
	   pressureRight = (gamma_ - 1)*(qSolutionRight(DIM + 1) - 0.5*ruSquaredRight*rhoInverseRight);

	   qAverage(DIM) = (qSolutionLeft(DIM+1) + pressureLeft)*tmp1 + (qSolutionRight(DIM+1) + pressureRight)*tmp2; //H_average

	   //Compute useful variables for constructing the flux function
	   double alphaAvg = 0.0;
	   double unAvg = 0.0;
	   for (std::size_t iD = 0; iD < DIM; iD++)
	   {
		   alphaAvg += (qAverage(iD)*qAverage(iD));
		   unAvg += qAverage(iD)*normal(iD);
	   }
	   alphaAvg *= 0.5;

	   const double a2Avg = std::abs((gamma_ -1)*(qAverage(DIM) - alphaAvg));
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
	   for (std::size_t iD = 0; iD < DIM; iD++)
	   {
		   abv4 += -qAverage(iD)*qDifference(iD+1);
		   abv5 += normal(iD)*qDifference(iD+1);
	   }
	   abv4 += qDifference(DIM+1);
	   abv4 *= (gamma_ -1);

	   const double abv6 = abv3*abv4*ova2Avg + abv2*abv5*ovaAvg;
	   const double abv7 = abv2*abv4*ovaAvg + abv3*abv5;

	   //Compute the Roe Riemann Flux function: h(u_L,u_R) = 0.5*(F(u_L) + F(u_R) - |A|(u_R - u_L))
	   LinearAlgebra::MiddleSizeVector flux(DIM + 2);

	    double pLR = pressureLeft + pressureRight;

	    double runR = 0.0;
	    double runL = 0.0;
	    for (std::size_t iD = 0; iD < DIM; iD++)
	    {
	    	runL += qSolutionLeft(iD+1)*normal(iD);
	    	runR += qSolutionRight(iD+1)*normal(iD);
	    }

	    double unL = runL*rhoInverseLeft;
	    double unR = runR*rhoInverseRight;

	    //continuity equation
	    flux(0) = (runL + runR - (lam3*qDifference(0) + abv6));

	    //momentum equations
	    for (std::size_t iD = 0; iD < DIM; iD++)
	    {
	    	flux(iD+1) = runL*qSolutionLeft(iD+1)*rhoInverseLeft + runR*qSolutionRight(iD+1)*rhoInverseRight + pLR*normal(iD) - (lam3*qDifference(iD+1) + qAverage(iD)*abv6 + normal(iD)*abv7);
	    }

	    //energy equation
	    flux(DIM+1) = (unL*(qSolutionLeft(DIM+1)+pressureLeft) + unR*(qSolutionRight(DIM+1)+pressureRight) - (lam3*qDifference(DIM+1) + qAverage(DIM)*abv6 + unAvg*abv7));

	    //Note: Twice the flux is computed above, hence the factor 0.5 in front of the equation
	    //Note: Correction is made to the flux since F*n was computed above, where n is the normal unit vector
	    //However the face integral function does not use a normalised vector.
        return 0.5*flux*area;
   }

   /// \brief Compute the integrand for the right hand side for the reference face corresponding to an external face.
    template<std::size_t DIM>
   LinearAlgebra::MiddleSizeVector Euler<DIM>::integrandRightHandSideOnRefFace(Base::PhysicalFace<DIM>& face, const double &time, const LinearAlgebra::MiddleSizeVector &solutionCoefficients)
   {
	   //Get the number of basis functions
	   std::size_t numOfBasisFunctionsLeft= face.getFace()->getPtrElementLeft()->getNrOfBasisFunctions(); //Get the number of basis functions on the left

	   LinearAlgebra::MiddleSizeVector& integrand = face.getResultVector();
	   LinearAlgebra::MiddleSizeVector qReconstructionLeft(numOfVariables_);
	   LinearAlgebra::MiddleSizeVector qReconstructionRight(numOfVariables_);

	    // Compute the numerical solution at the given point.
	    std::size_t jVB; // Index for both variable and basis function.
	    for (std::size_t jV = 0; jV < numOfVariables_; jV++)
	    {
	        qReconstructionLeft(jV) = 0;
	        for (std::size_t jB = 0; jB < numOfBasisFunctionsLeft; jB++)
	        {
	            jVB = face.getFace()->getPtrElementLeft()->convertToSingleIndex(jB, jV);
	            qReconstructionLeft(jV) += face.basisFunction(Base::Side::LEFT, jB) * solutionCoefficients(jVB);
	        }
	    }

	   // Boundary face is assumed to be a solid wall: set reflective solution on the other side
	   qReconstructionRight(0) = qReconstructionLeft(0);
	   for (std::size_t iD = 0; iD < DIM; iD++)
	   {
		   qReconstructionRight(iD+1) = -qReconstructionLeft(iD+1);
	   }
	   qReconstructionRight(DIM+1) = qReconstructionLeft(DIM+1);

	   // Compute normal vector, with size of the ref-to-phys face scale, pointing outward of the left element.
	   LinearAlgebra::SmallVector<DIM> normal = face.getNormalVector();


	   //Compute flux
	   //todo:RoeRiemannFluxFunction might benefit from features of PhysicalFace
	   LinearAlgebra::MiddleSizeVector flux = RoeRiemannFluxFunction(qReconstructionRight, qReconstructionLeft, normal);

	   //todo: Other implementation is simply to put the flux to zero. No mass, momentum and energy should get out. This saves computational time

	   // Compute integrand on the reference element.
	   std::size_t iVB; // Index for both variable and basis function.
	   for (std::size_t iB = 0; iB < numOfBasisFunctionsLeft; iB++)
	   {
	       for (std::size_t iV = 0; iV < numOfVariables_; iV++) // Index for direction
	       {
	           iVB = face.getFace()->getPtrElementLeft()->convertToSingleIndex(iB, iV);
	           integrand(iVB) = -flux(iV)*face.basisFunction(Base::Side::LEFT, iB);
	       }
	   }

	   return  integrand;
   }

   /// \brief Compute the integrand for the right hand side for the reference face corresponding to an internal face.
    template<std::size_t DIM>
   LinearAlgebra::MiddleSizeVector Euler<DIM>::integrandRightHandSideOnRefFace(Base::PhysicalFace<DIM> &face, const double &time, const Base::Side &iSide, const LinearAlgebra::MiddleSizeVector &solutionCoefficientsLeft, const LinearAlgebra::MiddleSizeVector &solutionCoefficientsRight)
   {
	   //Get the number of basis functions
	   std::size_t numOfTestBasisFunctions = face.getFace()->getPtrElement(iSide)->getNrOfBasisFunctions(); // Get the number of test basis functions on a given side, iSide
	    std::size_t numOfSolutionBasisFunctionsLeft = face.getFace()->getPtrElementLeft()->getNrOfBasisFunctions(); //Get the number of basis functions on the left
	    std::size_t numOfSolutionBasisFunctionsRight = face.getFace()->getPtrElementRight()->getNrOfBasisFunctions(); //Get the number of basis functions on the right side


	   LinearAlgebra::MiddleSizeVector& integrand = face.getResultVector(iSide);
	   LinearAlgebra::MiddleSizeVector qReconstructionLeft(numOfVariables_);
	   LinearAlgebra::MiddleSizeVector qReconstructionRight(numOfVariables_);

	    // Compute the numerical solution at the given point at the left and right side.
	    std::size_t jVB; // Index for both variable and basis function.
	    for (std::size_t jV = 0; jV < numOfVariables_; jV++)
	    {
	        qReconstructionLeft(jV) = 0;
	        qReconstructionRight(jV) = 0;
	        for (std::size_t jB = 0; jB < numOfSolutionBasisFunctionsLeft; jB++)
	        {
	            jVB = face.getFace()->getPtrElementLeft()->convertToSingleIndex(jB, jV);
	            qReconstructionLeft(jV) += face.basisFunction(Base::Side::LEFT, jB) * solutionCoefficientsLeft(jVB);
	        }
	        for (std::size_t jB = 0; jB < numOfSolutionBasisFunctionsRight; jB++)
	        {
	            jVB = face.getFace()->getPtrElementRight()->convertToSingleIndex(jB, jV);
	            qReconstructionRight(jV) += face.basisFunction(Base::Side::RIGHT, jB) * solutionCoefficientsRight(jVB);
	        }
	    }

	   // Compute normal vector, with size of the ref-to-phys face scale, pointing outward of the left element.
	   LinearAlgebra::SmallVector<DIM> normal = face.getNormalVector();

	   //Compute flux
	   LinearAlgebra::MiddleSizeVector flux;
       //todo:RoeRiemannFluxFunction might benefit from features of PhysicalFace
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
	           iVB = face.getFace()->getPtrElement(iSide)->convertToSingleIndex(iB, iV);
	           integrand(iVB) = -flux(iV)*face.basisFunction(iSide, iB);
	       }
	   }

	    if (integrand[0] != integrand[0])
	    {
	    	logger(ERROR,"integrand on ref face is nan");
	    }
	   return  integrand;
   }

   /// \brief Compute the right-hand side corresponding to a boundary face
    template<std::size_t DIM>
   LinearAlgebra::MiddleSizeVector Euler<DIM>::computeRightHandSideAtFace(Base::Face *ptrFace, LinearAlgebra::MiddleSizeVector &solutionCoefficients, const double time)
   {
	    // Define the integrand function for the right hand side for the reference face.
	    std::function<LinearAlgebra::MiddleSizeVector(Base::PhysicalFace<DIM>& face)> integrandFunction = [&](Base::PhysicalFace<DIM>& face) -> LinearAlgebra::MiddleSizeVector
	    {   return this->integrandRightHandSideOnRefFace(face, time, solutionCoefficients);};

	    return this->faceIntegrator_.integrate(ptrFace, integrandFunction);
   }

   /// \brief Compute the right-hand side corresponding to an internal face
    template<std::size_t DIM>
   LinearAlgebra::MiddleSizeVector Euler<DIM>::computeRightHandSideAtFace(Base::Face *ptrFace, const Base::Side side, LinearAlgebra::MiddleSizeVector &solutionCoefficientsLeft, LinearAlgebra::MiddleSizeVector &solutionCoefficientsRight, const double time)
   {
	    // Define the integrand function for the right hand side for the reference face.
	    std::function<LinearAlgebra::MiddleSizeVector(Base::PhysicalFace<DIM> &)> integrandFunction = [&](Base::PhysicalFace<DIM> &face) -> LinearAlgebra::MiddleSizeVector
	    {   return this->integrandRightHandSideOnRefFace(face, time, side, solutionCoefficientsLeft, solutionCoefficientsRight);};
	    return this->faceIntegrator_.integrate(ptrFace, integrandFunction);
   }


   /// *****************************************
   /// ***    		Various Functions        ***
   /// *****************************************

    template<std::size_t DIM>
    void Euler<DIM>::tasksBeforeSolving()
    {
        this->faceIntegrator_.setTransformation(std::shared_ptr<Base::CoordinateTransformation<DIM> >(new Base::DoNotScaleIntegrands<DIM>(new Base::H1ConformingTransformation<DIM>())));
        Base::HpgemAPISimplified<DIM>::tasksBeforeSolving();
    }


    template<std::size_t DIM>
   LinearAlgebra::MiddleSizeVector Euler<DIM>::getExactSolution(const PointPhysicalT &pPhys, const double &time, const std::size_t orderTimeDerivative)
   {
		LinearAlgebra::MiddleSizeVector exactSolution(numOfVariables_);

		double amplitude = 0.2;
		double frequency = 2.0*M_PI;
		double function = amplitude*std::cos(frequency*time);

		for (std::size_t iD = 0; iD < DIM; iD++)
		{
			function *= std::cos(frequency*pPhys[iD]);
		}

		exactSolution(0) = 1.5 + function;

		for (std::size_t iD = 0; iD < DIM; iD++)
		{
			exactSolution(iD+1) = function;
		}

		exactSolution(DIM + 1) = 30.0 + function;

	    return exactSolution;
   }

   /// \brief Compute the initial solution at a given point in space and time.
    template<std::size_t DIM>
   LinearAlgebra::MiddleSizeVector Euler<DIM>::getInitialSolution(const PointPhysicalT &pPhys, const double &startTime, const std::size_t orderTimeDerivative)
   {
       return getExactSolution(pPhys, startTime, orderTimeDerivative);
   }

   /// \brief Computes the error for output purposes
    template<std::size_t DIM>
   LinearAlgebra::MiddleSizeVector Euler<DIM>::Error(const double time)
   {
	   return this->computeMaxError(this->solutionTimeLevel_, time);
   }
















