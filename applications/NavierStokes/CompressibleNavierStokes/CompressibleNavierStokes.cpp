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

#include <cmath>

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

// todo: rewrite all the normal vector queries etc.

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
    }
    description.boundaryConditions_[0] = Base::BoundaryType::PERIODIC;
    description.boundaryConditions_[1] = Base::BoundaryType::PERIODIC;

    return description;
}

void CompressibleNavierStokes::setStabilityMassMatrix()
{
	//For a single element create the mass matrix: note this breaks down with p-refinement or limiters
	LinearAlgebra::MiddleSizeMatrix stabilityMassMatrix  = computeMassMatrixAtElement(meshes_[0]->getElementsList()[0]);
	viscousTerms_.setStabilityMassMatrix(stabilityMassMatrix);
}

double CompressibleNavierStokes::computePressure(const LinearAlgebra::MiddleSizeVector &state)
{
	//Compute pressure term
	double pressure;
	double inverseDensity = 1.0/state(0);
	double momentumSquared = 0.0;
	for(std::size_t iD = 0; iD < DIM_; iD++)
	{
		momentumSquared += state(iD+1)*state(iD+1); // (u^2 + v^2 + w^2)*rho^2
	}


	pressure = (gamma_ -1)*(state(DIM_+1)*nonDIM1_ - 0.5*inverseDensity*(momentumSquared)); // (gamma-1)*rho*(e*nondim1- (u^2 + v^2 + w^2)/2), where nondim1 is scaling such that it is dimensionless

	logger.assert(pressure > 0, "Negative pressure.");

	return pressure;
}


/// *************************************************
/// ***   Element integration support functions   ***
/// *************************************************

///  \brief Constructs the solution based on the solutionCoefficients.
LinearAlgebra::MiddleSizeVector CompressibleNavierStokes::computeStateOnElement(Base::PhysicalElement<DIM> &element, const LinearAlgebra::MiddleSizeVector &solutionCoefficients)
{
		std::size_t numberOfBasisFunctions =  element.getNumOfBasisFunctions();
		LinearAlgebra::MiddleSizeVector elementState(numOfVariables_);
		std::size_t iVB; // Index in solution coefficients for variable i and basisfunction j

		for(std::size_t iV = 0; iV < numOfVariables_; iV++)
		{
			elementState(iV) = 0.0;
			for (std::size_t iB = 0; iB < numberOfBasisFunctions; iB++)
			{
				iVB = element.convertToSingleIndex(iB,iV);
				elementState(iV) += solutionCoefficients(iVB)*element.basisFunction(iB);
			}
		}

		return elementState;
}

LinearAlgebra::MiddleSizeMatrix CompressibleNavierStokes::computeStateJacobianAtElement(Base::PhysicalElement<DIM> &element, const LinearAlgebra::MiddleSizeVector &solutionCoefficients)
{
		std::size_t numberOfBasisFunctions =  element.getNumOfBasisFunctions();
		std::size_t iVB; // Index in solution coefficients for variable i and basisfunction j

		LinearAlgebra::MiddleSizeMatrix solutionGradient(numOfVariables_,DIM_);
		LinearAlgebra::SmallVector<DIM> gradientBasisFunction;

		for (std::size_t iB = 0; iB < numberOfBasisFunctions; iB++) //Note, For the current purpose the derivative of rhoE is not required.
		{
			gradientBasisFunction = element.basisFunctionDeriv(iB);
			for (std::size_t iV = 0; iV < numOfVariables_; iV++)
			{
				iVB = element.convertToSingleIndex(iB,iV);
				for (std::size_t iD = 0; iD < DIM_; iD++)
				{
					solutionGradient(iV,iD) += solutionCoefficients(iVB)*gradientBasisFunction(iD);
				}
			}
		}

		return solutionGradient;
}

LinearAlgebra::MiddleSizeVector CompressibleNavierStokes::computePartialState(const LinearAlgebra::MiddleSizeVector &state)
{
	LinearAlgebra::MiddleSizeVector partialState(DIM_ + 2);
	double q1Inverse = 1.0/state(0);

	partialState(0) = state(0);

	for (std::size_t iD = 0; iD < DIM_ + 1; iD++)
	{
		partialState(iD+1) = state(iD+1)*q1Inverse;
	}
	return partialState;
}


/// *****************************************
/// ***   Element integration functions   ***
/// *****************************************

/// \brief computes the source at an element
LinearAlgebra::MiddleSizeVector CompressibleNavierStokes::integrandSourceAtElement(Base::PhysicalElement<DIM>& element, const LinearAlgebra::MiddleSizeVector &state, const double &pressureTerm, const double &time)
{
	std::size_t numOfBasisFunctions = element.getNumOfBasisFunctions();

	LinearAlgebra::MiddleSizeVector integrandSource(numOfVariables_ * numOfBasisFunctions);

	//Convert pRef to pPhys
	Geometry::PointPhysical<DIM> pPhys = element.getPointPhysical();

	//create datastructures
	double A = 1.225;
	double B = 239750.0;;
	double C = 1.0; //std::exp(-1000*time);
	double C_t = 0.0; //-1000*std::exp(-1000*time);
	double D = 50.0;
	double pi = M_PI;
	LinearAlgebra::MiddleSizeVector zVal(4);
	LinearAlgebra::MiddleSizeVector exactState(4);
	LinearAlgebra::MiddleSizeVector state_t(4);
	LinearAlgebra::MiddleSizeVector state_x(4);
	LinearAlgebra::MiddleSizeVector state_y(4);
	LinearAlgebra::MiddleSizeVector state_xx(4);
	LinearAlgebra::MiddleSizeVector state_yy(4);
	LinearAlgebra::MiddleSizeVector state_xy(4);

	LinearAlgebra::MiddleSizeVector partialState(4);
	LinearAlgebra::MiddleSizeVector partialState_x(4);
	LinearAlgebra::MiddleSizeVector partialState_y(4);
	LinearAlgebra::MiddleSizeVector partialState_xx(4);
	LinearAlgebra::MiddleSizeVector partialState_yy(4);
	LinearAlgebra::MiddleSizeVector partialState_xy(4);

	//todo: Add some amplitude scaling; how does it propogate through the derivation?
	//Compute Zval
	zVal(0) = std::cos(2*pi*pPhys[0])*std::cos(2*pi*pPhys[1]);
	zVal(1) = std::sin(2*pi*pPhys[0])*std::sin(2*pi*pPhys[1]);
	zVal(2) = std::cos(2*pi*pPhys[0])*std::sin(2*pi*pPhys[1]);
	zVal(3) = std::sin(2*pi*pPhys[0])*std::cos(2*pi*pPhys[1]);

	//compute state
	exactState(0) = A + C*zVal(0);
	exactState(1) = D + C*zVal(3);
	exactState(2) = D + C*zVal(2);
	exactState(3) = B + C*zVal(0);
	double inverseDensity = 1.0/exactState(0);

	//Compute state_t
	state_t(0) 	= C_t*zVal(0);
	state_t(1) 	= C_t*zVal(3);
	state_t(2) 	= C_t*zVal(2);
	state_t(3)  = C_t*zVal(0);

	//compute state_x
	state_x(0) = -2*pi*C*zVal(3);
	state_x(1) = 2*pi*C*zVal(0);
	state_x(2) = -2*pi*C*zVal(1);
	state_x(3) = -2*pi*C*zVal(3);

	//compute state_y
	state_y(0) = -2*pi*C*zVal(2);
	state_y(1) = -2*pi*C*zVal(1);
	state_y(2) = 2*pi*C*zVal(0);
	state_y(3) = -2*pi*C*zVal(2);

	//compute state_xx
	state_xx(0) = -4*pi*pi*C*zVal(0);
	state_xx(1) = -4*pi*pi*C*zVal(3);
	state_xx(2) = -4*pi*pi*C*zVal(2);
	state_xx(3) = -4*pi*pi*C*zVal(0);

	//compute state_yy
	state_yy(0) = -4*pi*pi*C*zVal(0);
	state_yy(1) = -4*pi*pi*C*zVal(3);
	state_yy(2) = -4*pi*pi*C*zVal(2);
	state_yy(3) = -4*pi*pi*C*zVal(0);

	//compute state_xy
	state_xy(0) = -4*pi*pi*C*zVal(1);
	state_xy(1) = -4*pi*pi*C*zVal(2);
	state_xy(2) = -4*pi*pi*C*zVal(3);
	state_xy(3) = -4*pi*pi*C*zVal(1);

	//compute partialState
	partialState(0) = exactState(0);
	partialState(1) = exactState(1)*inverseDensity;
	partialState(2) = exactState(2)*inverseDensity;
	partialState(3) = exactState(3)*inverseDensity;

	//compute partialState_x
	partialState_x(0) = state_x(0);
	partialState_x(1) = inverseDensity*(state_x(1)-partialState(1)*state_x(0));
	partialState_x(2) = inverseDensity*(state_x(2)-partialState(2)*state_x(0));
	partialState_x(3) = inverseDensity*(state_x(3)-partialState(3)*state_x(0));

	//compute partialState_y
	partialState_y(0) = state_y(0);
	partialState_y(1) = inverseDensity*(state_y(1)-partialState(1)*state_y(0));
	partialState_y(2) = inverseDensity*(state_y(2)-partialState(2)*state_y(0));
	partialState_y(3) = inverseDensity*(state_y(3)-partialState(3)*state_y(0));

	//compute partialState_xx
	partialState_xx(0) = state_xx(0);
	partialState_xx(1) = (-inverseDensity*partialState_x(1)
								+ inverseDensity*inverseDensity*partialState(1)*state_x(0)
								- inverseDensity*inverseDensity*state_x(1))*state_x(0)
								+ inverseDensity*state_xx(1)
								- inverseDensity*partialState(1)*state_xx(0);
	partialState_xx(2) = (-inverseDensity*partialState_x(2)
								+ inverseDensity*inverseDensity*partialState(2)*state_x(0)
								- inverseDensity*inverseDensity*state_x(2))*state_x(0)
								+ inverseDensity*state_xx(2)
								- inverseDensity*partialState(2)*state_xx(0);
	partialState_xx(3) = (-inverseDensity*partialState_x(3)
								+ inverseDensity*inverseDensity*partialState(3)*state_x(0)
								- inverseDensity*inverseDensity*state_x(3))*state_x(0)
								+ inverseDensity*state_xx(3)
								- inverseDensity*partialState(3)*state_xx(0);

	//compute partialState_yy
	partialState_yy(0) = state_yy(0);
	partialState_xx(1) = (-inverseDensity*partialState_y(1)
								+ inverseDensity*inverseDensity*partialState(1)*state_y(0)
								- inverseDensity*inverseDensity*state_y(1))*state_y(0)
								+ inverseDensity*state_yy(1)
								- inverseDensity*partialState(1)*state_yy(0);
	partialState_xx(2) = (-inverseDensity*partialState_y(2)
								+ inverseDensity*inverseDensity*partialState(2)*state_y(0)
								- inverseDensity*inverseDensity*state_y(2))*state_y(0)
								+ inverseDensity*state_yy(2)
								- inverseDensity*partialState(2)*state_yy(0);
	partialState_xx(3) = (-inverseDensity*partialState_y(3)
								+ inverseDensity*inverseDensity*partialState(3)*state_y(0)
								- inverseDensity*inverseDensity*state_y(3))*state_y(0)
								+ inverseDensity*state_yy(3)
								- inverseDensity*partialState(3)*state_yy(0);

	//compute partialState_xy
	partialState_xy(0) = state_xy(0);
	partialState_xy(1) = -inverseDensity*inverseDensity*state_x(0)*state_y(1)
								+inverseDensity*state_xy(1)
								-inverseDensity*partialState_x(1)*state_y(0)
								+inverseDensity*inverseDensity*partialState(1)*state_x(0)*state_y(0)
								-inverseDensity*partialState(1)*state_xy(0);
	partialState_xy(2) = -inverseDensity*inverseDensity*state_x(0)*state_y(2)
								+inverseDensity*state_xy(2)
								-inverseDensity*partialState_x(2)*state_y(0)
								+inverseDensity*inverseDensity*partialState(2)*state_x(0)*state_y(0)
								-inverseDensity*partialState(2)*state_xy(0);
	partialState_xy(1) = -inverseDensity*inverseDensity*state_x(0)*state_y(3)
								+inverseDensity*state_xy(3)
								-inverseDensity*partialState_x(3)*state_y(0)
								+inverseDensity*inverseDensity*partialState(3)*state_x(0)*state_y(0)
								-inverseDensity*partialState(3)*state_xy(0);
	//*****************************
	//*** 	Additional values	***
	//*****************************

	double conv_vv_x = 2*partialState(2)*state_x(2)
   		   	   -partialState(2)*partialState(2)*state_x(0);

	double conv_uu_y = 2*partialState(1)*state_y(1)
				-partialState(1)*partialState(1)*state_y(0);

	//*************************************
	//***	Compute closure variables	***
	//*************************************

	double p =  computePressure(exactState);

	//double e = partialState(3) - 0.5*(partialState(1)*partialState(1) + partialState(2)*partialState(2));
	double e_x = partialState_x(3) - partialState(1)*partialState_x(1) - partialState(2)*partialState_x(2);
	double e_y = partialState_y(3) - partialState(1)*partialState_y(1) - partialState(2)*partialState_y(2);
	double e_xx = partialState_xx(3)
						-partialState_x(1)*partialState_x(1) - partialState(1)*partialState_xx(1)
						-partialState_x(2)*partialState_x(2) - partialState(2)*partialState_xx(2);
	double e_yy = partialState_yy(3)
						-partialState_y(1)*partialState_y(1) - partialState(1)*partialState_yy(1)
						-partialState_y(2)*partialState_y(2) - partialState(2)*partialState_yy(2);

	//WARNING cv is hardcoded
	double cv = 1000/gamma_;
	double T = viscousTerms_.computeTemperature(exactState,p);
	double T_x = e_x/cv;
	double T_y = e_y/cv;
	double T_xx = e_xx/cv;
	double T_yy = e_yy/cv;

	//WARNING: Tref Lambda and Ts are hardcoded
	double Tref = 288.16;
	double Ts = 110;
	double lambda = -2/3;
	double mu = viscousTerms_.computeViscosity(T);
	double mu_x = mu*(3*Tref/(2*T) - 1/(T + Ts))*T_x;
	double mu_y = mu*(3*Tref/(2*T) - 1/(T + Ts))*T_y;

	double Txx = (2 + lambda)*mu*partialState_x(1) + lambda*mu*partialState_y(2);
	double Txy = mu*(partialState_y(1) + partialState_x(2));
	double Tyy = (2 + lambda)*mu*partialState_y(2) + lambda*mu*partialState_x(1);

	//WARNING: Kappa is hardcoded
	double kappa = gamma_*mu*0.71*cv;

	//**********************************
	//***	Compute equation terms	 ***
	//**********************************

	double conv_uu_x = 2*partialState(1)*state_x(1)
								-partialState(1)*partialState(1)*state_x(0); // conv_uu_x = d (rho * u * u ) / dx

	double conv_vu_y = partialState(1)*state_y(2)
								+partialState(2)*state_y(1)
								-partialState(1)*partialState(2)*state_y(0);

	double conv_uv_x = partialState(1)*state_x(2)
								+partialState(2)*state_x(1)
								-partialState(1)*partialState(2)*state_x(0);

	double conv_vv_y = 2*partialState(2)*state_y(2)
	   							-partialState(2)*partialState(2)*state_y(0);

	double p_x = (gamma_ - 1)*(state_x(3) - 0.5*conv_uu_x - 0.5*conv_vv_x);

	double p_y = (gamma_ - 1)*(state_y(3) - 0.5*conv_uu_y - 0.5*conv_vv_y);

	double Txx_x = (2 + lambda)*mu_x*partialState_x(1)
						+ (2 + lambda)*mu*partialState_xx(1)
						+ lambda*mu_x*partialState_y(2)
						+ lambda*mu*partialState_xy(2);

	double Tyx_y = mu_x*(partialState_y(1) + partialState_x(2))
						+ mu*(partialState_xy(1) + partialState_xx(2));

	double Txy_x = mu_y*(partialState_y(1) + partialState_x(2))
						+ mu*(partialState_yy(1) + partialState_xy(2));

	double Tyy_y = (2 + lambda)*mu_y*partialState_y(2)
						+ (2 + lambda)*mu*partialState_yy(2)
						+ lambda*mu_y*partialState_x(1)
						+ lambda*mu*partialState_xy(1);

	double conv_Hu_x = exactState(3)*partialState_x(1)
								+ partialState(1)*state_x(3)
								+ p*partialState_x(1)
								+ partialState(1)*p_x;

	double conv_Hv_y = exactState(3)*partialState_y(2)
							+ partialState(2)*state_y(3)
							+ p*partialState_y(2)
							+ partialState(2)*p_y;

	double Ax_x = partialState(1)*Txx_x
						+ Txx*partialState_x(1)
						+ partialState(2)*Txy_x
						+ Txy*partialState_x(2)
						+ kappa*T_xx;

	double Ay_y = partialState(1)*Tyx_y
						+ Txy*partialState_y(1)
						+ partialState(2)*Tyy_y
						+ Tyy*partialState_y(2)
						+ kappa*T_yy;



	//*********************************
	//***	Compute Source term		***
	//*********************************

	double sDensity = state_t(0) + state_x(1) + state_y(2);
	double sMomentumX = state_t(1) + conv_uu_x + conv_vu_y + p_x - Txx_x - Tyx_y;
	double sMomentumY = state_t(2) + conv_uv_x + conv_vv_y + p_y - Txy_x - Tyy_y;
	double sEnergy = state_t(3) + conv_Hu_x + conv_Hv_y - Ax_x -Ay_y;

	//*************************************
	//***	Compute source integrand	***
	//*************************************

	std::size_t iVB;
	for(std::size_t iB = 0; iB < numOfBasisFunctions; iB++)
	{
		// Density
		iVB = element.convertToSingleIndex(iB,0);
		integrandSource(iVB) = sDensity*element.basisFunction(iB);

		// Momentum x
		iVB = element.convertToSingleIndex(iB,1);
		integrandSource(iVB) = sMomentumX*element.basisFunction(iB);

		// Momentum Y
		iVB = element.convertToSingleIndex(iB,2);
		integrandSource(iVB) = sMomentumY*element.basisFunction(iB);

		// Energy
		iVB = element.convertToSingleIndex(iB,3);
		integrandSource(iVB) = sEnergy*element.basisFunction(iB);
	}



	return integrandSource;
}

LinearAlgebra::MiddleSizeVector CompressibleNavierStokes::integrandRightHandSideOnElement(Base::PhysicalElement<DIM>& element, const double &time, const LinearAlgebra::MiddleSizeVector &stateCoefficients)
{

	//reconstruct the solution, partial state jacobian and pressure at pRef
	//todo: check if partialState and its Jacobian can somehow be computed efficiently together
	const LinearAlgebra::MiddleSizeVector state = computeStateOnElement(element, stateCoefficients);
	const LinearAlgebra::MiddleSizeMatrix stateJacobian = computeStateJacobianAtElement(element, stateCoefficients);
	const LinearAlgebra::MiddleSizeVector partialState = computePartialState(state);
	const double pressure = computePressure(state);

	//Compute inviscid terms
	LinearAlgebra::MiddleSizeVector integrandInviscid = inviscidTerms_.integrandAtElement(element, time, pressure, state);

	//Compute viscous terms
	LinearAlgebra::MiddleSizeVector integrandViscous = viscousTerms_.integrandAtElement(element,  state, stateJacobian, pressure, partialState);

	//Compute source terms
	LinearAlgebra::MiddleSizeVector integrandSource = integrandSourceAtElement(element, state, pressure, time);

	//DEBUG
	//std::cout << "Element_Inviscid: " << integrandInviscid << std::endl;
/*	std::size_t iVBx;
	std::size_t iVBy;
	LinearAlgebra::MiddleSizeVector integrandDif(element.getNumOfBasisFunctions());
	for(std::size_t iB = 0; iB < element.getNumOfBasisFunctions(); iB++)
		{
			// Momentum x
			iVBx = element.convertToSingleIndex(iB,1);
			iVBy = element.convertToSingleIndex(iB,2);
			integrandDif(iB) = integrandSource(iVBx) - integrandSource(iVBy);
		}
	std::cout << integrandDif << std::endl;*/
	//std::cout << "Element_Viscous: " << integrandViscous << std::endl;
	//std::cout << "Element_Source: " << integrandSource << std::endl;
    return (integrandInviscid + integrandViscous + integrandSource);

}

LinearAlgebra::MiddleSizeVector CompressibleNavierStokes::computeRightHandSideAtElement(Base::Element *ptrElement, LinearAlgebra::MiddleSizeVector &stateCoefficients, const double time)
{
	std::function<LinearAlgebra::MiddleSizeVector(Base::PhysicalElement<DIM> &)> integrandFunction = [=](Base::PhysicalElement<DIM>& element) -> LinearAlgebra::MiddleSizeVector
	    {   return this->integrandRightHandSideOnElement(element, time, stateCoefficients);};

    return elementIntegrator_.integrate(ptrElement, integrandFunction, ptrElement->getGaussQuadratureRule());

}


/// *************************************************
/// ***    face integration support functions     ***
/// *************************************************

LinearAlgebra::MiddleSizeVector CompressibleNavierStokes::computeStateOnFace(Base::PhysicalFace<DIM> &face, const Base::Side &iSide, const LinearAlgebra::MiddleSizeVector &stateCoefficients) const
{
	std::size_t numOfBasisFunctions =  face.getPhysicalElement(iSide).getNumOfBasisFunctions();
	LinearAlgebra::MiddleSizeVector elementState(numOfVariables_);
	std::size_t iVB; // Index in solution coefficients for variable i and basisfunction j

	for(std::size_t iV = 0; iV < numOfVariables_; iV++)
	{
		elementState(iV) = 0.0;
		for (std::size_t iB = 0; iB < numOfBasisFunctions; iB++)
		{
			iVB = face.getPhysicalElement(iSide).convertToSingleIndex(iB,iV);
			elementState(iV) += stateCoefficients(iVB)*face.basisFunction(iSide, iB); //basisFunction returns physical value
		}
	}

	return elementState;
}


LinearAlgebra::MiddleSizeMatrix CompressibleNavierStokes::computeStateJacobianAtFace(Base::PhysicalFace<DIM> &face, const Base::Side &iSide, const LinearAlgebra::MiddleSizeVector &stateCoefficients)
{
	std::size_t numberOfBasisFunctions =  face.getPhysicalElement(iSide).getNumOfBasisFunctions();
	std::size_t iVB; // Index in solution coefficients for variable i and basisfunction j

	LinearAlgebra::MiddleSizeMatrix stateJacobian(numOfVariables_,DIM_);
	LinearAlgebra::SmallVector<DIM> gradientBasisFunction;

	for (std::size_t iB = 0; iB < numberOfBasisFunctions; iB++)
	{
		gradientBasisFunction = face.basisFunctionDeriv(iSide,iB);
		for (std::size_t iV = 0; iV < numOfVariables_; iV++)
		{
			iVB = face.getPhysicalElement(iSide).convertToSingleIndex(iB,iV);
			for (std::size_t iD = 0; iD < DIM_; iD++)
			{
				stateJacobian(iV,iD) += stateCoefficients(iVB)*gradientBasisFunction(iD);
			}
		}
	}

	return stateJacobian;
}


/// **************************************************
/// ***    external face integration functions     ***
/// **************************************************

/// \brief Compute the integrand for the right hand side for the reference face corresponding to an external face.
//This function is written for the Couette type flow, a plate at top and a plate at bottom
//The state is reconstructed based on the plate's temperature and movement speed, and therefore this is a dirichlet type of BC
//It is handled like an internal element.
LinearAlgebra::MiddleSizeVector CompressibleNavierStokes::integrandRightHandSideOnFace(Base::PhysicalFace<DIM> &face, const double &time, const LinearAlgebra::MiddleSizeVector &stateCoefficients)
{
	//Compute the internal state
	const LinearAlgebra::MiddleSizeVector stateInternal = computeStateOnFace(face, Base::Side::LEFT, stateCoefficients);
	const LinearAlgebra::SmallVector<DIM> normalInternal = face.getNormalVector();
	//const double pressureInternal = computePressure(stateInternal);
	const LinearAlgebra::MiddleSizeVector partialStateInternal = computePartialState(stateInternal);
	const LinearAlgebra::MiddleSizeMatrix stateJacobianInternal = computeStateJacobianAtFace(face, Base::Side::LEFT, stateCoefficients);

	//Determine if this is the top or bottom boundary and set boundary state
	LinearAlgebra::MiddleSizeVector stateBoundary(DIM_+2);
	double temperatureBoundary;
	double velocityBoundary;
	double exponent;
	const Geometry::PointPhysical<DIM> pPhys = face.getPointPhysical();
	if(pPhys[1] > 0.8)
	{
		exponent = time*time/(Tc_*Tc_);
		velocityBoundary = uPlateTop_ - std::exp(-exponent);							// Spatial blending factor for the transient phase
		temperatureBoundary = tPlateTop_;
		stateBoundary(0) = (gamma_ - 1.0)*stateInternal(DIM_ + 1)/tPlateTop_/Rs_;	 	// density
		stateBoundary(1) = stateBoundary(0)*velocityBoundary; 								// velocity u on plate is uPlateTop_
		stateBoundary(2) = 0; 															//velocity v on plate is zero
		stateBoundary(DIM_+1) = stateInternal(DIM_ + 1);								// energy is copied from internal
	}
	else
	{
		exponent = time*time/(Tc_*Tc_);
		velocityBoundary = uPlateBottom_ - std::exp(-exponent);
		temperatureBoundary = tPlateBottom_;
		stateBoundary(0) = (gamma_ - 1.0)*stateInternal(DIM_ + 1)/tPlateBottom_/Rs_;	// density
		stateBoundary(1) = stateBoundary(0)*uPlateBottom_; 								// velocity u on plate is uPlateBottom_
		stateBoundary(2) = 0; 															// velocity v on plate is zero
		stateBoundary(DIM_+1) = stateInternal(DIM_ + 1);								// energy is copied from internal
	}
	const LinearAlgebra::MiddleSizeVector partialStateBoundary = computePartialState(stateBoundary);

	//Compute inviscid terms
	//todo:FIX ERROR: unitNormal goes in here
	LinearAlgebra::MiddleSizeVector integrandInviscid = inviscidTerms_.integrandAtFace(face, time, Base::Side::LEFT, stateInternal, stateBoundary, normalInternal);

	//Compute viscous terms
	LinearAlgebra::MiddleSizeVector integrandViscous = integrandInviscid; //viscousTerms_.integrandViscousAtFace(face, Base::Side::LEFT, stateInternal, stateBoundary, pressureInternal, partialStateInternal, normalInternal);


	//Compute support variable terms
	LinearAlgebra::MiddleSizeVector integrandAuxilliary = integrandInviscid; //viscousTerms_.integrandAuxilliaryAtFace(face, temperatureBoundary, stateInternal, stateBoundary, partialStateBoundary, normalInternal, stateJacobianInternal);

	return  integrandInviscid + integrandViscous + integrandAuxilliary;
}

/// \brief Compute the right-hand side corresponding to a boundary face
LinearAlgebra::MiddleSizeVector CompressibleNavierStokes::computeRightHandSideAtFace(Base::Face *ptrFace, LinearAlgebra::MiddleSizeVector &solutionCoefficients, const double time)
{
	    // Define the integrand function for the right hand side for the reference face.
	    std::function<LinearAlgebra::MiddleSizeVector(Base::PhysicalFace<DIM>&)> integrandFunction = [=](Base::PhysicalFace<DIM>& face) -> LinearAlgebra::MiddleSizeVector
	    {   return this->integrandRightHandSideOnFace(face, time, solutionCoefficients);};

	    return faceIntegrator_.integrate(ptrFace, integrandFunction);
}


/// **************************************************
/// ***    internal face integration functions     ***
/// **************************************************

/// \brief Compute the integrand for the right hand side for the face corresponding to an internal face.
   LinearAlgebra::MiddleSizeVector CompressibleNavierStokes::integrandRightHandSideOnFace(Base::PhysicalFace<DIM>& face, const double &time, const Base::Side &iSide, const LinearAlgebra::MiddleSizeVector &stateCoefficientsLeft, const LinearAlgebra::MiddleSizeVector &stateCoefficientsRight)
   {
       //reconstruct the solution, partial state Jacobian and pressure at pRef, left and right of the interface
       const LinearAlgebra::MiddleSizeVector stateLeft = computeStateOnFace(face, Base::Side::LEFT, stateCoefficientsLeft);
       const LinearAlgebra::MiddleSizeVector stateRight = computeStateOnFace(face, Base::Side::RIGHT, stateCoefficientsRight);
       const LinearAlgebra::MiddleSizeMatrix stateJacobianLeft = computeStateJacobianAtFace(face, Base::Side::LEFT, stateCoefficientsLeft);
       const LinearAlgebra::MiddleSizeMatrix stateJacobianRight = computeStateJacobianAtFace(face, Base::Side::RIGHT, stateCoefficientsRight);

		//Determine internal and external solutions
		LinearAlgebra::MiddleSizeVector stateInternal;
		LinearAlgebra::MiddleSizeVector stateExternal;
		LinearAlgebra::MiddleSizeMatrix stateJacobianInternal;
		LinearAlgebra::MiddleSizeMatrix stateJacobianExternal;
 		LinearAlgebra::SmallVector<DIM> unitNormalInternal;
		LinearAlgebra::SmallVector<DIM> unitNormalExternal;

		if (iSide == Base::Side::RIGHT)
		{
			stateInternal = stateRight;
			stateExternal = stateLeft;
			stateJacobianInternal = stateJacobianRight;
			stateJacobianExternal = stateJacobianLeft;
			unitNormalInternal = -face.getUnitNormalVector();
		}
		else
		{
			stateInternal = stateLeft;
			stateExternal = stateRight;
			stateJacobianInternal = stateJacobianLeft;
			stateJacobianExternal = stateJacobianRight;
			unitNormalInternal = face.getUnitNormalVector();
		}

		const double pressureInternal = computePressure(stateInternal);
		const double pressureExternal = computePressure(stateExternal);
		const LinearAlgebra::MiddleSizeVector partialStateInternal = computePartialState(stateInternal);
		const LinearAlgebra::MiddleSizeVector partialStateExternal = computePartialState(stateExternal);

		//Compute inviscid terms
		LinearAlgebra::MiddleSizeVector integrandInviscid = inviscidTerms_.integrandAtFace(face, time, iSide, stateInternal, stateExternal, unitNormalInternal);

		//Compute viscous terms
		LinearAlgebra::MiddleSizeVector integrandViscous = viscousTerms_.integrandViscousAtFace(face, iSide, stateInternal, stateExternal, pressureInternal, partialStateInternal, unitNormalInternal);

		//Compute support variable terms
		//todo: write out the integral as summ, see if things cancel
		//todo: Integrate over whole face in one go? More efficient?
		LinearAlgebra::MiddleSizeVector integrandAuxilliary = viscousTerms_.integrandAuxilliaryAtFace(face, iSide, stateInternal, stateExternal, pressureInternal, pressureExternal, partialStateInternal, partialStateExternal, stateJacobianInternal, stateJacobianExternal, unitNormalInternal);


		//DEBUG
		/*	std::size_t iVBx;
			std::size_t iVBy;
			LinearAlgebra::MiddleSizeVector integrandDif(element.getNumOfBasisFunctions());
			for(std::size_t iB = 0; iB < element.getNumOfBasisFunctions(); iB++)
				{
					// Momentum x
					iVBx = element.convertToSingleIndex(iB,1);
					iVBy = element.convertToSingleIndex(iB,2);
					integrandDif(iB) = integrandSource(iVBx) - integrandSource(iVBy);
				}
			std::cout << integrandDif << std::endl;*/

		//std::cout << integrandAuxilliary << std::endl;
		//std::cout << integrandViscous << std::endl;
		return  integrandInviscid + integrandViscous + integrandAuxilliary;

   }

   /// \brief Compute the right-hand side corresponding to an internal face
   LinearAlgebra::MiddleSizeVector CompressibleNavierStokes::computeRightHandSideAtFace(Base::Face *ptrFace, const Base::Side side, LinearAlgebra::MiddleSizeVector &solutionCoefficientsLeft, LinearAlgebra::MiddleSizeVector &solutionCoefficientsRight, const double time)
   {
	    // Define the integrand function for the right hand side for the face.
	    std::function<LinearAlgebra::MiddleSizeVector(Base::PhysicalFace<DIM>&)> integrandFunction = [=](Base::PhysicalFace<DIM>& face) -> LinearAlgebra::MiddleSizeVector
	    {   return this->integrandRightHandSideOnFace(face, time, side, solutionCoefficientsLeft, solutionCoefficientsRight);};
	    return faceIntegrator_.integrate(ptrFace, integrandFunction);
   }



   /// *****************************************
   /// ***    		Various Functions        ***
   /// *****************************************

   LinearAlgebra::MiddleSizeVector CompressibleNavierStokes::getExactSolution(const PointPhysicalT &pPhys, const double &time, const std::size_t orderTimeDerivative)
   {

		LinearAlgebra::MiddleSizeVector exactSolution(numOfVariables_);
/*

		exactSolution(0) = 1.225; // + function;

		exactSolution(1) = 2.0;//*function;

		exactSolution(2) = 2.0;//*function;

		exactSolution(DIM_ + 1) = exactSolution(0)*cp_/gamma_*288.0 + 0.5/exactSolution(0)*(exactSolution(1)*exactSolution(1) + exactSolution(2)*exactSolution(2)); // + 0.5*exactSolution(1)*exactSolution(1)/exactSolution(0); // + function;
*/



/*		double Val = std::sqrt((pPhys[0]-0.5)*(pPhys[0]-0.5) + (pPhys[1]-0.5)*(pPhys[1]-0.5));

		if (Val < 0.15)
		{
			exactSolution(0) = 1.5;
		}
		else
		{
			exactSolution(0) = 1;
		}


		exactSolution(1) = 0.0;


		exactSolution(DIM + 1) = exactSolution(0)*30.0;*/

		double A = 1.225;
		double B = 239750.0;
		double C = 1.0; //std::exp(-1000*time);
		double D = 50.0;
		double pi = M_PI;
		LinearAlgebra::MiddleSizeVector zVal(4);

		//todo: Add some amplitude scaling; how does it propagate through the derivation?
		//Compute Zval
		zVal(0) = std::cos(2*pi*pPhys[0])*std::cos(2*pi*pPhys[1]);
		zVal(1) = std::sin(2*pi*pPhys[0])*std::sin(2*pi*pPhys[1]);
		zVal(2) = std::cos(2*pi*pPhys[0])*std::sin(2*pi*pPhys[1]);
		zVal(3) = std::sin(2*pi*pPhys[0])*std::cos(2*pi*pPhys[1]);

		//compute state
		exactSolution(0) = A + C*zVal(0);
		exactSolution(1) = D + C*zVal(3);
		exactSolution(2) = D + C*zVal(2);
		exactSolution(3) = B + C*zVal(0);

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

   /// \brief Show the progress of the time integration.
   void CompressibleNavierStokes::showProgress(const double time, const std::size_t timeStepID)
   {
	   std::cout << "Time is " << time << std::endl;
	   std::cout << "timeStepID is " << timeStepID << std::endl;
   }




































