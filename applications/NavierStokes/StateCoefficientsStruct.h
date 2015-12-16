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

#ifndef STATECOEFFICIENTSSTRUCT_H_
#define STATECOEFFICIENTSSTRUCT_H_

#include "Base/PhysicalElement.h"
#include "Base/PhysicalFace.h"
#include "LinearAlgebra/MiddleSizeVector.h"
#include "LinearAlgebra/MiddleSizeMatrix.h"
#include "StateCoefficientsFunctions.h"
#include <vector>
#include <limits>

//todo: Optimise this code

template<std::size_t DIM, std::size_t NUMBER_OF_VARIABLES>
class StateCoefficientsStruct {
public:

	/// \brief Constructor to obtain constitutive relations
	StateCoefficientsStruct()
	{
	}

	/// \brief Constructor for solution functions at an element
	StateCoefficientsStruct(
			Base::PhysicalElement<DIM> &element,
			const LinearAlgebra::MiddleSizeVector &stateCoefficients,
			const double time) :
				stateCoefficients_(stateCoefficients),
				state_(computeStateOnElement<DIM,NUMBER_OF_VARIABLES>(element, stateCoefficients_)),
				stateJacobian_(computeStateJacobianAtElement<DIM,NUMBER_OF_VARIABLES>(element, stateCoefficients_)),
				partialState_(computePartialState<DIM,NUMBER_OF_VARIABLES>(state_))
	{
	}

	/// \brief Constructor for solution functions at a face
	StateCoefficientsStruct(
			Base::PhysicalFace<DIM> &face,
			const LinearAlgebra::MiddleSizeVector &stateCoefficients,
			const Base::Side side,
			const double time) :
				stateCoefficients_(stateCoefficients),
				state_(computeStateOnFace<DIM,NUMBER_OF_VARIABLES>(face, side, stateCoefficients_)),
				stateJacobian_(computeStateJacobianAtFace<DIM,NUMBER_OF_VARIABLES>(face, side, stateCoefficients_)),
				partialState_(computePartialState<DIM,NUMBER_OF_VARIABLES>(state_))
	{
	}

	virtual ~StateCoefficientsStruct()
	{
	}

	/// ***************************************
	/// ***      Constutive relations       ***
	/// ***************************************

	virtual double computePressure(const LinearAlgebra::MiddleSizeVector &state) const
	{
		logger(ERROR, "I want to compute the pressure, but I am not yet implemented!");
		double result = 15;
		return result;
	}

	virtual double computeSpeedOfSound(const LinearAlgebra::MiddleSizeVector &state, const double pressure) const
	{
		logger(ERROR, "I want to compute the speed of sound, but I am not yet implemented!");
		double result = 15;
		return result;
	}

	virtual LinearAlgebra::MiddleSizeMatrix computeHyperbolicMatrix(const LinearAlgebra::MiddleSizeVector &state, const double pressure)
	{
		logger(ERROR, "I want to compute the hyperbolic flux matrix, but I am not yet implemented!");
		LinearAlgebra::MiddleSizeMatrix matrix;
		return matrix;
	}

	virtual double computeViscosity(const LinearAlgebra::MiddleSizeVector &state, const LinearAlgebra::MiddleSizeVector &partialState, const LinearAlgebra::MiddleSizeMatrix stateJacobian, const double pressure)
	{
		logger(ERROR, "I want to compute the viscosity, but I am not yet implemented!");
		double viscosity = 30;
		return viscosity;
	}

	virtual std::vector<LinearAlgebra::MiddleSizeMatrix> computeEllipticTensor(const LinearAlgebra::MiddleSizeVector partialState, const double viscosity)
	{
		logger(ERROR, "I want to compute the elliptic tensor, but I am not yet implemented!");
		std::vector<LinearAlgebra::MiddleSizeMatrix> ellipticTensor;
		return ellipticTensor;
	}

	/// ***************************************
	/// ***      Additional Functions       ***
	/// ***************************************

	//Note: This is a naive implementation of the elliptic tensor matrix times matrix computation. You can specify your own, faster one in the class.
	virtual LinearAlgebra::MiddleSizeMatrix computeEllipticTensorMatrixContractionFast(
			const LinearAlgebra::MiddleSizeMatrix &matrix) const
	{
		return computeEllipticTensorMatrixContraction<DIM,NUMBER_OF_VARIABLES>(ellipticTensor_, matrix);
	}

    /// ********************************
    /// ***      Get Functions       ***
    /// ********************************

	LinearAlgebra::MiddleSizeVector getStateCoefficients() const
	{
		return stateCoefficients_;
	}

	LinearAlgebra::MiddleSizeVector getState() const
	{
		return state_;
	}

	LinearAlgebra::MiddleSizeVector getPartialState() const
	{
		return partialState_;
	}

	LinearAlgebra::MiddleSizeMatrix getStateJacobian() const
	{
		return stateJacobian_;
	}

	double getPressure() const
	{
		return pressure_;
	}

	double getSpeedOfSound() const
	{
		return speedOfSound_;
	}

	double getViscosity() const
	{
		return viscosity_;
	}

	LinearAlgebra::MiddleSizeMatrix getHyperbolicMatrix() const
	{
		return hyperbolicMatrix_;
	}

	std::vector<LinearAlgebra::MiddleSizeMatrix> getEllipticTensor() const
	{
		return ellipticTensor_;
	}


protected:
	/// \brief This contains the solutionCoefficients where all other structures are based on
	LinearAlgebra::MiddleSizeVector stateCoefficients_;

	/// \brief This is the state at the given point of the element
	LinearAlgebra::MiddleSizeVector state_;

	/// \brief The statejacobian is the derivative of the states with respect to the coordinate system. i.e. the first row is d/dx rho(x,y), d/dy rho(x,y)
	LinearAlgebra::MiddleSizeMatrix stateJacobian_;

	/// \brief This is the partial state at the given point of the element. By partial state all state variables except the density, is divided by the density (i.e. rho, u , v)
	LinearAlgebra::MiddleSizeVector partialState_;

	/// \brief This is the pressure calculated from the constitutive relations, which is required to calculate the ellipticTensor
	double pressure_ = std::numeric_limits<double>::quiet_NaN();

	/// \brief This is the speed of sound for a given state
	double speedOfSound_ = std::numeric_limits<double>::quiet_NaN();

	/// \brief This is the viscosity calculated from the constitutive relations, which is required to calculate the elliptic tensor
	double viscosity_ = std::numeric_limits<double>::quiet_NaN();

	/// \brief The hyperbolic matrix is used for the hyperbolic part of the Navier-Stokes equations.
	LinearAlgebra::MiddleSizeMatrix hyperbolicMatrix_;

	/// \brief The elliptic tensor is used for the elliptic part of the Navier-Stokes equations.
	std::vector<LinearAlgebra::MiddleSizeMatrix> ellipticTensor_;
};


#endif /* STATECOEFFICIENTSSTRUCT_H_ */
