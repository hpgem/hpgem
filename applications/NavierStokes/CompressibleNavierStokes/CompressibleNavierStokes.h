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

#ifndef COMPRESSIBLENAVIERSTOKES_H_
#define COMPRESSIBLENAVIERSTOKES_H_

#include "Base/HpgemAPISimplified.h"

class CompressibleNavierStokes;

#include "Inviscid.h"
#include "Viscous.h"

class CompressibleNavierStokes : public Base::HpgemAPISimplified
{
public:
	CompressibleNavierStokes
	(
		const std::size_t dimension,
		const std::size_t numOfVariables,
		const double endTime,
		const std::size_t polynomialOrder,
		const Base::ButcherTableau * const ptrButcherTableau
	);

    /// \brief Create a domain
    Base::RectangularMeshDescriptor createMeshDescription(const std::size_t numOfElementPerDirection) override;

    void setStabilityMassMatrix();

    // Computes pressure for a given solution q
    double computePressure(const LinearAlgebra::NumericalVector &qSolution);


    /// *****************************************
    /// ***   Element integration functions   ***
    /// *****************************************

    /// Compute source function at an element
    LinearAlgebra::NumericalVector integrandSourceAtElement(const Base::Element *ptrElement, const LinearAlgebra::NumericalVector qSolution, const double pressureTerm, const double &time, const Geometry::PointReference &pRef);

    /// Compute solution at an element
    LinearAlgebra::NumericalVector computeSolutionOnElement(const Base::Element *ptrElement, const LinearAlgebra::NumericalVector &solutionCoefficients, const Geometry::PointReference &pRef);

    /// Compute solution derivatives at an element. Energy state is not included.
    LinearAlgebra::Matrix computeSolutionJacobianAtElement(const Base::Element *ptrElement, const LinearAlgebra::NumericalVector &solutionCoefficients, const Geometry::PointReference &pRef);

    /// Compute the Jacobian of the velocities
    LinearAlgebra::Matrix computePartialStateJacobian(const LinearAlgebra::Matrix qSolutionGradient, const LinearAlgebra::NumericalVector qSolution);

    /// Computes the partial States: all states except the density are divided by the density
    LinearAlgebra::NumericalVector computePartialState(const LinearAlgebra::NumericalVector qSolution);

    /// Compute integrand of righthandside on an element
    LinearAlgebra::NumericalVector integrandRightHandSideOnRefElement(const Base::Element *ptrElement, const double &time, const Geometry::PointReference &pRef, const LinearAlgebra::NumericalVector &solutionCoefficients);

    /// \brief Compute the right hand side on an element
    LinearAlgebra::NumericalVector computeRightHandSideAtElement(Base::Element *ptrElement,	LinearAlgebra::NumericalVector &solutionCoefficients, const double time) override;

    /// *****************************************
    /// ***    face integration functions     ***
    /// *****************************************

    /// Compute solution at a face
    LinearAlgebra::NumericalVector computeSolutionOnFace(const Base::Face *ptrFace, const Base::Side &iSide, const LinearAlgebra::NumericalVector &solutionCoefficients, const Geometry::PointReference &pRef) const;

    /// Compute solution Jacobian on face
    LinearAlgebra::Matrix computeSolutionJacobianAtFace(const Base::Face *ptrFace, const Base::Side &iSide, const LinearAlgebra::NumericalVector &solutionCoefficients, const Geometry::PointReference &pRef);

    /// \brief Compute the integrand for the right hand side for the reference face corresponding to a boundary face.
    LinearAlgebra::NumericalVector integrandRightHandSideOnRefFace(const Base::Face *ptrFace, const double &time, const Geometry::PointReference &pRef, const LinearAlgebra::NumericalVector &solutionCoefficients);

    /// \brief Compute the integrand for the right hand side for the reference face corresponding to an internal face.
    LinearAlgebra::NumericalVector integrandRightHandSideOnRefFace(const Base::Face *ptrFace, const double &time, const Geometry::PointReference &pRef, const Base::Side &iSide, const LinearAlgebra::NumericalVector &solutionCoefficientsLeft, const LinearAlgebra::NumericalVector &solutionCoefficientsRight);

    /// \brief Compute the right-hand side corresponding to a boundary face
    LinearAlgebra::NumericalVector computeRightHandSideAtFace(Base::Face *ptrFace, LinearAlgebra::NumericalVector &solutionCoefficients, const double time) override final;

    /// \brief Compute the right-hand side corresponding to an internal face
    LinearAlgebra::NumericalVector computeRightHandSideAtFace(Base::Face *ptrFace, const Base::Side side, LinearAlgebra::NumericalVector &solutionCoefficientsLeft, LinearAlgebra::NumericalVector &solutionCoefficientsRight, const double time) override final;


    /// *****************************************
    /// ***    		Various Functions         ***
    /// *****************************************

    LinearAlgebra::NumericalVector getExactSolution(const PointPhysicalT &pPhys, const double &time, const std::size_t orderTimeDerivative) override final;

    /// \brief Compute the initial solution at a given point in space and time.
	LinearAlgebra::NumericalVector getInitialSolution(const PointPhysicalT &pPhys, const double &startTime, const std::size_t orderTimeDerivative = 0) override final;

	LinearAlgebra::NumericalVector Error(const double time);

private:
    /// Dimension of the domain
    const std::size_t DIM_;

    /// Specific heat ratio
    const double gamma_ = 1.4;

	/// Number of variables
	const std::size_t numOfVariables_;

	/// Inviscid class, treating the inviscid part of the NS equations
	Inviscid inviscidTerms_;

	/// Viscous class, treating the viscosity part of the NS equations
	Viscous viscousTerms_;

	friend class Inviscid;
	friend class Viscous;


};

#endif /* COMPRESSIBLENAVIERSTOKES_H_ */
