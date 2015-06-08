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

#ifndef EULER_H_
#define EULER_H_

#include "Base/HpgemAPISimplified.h"

template<std::size_t DIM>
class Euler : public Base::HpgemAPISimplified<DIM>
{
public:

    using typename Base::HpgemAPIBase<DIM>::PointPhysicalT;
    using typename Base::HpgemAPIBase<DIM>::PointReferenceT;
    using typename Base::HpgemAPIBase<DIM>::PointReferenceOnFaceT;

	Euler
	(
	 const std::size_t numOfVariables,
	 const double endTime,
	 const std::size_t polynomialOrder,
	 const Base::ButcherTableau * const ptrButcherTableau
	);

    /// \brief Create a domain
    Base::RectangularMeshDescriptor<DIM> createMeshDescription(const std::size_t numOfElementPerDirection) override;


    /// *****************************************
    /// ***   Element integration functions   ***
    /// *****************************************

    /// Compute source function at an element
    LinearAlgebra::MiddleSizeVector integrandSourceAtElement(Base::PhysicalElement<DIM> &ptrElement, const LinearAlgebra::MiddleSizeVector qSolution, const double pressureTerm, const double &time);

    /// Compute solution at an element
    LinearAlgebra::MiddleSizeVector computeSolutionAtElement(const Base::Element *ptrElement, const LinearAlgebra::MiddleSizeVector &solutionCoefficients, const PointReferenceT &pRef);

    /// Compute integrand of righthandside on an element
    LinearAlgebra::MiddleSizeVector integrandRightHandSideOnRefElement(Base::PhysicalElement<DIM> &ptrElement, const double &time, const LinearAlgebra::MiddleSizeVector &solutionCoefficients);

    /// \brief Compute the right hand side on an element
    LinearAlgebra::MiddleSizeVector computeRightHandSideAtElement(Base::Element *ptrElement,	LinearAlgebra::MiddleSizeVector &solutionCoefficients, const double time) override;

    /// *****************************************
    /// ***    face integration functions     ***
    /// *****************************************

    /// \brief Compute the Roe Riemann Flux.
    LinearAlgebra::MiddleSizeVector RoeRiemannFluxFunction(const LinearAlgebra::MiddleSizeVector &qReconstructionLeft, const LinearAlgebra::MiddleSizeVector &qReconstructionRight, LinearAlgebra::SmallVector<DIM> &normal);

    /// \brief Compute the integrand for the right hand side for the reference face corresponding to a boundary face.
    LinearAlgebra::MiddleSizeVector integrandRightHandSideOnRefFace(Base::PhysicalFace<DIM>& face, const double &time, const LinearAlgebra::MiddleSizeVector &solutionCoefficients);

    /// \brief Compute the integrand for the right hand side for the reference face corresponding to an internal face.
    LinearAlgebra::MiddleSizeVector integrandRightHandSideOnRefFace(Base::PhysicalFace<DIM>& face, const double &time, const Base::Side &iSide, const LinearAlgebra::MiddleSizeVector &solutionCoefficientsLeft, const LinearAlgebra::MiddleSizeVector &solutionCoefficientsRight);

    /// \brief Compute the right-hand side corresponding to a boundary face
    LinearAlgebra::MiddleSizeVector computeRightHandSideAtFace(Base::Face *ptrFace, LinearAlgebra::MiddleSizeVector &solutionCoefficients, const double time) override final;

    /// \brief Compute the right-hand side corresponding to an internal face
    LinearAlgebra::MiddleSizeVector computeRightHandSideAtFace(Base::Face *ptrFace, const Base::Side side, LinearAlgebra::MiddleSizeVector &solutionCoefficientsLeft, LinearAlgebra::MiddleSizeVector &solutionCoefficientsRight, const double time) override final;


    /// *****************************************
    /// ***    		Various Functions         ***
    /// *****************************************

    void tasksBeforeSolving() override final;

    LinearAlgebra::MiddleSizeVector getExactSolution(const PointPhysicalT &pPhys, const double &time, const std::size_t orderTimeDerivative) override final;

    /// \brief Compute the initial solution at a given point in space and time.
	LinearAlgebra::MiddleSizeVector getInitialSolution(const PointPhysicalT &pPhys, const double &startTime, const std::size_t orderTimeDerivative = 0) override final;

	LinearAlgebra::MiddleSizeVector Error(const double time);


private:

    /// Specific heat ratio
    const double gamma_ = 1.4;

	/// Number of variables
	const std::size_t numOfVariables_;
};

#include "Euler_Impl.h"

#endif /* EULER_H_ */
