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

#ifndef SAVAGEHUTTERBASE_H
#define	SAVAGEHUTTERBASE_H

#include "SlopeLimiters/SlopeLimiter.h"
#include "HeightLimiters/HeightLimiter.h"
#include "Base/HpgemAPISimplified.h"
#include "SlopeLimiters/EmptySlopeLimiter.h"
#include "HeightLimiters/EmptyHeightLimiter.h"


template <std::size_t DIM>
class SavageHutterBase : public Base::HpgemAPISimplified<DIM>
{
public:
    using MiddleSizeVector = LinearAlgebra::MiddleSizeVector;
    
    SavageHutterBase(std::size_t numberOfVariables, std::size_t polyOrder);

    virtual ~SavageHutterBase()
    {
        delete slopeLimiter_;
        delete heightLimiter_;
    }    
    
    ///\brief Show the progress of the time integration.
    void showProgress(const double time, const std::size_t timeStepID) override final;

protected:

    /// \brief Purely virtual function to compute the integrand for the right hand side for the reference element.
    virtual const MiddleSizeVector integrandRightHandSideOnElement
    (
        Base::PhysicalElement<DIM> &element,
        const double &time,
        const MiddleSizeVector &solutionCoefficients
        ) = 0;

    /// \brief Purely virtual function to compute the integrand for the right hand side for the reference face corresponding to a boundary face.
    virtual const MiddleSizeVector integrandRightHandSideOnRefFace
    (
        Base::PhysicalFace<DIM> &face,
        const MiddleSizeVector &solutionCoefficients,
        const double time
        ) = 0;

    /// \brief Purely virtual function to compute the integrand for the right hand side for the reference face corresponding to an internal face.
    virtual const MiddleSizeVector integrandRightHandSideOnRefFace
    (
        Base::PhysicalFace<DIM> &face,
        const Base::Side &iSide,
        const MiddleSizeVector &solutionCoefficientsLeft,
        const MiddleSizeVector &solutionCoefficientsRight
            ) = 0;



    virtual std::pair<LinearAlgebra::MiddleSizeVector,LinearAlgebra::MiddleSizeVector> integrandsAtFace(
		Base::PhysicalFace<DIM> &face,
		const double &time,
		const LinearAlgebra::MiddleSizeVector &stateLeft,
		const LinearAlgebra::MiddleSizeVector &stateRight) = 0;

    void tasksBeforeSolving() override
    {
        slopeLimiter_ = createSlopeLimiter();
        heightLimiter_ = createHeightLimiter();
    }

    void tasksAfterTimeStep()
    {
        limitSolutionOuterLoop();
    }
    
    ///\brief Creates an empty slope limiter.
    ///\details This function creates an empty slope limiter, but can be overwritten
    ///by an application in order to create the slope limiter you want. This is not 
    ///available yet in 2D.
    virtual SlopeLimiter* createSlopeLimiter()
    {
        return new EmptySlopeLimiter;
    }
    
    ///\brief Creates an empty height limiter.
    ///\details This function creates an empty height limiter, but can be overwritten
    ///by an application in order to create the non-negativity limiter you want.
    virtual HeightLimiter* createHeightLimiter()
    {
        return new EmptyHeightLimiter;
    }
    
    /// Number of variables
    const std::size_t numberOfVariables_;

    SlopeLimiter* slopeLimiter_;

    HeightLimiter* heightLimiter_;

    ///If the minimum height in an element is below this number, the element is considered to be dry.
    double dryLimit_;
    
    double time_;
    
    double epsilon_;
    
    double chuteAngle_;

    MiddleSizeVector inflowBC_;

private:
    
    virtual void setInflowBC(double time) {  }
    
    LinearAlgebra::MiddleSizeVector computeRightHandSideAtElement(Base::Element *ptrElement, 
    LinearAlgebra::MiddleSizeVector &solutionCoefficients, 
    const double time) override final;

    /// \brief Compute the right-hand side corresponding to an internal face
    LinearAlgebra::MiddleSizeVector computeRightHandSideAtFace(Base::Face *ptrFace,
        const Base::Side side,
        LinearAlgebra::MiddleSizeVector &solutionCoefficientsLeft,
        LinearAlgebra::MiddleSizeVector &solutionCoefficientsRight,
        const double time) override final;

    LinearAlgebra::MiddleSizeVector computeRightHandSideAtFace
    (
        Base::Face *ptrFace,
        LinearAlgebra::MiddleSizeVector &solutionCoefficients,
        const double time
        ) override final;

    std::pair<LinearAlgebra::MiddleSizeVector,LinearAlgebra::MiddleSizeVector> computeBothRightHandSidesAtFace
        (
         Base::Face *ptrFace,
         LinearAlgebra::MiddleSizeVector &solutionCoefficientsLeft,
         LinearAlgebra::MiddleSizeVector &solutionCoefficientsRight,
         const double time
         ) override final;

    void computeOneTimeStep(double &time, const double dt) override final;
    void limitSolutionOuterLoop();
    void limitSolutionInnerLoop();

    ///Compute the minimum of the height in the given element
    const double getMinimumHeight(const Base::Element *element);

};

#include "SavageHutterBase_Impl.h"
#endif	/* SAVAGEHUTTERBASE_H */

