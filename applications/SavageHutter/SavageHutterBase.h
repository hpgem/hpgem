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

#include "RightHandSideComputer.h"
#include "SlopeLimiters/SlopeLimiter.h"
#include "HeightLimiters/HeightLimiter.h"
#include "Base/HpgemAPISimplified.h"

/// \param[in] numberOfVariables Number of variables in the PDE
/// \param[in] polynomialOrder Polynomial order of the basis functions
/// \param[in] useMatrixStorage Boolean to indicate if element and face matrices for the PDE should be stored
/// \param[in] ptrButcherTableau Pointer to a Butcher Tableau used to do the time integration with a Runge-Kutta scheme. By default this is a RK4 scheme.
struct SHConstructorStruct
{
    std::size_t numOfVariables;
    std::size_t polyOrder;
    std::size_t numElements;
    Base::MeshType meshType;
    TimeIntegration::ButcherTableau * ptrButcherTableau;
};

class SavageHutterBase : public Base::HpgemAPISimplified<DIM>
{
public:
    SavageHutterBase(const SHConstructorStruct & inputValues);

    virtual ~SavageHutterBase()
    {
        delete rhsComputer_;
        delete slopeLimiter_;
        delete heightLimiter_;
    }
    
    std::vector<std::pair<double, LinearAlgebra::MiddleSizeVector>> widthAverage();

protected:

    /// \brief Create a domain
    Base::RectangularMeshDescriptor<DIM> createMeshDescription(const std::size_t numOfElementPerDirection) override final;
    
    /// Number of variables
    const std::size_t numOfVariables_;

    RightHandSideComputer* rhsComputer_;

    SlopeLimiter* slopeLimiter_;

    HeightLimiter* heightLimiter_;

    ///If the minimum height in an element is below this number, the element is considered to be dry.
    const double dryLimit_;
    
    
    double time_;

private:
    virtual SlopeLimiter * createSlopeLimiter(const SHConstructorStruct &inputValues) = 0;
    virtual HeightLimiter * createHeightLimiter(const SHConstructorStruct &inputValues) = 0;
    virtual RightHandSideComputer * createRightHandSideComputer(const SHConstructorStruct &inputValues) = 0;
    
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

    void computeOneTimeStep(double &time, const double dt) override final;
    void limitSolutionOuterLoop();
    void limitSolutionInnerLoop();

    ///Compute the minimum of the height in the given element
    double getMinimumHeight(const Base::Element *element);

};

#endif	/* SAVAGEHUTTERBASE_H */

