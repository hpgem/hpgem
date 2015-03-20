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

#ifndef SHORTTERMSTORAGEELEMENTH1_HPP_
#define SHORTTERMSTORAGEELEMENTH1_HPP_

#include "Base/ShortTermStorageElementBase.h"

namespace Base
{
    
    /**
     * H1 conforming specialization does not alter function values and will transform derivatives with a factor J^-T
     * does not store curls or directional derivatives
     */
    class ShortTermStorageElementH1 : public ShortTermStorageElementBase
    {
    public:
        
        ShortTermStorageElementH1(std::size_t dimension)
                : ShortTermStorageElementBase(dimension)
        {
        }
        
        void computeData() override;

        double basisFunction(std::size_t i, const PointReferenceT& p) override final;
        double basisFunction(std::size_t i, const PointReferenceT& p) const override final;

        void basisFunction(std::size_t i, const PointReferenceT& p, LinearAlgebra::NumericalVector& ret) override final;
        void basisFunction(std::size_t i, const PointReferenceT& p, LinearAlgebra::NumericalVector& ret) const override final;

        LinearAlgebra::NumericalVector basisFunctionDeriv(std::size_t i, const PointReferenceT& p, const Element* = nullptr) override final;
        LinearAlgebra::NumericalVector basisFunctionDeriv(std::size_t i, const PointReferenceT& p, const Element* = nullptr) const override final;

        ///special case derivative: compute individual components, then mix and match as desired !warning! this routine assumes the user wants to construct a specialized transformation and will not premultiply by the Jacobian
        double basisFunctionDeriv(std::size_t i, std::size_t jDir, const PointReferenceT& p);
        double basisFunctionDeriv(std::size_t i, std::size_t jDir, const PointReferenceT& p) const override final;

    private:
        
        std::vector<LinearAlgebra::NumericalVector> basisFunctionIndividualDerivatives_;
    };
}

#endif /* SHORTTERMSTORAGEELEMENTH1_HPP_ */
