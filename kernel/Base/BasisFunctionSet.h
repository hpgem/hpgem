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

#ifndef BasisFunctionSet_h
#define BasisFunctionSet_h

#include <vector>
#include "Logger.h"

namespace LinearAlgebra
{
    class NumericalVector;
}

namespace Geometry
{
    class PointReference;
}

namespace Base
{
    class BaseBasisFunction;
    
    class BasisFunctionSet
    {
    public:
        using BaseBasisFunctionT = BaseBasisFunction;
        using BaseBasisFunctions = std::vector<BaseBasisFunctionT*>; //check again
        using PointReferenceT = Geometry::PointReference;

    public:
        BasisFunctionSet(std::size_t order);

        virtual ~BasisFunctionSet();

        virtual std::size_t size() const;

        virtual std::size_t getOrder() const;

        virtual void addBasisFunction(BaseBasisFunctionT* bf);

        virtual double eval(std::size_t i, const PointReferenceT& p) const;

        ///\brief returns the value of the i-th basisfunction at point p in ret
        virtual void eval(std::size_t i, const PointReferenceT& p, LinearAlgebra::NumericalVector& ret) const;

        virtual double evalDeriv(std::size_t i, std::size_t jDir, const PointReferenceT& p) const;

        ///\brief returns the curl of the i-th basisfunction at point p in ret
        virtual LinearAlgebra::NumericalVector evalCurl(std::size_t i, const PointReferenceT& p) const;

        virtual const BaseBasisFunction* operator[](int i) const
        {
            logger.assert(i<size(), "Asked for basis function %, but there are only % basis functions", i, size());
            return vecOfBasisFcn_[i];
        }
        
    private:
        BasisFunctionSet();
        BasisFunctionSet(const BasisFunctionSet& other);

    private:
        std::size_t order_;
        BaseBasisFunctions vecOfBasisFcn_;
    };
}

#endif
