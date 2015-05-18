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
    class MiddleSizeVector;
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
        using BaseBasisFunctions = std::vector<BaseBasisFunction*>; //check again
        using PointReferenceT = Geometry::PointReference;

        explicit BasisFunctionSet(std::size_t order);      
        
        //BasisFunctionSets should not be copied, therefore the copy constructor is deleted.
        BasisFunctionSet(const BasisFunctionSet& other) = delete;

        virtual ~BasisFunctionSet();

        std::size_t size() const;

        std::size_t getOrder() const;

        void addBasisFunction(BaseBasisFunction* bf);

        double eval(std::size_t i, const PointReferenceT& p) const;

        ///\brief returns the value of the i-th basisfunction at point p in ret
        void eval(std::size_t i, const PointReferenceT& p, LinearAlgebra::MiddleSizeVector& ret) const;

        double evalDeriv(std::size_t i, std::size_t jDir, const PointReferenceT& p) const;

        ///\brief returns the curl of the i-th basisfunction at point p in ret
        LinearAlgebra::MiddleSizeVector evalCurl(std::size_t i, const PointReferenceT& p) const;

        const BaseBasisFunction* operator[](std::size_t i) const
        {
            logger.assert(i<size(), "Asked for basis function %, but there are only % basis functions", i, size());
            return vecOfBasisFcn_[i];
        }

        ///iterators (for range-based for loop)
        BaseBasisFunctions::const_iterator begin() const
        {
            return vecOfBasisFcn_.begin();
        }

        BaseBasisFunctions::iterator begin()
        {
            return vecOfBasisFcn_.begin();
        }

        BaseBasisFunctions::const_iterator end() const
        {
            return vecOfBasisFcn_.end();
        }

        BaseBasisFunctions::iterator end()
        {
            return vecOfBasisFcn_.end();
        }

    private:
        std::size_t order_;
        BaseBasisFunctions vecOfBasisFcn_;
    };
}

#endif
