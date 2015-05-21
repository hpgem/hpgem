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

#ifndef PointReference_h
#define PointReference_h

#include "Point.h"
#include "PointReferenceBase.h"
#include "PointReferenceFactory.h"
#include "Base/BaseBasisFunction.h"

#include <unordered_map>

namespace Geometry
{
    template<std::size_t DIM>
    class PointReferenceFactory;

    template<std::size_t DIM>
    class PointReference : public Point<DIM>, public PointReferenceBase
    {
    public:
        void removeBasisFunctionData(const Base::BaseBasisFunction* function)
        {
            basisfunctionValues_.erase(function);
            basisfunctionDerivatives_.erase(function);
        }

        double getBasisFunctionValue(const Base::BaseBasisFunction* function) const;
        const LinearAlgebra::SmallVector<DIM>& getBasisFunctionDerivative(const Base::BaseBasisFunction* function) const;
        //do not trust any other class to not create duplicates
        friend PointReferenceFactory<DIM>;
    private:
        
        PointReference()
                : Point<DIM>()
        {
        }
        
        //do not copy a pointReference, its memory address is used to quickly collect precomputed values of basis functions
        PointReference(const PointReference& p) = delete;
        PointReference(PointReference&& p) = delete;

        explicit PointReference(const Point<DIM>& p)
                : Point<DIM>(p)
        {
        }
        
        PointReference(std::initializer_list<double> data)
                : Point<DIM>(data)
        {
        }

        PointReference(double coords[])
                : Point<DIM>(coords)
        {
        }
        
        explicit PointReference(const LinearAlgebra::SmallVector<DIM>& coord)
                : Point<DIM>(coord)
        {
        }

        PointReference& operator =(const PointReference& rhs) = delete;
        PointReference& operator =(PointReference&& rhs) = delete;

        std::unordered_map<const Base::BaseBasisFunction*, double > basisfunctionValues_;
        std::unordered_map<const Base::BaseBasisFunction*, LinearAlgebra::SmallVector<DIM> > basisfunctionDerivatives_;
        
    };

    //PointReference operator*(double left, const PointReference& right);
    

}

#include "PointReference.cpp"

#endif /* POINTREFERENCE_HPP_ */
