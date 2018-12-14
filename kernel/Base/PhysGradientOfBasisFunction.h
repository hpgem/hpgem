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

#ifndef ____PhysGradientOfBasisFunction__
#define ____PhysGradientOfBasisFunction__

#include "Logger.h"

namespace LinearAlgebra
{
    template<std::size_t DIM>
    class SmallVector;
}

namespace Geometry
{
    template<std::size_t DIM>
    class PointReference;
}

namespace Base
{
    class Element;
    class BaseBasisFunction;
}

namespace Utilities
{
    /*! For a basis function we also need its physical space gradient as opposed
     *  to the one in reference space (which can be harvested by evaluating the
     *  derivatives of basis functions). Hence this class computes the
     *  reference space gradient, transforms it with the Jacobian of the mapping
     *  and thus yields the physical space gradient.
     *  \deprecated functionality is specific for H1 conforming basisfunctions*/
    struct PhysGradientOfBasisFunction
    {

        PhysGradientOfBasisFunction(const Base::Element* e, const Base::BaseBasisFunction* function)
                : myElement_(e), myFunction_(function)
        {
            logger.assert_debug(e != nullptr, "Invalid element passed");
            logger.assert_debug(function != nullptr, "Invalid function passed");
        }
        
        //! Evaluation operator, also compatible with integration routines.
        template<std::size_t DIM>
        LinearAlgebra::SmallVector<DIM> operator ()(const Geometry::PointReference<DIM>& p) const;

    private:
        const Base::Element* myElement_;
        const Base::BaseBasisFunction* myFunction_;
    };

} // namespace

#include "PhysGradientOfBasisFunction_Impl.h"

#endif /* defined(____PhysGradientOfBasisFunction__) */
