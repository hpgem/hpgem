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

#include "Geometry/PointReference.h"

namespace Geometry
{
    double PointReference::getBasisFunctionValue(const Base::BaseBasisFunction* function) const
    {

        logger.assert(function!=nullptr, "Invalid basis function passed");
        try
        {
            return basisfunctionValues_.at(function);
        }
        catch (std::out_of_range&)
        {
            const_cast<std::unordered_map<const Base::BaseBasisFunction*, double>&>(basisfunctionValues_)[function] = function->eval(*this);
            return basisfunctionValues_.at(function);
        }
    }

    const LinearAlgebra::NumericalVector& PointReference::getBasisFunctionDerivative(const Base::BaseBasisFunction* function) const
    {

        logger.assert(function!=nullptr, "Invalid basis function passed");
        try
        {
            return basisfunctionDerivatives_.at(function);
        }
        catch (std::out_of_range&)
        {
            const_cast<std::unordered_map<const Base::BaseBasisFunction*, LinearAlgebra::NumericalVector>&>(basisfunctionDerivatives_)[function] = function->evalDeriv(*this);
            return basisfunctionDerivatives_.at(function);
        }
    }


}

