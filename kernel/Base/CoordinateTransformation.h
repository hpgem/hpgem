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

#ifndef COORDINATETRANSFORMATION_H_
#define COORDINATETRANSFORMATION_H_

#include <cstdlib>
#include "LinearAlgebra/SmallVector.h"
#include "PhysicalElement.h"
#include "PhysicalFace.h"

namespace Base
{
    ///base class for coordinate transformations. All functions are implemented to generate errors so specializations only need to override the functions they need
    template<std::size_t DIM>
    class CoordinateTransformation
    {
    public:
        virtual ~CoordinateTransformation() = default;

        virtual double transform(double referenceData, PhysicalElement<DIM>& element) const
        {
            logger(ERROR, "Transforming scalar data is not supported, please set a different transformation");
            return 0.;
        }

        virtual LinearAlgebra::SmallVector<DIM> transform(LinearAlgebra::SmallVector<DIM> referenceData, PhysicalElement<DIM>& element) const
        {
            logger(ERROR, "Transforming vector data is not supported, please set a different transformation");
            return LinearAlgebra::SmallVector<DIM>();
        }

        virtual LinearAlgebra::SmallVector<DIM> transformDeriv(LinearAlgebra::SmallVector<DIM> referenceData, PhysicalElement<DIM>& element) const
        {
            logger(ERROR, "Transforming derivative data is not supported, please set a different transformation");
            return LinearAlgebra::SmallVector<DIM>();
        }

        virtual LinearAlgebra::SmallVector<DIM> transformCurl(LinearAlgebra::SmallVector<DIM> referenceData, PhysicalElement<DIM>& element) const
        {
            logger(ERROR, "Transforming curl data is not supported, please set a different transformation");
            return LinearAlgebra::SmallVector<DIM>();
        }

        virtual double getIntegrandScaleFactor(PhysicalElement<DIM>& element) const
        {
            logger(ERROR, "Scaling integrands is not supported, please set a different transformation");
            return 0.;
        }

        virtual double getIntegrandScaleFactor(PhysicalFace<DIM>& face) const
        {
            logger(ERROR, "Scaling integrands is not supported, please set a different transformation");
            return 0.;
        }
    };
}


#endif /* COORDINATETRANSFORMATION_H_ */
