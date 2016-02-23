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

#include "DoNotScaleIntegrands.h"

namespace Base
{

    template<std::size_t DIM>
    double DoNotScaleIntegrands<DIM>::transform(double referenceData,
                                                PhysicalElement<DIM> &element) const
    {
        return underlying_->transform(referenceData, element);
    }

    template<std::size_t DIM>
    LinearAlgebra::SmallVector<DIM> DoNotScaleIntegrands<DIM>::transform(LinearAlgebra::SmallVector<DIM> referenceData,
                                                                         PhysicalElement<DIM> &element) const
    {
        return underlying_->transform(referenceData, element);
    }

    template<std::size_t DIM>
    LinearAlgebra::SmallVector<DIM> DoNotScaleIntegrands<DIM>::transformDeriv(
            LinearAlgebra::SmallVector<DIM> referenceData, PhysicalElement<DIM> &element) const
    {
        return underlying_->transformDeriv(referenceData, element);
    }

    template<std::size_t DIM>
    double DoNotScaleIntegrands<DIM>::transformDiv(double referenceData,
                                                   PhysicalElement<DIM> &element) const
    {
        return underlying_->transformDiv(referenceData, element);
    }

    template<std::size_t DIM>
    LinearAlgebra::SmallVector<DIM> DoNotScaleIntegrands<DIM>::transformCurl(
            LinearAlgebra::SmallVector<DIM> referenceData, PhysicalElement<DIM> &element) const
    {
        return underlying_->transformCurl(referenceData, element);
    }

    template<std::size_t DIM>
    double DoNotScaleIntegrands<DIM>::getIntegrandScaleFactor(PhysicalElement<DIM> &element) const
    {
        return 1.;
    }

    template<std::size_t DIM>
    double DoNotScaleIntegrands<DIM>::getIntegrandScaleFactor(PhysicalFace<DIM> &face) const
    {
        return 1.;
    }

    template<std::size_t DIM>
    const CoordinateTransformation<DIM> *DoNotScaleIntegrands<DIM>::getUnderlying() const
    {
        return underlying_;
    }
}


BOOST_CLASS_EXPORT_IMPLEMENT(Base::DoNotScaleIntegrands<1>);
BOOST_CLASS_EXPORT_IMPLEMENT(Base::DoNotScaleIntegrands<2>);
BOOST_CLASS_EXPORT_IMPLEMENT(Base::DoNotScaleIntegrands<3>);
BOOST_CLASS_EXPORT_IMPLEMENT(Base::DoNotScaleIntegrands<4>);