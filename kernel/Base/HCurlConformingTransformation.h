/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2014, University of Twente
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.

 3. Neither the name of the copyright holder nor the names of its contributors
 may be used to endorse or promote products derived from this software without
 specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef HPGEM_KERNEL_HCURLCONFORMINGTRANSFORMATION_H
#define HPGEM_KERNEL_HCURLCONFORMINGTRANSFORMATION_H

#include <cstdlib>
#include "LinearAlgebra/SmallVector.h"
#include "CoordinateTransformation.h"

namespace hpgem {

namespace Base {
/// transforms vector functions and their curl in a conforming way
template <std::size_t DIM>
class HCurlConformingTransformation : public CoordinateTransformation<DIM> {
   public:
    /// transform functions as if they are the gradient of some other function.
    /// This will exactly map the kernel of the physical curl-operator to the
    /// kernel of the reference curl-operator.
    LinearAlgebra::SmallVector<DIM> transform(
        LinearAlgebra::SmallVector<DIM> referenceData,
        const CoordinateTransformationData<DIM>& data) const final {
        data.getTransposeJacobian().solve(referenceData);
        return referenceData;
    }

    /// transform the curl by using the chain rule
    LinearAlgebra::SmallVector<DIM> transformCurl(
        LinearAlgebra::SmallVector<DIM> referenceData,
        const CoordinateTransformationData<DIM>& data) const final {
        return data.getJacobian() * referenceData / data.getJacobianDet();
    }
};

template <>
inline LinearAlgebra::SmallVector<2>
    HCurlConformingTransformation<2>::transformCurl(
        LinearAlgebra::SmallVector<2> referenceData,
        const CoordinateTransformationData<2>& data) const {
    return referenceData / data.getJacobianDet();
}
}  // namespace Base

}  // namespace hpgem

#endif  // HPGEM_KERNEL_HCURLCONFORMINGTRANSFORMATION_H
