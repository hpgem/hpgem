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

#ifndef HDIVCONFORMINGTRANSFORMATION_H_
#define HDIVCONFORMINGTRANSFORMATION_H_

#include <cstdlib>
#include "LinearAlgebra/SmallVector.h"
#include "CoordinateTransformation.h"

namespace Base {
/// transforms vector functions and their divergence in a conforming way
template <std::size_t DIM>
class HDivConformingTransformation : public CoordinateTransformation<DIM> {
   public:
    /// transform functions as if they are the curl of some other function. This
    /// will exactly map the kernel of the physical div-operator to the kernel of
    /// the reference div-operator.
    LinearAlgebra::SmallVector<DIM> transform(
        LinearAlgebra::SmallVector<DIM> referenceData,
        PhysicalElement<DIM>& element) const override final {
        return element.getJacobian() * referenceData / element.getJacobianDet();
    }

    /// transform the div by using the chain rule
    double transformDiv(double referenceData,
                        PhysicalElement<DIM>& element) const override final {
        return referenceData / element.getJacobianDet();
    }

    /// integrands for elements are multiplied by the absolute value of the
    /// determinant of the Jacobian to correct for the difference in volume
    double getIntegrandScaleFactor(
        PhysicalElement<DIM>& element) const override final {
        return element.getJacobianAbsDet();
    }

    /// integrands for faces are multiplied by the norm of the outward normal
    /// vector to correct for the difference in area
    double getIntegrandScaleFactor(
        PhysicalFace<DIM>& face) const override final {
        return face.getRelativeSurfaceArea();
    }
};
}  // namespace Base

#endif /* HDIVCONFORMINGTRANSFORMATION_H_ */
