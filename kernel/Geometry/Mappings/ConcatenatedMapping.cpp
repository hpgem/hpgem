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

#include "ConcatenatedMapping.h"

#include "Geometry/PointReference.h"
#include "Geometry/Jacobian.h"

namespace Geometry {

PointReference<1> ConcatenatedMapping::transform(
    const PointReference<0>& pIn) const {
    return map2_.transform(map1_.transform(pIn));
}

PointReference<2> ConcatenatedMapping::transform(
    const PointReference<1>& pIn) const {
    return map2_.transform(map1_.transform(pIn));
}

PointReference<3> ConcatenatedMapping::transform(
    const PointReference<2>& pIn) const {
    return map2_.transform(map1_.transform(pIn));
}

PointReference<4> ConcatenatedMapping::transform(
    const PointReference<3>& pIn) const {
    return map2_.transform(map1_.transform(pIn));
}

Jacobian<0, 1> ConcatenatedMapping::calcJacobian(
    const PointReference<0>& p) const {
    // degenerate case
    return Jacobian<0, 1>();
}

Jacobian<1, 2> ConcatenatedMapping::calcJacobian(
    const PointReference<1>& p) const {
    Jacobian<1, 1> j1 = map1_.calcJacobian(p);
    Jacobian<1, 2> j2 = map2_.calcJacobian(map1_.transform(p));
    return j2.multiplyJacobiansInto(j1);
}

Jacobian<2, 3> ConcatenatedMapping::calcJacobian(
    const PointReference<2>& p) const {
    Jacobian<2, 2> j1 = map1_.calcJacobian(p);
    Jacobian<2, 3> j2 = map2_.calcJacobian(map1_.transform(p));
    return j2.multiplyJacobiansInto(j1);
}

Jacobian<3, 4> ConcatenatedMapping::calcJacobian(
    const PointReference<3>& p) const {
    Jacobian<3, 3> j1 = map1_.calcJacobian(p);
    Jacobian<3, 4> j2 = map2_.calcJacobian(map1_.transform(p));
    return j2.multiplyJacobiansInto(j1);
}

std::size_t ConcatenatedMapping::getTargetDimension() const {
    return map2_.getTargetDimension();
}

}  // namespace Geometry
