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

#include "MappingToRefPointToLine.h"
#include "Geometry/PointReference.h"

namespace hpgem {

namespace Geometry {
/*
 * The reference line:
 *
 * (-1) 0-------1 (+1)
 *
 * Linear map a point into a line. There are only two possible mappings:
 *
 *      index 0: () -> -1.0
 *      index 1: () -> 1.0
 *
 *
 */

// ~~~~~~~~~~~~~~~==============================================================================
// ~~~ index 0
// ~~~==============================================================================
// ~~~~~~~~~~~~~~~==============================================================================
const MappingToRefPointToLine0& MappingToRefPointToLine0::Instance() {
    static const MappingToRefPointToLine0 theInstance;
    return theInstance;
}

PointReference<1> MappingToRefPointToLine0::transform(
    const Geometry::PointReference<0>& p1) const {
    logger.assert_debug(p1.size() == 0,
                        "Reference point has the wrong dimension");
    return {-1.};
}

Jacobian<0, 1> MappingToRefPointToLine0::calcJacobian(
    const Geometry::PointReference<0>& p1) const {
    logger.assert_debug(p1.size() == 0,
                        "Reference point has the wrong dimension");
    return Jacobian<0, 1>();
}

MappingToRefPointToLine0::MappingToRefPointToLine0()
    : BoundaryFaceMapping(-1.0) {}

// ~~~~~~~~~~~~~~~==============================================================================
// ~~~ index 1
// ~~~==============================================================================
// ~~~~~~~~~~~~~~~==============================================================================

const MappingToRefPointToLine1& MappingToRefPointToLine1::Instance() {
    static const MappingToRefPointToLine1 theInstance;
    return theInstance;
}

PointReference<1> MappingToRefPointToLine1::transform(
    const Geometry::PointReference<0>& p1) const {
    logger.assert_debug(p1.size() == 0,
                        "Reference point has the wrong dimension");
    return {1.};
}

Jacobian<0, 1> MappingToRefPointToLine1::calcJacobian(
    const Geometry::PointReference<0>& p1) const {
    logger.assert_debug(p1.size() == 0,
                        "Reference point has the wrong dimension");
    return Jacobian<0, 1>();
}

MappingToRefPointToLine1::MappingToRefPointToLine1()
    : BoundaryFaceMapping(1.0) {}

}  // namespace Geometry

}  // namespace hpgem
