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

#ifndef HPGEM_KERNEL_MAPPINGTOPHYSPYRAMID_H
#define HPGEM_KERNEL_MAPPINGTOPHYSPYRAMID_H

#include "MappingReferenceToPhysical.h"
#include <vector>

namespace hpgem {

namespace Geometry {
/*!
 * "In geometry, a pyramid is a polyhedron formed by connecting a polygonal base
 * and a point, called the apex. Each base edge and apex form a triangle."
 * -Wikipedia.
 *
 * This class defines the mappings between reference and physical pyramids.
 */

class MappingToPhysPyramid : public MappingReferenceToPhysicalDim<3> {
   public:
    MappingToPhysPyramid(const PhysicalGeometry<3>* const physicalGeometry);

    MappingToPhysPyramid(const MappingToPhysPyramid& other) = default;

    PointPhysical<3> transform(const PointReference<3>&) const final;

    PointReference<3> inverseTransform(const PointPhysical<3>&) const final;

    Jacobian<3, 3> calcJacobian(const PointReference<3>&) const final;

    void reinit() final;

    MappingToPhysPyramid* copy() const final {
        return new MappingToPhysPyramid(*this);
    }

    bool isValidPoint(const PointReference<3>&) const;
};
}  // namespace Geometry
}  // namespace hpgem

#endif  // HPGEM_KERNEL_MAPPINGTOPHYSPYRAMID_H
