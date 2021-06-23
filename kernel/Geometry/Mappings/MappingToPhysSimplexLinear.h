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

#ifndef HPGEM_KERNEL_MAPPINGTOPHYSSIMPLEXLINEAR_H
#define HPGEM_KERNEL_MAPPINGTOPHYSSIMPLEXLINEAR_H

#include "MappingReferenceToPhysical.h"

// maybe it is a better idea to specialize this template in a source file even
// though it can be coded up generally in a header
#include "Geometry/PointPhysical.h"
#include "Geometry/PhysicalGeometry.h"
namespace hpgem {
namespace Geometry {
/*!
 * "In geometry, a simplex (plural simplexes or simplices) is a generalization
 * of the notion of a triangle or tetrahedron to arbitrary dimension."
 * -Wikipedia.
 *
 * This class defines the linear mappings between simplexes, namely from the
 * reference domain to the physical domain.
 * See the comments in the Physical<Simplex>.cpp files to know the order of the
 * vertex of each simplex, an order which is kept by the mappings.
 * No specialization is needed because the mapping is general for geometries of
 * vertex number (vn) one greater than dimension (d), that is, vn = d+1; Do not
 * use this mapping for DIM=1, use MappingToPhysHypercubeLinear<1> instead
 */

template <std::size_t DIM>
class MappingToPhysSimplexLinear : public MappingReferenceToPhysical<DIM> {
   public:
    MappingToPhysSimplexLinear(const PhysicalGeometry<DIM>* const& pG)
        : MappingReferenceToPhysical<DIM>(pG) {
        reinit();
    }

    MappingToPhysSimplexLinear(const MappingToPhysSimplexLinear<DIM>& other) =
        default;

    PointPhysical<DIM> transform(const PointReference<DIM>&) const final;
    PointReference<DIM> inverseTransform(const PointPhysical<DIM>&) const final;
    Jacobian<DIM, DIM> calcJacobian(const PointReference<DIM>&) const final;
    void reinit() final;

    MappingToPhysSimplexLinear<DIM>* copy() const final {
        return new MappingToPhysSimplexLinear<DIM>(*this);
    }

   private:
    ///\todo: Implement this function.
    bool isValidPoint(const PointReference<DIM>&) const;
};
}  // namespace Geometry
}  // namespace hpgem
#include "MappingToPhysSimplexLinear_Impl.h"

#endif  // HPGEM_KERNEL_MAPPINGTOPHYSSIMPLEXLINEAR_H
