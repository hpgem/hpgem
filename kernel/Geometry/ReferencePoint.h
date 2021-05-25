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

#ifndef HPGEM_KERNEL_REFERENCEPOINT_H
#define HPGEM_KERNEL_REFERENCEPOINT_H

#include <iostream>

#include "ReferenceGeometry.h"
#include <vector>
#include "Logger.h"

namespace hpgem {

namespace Geometry {
/// \class ReferencePoint
/// \brief Reference geometry of dimensions 0.
/// \details
/// Most of the codimension functions are not well defined for dimension 0, so
/// be careful.
class ReferencePoint : public ReferenceGeometry {
   public:
    static ReferencePoint& Instance() {
        static ReferencePoint theInstance;
        return theInstance;
    }

    ReferencePoint(const ReferencePoint&) = delete;

    /// \brief Return true.
    bool isInternalPoint(const PointReference<0>& p) const final;

    // ================================== Codimension 0
    // ========================================

    /// \brief Returns 0.
    std::size_t getCodim0MappingIndex(
        const std::vector<std::size_t>&,
        const std::vector<std::size_t>&) const final;

    /// \brief Returns 0.
    const MappingReferenceToReference<0>* getCodim0MappingPtr(
        const std::size_t a) const final;

    const PointReferenceBase& getCenter() const final { return center_; }

    std::size_t getNumberOfNodes() const final { return 1; }

    const PointReferenceBase& getReferenceNodeCoordinate(
        const std::size_t& i) const final {
        logger.assert_debug(i < getNumberOfNodes(),
                            "Asked for node %, but there are only % nodes", i,
                            getNumberOfNodes());
        return points_[i];
    }


   private:
    /// \brief List of valid quadrature rules for this reference geometry
    std::vector<QuadratureRules::GaussQuadratureRule*> lstGaussQuadratureRules_;

    //! Codimension 0 mappings, from a line to a line. Used to rotate the face
    //! when the left and right elements dont think it has the same orientation
    /// Provided for consistency with other dimensions
    const MappingReferenceToReference<0>* mappingsPointToPoint_;

    std::vector<PointReference<0> > points_;

    PointReference<0> center_;

    ReferencePoint();
};
}  // namespace Geometry
}  // namespace hpgem

#endif  // HPGEM_KERNEL_REFERENCEPOINT_H
