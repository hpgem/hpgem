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

#ifndef HPGEM_KERNEL_REFERENCELINE_H
#define HPGEM_KERNEL_REFERENCELINE_H

#include "ReferenceGeometry.h"

#include <iostream>
#include <vector>

namespace Geometry {
/* Behold the reference line:
 *
 * (-1) 0---1---1 (+1)
 *
 */
class ReferenceLine : public ReferenceGeometry {
   public:
    static ReferenceLine& Instance() {
        static ReferenceLine theInstance;
        return theInstance;
    }

    ReferenceLine(const ReferenceLine&) = delete;

    //! (see ReferenceGeometry.h)
    bool isInternalPoint(const PointReference<1>&) const final;

    //! Output routine.
    friend std::ostream& operator<<(std::ostream& os,
                                    const ReferenceLine& point);

    const PointReferenceBase& getCenter() const final { return center_; }

    std::size_t getNumberOfNodes() const final { return 2; }

    const PointReferenceBase& getReferenceNodeCoordinate(
        const std::size_t& i) const final {
        logger.assert_debug(i < getNumberOfNodes(),
                            "Asked for node %, but there are only % nodes", i,
                            getNumberOfNodes());
        return points_[i];
    }

    // ================================== Codimension 0
    // ========================================

    //! (see MappingCodimensions.h)
    std::size_t getCodim0MappingIndex(
        const std::vector<std::size_t>&,
        const std::vector<std::size_t>&) const final;

    //! (see MappingCodimensions.h)
    const MappingReferenceToReference<0>* getCodim0MappingPtr(
        const std::size_t) const final;

    using MappingCodimensions::getCodim0MappingPtr;

    // ================================== Codimension 1
    // ========================================

    //! (see MappingCodimensions.h)
    std::size_t getNumberOfCodim1Entities() const final { return 2; }

    //! (see MappingCodimensions.h)
    std::vector<std::size_t> getCodim1EntityLocalIndices(
        const std::size_t) const final;

    //! (see MappingCodimensions.h)
    const MappingReferenceToReference<1>* getCodim1MappingPtr(
        const std::size_t) const final;

    //! (see MappingCodimensions.h)
    const ReferenceGeometry* getCodim1ReferenceGeometry(
        const std::size_t) const final;

   private:
    ReferenceLine();

    //! Local node indexes contains the numbering of the vertex of the shape,
    //! ordered by faces. See top comment for the corresponding numbering. (In a
    //! line, the 'faces' are nodes, and nodes are just a coordinate).
    static std::size_t localNodeIndexes_[2][1];

    //! Codimension 0 mappings, from a line to a line. Used to rotate the face
    //! when the left and right elements dont think it has the same orientation
    const MappingReferenceToReference<0>* mappingsLineToLine_[2];

    //! Codimension 1 mappings, from a point to a line. This is the 1D
    //! face->element map
    const MappingReferenceToReference<1>* mappingsPointToLine_[2];

    //! Pointer to the Codimension 1 reference geometry, in this case, to
    //! ReferencePoint.
    ReferenceGeometry* const referenceGeometryCodim1Ptr_;

    //! List of valid quadrature rules for this reference geometry
    std::vector<QuadratureRules::GaussQuadratureRule*> lstGaussQuadratureRules_;

    std::vector<PointReference<1> > points_;

    PointReference<1> center_;
};
}  // namespace Geometry
#endif  // HPGEM_KERNEL_REFERENCELINE_H
