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

#ifndef HPGEM_KERNEL_REFINEMENTMAPSFORPOINT_H
#define HPGEM_KERNEL_REFINEMENTMAPSFORPOINT_H

#include "RefinementMapping.h"
#include "Geometry/PointReference.h"
#include "Geometry/ReferencePoint.h"

namespace Geometry {
/**
 * Stores all the refinement mappings for a point.
 * This will likely only be the identity map.
 * The legacy implementation used to have an additional row and column in its
 * transformation matrix The final row contains a 1 on the diagonal and 0's
 * elsewhere, the final column is mapping-depended and will be provided in
 * comments with the appropriate function. This column seems to be the constant
 * component of the transformation
 */

class RefinementMapForPoint0 : public RefinementMapping {
   public:
    static const RefinementMapForPoint0* instance() {
        static RefinementMapForPoint0 theInstance;
        return &theInstance;
    }

    std::string getName() const final { return "Identity map"; }

    PointReference<0> refinementTransform(
        std::size_t subElementIndex, const PointReference<0>& p) const final {
        logger.assert_debug(
            subElementIndex == 0,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        return p;
    }

    LinearAlgebra::SmallMatrix<0, 0> getRefinementMappingMatrixL(
        std::size_t subElementIndex, const PointReference<0>& p) const final {
        logger.assert_debug(
            subElementIndex == 0,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        // the bonus column is most likely {1}
        return LinearAlgebra::SmallMatrix<0, 0>();
    }

    LinearAlgebra::SmallMatrix<0, 0> getRefinementMappingMatrixR(
        std::size_t subElementIndex, const PointReference<0>& p) const final {
        logger.assert_debug(
            subElementIndex == 0,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        // the bonus column is most likely {1}
        return LinearAlgebra::SmallMatrix<0, 0>();
    }

    std::size_t getNumberOfNewNodes() const final { return 0; }

    std::size_t getNumberOfSubElements() const final { return 1; }

    Geometry::ReferenceGeometry* getBigElementReferenceGeometry() const final {
        return &Geometry::ReferencePoint::Instance();
    }

    Geometry::ReferenceGeometry* getSubElementReferenceGeometry(
        std::size_t subElement) const final {
        logger.assert_debug(
            subElement < getNumberOfSubElements(),
            "asked for subElement %, but the % has only % subElements",
            subElement, getName(), getNumberOfSubElements());
        return &Geometry::ReferencePoint::Instance();
    }

    std::vector<PointReference<0>> getNewNodeLocations(
        const PointReference<0>&) const final {
        return {};
    }

    std::vector<std::size_t> getSubElementLocalNodeIndices(
        std::size_t subElementIndex) const final {
        logger.assert_debug(
            subElementIndex == 0,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        return std::vector<std::size_t>(1, 0);
    }

    std::vector<const RefinementMapping*> getCodim1RefinementMaps()
        const final {
        return std::vector<const RefinementMapping*>();
    }

    std::vector<std::size_t> getCodim1LocalNodeIndices(
        std::size_t localFaceNumber) const final {
        logger(ERROR, "Asked for face number %, but there are only 0 faces",
               localFaceNumber);
        return {};
    }

    std::tuple<std::size_t, std::size_t> getSubElementAndLocalFaceIndex(
        std::size_t face, std::size_t subFaceIndex) const final {
        logger(ERROR, "Asked for face number %, but there are only 0 faces",
               face);
        return {};
    }

   private:
    RefinementMapForPoint0() = default;
    RefinementMapForPoint0(const RefinementMapForPoint0&) = delete;
};
}  // namespace Geometry

#endif // HPGEM_KERNEL_REFINEMENTMAPSFORPOINT_H
