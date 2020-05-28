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

#ifndef KERNEL_GEOMETRY_MAPPINGS_REFINEMENTMAPSFORLINE_H_
#define KERNEL_GEOMETRY_MAPPINGS_REFINEMENTMAPSFORLINE_H_

#include "RefinementMapping.h"
#include "RefinementMapsForPoint.h"
#include "Geometry/PointReference.h"
#include "Geometry/ReferenceLine.h"

namespace Geometry {
/**
 * Stores all the refinement mappings for a line.
 *
 * o-----0-----o          o--0--o--1--o
 * 0           1          0     2     1
 *   (identity)       (split in the middle)
 *      new               legacy index 0
 *    index 0                 index 1
 *
 * The legacy implementation used to have an additional row and column in its
 * transformation matrix The final row contains a 1 on the diagonal and 0's
 * elsewhere, the final column is mapping-depended and will be provided in
 * comments with the appropriate function. This column seems to be the constant
 * component of the transformation
 */

class RefinementMapForLine0 : public RefinementMapping {
   public:
    static const RefinementMapForLine0* instance() {
        static RefinementMapForLine0 theInstance;
        return &theInstance;
    }

    std::string getName() const final { return "Identity map"; }

    PointReference<1> refinementTransform(
        std::size_t subElementIndex, const PointReference<1>& p) const final {
        logger.assert_debug(
            subElementIndex == 0,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        return p;
    }

    LinearAlgebra::SmallMatrix<1, 1> getRefinementMappingMatrixL(
        std::size_t subElementIndex, const PointReference<1>& p) const final {
        logger.assert_debug(
            subElementIndex == 0,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        // the bonus column is most likely {0, 1}
        return LinearAlgebra::SmallMatrix<1, 1>{1.};
    }

    LinearAlgebra::SmallMatrix<1, 1> getRefinementMappingMatrixR(
        std::size_t subElementIndex, const PointReference<1>& p) const final {
        logger.assert_debug(
            subElementIndex == 0,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        // the bonus column is most likely {0, 1}
        return LinearAlgebra::SmallMatrix<1, 1>{1.};
    }

    std::size_t getNumberOfNewNodes() const final { return 0; }

    std::size_t getNumberOfSubElements() const final { return 1; }

    std::vector<std::size_t> getSubElementLocalNodeIndices(
        std::size_t subElementIndex) const final {
        logger.assert_debug(
            subElementIndex == 0,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        return std::vector<std::size_t>{0, 1};
    }

    Geometry::ReferenceGeometry* getBigElementReferenceGeometry() const final {
        return &Geometry::ReferenceLine::Instance();
    }

    Geometry::ReferenceGeometry* getSubElementReferenceGeometry(
        std::size_t subElement) const final {
        logger.assert_debug(
            subElement < getNumberOfSubElements(),
            "asked for subElement %, but the % has only % subElements",
            subElement, getName(), getNumberOfSubElements());
        return &Geometry::ReferenceLine::Instance();
    }

    std::vector<PointReference<1>> getNewNodeLocations(
        const PointReference<1>&) const final {
        return {};
    }

    std::vector<const RefinementMapping*> getCodim1RefinementMaps()
        const final {
        return std::vector<const RefinementMapping*>{
            RefinementMapForPoint0::instance(),
            RefinementMapForPoint0::instance()};
    }

    std::vector<std::size_t> getCodim1LocalNodeIndices(
        std::size_t localFaceNumber) const final {
        switch (localFaceNumber) {
            case 0:
                return {0};
            case 1:
                return {1};
            default:
                logger(ERROR, "asked for Face %, but there are only % Faces",
                       localFaceNumber,
                       getBigElementReferenceGeometry()
                           ->getNumberOfCodim1Entities());
                return {};
        }
    }

    std::tuple<std::size_t, std::size_t> getSubElementAndLocalFaceIndex(
        std::size_t face, std::size_t subFaceIndex) const final {
        logger.assert_debug(
            face <
                getBigElementReferenceGeometry()->getNumberOfCodim1Entities(),
            "asked for Face %, but the % has only % faces", face, getName(),
            getBigElementReferenceGeometry()->getNumberOfCodim1Entities());
        logger.assert_debug(
            subFaceIndex <
                getCodim1RefinementMaps()[face]->getNumberOfSubElements(),
            "asked for subFace %, but the % has only % subFaces", subFaceIndex,
            getCodim1RefinementMaps()[face]->getName(),
            getCodim1RefinementMaps()[face]->getNumberOfSubElements());
        return std::make_tuple(0, face);
    }

   private:
    RefinementMapForLine0() = default;
    RefinementMapForLine0(const RefinementMapForLine0&) = delete;
};

class RefinementMapForLine1 : public RefinementMapping {
   public:
    static const RefinementMapForLine1* instance() {
        static RefinementMapForLine1 theInstance;
        return &theInstance;
    }

    std::string getName() const final { return "split"; }

    PointReference<1> refinementTransform(
        std::size_t subElementIndex, const PointReference<1>& p) const final {
        logger.assert_debug(
            subElementIndex < 2,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        return {p[0] / 2. + subElementIndex - 0.5};
    }

    LinearAlgebra::SmallMatrix<1, 1> getRefinementMappingMatrixL(
        std::size_t subElementIndex, const PointReference<1>& p) const final {
        logger.assert_debug(
            subElementIndex < 2,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        // the bonus column is probably {-0.5 + subElementIndex, 1}
        return LinearAlgebra::SmallMatrix<1, 1>{0.5};
    }

    LinearAlgebra::SmallMatrix<1, 1> getRefinementMappingMatrixR(
        std::size_t subElementIndex, const PointReference<1>& p) const final {
        logger.assert_debug(
            subElementIndex < 2,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        // the bonus column is probably {1 - 2 * subElementIndex, 1}
        return LinearAlgebra::SmallMatrix<1, 1>{2.};
    }

    std::size_t getNumberOfNewNodes() const final { return 1; }

    std::size_t getNumberOfSubElements() const final { return 2; }

    std::vector<std::size_t> getSubElementLocalNodeIndices(
        std::size_t subElementIndex) const final {
        logger.assert_debug(
            subElementIndex < 2,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        if (subElementIndex == 0) {
            return std::vector<std::size_t>{0, 2};
        }
        return std::vector<std::size_t>{2, 1};
    }

    Geometry::ReferenceGeometry* getBigElementReferenceGeometry() const final {
        return &Geometry::ReferenceLine::Instance();
    }

    Geometry::ReferenceGeometry* getSubElementReferenceGeometry(
        std::size_t subElement) const final {
        logger.assert_debug(
            subElement < getNumberOfSubElements(),
            "asked for subElement %, but the % has only % subElements",
            subElement, getName(), getNumberOfSubElements());
        return &Geometry::ReferenceLine::Instance();
    }

    std::vector<PointReference<1>> getNewNodeLocations(
        const PointReference<1>&) const final {
        return {{0.}};
    }

    std::vector<const RefinementMapping*> getCodim1RefinementMaps()
        const final {
        return std::vector<const RefinementMapping*>{
            RefinementMapForPoint0::instance(),
            RefinementMapForPoint0::instance()};
    }

    std::vector<std::size_t> getCodim1LocalNodeIndices(
        std::size_t localFaceNumber) const final {
        switch (localFaceNumber) {
            case 0:
                return {0};
            case 1:
                return {1};
            default:
                logger(ERROR, "asked for Face %, but there are only % Faces",
                       localFaceNumber,
                       getBigElementReferenceGeometry()
                           ->getNumberOfCodim1Entities());
                return {};
        }
    }

    std::tuple<std::size_t, std::size_t> getSubElementAndLocalFaceIndex(
        std::size_t face, std::size_t subFaceIndex) const final {
        logger.assert_debug(
            face <
                getBigElementReferenceGeometry()->getNumberOfCodim1Entities(),
            "asked for Face %, but the % has only % faces", face, getName(),
            getBigElementReferenceGeometry()->getNumberOfCodim1Entities());
        logger.assert_debug(
            subFaceIndex <
                getCodim1RefinementMaps()[face]->getNumberOfSubElements(),
            "asked for subFace %, but the % has only % subFaces", subFaceIndex,
            getCodim1RefinementMaps()[face]->getName(),
            getCodim1RefinementMaps()[face]->getNumberOfSubElements());
        return std::make_tuple(face, face);
    }

   private:
    RefinementMapForLine1() = default;
    RefinementMapForLine1(const RefinementMapForLine1&) = delete;
};
}  // namespace Geometry

#endif /* KERNEL_GEOMETRY_MAPPINGS_REFINEMENTMAPSFORLINE_H_ */
