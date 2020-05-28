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

#ifndef KERNEL_GEOMETRY_MAPPINGS_REFINEMENTMAPSFORSQUARE_H_
#define KERNEL_GEOMETRY_MAPPINGS_REFINEMENTMAPSFORSQUARE_H_

#include "RefinementMapping.h"
#include "RefinementMapsForLine.h"
#include "Geometry/PointReference.h"
#include "Geometry/ReferenceSquare.h"

namespace Geometry {
/**
 * Stores all the refinement mappings for a square.
 *
 *                                  5
 * 2 o-------3-------o 3  2 o---3---o---5---o 3
 *   |               |      |       |       |
 *   |               |      |       |       |
 *   |               |      |       |       |
 *   1       0       2      1   0   6   1   2
 *   |               |      |       |       |
 *   |               |      |       |       |
 *   |               |      |       |       |
 * 0 o-------0-------o 1  0 o---0---o---4---o 1
 *                                  4
 *    (identity)            (vertical split)
 *       new                 legacy index 0
 *     index 0                  index 1
 *
 *                                  7
 * 2 o-------3-------o 3  2 o---4---o---5---o 3
 *   |               |      |       |       |
 *   4       1       5      9   2   a   3   b
 *   |               |      |       |       |
 * 4 o-------6-------o 5  5 o---2---o---3---o 6
 *   |               |      |       |       |
 *   1       0       2      6   0   7   1   8
 *   |               |      |       |       |
 * 0 o-------0-------o 1  0 o---0---o---1---o 1
 *                                  4
 *  (horizontal split)        (4-way split)
 *    legacy index 1          legacy index 2
 *       index 2                  index 3
 *                         centroid is node 8
 *
 * The legacy implementation used to have an additional row and column in its
 * transformation matrix The final row contains a 1 on the diagonal and 0's
 * elsewhere, the final column is mapping-depended and will be provided in
 * comments with the appropriate function. This column seems to be the constant
 * component of the transformation
 */

class RefinementMapForSquare0 : public RefinementMapping {
   public:
    static const RefinementMapForSquare0* instance() {
        static RefinementMapForSquare0 theInstance;
        return &theInstance;
    }

    std::string getName() const override final { return "Identity map"; }

    PointReference<2> refinementTransform(
        std::size_t subElementIndex,
        const PointReference<2>& p) const override final {
        logger.assert_debug(
            subElementIndex == 0,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        return p;
    }

    LinearAlgebra::SmallMatrix<2, 2> getRefinementMappingMatrixL(
        std::size_t subElementIndex,
        const PointReference<2>& p) const override final {
        logger.assert_debug(
            subElementIndex == 0,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        // the bonus column is most likely {0, 0, 1}
        return LinearAlgebra::SmallMatrix<2, 2>{{{{1., 0.}}, {{0., 1.}}}};
    }

    LinearAlgebra::SmallMatrix<2, 2> getRefinementMappingMatrixR(
        std::size_t subElementIndex,
        const PointReference<2>& p) const override final {
        logger.assert_debug(
            subElementIndex == 0,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        // the bonus column is most likely {0, 0, 1}
        return LinearAlgebra::SmallMatrix<2, 2>{{{{1., 0.}}, {{0., 1.}}}};
    }

    std::size_t getNumberOfNewNodes() const override final { return 0; }

    std::size_t getNumberOfSubElements() const override final { return 1; }

    std::vector<std::size_t> getSubElementLocalNodeIndices(
        std::size_t subElementIndex) const override final {
        logger.assert_debug(
            subElementIndex == 0,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        return std::vector<std::size_t>{0, 1, 2, 3};
    }

    Geometry::ReferenceGeometry* getBigElementReferenceGeometry()
        const override final {
        return &Geometry::ReferenceSquare::Instance();
    }

    Geometry::ReferenceGeometry* getSubElementReferenceGeometry(
        std::size_t subElement) const override final {
        logger.assert_debug(
            subElement < getNumberOfSubElements(),
            "asked for subElement %, but the % has only % subElements",
            subElement, getName(), getNumberOfSubElements());
        return &Geometry::ReferenceSquare::Instance();
    }

    std::vector<Geometry::PointReference<2>> getNewNodeLocations(
        const PointReference<2>&) const override final {
        return {};
    }

    std::vector<const RefinementMapping*> getCodim1RefinementMaps()
        const override final {
        return std::vector<const RefinementMapping*>{
            RefinementMapForLine0::instance(),
            RefinementMapForLine0::instance(),
            RefinementMapForLine0::instance(),
            RefinementMapForLine0::instance()};
    }

    std::vector<std::size_t> getCodim1LocalNodeIndices(
        std::size_t localFaceNumber) const override final {
        switch (localFaceNumber) {
            case 0:
                return {0, 1};
            case 1:
                return {0, 2};
            case 2:
                return {1, 3};
            case 3:
                return {2, 3};
            default:
                logger(ERROR, "asked for Face %, but there are only % Faces",
                       localFaceNumber,
                       getBigElementReferenceGeometry()
                           ->getNumberOfCodim1Entities());
                return {};
        }
    }

    std::tuple<std::size_t, std::size_t> getSubElementAndLocalFaceIndex(
        std::size_t face, std::size_t subFaceIndex) const override final {
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
    RefinementMapForSquare0() = default;
    RefinementMapForSquare0(const RefinementMapForSquare0&) = delete;
};

class RefinementMapForSquare1 : public RefinementMapping {
   public:
    static const RefinementMapForSquare1* instance() {
        static RefinementMapForSquare1 theInstance;
        return &theInstance;
    }

    std::string getName() const override final { return "vertical split"; }

    PointReference<2> refinementTransform(
        std::size_t subElementIndex,
        const PointReference<2>& p) const override final {
        logger.assert_debug(
            subElementIndex < 2,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        return {p[0] / 2. + subElementIndex - 0.5, p[1]};
    }

    LinearAlgebra::SmallMatrix<2, 2> getRefinementMappingMatrixL(
        std::size_t subElementIndex,
        const PointReference<2>& p) const override final {
        logger.assert_debug(
            subElementIndex < 2,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        // the bonus column is {-0.5 + subElementIndex, 0, 1}
        return LinearAlgebra::SmallMatrix<2, 2>{{{{0.5, 0.}}, {{0., 1.}}}};
    }

    LinearAlgebra::SmallMatrix<2, 2> getRefinementMappingMatrixR(
        std::size_t subElementIndex,
        const PointReference<2>& p) const override final {
        logger.assert_debug(
            subElementIndex < 2,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        // the bonus column is {1 - 2 * subelementIndex, 0, 1}
        return LinearAlgebra::SmallMatrix<2, 2>{{{{2., 0.}}, {{0., 1.}}}};
    }

    std::size_t getNumberOfNewNodes() const override final { return 2; }

    std::size_t getNumberOfSubElements() const override final { return 2; }

    std::vector<std::size_t> getSubElementLocalNodeIndices(
        std::size_t subElementIndex) const override final {
        logger.assert_debug(
            subElementIndex < 2,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        if (subElementIndex == 0) {
            return std::vector<std::size_t>{0, 4, 2, 5};
        } 
            return std::vector<std::size_t>{4, 1, 5, 3};
        
    }

    Geometry::ReferenceGeometry* getBigElementReferenceGeometry()
        const override final {
        return &Geometry::ReferenceSquare::Instance();
    }

    Geometry::ReferenceGeometry* getSubElementReferenceGeometry(
        std::size_t subElement) const override final {
        logger.assert_debug(
            subElement < getNumberOfSubElements(),
            "asked for subElement %, but the % has only % subElements",
            subElement, getName(), getNumberOfSubElements());
        return &Geometry::ReferenceSquare::Instance();
    }

    std::vector<Geometry::PointReference<2>> getNewNodeLocations(
        const PointReference<2>&) const override final {
        return {{0., -1.}, {0., 1.}};
    }

    std::vector<const RefinementMapping*> getCodim1RefinementMaps()
        const override final {
        return std::vector<const RefinementMapping*>{
            RefinementMapForLine1::instance(),
            RefinementMapForLine0::instance(),
            RefinementMapForLine0::instance(),
            RefinementMapForLine1::instance()};
    }

    std::vector<std::size_t> getCodim1LocalNodeIndices(
        std::size_t localFaceNumber) const override final {
        switch (localFaceNumber) {
            case 0:
                return {0, 1, 4};
            case 1:
                return {0, 2};
            case 2:
                return {1, 3};
            case 3:
                return {2, 3, 5};
            default:
                logger(ERROR, "asked for Face %, but there are only % Faces",
                       localFaceNumber,
                       getBigElementReferenceGeometry()
                           ->getNumberOfCodim1Entities());
                return {};
        }
    }

    std::tuple<std::size_t, std::size_t> getSubElementAndLocalFaceIndex(
        std::size_t face, std::size_t subFaceIndex) const override final {
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
        switch (face) {
            case 0:
            case 3:
                return std::make_tuple(subFaceIndex, face);
            case 1:
                return std::make_tuple(0, face);
            case 2:
                return std::make_tuple(1, face);
            default:
                logger(ERROR, "something broke, please check the assertions");
        }
        return {};
    }

   private:
    RefinementMapForSquare1() = default;
    RefinementMapForSquare1(const RefinementMapForSquare1&) = delete;
};

class RefinementMapForSquare2 : public RefinementMapping {
   public:
    static const RefinementMapForSquare2* instance() {
        static RefinementMapForSquare2 theInstance;
        return &theInstance;
    }

    std::string getName() const override final { return "horizontal split"; }

    PointReference<2> refinementTransform(
        std::size_t subElementIndex,
        const PointReference<2>& p) const override final {
        logger.assert_debug(
            subElementIndex < 2,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        return {p[0], p[1] / 2. + subElementIndex - 0.5};
    }

    LinearAlgebra::SmallMatrix<2, 2> getRefinementMappingMatrixL(
        std::size_t subElementIndex,
        const PointReference<2>& p) const override final {
        logger.assert_debug(
            subElementIndex < 2,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        // the bonus column is {0, -0.5 + subElementIndex, 1}
        return LinearAlgebra::SmallMatrix<2, 2>{{{{1., 0.}}, {{0., 0.5}}}};
    }

    LinearAlgebra::SmallMatrix<2, 2> getRefinementMappingMatrixR(
        std::size_t subElementIndex,
        const PointReference<2>& p) const override final {
        logger.assert_debug(
            subElementIndex < 2,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        // the bonus column is {0, 1 - 2 * subelementIndex, 1}
        return LinearAlgebra::SmallMatrix<2, 2>{{{{1., 0.}}, {{0., 2.}}}};
    }

    std::size_t getNumberOfNewNodes() const override final { return 2; }

    std::size_t getNumberOfSubElements() const override final { return 2; }

    std::vector<std::size_t> getSubElementLocalNodeIndices(
        std::size_t subElementIndex) const override final {
        logger.assert_debug(
            subElementIndex < 2,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        if (subElementIndex == 0) {
            return std::vector<std::size_t>{0, 1, 4, 5};
        } 
            return std::vector<std::size_t>{4, 5, 2, 3};
        
    }

    Geometry::ReferenceGeometry* getBigElementReferenceGeometry()
        const override final {
        return &Geometry::ReferenceSquare::Instance();
    }

    Geometry::ReferenceGeometry* getSubElementReferenceGeometry(
        std::size_t subElement) const override final {
        logger.assert_debug(
            subElement < getNumberOfSubElements(),
            "asked for subElement %, but the % has only % subElements",
            subElement, getName(), getNumberOfSubElements());
        return &Geometry::ReferenceSquare::Instance();
    }

    std::vector<Geometry::PointReference<2>> getNewNodeLocations(
        const PointReference<2>&) const override final {
        return {{-1., 0.}, {1., 0.}};
    }

    std::vector<const RefinementMapping*> getCodim1RefinementMaps()
        const override final {
        return std::vector<const RefinementMapping*>{
            RefinementMapForLine0::instance(),
            RefinementMapForLine1::instance(),
            RefinementMapForLine1::instance(),
            RefinementMapForLine0::instance()};
    }

    std::vector<std::size_t> getCodim1LocalNodeIndices(
        std::size_t localFaceNumber) const override final {
        switch (localFaceNumber) {
            case 0:
                return {0, 1};
            case 1:
                return {0, 2, 4};
            case 2:
                return {1, 3, 5};
            case 3:
                return {2, 3};
            default:
                logger(ERROR, "asked for Face %, but there are only % Faces",
                       localFaceNumber,
                       getBigElementReferenceGeometry()
                           ->getNumberOfCodim1Entities());
                return {};
        }
    }

    std::tuple<std::size_t, std::size_t> getSubElementAndLocalFaceIndex(
        std::size_t face, std::size_t subFaceIndex) const override final {
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
        switch (face) {
            case 1:
            case 2:
                return std::make_tuple(subFaceIndex, face);
            case 0:
                return std::make_tuple(0, face);
            case 3:
                return std::make_tuple(1, face);
            default:
                logger(ERROR, "something broke, please check the assertions");
        }
        return {};
    }

   private:
    RefinementMapForSquare2() = default;
    RefinementMapForSquare2(const RefinementMapForSquare2&) = delete;
};

class RefinementMapForSquare3 : public RefinementMapping {
   public:
    static const RefinementMapForSquare3* instance() {
        static RefinementMapForSquare3 theInstance;
        return &theInstance;
    }

    std::string getName() const override final { return "split all faces"; }

    PointReference<2> refinementTransform(
        std::size_t subElementIndex,
        const PointReference<2>& p) const override final {
        switch (subElementIndex) {
            case 0:
                return {-0.5 + 0.5 * p[0], -0.5 + 0.5 * p[1]};
            case 1:
                return {0.5 + 0.5 * p[0], -0.5 + 0.5 * p[1]};
            case 2:
                return {-0.5 + 0.5 * p[0], 0.5 + 0.5 * p[1]};
            case 3:
                return {0.5 + 0.5 * p[0], 0.5 + 0.5 * p[1]};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return p;
        }
    }

    LinearAlgebra::SmallMatrix<2, 2> getRefinementMappingMatrixL(
        std::size_t subElementIndex,
        const PointReference<2>& p) const override final {
        switch (subElementIndex) {
            case 0:
                // the bonus column is {-0.5, -0.5, 1}
                return LinearAlgebra::SmallMatrix<2, 2>{
                    {{{0.5, 0.}}, {{0., 0.5}}}};
            case 1:
                // the bonus column is {0.5, -0.5, 1}
                return LinearAlgebra::SmallMatrix<2, 2>{
                    {{{0.5, 0.}}, {{0., 0.5}}}};
            case 2:
                // the bonus column is {-0.5, 0.5, 1}
                return LinearAlgebra::SmallMatrix<2, 2>{
                    {{{0.5, 0.}}, {{0., 0.5}}}};
            case 3:
                // the bonus column is {0.5, 0.5, 1}
                return LinearAlgebra::SmallMatrix<2, 2>{
                    {{{0.5, 0.}}, {{0., 0.5}}}};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return LinearAlgebra::SmallMatrix<2, 2>();
        }
    }

    LinearAlgebra::SmallMatrix<2, 2> getRefinementMappingMatrixR(
        std::size_t subElementIndex,
        const PointReference<2>& p) const override final {
        switch (subElementIndex) {
            case 0:
                // the bonus column is {1, 1, 1}
                return LinearAlgebra::SmallMatrix<2, 2>{
                    {{{2., 0.}}, {{0., 2.0}}}};
            case 1:
                // the bonus column is {-1, 1, 1}
                return LinearAlgebra::SmallMatrix<2, 2>{
                    {{{2., 0.}}, {{0., 2.0}}}};
            case 2:
                // the bonus column is {1, -1, 1}
                return LinearAlgebra::SmallMatrix<2, 2>{
                    {{{2., 0.}}, {{0., 2.0}}}};
            case 3:
                // the bonus column is {-1, -1, 1}
                return LinearAlgebra::SmallMatrix<2, 2>{
                    {{{2., 0.}}, {{0., 2.0}}}};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return LinearAlgebra::MiddleSizeMatrix();
        }
    }

    std::size_t getNumberOfNewNodes() const override final { return 5; }

    std::size_t getNumberOfSubElements() const override final { return 4; }

    std::vector<std::size_t> getSubElementLocalNodeIndices(
        std::size_t subElementIndex) const override final {
        switch (subElementIndex) {
            case 0:
                return std::vector<std::size_t>{0, 4, 5, 8};
            case 1:
                return std::vector<std::size_t>{4, 1, 8, 6};
            case 2:
                return std::vector<std::size_t>{5, 8, 2, 7};
            case 3:
                return std::vector<std::size_t>{8, 6, 7, 3};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return std::vector<std::size_t>();
        }
    }

    Geometry::ReferenceGeometry* getBigElementReferenceGeometry()
        const override final {
        return &Geometry::ReferenceSquare::Instance();
    }

    Geometry::ReferenceGeometry* getSubElementReferenceGeometry(
        std::size_t subElement) const override final {
        logger.assert_debug(
            subElement < getNumberOfSubElements(),
            "asked for subElement %, but the % has only % subElements",
            subElement, getName(), getNumberOfSubElements());
        return &Geometry::ReferenceSquare::Instance();
    }

    std::vector<Geometry::PointReference<2>> getNewNodeLocations(
        const PointReference<2>&) const override final {
        return {{0., -1.}, {-1., 0.}, {1., 0}, {0., 1.}, {0., 0.}};
    }

    std::vector<const RefinementMapping*> getCodim1RefinementMaps()
        const override final {
        return std::vector<const RefinementMapping*>{
            RefinementMapForLine1::instance(),
            RefinementMapForLine1::instance(),
            RefinementMapForLine1::instance(),
            RefinementMapForLine1::instance()};
    }

    std::vector<std::size_t> getCodim1LocalNodeIndices(
        std::size_t localFaceNumber) const override final {
        switch (localFaceNumber) {
            case 0:
                return {0, 1, 4};
            case 1:
                return {0, 2, 5};
            case 2:
                return {1, 3, 6};
            case 3:
                return {2, 3, 7};
            default:
                logger(ERROR, "asked for Face %, but there are only % Faces",
                       localFaceNumber,
                       getBigElementReferenceGeometry()
                           ->getNumberOfCodim1Entities());
                return {};
        }
    }

    std::tuple<std::size_t, std::size_t> getSubElementAndLocalFaceIndex(
        std::size_t face, std::size_t subFaceIndex) const override final {
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
        switch (face) {
            case 0:
                return std::make_tuple(subFaceIndex, face);
            case 1:
                return std::make_tuple(subFaceIndex * 2, face);
            case 2:
                return std::make_tuple(subFaceIndex * 2 + 1, face);
            case 3:
                return std::make_tuple(subFaceIndex + 2, face);
            default:
                logger(ERROR, "something broke, please check the assertions");
        }
        return {};
    }

   private:
    RefinementMapForSquare3() = default;
    RefinementMapForSquare3(const RefinementMapForSquare3&) = delete;
};
}  // namespace Geometry

#endif /* KERNEL_GEOMETRY_MAPPINGS_REFINEMENTMAPSFORSQUARE_H_ */
