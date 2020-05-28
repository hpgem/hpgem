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

#ifndef KERNEL_GEOMETRY_MAPPINGS_REFINEMENTMAPSFORTRIANGLE_H_
#define KERNEL_GEOMETRY_MAPPINGS_REFINEMENTMAPSFORTRIANGLE_H_

#include "RefinementMapping.h"
#include "RefinementMapsForLine.h"
#include "Geometry/PointReference.h"
#include "Geometry/ReferenceSquare.h"
#include "Geometry/ReferenceTriangle.h"

namespace Geometry {
/**
 * Stores all the refinement mappings for a triangle.
 *
 *  2 o                       2 o
 *    | \                       |\\
 *    |   \                     | \ \
 *    |     \                   |  \  \
 *    1       2                 1   4   2
 *    |         \               |    \    \
 *    |   0       \             |  0  \  1  \
 *    |             \           |      \      \
 *  0 o-------0-------o 1     0 o---0---o---3---o 1
 *                                      3
 *      (identity)           (split in the middle of 0)
 *         new                    legacy index 0
 *       index 0                     index 1
 *
 *  2 o                       2 o
 *    | \                       | \
 *    3   \                     |   2
 *    |  1  \                   |     \  3
 *  3 o _     2                 1  0    o
 *    |   -     \               |     /   \
 *    1     4 _   \             |   4   1   3
 *    |  0      - _ \           | /           \
 *  0 o-------0-------o 1     0 o-------0-------o 1
 *
 *(split in the middle of 1) (split in the middle of 2)
 *     legacy index 1             legacy index 2
 *        index 2                    index 3
 *
 *  2 o                       2 o
 *    | \                       | \
 *    3   4                     3   4
 *    | 2   \  5                |  2  \
 *  4 o---8---o               4 o_      o 5
 *    | \   3 | \               |  8_  7  \
 *    2   6   7   5             2     o6    5
 *    | 0   \ | 1   \           | 0    6  1   \
 *  0 o---0---o---1---o 1     0 o---0---o---1---o 1
 *            3                         3
 *   (split into 4 parts)      (split into 3 parts)
 *     legacy index 3                  new
 *        index 4                    index 5
 *
 *  2 o                       2 o
 *    | \                       | \
 *    |   2                     4   5
 *    |     \  4                | 1   \  4
 *    1       o               3 o---3---o
 *    |   0   | \               |         \
 *    |       4   5             1     0     2
 *    |       | 1   \           |             \
 *  0 o---0---o---3---o 1     0 o-------0-------o 1
 *            3                         3
 *     (not split in 1)          (not split in 0)
 *           new                       new
 *         index 6                   index 7
 *
 *  2 o
 *    | \
 *    4   \
 *    |     \
 *  4 o       5
 *    | \   1   \
 *    1   2       \
 *    | 0   \       \
 *  0 o---0---o---3---o 1
 *            3
 *     (not split in 1)
 *           new
 *         index 8
 *
 * The legacy implementation used to have an additional row and column in its
 *transformation matrix The final row contains a 1 on the diagonal and 0's
 *elsewhere, the final column is mapping-depended and will be provided in
 *comments with the appropriate function. This column seems to be the constant
 * component of the transformation
 */

class RefinementMapForTriangle0 : public RefinementMapping {
   public:
    static const RefinementMapForTriangle0* instance() {
        static RefinementMapForTriangle0 theInstance;
        return &theInstance;
    }

    std::string getName() const final { return "Identity map"; }

    PointReference<2> refinementTransform(
        std::size_t subElementIndex,
        const PointReference<2>& p) const final {
        logger.assert_debug(
            subElementIndex == 0,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        return p;
    }

    LinearAlgebra::SmallMatrix<2, 2> getRefinementMappingMatrixL(
        std::size_t subElementIndex,
        const PointReference<2>& p) const final {
        logger.assert_debug(
            subElementIndex == 0,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        // the bonus column is most likely {0, 0, 1}
        return LinearAlgebra::SmallMatrix<2, 2>{{{{1., 0.}}, {{0., 1.}}}};
    }

    LinearAlgebra::SmallMatrix<2, 2> getRefinementMappingMatrixR(
        std::size_t subElementIndex,
        const PointReference<2>& p) const final {
        logger.assert_debug(
            subElementIndex == 0,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        // the bonus column is most likely {0, 0, 1}
        return LinearAlgebra::SmallMatrix<2, 2>{{{{1., 0.}}, {{0., 1.}}}};
    }

    std::size_t getNumberOfNewNodes() const final { return 0; }

    std::size_t getNumberOfSubElements() const final { return 1; }

    std::vector<std::size_t> getSubElementLocalNodeIndices(
        std::size_t subElementIndex) const final {
        logger.assert_debug(
            subElementIndex == 0,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        return std::vector<std::size_t>{0, 1, 2};
    }

    Geometry::ReferenceGeometry* getBigElementReferenceGeometry()
        const final {
        return &Geometry::ReferenceTriangle::Instance();
    }

    Geometry::ReferenceGeometry* getSubElementReferenceGeometry(
        std::size_t subElement) const final {
        logger.assert_debug(
            subElement < getNumberOfSubElements(),
            "asked for subElement %, but the % has only % subElements",
            subElement, getName(), getNumberOfSubElements());
        return &Geometry::ReferenceTriangle::Instance();
    }

    std::vector<Geometry::PointReference<2>> getNewNodeLocations(
        const PointReference<2>&) const final {
        return {};
    }

    std::vector<const RefinementMapping*> getCodim1RefinementMaps()
        const final {
        return std::vector<const RefinementMapping*>{
            RefinementMapForLine0::instance(),
            RefinementMapForLine0::instance(),
            RefinementMapForLine0::instance()};
    }

    std::vector<std::size_t> getCodim1LocalNodeIndices(
        std::size_t localFaceNumber) const final {
        switch (localFaceNumber) {
            case 0:
                return {0, 1};
            case 1:
                return {0, 2};
            case 2:
                return {1, 2};
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
    RefinementMapForTriangle0() = default;
    RefinementMapForTriangle0(const RefinementMapForTriangle0&) = delete;
};

class RefinementMapForTriangle1 : public RefinementMapping {
   public:
    static const RefinementMapForTriangle1* instance() {
        static RefinementMapForTriangle1 theInstance;
        return &theInstance;
    }

    std::string getName() const final { return "split face 0"; }

    PointReference<2> refinementTransform(
        std::size_t subElementIndex,
        const PointReference<2>& p) const final {
        logger.assert_debug(
            subElementIndex < 2,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        return {(p[0] + subElementIndex) / 2., p[1]};
    }

    LinearAlgebra::SmallMatrix<2, 2> getRefinementMappingMatrixL(
        std::size_t subElementIndex,
        const PointReference<2>& p) const final {
        logger.assert_debug(
            subElementIndex < 2,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        // the bonus column is probably {subElementIndex / 2, 0, 1}
        return LinearAlgebra::SmallMatrix<2, 2>{{{{0.5, 0.}}, {{0., 1.}}}};
    }

    LinearAlgebra::SmallMatrix<2, 2> getRefinementMappingMatrixR(
        std::size_t subElementIndex,
        const PointReference<2>& p) const final {
        logger.assert_debug(
            subElementIndex < 2,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        // the bonus column is probably {-subelementIndex, 0, 1}
        return LinearAlgebra::SmallMatrix<2, 2>{{{{2., 0.}}, {{0., 1.}}}};
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
            return std::vector<std::size_t>{0, 3, 2};
        } 
            return std::vector<std::size_t>{3, 1, 2};
        
    }

    Geometry::ReferenceGeometry* getBigElementReferenceGeometry()
        const final {
        return &Geometry::ReferenceTriangle::Instance();
    }

    Geometry::ReferenceGeometry* getSubElementReferenceGeometry(
        std::size_t subElement) const final {
        logger.assert_debug(
            subElement < getNumberOfSubElements(),
            "asked for subElement %, but the % has only % subElements",
            subElement, getName(), getNumberOfSubElements());
        return &Geometry::ReferenceTriangle::Instance();
    }

    std::vector<Geometry::PointReference<2>> getNewNodeLocations(
        const PointReference<2>&) const final {
        return {{0.5, 0.}};
    }

    std::vector<const RefinementMapping*> getCodim1RefinementMaps()
        const final {
        return std::vector<const RefinementMapping*>{
            RefinementMapForLine1::instance(),
            RefinementMapForLine0::instance(),
            RefinementMapForLine0::instance()};
    }

    std::vector<std::size_t> getCodim1LocalNodeIndices(
        std::size_t localFaceNumber) const final {
        switch (localFaceNumber) {
            case 0:
                return {0, 1, 3};
            case 1:
                return {0, 2};
            case 2:
                return {1, 2};
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
        switch (face) {
            case 0:
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
    RefinementMapForTriangle1() = default;
    RefinementMapForTriangle1(const RefinementMapForTriangle1&) = delete;
};

class RefinementMapForTriangle2 : public RefinementMapping {
   public:
    static const RefinementMapForTriangle2* instance() {
        static RefinementMapForTriangle2 theInstance;
        return &theInstance;
    }

    std::string getName() const final { return "split face 1"; }

    PointReference<2> refinementTransform(
        std::size_t subElementIndex,
        const PointReference<2>& p) const final {
        logger.assert_debug(
            subElementIndex < 2,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        return {p[0], (p[1] + subElementIndex) / 2.};
    }

    LinearAlgebra::SmallMatrix<2, 2> getRefinementMappingMatrixL(
        std::size_t subElementIndex,
        const PointReference<2>& p) const final {
        logger.assert_debug(
            subElementIndex < 2,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        // the bonus column is probably {0, subElementIndex / 2, 1}
        return LinearAlgebra::SmallMatrix<2, 2>{{{{1., 0.}}, {{0., 0.5}}}};
    }

    LinearAlgebra::SmallMatrix<2, 2> getRefinementMappingMatrixR(
        std::size_t subElementIndex,
        const PointReference<2>& p) const final {
        logger.assert_debug(
            subElementIndex < 2,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        // the bonus column is probably {0, -subelementIndex, 1}
        return LinearAlgebra::SmallMatrix<2, 2>{{{{1., 0.}}, {{0., 2.}}}};
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
            return std::vector<std::size_t>{0, 1, 3};
        } 
            return std::vector<std::size_t>{3, 1, 2};
        
    }

    Geometry::ReferenceGeometry* getBigElementReferenceGeometry()
        const final {
        return &Geometry::ReferenceTriangle::Instance();
    }

    Geometry::ReferenceGeometry* getSubElementReferenceGeometry(
        std::size_t subElement) const final {
        logger.assert_debug(
            subElement < getNumberOfSubElements(),
            "asked for subElement %, but the % has only % subElements",
            subElement, getName(), getNumberOfSubElements());
        return &Geometry::ReferenceTriangle::Instance();
    }

    std::vector<Geometry::PointReference<2>> getNewNodeLocations(
        const PointReference<2>&) const final {
        return {{0., 0.5}};
    }

    std::vector<const RefinementMapping*> getCodim1RefinementMaps()
        const final {
        return std::vector<const RefinementMapping*>{
            RefinementMapForLine0::instance(),
            RefinementMapForLine1::instance(),
            RefinementMapForLine0::instance()};
    }

    std::vector<std::size_t> getCodim1LocalNodeIndices(
        std::size_t localFaceNumber) const final {
        switch (localFaceNumber) {
            case 0:
                return {0, 1};
            case 1:
                return {0, 2, 3};
            case 2:
                return {1, 2};
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
        switch (face) {
            case 0:
                return std::make_tuple(0, face);
            case 1:
                return std::make_tuple(subFaceIndex, face);
            case 2:
                return std::make_tuple(1, face);
            default:
                logger(ERROR, "something broke, please check the assertions");
        }
        return {};
    }

   private:
    RefinementMapForTriangle2() = default;
    RefinementMapForTriangle2(const RefinementMapForTriangle2&) = delete;
};

class RefinementMapForTriangle3 : public RefinementMapping {
   public:
    static const RefinementMapForTriangle3* instance() {
        static RefinementMapForTriangle3 theInstance;
        return &theInstance;
    }

    std::string getName() const final { return "split face 2"; }

    PointReference<2> refinementTransform(
        std::size_t subElementIndex,
        const PointReference<2>& p) const final {
        logger.assert_debug(
            subElementIndex < 2,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        if (subElementIndex == 0) {
            return {0.5 * p[0], 0.5 * p[0] + p[1]};
        } 
            return {p[0] + 0.5 * p[1], 0.5 * p[1]};
        
    }

    LinearAlgebra::SmallMatrix<2, 2> getRefinementMappingMatrixL(
        std::size_t subElementIndex,
        const PointReference<2>& p) const final {
        logger.assert_debug(
            subElementIndex < 2,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        if (subElementIndex == 0) {
            // the bonus column might be {0, 0, 1}
            return LinearAlgebra::SmallMatrix<2, 2>{{{{0.5, 0.5}}, {{0., 1.}}}};
        } 
            // the bonus column might be {0, 0, 1}
            return LinearAlgebra::SmallMatrix<2, 2>{{{{1., 0.}}, {{0.5, 0.5}}}};
        
    }

    LinearAlgebra::SmallMatrix<2, 2> getRefinementMappingMatrixR(
        std::size_t subElementIndex,
        const PointReference<2>& p) const final {
        logger.assert_debug(
            subElementIndex < 2,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        if (subElementIndex == 0) {
            // the bonus column might be {0, 0, 1}
            return LinearAlgebra::SmallMatrix<2, 2>{{{{2., -1.}}, {{0., 1.}}}};
        } 
            // the bonus column might be {0, 0, 1}
            return LinearAlgebra::SmallMatrix<2, 2>{{{{1., 0.}}, {{-1., 2.}}}};
        
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
            return std::vector<std::size_t>{0, 1, 3};
        } 
            return std::vector<std::size_t>{0, 3, 2};
        
    }

    Geometry::ReferenceGeometry* getBigElementReferenceGeometry()
        const final {
        return &Geometry::ReferenceTriangle::Instance();
    }

    Geometry::ReferenceGeometry* getSubElementReferenceGeometry(
        std::size_t subElement) const final {
        logger.assert_debug(
            subElement < getNumberOfSubElements(),
            "asked for subElement %, but the % has only % subElements",
            subElement, getName(), getNumberOfSubElements());
        return &Geometry::ReferenceTriangle::Instance();
    }

    std::vector<Geometry::PointReference<2>> getNewNodeLocations(
        const PointReference<2>&) const final {
        return {{0.5, 0.5}};
    }

    std::vector<const RefinementMapping*> getCodim1RefinementMaps()
        const final {
        return std::vector<const RefinementMapping*>{
            RefinementMapForLine0::instance(),
            RefinementMapForLine0::instance(),
            RefinementMapForLine1::instance()};
    }

    std::vector<std::size_t> getCodim1LocalNodeIndices(
        std::size_t localFaceNumber) const final {
        switch (localFaceNumber) {
            case 0:
                return {0, 1};
            case 1:
                return {0, 2};
            case 2:
                return {1, 2, 3};
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
        switch (face) {
            case 0:
                return std::make_tuple(0, face);
            case 1:
                return std::make_tuple(1, face);
            case 2:
                return std::make_tuple(subFaceIndex, face);
            default:
                logger(ERROR, "something broke, please check the assertions");
        }
        return {};
    }

   private:
    RefinementMapForTriangle3() = default;
    RefinementMapForTriangle3(const RefinementMapForTriangle3&) = delete;
};

class RefinementMapForTriangle4 : public RefinementMapping {
   public:
    static const RefinementMapForTriangle4* instance() {
        static RefinementMapForTriangle4 theInstance;
        return &theInstance;
    }

    std::string getName() const final { return "split all faces"; }

    PointReference<2> refinementTransform(
        std::size_t subElementIndex,
        const PointReference<2>& p) const final {
        switch (subElementIndex) {
            case 0:
                return {0.5 * p[0], 0.5 * p[1]};
            case 1:
                return {0.5 + 0.5 * p[0], 0.5 * p[1]};
            case 2:
                return {0.5 * p[0], 0.5 + 0.5 * p[1]};
            case 3:
                return {0.5 - 0.5 * p[0], 0.5 - 0.5 * p[1]};
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
        const PointReference<2>& p) const final {
        switch (subElementIndex) {
            case 0:
                // the bonus column might be {0, 0, 1}
                return LinearAlgebra::SmallMatrix<2, 2>{
                    {{{0.5, 0.}}, {{0., 0.5}}}};
            case 1:
                // the bonus column might be {0.5, 0, 1}
                return LinearAlgebra::SmallMatrix<2, 2>{
                    {{{0.5, 0.}}, {{0., 0.5}}}};
            case 2:
                // the bonus column might be {0, 0.5, 1}
                return LinearAlgebra::SmallMatrix<2, 2>{
                    {{{0.5, 0.}}, {{0., 0.5}}}};
            case 3:
                // the bonus column might be {0.5, 0.5, 1}
                return LinearAlgebra::SmallMatrix<2, 2>{
                    {{{-0.5, 0.}}, {{0., -0.5}}}};
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
        const PointReference<2>& p) const final {
        // matrix was not present in the legacy code and I can't make a
        // reasonable guess about the bonus column
        switch (subElementIndex) {
            case 0:
                // the bonus column might be {0, 0, 1}
                return LinearAlgebra::SmallMatrix<2, 2>{
                    {{{2., 0.}}, {{0., 2.0}}}};
            case 1:
                // the bonus column might be {-1, 0, 1}
                return LinearAlgebra::SmallMatrix<2, 2>{
                    {{{2., 0.}}, {{0., 2.0}}}};
            case 2:
                // the bonus column might be {0, -1, 1}
                return LinearAlgebra::SmallMatrix<2, 2>{
                    {{{2., 0.}}, {{0., 2.0}}}};
            case 3:
                // the bonus column might be {-1, -1, 1}
                return LinearAlgebra::SmallMatrix<2, 2>{
                    {{{-2., 0.}}, {{0., -2.0}}}};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return LinearAlgebra::SmallMatrix<2, 2>();
        }
    }

    std::size_t getNumberOfNewNodes() const final { return 3; }

    std::size_t getNumberOfSubElements() const final { return 4; }

    std::vector<std::size_t> getSubElementLocalNodeIndices(
        std::size_t subElementIndex) const final {
        switch (subElementIndex) {
            case 0:
                return std::vector<std::size_t>{0, 3, 4};
            case 1:
                return std::vector<std::size_t>{3, 1, 5};
            case 2:
                return std::vector<std::size_t>{4, 5, 2};
            case 3:
                return std::vector<std::size_t>{5, 4, 3};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return std::vector<std::size_t>();
        }
    }

    Geometry::ReferenceGeometry* getBigElementReferenceGeometry()
        const final {
        return &Geometry::ReferenceTriangle::Instance();
    }

    Geometry::ReferenceGeometry* getSubElementReferenceGeometry(
        std::size_t subElement) const final {
        logger.assert_debug(
            subElement < getNumberOfSubElements(),
            "asked for subElement %, but the % has only % subElements",
            subElement, getName(), getNumberOfSubElements());
        return &Geometry::ReferenceTriangle::Instance();
    }

    std::vector<Geometry::PointReference<2>> getNewNodeLocations(
        const PointReference<2>&) const final {
        return {{0.5, 0.}, {0., 0.5}, {0.5, 0.5}};
    }

    std::vector<const RefinementMapping*> getCodim1RefinementMaps()
        const final {
        return std::vector<const RefinementMapping*>{
            RefinementMapForLine1::instance(),
            RefinementMapForLine1::instance(),
            RefinementMapForLine1::instance()};
    }

    std::vector<std::size_t> getCodim1LocalNodeIndices(
        std::size_t localFaceNumber) const final {
        switch (localFaceNumber) {
            case 0:
                return {0, 1, 3};
            case 1:
                return {0, 2, 4};
            case 2:
                return {1, 2, 5};
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
        switch (face) {
            case 0:
                return std::make_tuple(subFaceIndex, face);
            case 1:
                return std::make_tuple(subFaceIndex * 2, face);
            case 2:
                return std::make_tuple(subFaceIndex + 1, face);
            default:
                logger(ERROR, "something broke, please check the assertions");
        }
        return {};
    }

   private:
    RefinementMapForTriangle4() = default;
    RefinementMapForTriangle4(const RefinementMapForTriangle4&) = delete;
};

class RefinementMapForTriangle5 : public RefinementMapping {
   public:
    static const RefinementMapForTriangle5* instance() {
        static RefinementMapForTriangle5 theInstance;
        return &theInstance;
    }

    std::string getName() const final {
        return "split to 3 quadrilaterals";
    }

    PointReference<2> refinementTransform(
        std::size_t subElementIndex,
        const PointReference<2>& p) const final {
        switch (subElementIndex) {
            case 0:
                return {(p[1] + 1.) * (5. + p[0]) / 24.,
                        (1. - p[0]) * (5. - p[1]) / 24.};
            case 1:
                return {((p[0] + 2.) * (2. - p[1]) + 3.) / 12.,
                        (p[1] + 1.) * (5. + p[0]) / 24.};
            case 2:
                return {(1. - p[0]) * (5. - p[1]) / 24.,
                        ((p[0] + 2.) * (2. - p[1]) + 3.) / 12.};
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
        const PointReference<2>& p) const final {
        switch (subElementIndex) {
            case 0:
                // the bonus column should be {5./24., 5./24, 1}
                return LinearAlgebra::SmallMatrix<2, 2>{
                    {{{(p[1] + 1.) / 24., (p[1] - 5.) / 24.}},
                     {{(p[0] + 5.) / 24., (p[0] - 1.) / 24.}}}};
            case 1:
                // the bonus column should be {14./24., 5./24., 1}
                return LinearAlgebra::SmallMatrix<2, 2>{
                    {{{(2. - p[1]) / 12., (p[1] + 1.) / 24.}},
                     {{-(p[0] + 2.) / 12., (p[0] + 5.) / 24.}}}};
            case 2:
                // the bonus column should be {5./24., 14./24., 1}
                return LinearAlgebra::SmallMatrix<2, 2>{
                    {{{(p[1] - 5.) / 24., (2. - p[1]) / 12.}},
                     {{(p[0] - 1.) / 24., -(p[0] + 2.) / 12.}}}};
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
        const PointReference<2>& p) const final {
        switch (subElementIndex) {
            case 0:
                return LinearAlgebra::SmallMatrix<2, 2>{
                           {{{p[0] - 1., 5. - p[1]}},
                            {{-(p[0] + 5.), (p[1] + 1.)}}}} *
                       4. / (p[0] - p[1] + 4.);
            case 1:
                return LinearAlgebra::SmallMatrix<2, 2>{
                           {{{p[0] + 5., -(p[1] + 1.)}},
                            {{2. * p[0] + 4., 4. - 2. * p[1]}}}} *
                       4. / (p[0] - p[1] + 4.);
            case 2:
                return LinearAlgebra::SmallMatrix<2, 2>{
                           {{{-(2. * p[0] + 4.), 2. * p[1] - 4.}},
                            {{1. - p[0], p[1] - 5.}}}} *
                       4. / (p[0] - p[1] + 4.);
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return LinearAlgebra::SmallMatrix<2, 2>();
        }
    }

    std::size_t getNumberOfNewNodes() const final { return 4; }

    std::size_t getNumberOfSubElements() const final { return 3; }

    std::vector<std::size_t> getSubElementLocalNodeIndices(
        std::size_t subElementIndex) const final {
        switch (subElementIndex) {
            case 0:
                return std::vector<std::size_t>{4, 0, 6, 3};
            case 1:
                return std::vector<std::size_t>{3, 1, 6, 5};
            case 2:
                return std::vector<std::size_t>{5, 2, 6, 4};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return std::vector<std::size_t>();
        }
    }

    Geometry::ReferenceGeometry* getBigElementReferenceGeometry()
        const final {
        return &Geometry::ReferenceTriangle::Instance();
    }

    Geometry::ReferenceGeometry* getSubElementReferenceGeometry(
        std::size_t subElement) const final {
        logger.assert_debug(
            subElement < getNumberOfSubElements(),
            "asked for subElement %, but the % has only % subElements",
            subElement, getName(), getNumberOfSubElements());
        return &Geometry::ReferenceSquare::Instance();
    }

    std::vector<Geometry::PointReference<2>> getNewNodeLocations(
        const PointReference<2>&) const final {
        return {{0.5, 0.}, {0., 0.5}, {0.5, 0.5}, {1. / 3., 1. / 3.}};
    }

    std::vector<const RefinementMapping*> getCodim1RefinementMaps()
        const final {
        return std::vector<const RefinementMapping*>{
            RefinementMapForLine1::instance(),
            RefinementMapForLine1::instance(),
            RefinementMapForLine1::instance()};
    }

    std::vector<std::size_t> getCodim1LocalNodeIndices(
        std::size_t localFaceNumber) const final {
        switch (localFaceNumber) {
            case 0:
                return {0, 1, 3};
            case 1:
                return {0, 2, 4};
            case 2:
                return {1, 2, 5};
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
        switch (face) {
            case 0:
                return std::make_tuple(subFaceIndex, 2 - 2 * subFaceIndex);
            case 1:
                return std::make_tuple(subFaceIndex * 2, 2 * subFaceIndex);
            case 2:
                return std::make_tuple(subFaceIndex + 1, 2 - 2 * subFaceIndex);
            default:
                logger(ERROR, "something broke, please check the assertions");
        }
        return {};
    }

   private:
    RefinementMapForTriangle5() = default;
    RefinementMapForTriangle5(const RefinementMapForTriangle5&) = delete;
};

class RefinementMapForTriangle6 : public RefinementMapping {
   public:
    static const RefinementMapForTriangle6* instance() {
        static RefinementMapForTriangle6 theInstance;
        return &theInstance;
    }

    std::string getName() const final {
        return "split to a triangle and a quadrilateral";
    }

    PointReference<2> refinementTransform(
        std::size_t subElementIndex,
        const PointReference<2>& p) const final {
        switch (subElementIndex) {
            case 0:
                return {(1. + p[0]) / 4., (p[1] + 1.) * (3. - p[0]) / 8.};
            case 1:
                return {(1. + p[0]) / 2., p[1] / 2.};
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
        const PointReference<2>& p) const final {
        switch (subElementIndex) {
            case 0:
                // the bonus column should be {.25, .375, 1}
                return LinearAlgebra::SmallMatrix<2, 2>{
                    {{{.25, -(p[1] + 1.) / 8.}}, {{0., (3. - p[0]) / 8.}}}};
            case 1:
                // the bonus column should be {.5, 0, 1}
                return LinearAlgebra::SmallMatrix<2, 2>{
                    {{{.5, 0.}}, {{0., .5}}}};
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
        const PointReference<2>& p) const final {
        switch (subElementIndex) {
            case 0:
                return LinearAlgebra::SmallMatrix<2, 2>{
                           {{{3. - p[0], p[1] + 1.}}, {{0., 2.}}}} *
                       4. / (3. - p[0]);
            case 1:
                return LinearAlgebra::SmallMatrix<2, 2>{
                    {{{2., 0.}}, {{0., 2.}}}};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return LinearAlgebra::SmallMatrix<2, 2>();
        }
    }

    std::size_t getNumberOfNewNodes() const final { return 2; }

    std::size_t getNumberOfSubElements() const final { return 2; }

    std::vector<std::size_t> getSubElementLocalNodeIndices(
        std::size_t subElementIndex) const final {
        switch (subElementIndex) {
            case 0:
                return std::vector<std::size_t>{0, 3, 2, 4};
            case 1:
                return std::vector<std::size_t>{3, 1, 4};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return std::vector<std::size_t>();
        }
    }

    Geometry::ReferenceGeometry* getBigElementReferenceGeometry()
        const final {
        return &Geometry::ReferenceTriangle::Instance();
    }

    Geometry::ReferenceGeometry* getSubElementReferenceGeometry(
        std::size_t subElement) const final {
        switch (subElement) {
            case 0:
                return &Geometry::ReferenceSquare::Instance();
            case 1:
                return &Geometry::ReferenceTriangle::Instance();
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElement, getName(), getNumberOfSubElements());
                return nullptr;
        }
    }

    std::vector<Geometry::PointReference<2>> getNewNodeLocations(
        const PointReference<2>&) const final {
        return {{0.5, 0.}, {0.5, 0.5}};
    }

    std::vector<const RefinementMapping*> getCodim1RefinementMaps()
        const final {
        return std::vector<const RefinementMapping*>{
            RefinementMapForLine1::instance(),
            RefinementMapForLine0::instance(),
            RefinementMapForLine1::instance()};
    }

    std::vector<std::size_t> getCodim1LocalNodeIndices(
        std::size_t localFaceNumber) const final {
        switch (localFaceNumber) {
            case 0:
                return {0, 1, 3};
            case 1:
                return {0, 2};
            case 2:
                return {1, 2, 4};
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
        switch (face) {
            case 0:
                return std::make_tuple(subFaceIndex, face);
            case 1:
                return std::make_tuple(0, face);
            case 2:
                return std::make_tuple(subFaceIndex, subFaceIndex + 2);
            default:
                logger(ERROR, "something broke, please check the assertions");
        }
        return {};
    }

   private:
    RefinementMapForTriangle6() = default;
    RefinementMapForTriangle6(const RefinementMapForTriangle6&) = delete;
};

class RefinementMapForTriangle7 : public RefinementMapping {
   public:
    static const RefinementMapForTriangle7* instance() {
        static RefinementMapForTriangle7 theInstance;
        return &theInstance;
    }

    std::string getName() const final {
        return "split to a triangle and a quadrilateral";
    }

    PointReference<2> refinementTransform(
        std::size_t subElementIndex,
        const PointReference<2>& p) const final {
        switch (subElementIndex) {
            case 0:
                return {(p[0] + 1.) * (3. - p[1]) / 8., (1. + p[1]) / 4.};
            case 1:
                return {p[0] / 2., (1 + p[1]) / 2.};
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
        const PointReference<2>& p) const final {
        switch (subElementIndex) {
            case 0:
                // the bonus column should be {.375, .25, 1}
                return LinearAlgebra::SmallMatrix<2, 2>{
                    {{{(3. - p[1]) / 8., 0.}}, {{-(p[0] + 1.) / 8., .25}}}};
            case 1:
                // the bonus column should be {0, .5, 1}
                return LinearAlgebra::SmallMatrix<2, 2>{
                    {{{.5, 0.}}, {{0., .5}}}};
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
        const PointReference<2>& p) const final {
        switch (subElementIndex) {
            case 0:
                return LinearAlgebra::SmallMatrix<2, 2>{
                           {{{2., 0.}}, {{p[0] + 1., 3. - p[1]}}}} *
                       4. / (3. - p[1]);
            case 1:
                return LinearAlgebra::SmallMatrix<2, 2>{
                    {{{2., 0.}}, {{0., 2.}}}};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return LinearAlgebra::SmallMatrix<2, 2>();
        }
    }

    std::size_t getNumberOfNewNodes() const final { return 2; }

    std::size_t getNumberOfSubElements() const final { return 2; }

    std::vector<std::size_t> getSubElementLocalNodeIndices(
        std::size_t subElementIndex) const final {
        switch (subElementIndex) {
            case 0:
                return std::vector<std::size_t>{0, 1, 3, 4};
            case 1:
                return std::vector<std::size_t>{3, 4, 2};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return std::vector<std::size_t>();
        }
    }

    Geometry::ReferenceGeometry* getBigElementReferenceGeometry()
        const final {
        return &Geometry::ReferenceTriangle::Instance();
    }

    Geometry::ReferenceGeometry* getSubElementReferenceGeometry(
        std::size_t subElement) const final {
        switch (subElement) {
            case 0:
                return &Geometry::ReferenceSquare::Instance();
            case 1:
                return &Geometry::ReferenceTriangle::Instance();
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElement, getName(), getNumberOfSubElements());
                return nullptr;
        }
    }

    std::vector<Geometry::PointReference<2>> getNewNodeLocations(
        const PointReference<2>&) const final {
        return {{0., 0.5}, {0.5, 0.5}};
    }

    std::vector<const RefinementMapping*> getCodim1RefinementMaps()
        const final {
        return std::vector<const RefinementMapping*>{
            RefinementMapForLine0::instance(),
            RefinementMapForLine1::instance(),
            RefinementMapForLine1::instance()};
    }

    std::vector<std::size_t> getCodim1LocalNodeIndices(
        std::size_t localFaceNumber) const final {
        switch (localFaceNumber) {
            case 0:
                return {0, 1};
            case 1:
                return {0, 2, 3};
            case 2:
                return {1, 2, 4};
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
        switch (face) {
            case 0:
                return std::make_tuple(0, face);
            case 1:
                return std::make_tuple(subFaceIndex, face);
            case 2:
                return std::make_tuple(subFaceIndex, face);
            default:
                logger(ERROR, "something broke, please check the assertions");
        }
        return {};
    }

   private:
    RefinementMapForTriangle7() = default;
    RefinementMapForTriangle7(const RefinementMapForTriangle7&) = delete;
};

class RefinementMapForTriangle8 : public RefinementMapping {
   public:
    static const RefinementMapForTriangle8* instance() {
        static RefinementMapForTriangle8 theInstance;
        return &theInstance;
    }

    std::string getName() const final {
        return "split to a triangle and a quadrilateral";
    }

    PointReference<2> refinementTransform(
        std::size_t subElementIndex,
        const PointReference<2>& p) const final {
        switch (subElementIndex) {
            case 0:
                return {p[0] / 2., p[1] / 2.};
            case 1:
                return {(1. - p[1]) * (3. + p[0]) / 8.,
                        (1. + p[1]) * (3 + p[0]) / 8.};
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
        const PointReference<2>& p) const final {
        switch (subElementIndex) {
            case 0:
                // the bonus column should be {0, 0, 1}
                return LinearAlgebra::SmallMatrix<2, 2>{
                    {{{.5, 0.}}, {{0., .5}}}};
            case 1:
                // the bonus column should be {.375, .375, 1}
                return LinearAlgebra::SmallMatrix<2, 2>{
                    {{{(1. - p[1]) / 8., (1. + p[1]) / 8.}},
                     {{-(3. + p[0]) / 8., (3. + p[0]) / 8.}}}};
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
        const PointReference<2>& p) const final {
        switch (subElementIndex) {
            case 0:
                return LinearAlgebra::SmallMatrix<2, 2>{
                    {{{2., 0.}}, {{0., 2.}}}};
            case 1:
                return LinearAlgebra::SmallMatrix<2, 2>{
                           {{{3 + p[0], -(1. + p[1])}},
                            {{3. + p[0], 1. - p[1]}}}} *
                       4. / (p[0] + 3.);
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return LinearAlgebra::SmallMatrix<2, 2>();
        }
    }

    std::size_t getNumberOfNewNodes() const final { return 2; }

    std::size_t getNumberOfSubElements() const final { return 2; }

    std::vector<std::size_t> getSubElementLocalNodeIndices(
        std::size_t subElementIndex) const final {
        switch (subElementIndex) {
            case 0:
                return std::vector<std::size_t>{0, 3, 4};
            case 1:
                return std::vector<std::size_t>{3, 1, 4, 2};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return std::vector<std::size_t>();
        }
    }

    Geometry::ReferenceGeometry* getBigElementReferenceGeometry()
        const final {
        return &Geometry::ReferenceTriangle::Instance();
    }

    Geometry::ReferenceGeometry* getSubElementReferenceGeometry(
        std::size_t subElement) const final {
        switch (subElement) {
            case 0:
                return &Geometry::ReferenceTriangle::Instance();
            case 1:
                return &Geometry::ReferenceSquare::Instance();
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElement, getName(), getNumberOfSubElements());
                return nullptr;
        }
    }

    std::vector<Geometry::PointReference<2>> getNewNodeLocations(
        const PointReference<2>&) const final {
        return {{0.5, 0.}, {0., 0.5}};
    }

    std::vector<const RefinementMapping*> getCodim1RefinementMaps()
        const final {
        return std::vector<const RefinementMapping*>{
            RefinementMapForLine0::instance(),
            RefinementMapForLine1::instance(),
            RefinementMapForLine1::instance()};
    }

    std::vector<std::size_t> getCodim1LocalNodeIndices(
        std::size_t localFaceNumber) const final {
        switch (localFaceNumber) {
            case 0:
                return {0, 1, 3};
            case 1:
                return {0, 2, 4};
            case 2:
                return {1, 2};
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
        switch (face) {
            case 0:
                return std::make_tuple(subFaceIndex, face);
            case 1:
                return std::make_tuple(subFaceIndex, subFaceIndex * 2 + 1);
            case 2:
                return std::make_tuple(1, face);
            default:
                logger(ERROR, "something broke, please check the assertions");
        }
        return {};
    }

   private:
    RefinementMapForTriangle8() = default;
    RefinementMapForTriangle8(const RefinementMapForTriangle8&) = delete;
};
}  // namespace Geometry

#endif /* KERNEL_GEOMETRY_MAPPINGS_REFINEMENTMAPSFORTRIANGLE_H_ */
