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

#ifndef KERNEL_GEOMETRY_MAPPINGS_REFINEMENTMAPSFORTRIANGULARPRISM_H_
#define KERNEL_GEOMETRY_MAPPINGS_REFINEMENTMAPSFORTRIANGULARPRISM_H_

#include "RefinementMapping.h"
#include "RefinementMapsForTriangle.h"
#include "RefinementMapsForSquare.h"
#include "Geometry/PointReference.h"
#include "Geometry/ReferenceCube.h"
#include "Geometry/ReferenceTriangularPrism.h"

namespace Geometry {
/* Stores all the refinement mappings for a pyramid. Note that the ascii art is
 * provided in a different perspective than the other 3D geometries Numbering
 * for faces and edges has been omitted to prevent excessive cluttering
 *
 *          5 o                             5 o                             5 o
 *           /| \                            /|\\                            /|
 * \
 *         /  |   \                        /  | \ \                        /  |
 * \
 *        /   |     \                     /   |  \  \                     /   |
 * \
 *       /    |       \                  /    |   \   \                  /  7 o
 * \
 *     /      |         \              /      |    \    \              /     /|
 * - _     \
 *  2 o       |           \         2 o       |     \     \         2 o  1 /  |
 * - _   \
 *    | \     |             \         |\\     |      \7     \         | \ /   |
 * - _ \ |   \ 3 o---------------o 4     | \ \ 3 o-------o-------o 4     |  /\ 3
 * o---------------o 4 |     \/               /        |  \  \/       /       /
 * |/    \/               / |    /  \   0        /          |   \/  \    / / 6 o
 * /  \            / |   /     \         /           |   /\    \ /   1   / | -_/
 * \         / |  /        \      /            |  /  \    /\      /            |
 * /  - _   \      /
 *    |/            \  /              |/  0  \ /    \  /              |/  0 - _
 * \  / 0 o---------------o 1           0 o-------o-------o 1           0
 * o---------------o 1
 *                                            6
 *         (identity)              (split in the middle of 4)      (split in the
 * middle of 3) new                       legacy index 1                  legacy
 * index 0 index 0                        index 1                         index
 * 2
 *
 *          5 o                             5 o                             5 o
 *           /| \                            /| \                            /|
 * \
 *         /  |   \                        /  |   \                        /  |
 * \
 *        /   |     \                     /   |     \                   8 o   |
 * \
 *       /    |       o 7                / 10 o-------o 11               /|\  |
 * \
 *     /      |     //  \              /     /| \    /| \              /  |  \|
 * \
 *  2 o       |   //      \         2 o    /  |   \/  |   \         2 o   |   |\
 * \
 *    | \     | / /         \         | \ /   |   / \ |     \         | \ |   |
 * \          \ |   \ 3 o--/------------o 4     |2 /\ 3 o--/----o-------o 4 | \
 * 3 o----\----------o 4 | 0   \/ /             /        |/    \/ /     /9 / |
 * | \/       \   1   / |    /  o 6          /        7 o-------o 8  /       /
 * |   |/  \        \   / |   / /   \         /           | \ / 3 | \ /       /
 * | 6 o-----\---------o 7 |  //       \      /            |  /\   |  /\      /
 * |  /        \      /
 *    |//     1     \  /              |/ 0  \ |/  1 \  /              |/   0 \ /
 *  0 o---------------o 1           0 o-------o-------o 1           0
 * o---------------o 1
 *                                            6
 *  (split in the middle of 5)        (split into 4 parts)           (split in
 * the third axis) legacy index 2*                 legacy index 4 legacy index 3
 *          index 3                         index 4 index 5
 *
 *   *legacy implementation had the sub-elements flipped
 *
 *          5 o                             5 o                             5 o
 *           /| \                            /| \                            /|
 * \
 *         /  |   \                        /  |   \                        /  |
 * \
 *        /   |     \                     /   |     \                     /   |
 * \ / 11 o _     o 12               /    |       o 9                /  8
 * o-------o 9
 *     /  2  /|  - _ /  \              /      |      /| \              /     /|
 * /  \
 *  2 o    /  |     o 13  \         2 o    0  |    /  |   \         2 o  1 /  |
 * /      \
 *    | \ /   |   // \10    \         | \     |   /   |     \         | \ /   |
 * /         \ |  /\ 3 o--/----o-------o 4     |   \ 3 o--/----o-------o 4     |
 * /\ 3 o--/------------o 4
 *    |/    \/ //    /       /        |     \/ /     /8      /        |/    \/ /
 * / 7 o _  /  o 8  /       /          |    /  o 7  /       /        6 o-------o
 * 7          / |  -/_ /  \ /   1   /           |   /   | \ /   1   / |   / \ /
 *    |  /  o 9  /\      /            |  /    |  /\      /            |  /   0
 * \      /
 *    |/  0  \ /    \  /              |/      |/    \  /              |/ \  / 0
 * o-------o-------o 1           0 o-------o-------o 1           0
 * o---------------o 1 6                               6 (split into 3 parts)
 * (not split in the middle of 3)   (not split in the middle of 4) legacy index
 * 5                      legacy index 6                   legacy index 7 index
 * 6                             index 7                          index 8
 *
 *          5 o
 *           /| \
 *         /  |   \
 *        /   |     \
 *       /  9 o       \
 *     /     /| \       \
 *  2 o    /  |   \       \
 *    | \ /   |     \       \
 *    |  /\ 3 o-------o-------o 4
 *    |/    \/       /8      /
 *  7 o    /  \    /   1   /
 *    | \ /     \ /       /
 *    |  /\      /\      /
 *    |/ 0  \  /    \  /
 *  0 o-------o-------o 1
 *            6
 * (not split in the middle of 5)
 *             new
 *           index 9
 *
 * The legacy implementation used to have an additional row and column in its
 * transformation matrix The final row contains a 1 on the diagonal and 0's
 * elsewhere, the final column is mapping-depended and will be provided in
 * comments with the appropriate function. This column seems to be the constant
 * component of the transformation
 *
 *
 *
 * The legacy implementation of refinementTransform also provides the indices
 * 20-25, these appear to be refinement maps for cubes so this information will
 * be presented with the refinement maps for cubes
 */

class RefinementMapForTriangularPrism0 : public RefinementMapping {
   public:
    static const RefinementMapForTriangularPrism0* instance() {
        static RefinementMapForTriangularPrism0 theInstance;
        return &theInstance;
    }

    std::string getName() const override final { return "Identity map"; }

    PointReference<3> refinementTransform(
        std::size_t subElementIndex,
        const PointReference<3>& p) const override final {
        logger.assert_debug(
            subElementIndex == 0,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        return p;
    }

    LinearAlgebra::SmallMatrix<3, 3> getRefinementMappingMatrixL(
        std::size_t subElementIndex,
        const PointReference<3>& p) const override final {
        logger.assert_debug(
            subElementIndex == 0,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        // the bonus column is most likely {0, 0, 0, 1}
        return LinearAlgebra::SmallMatrix<3, 3>{
            {{{1., 0., 0.}}, {{0., 1., 0.}}, {{0., 0., 1.}}}};
    }

    LinearAlgebra::SmallMatrix<3, 3> getRefinementMappingMatrixR(
        std::size_t subElementIndex,
        const PointReference<3>& p) const override final {
        logger.assert_debug(
            subElementIndex == 0,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        // the bonus column is most likely {0, 0, 0, 1}
        return LinearAlgebra::SmallMatrix<3, 3>{
            {{{1., 0., 0.}}, {{0., 1., 0.}}, {{0., 0., 1.}}}};
    }

    std::size_t getNumberOfNewNodes() const override final { return 0; }

    std::size_t getNumberOfSubElements() const override final { return 1; }

    std::vector<std::size_t> getSubElementLocalNodeIndices(
        std::size_t subElementIndex) const override final {
        logger.assert_debug(
            subElementIndex == 0,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        return std::vector<std::size_t>{0, 1, 2, 3, 4, 5};
    }

    Geometry::ReferenceGeometry* getBigElementReferenceGeometry()
        const override final {
        return &Geometry::ReferenceTriangularPrism::Instance();
    }

    Geometry::ReferenceGeometry* getSubElementReferenceGeometry(
        std::size_t subElement) const override final {
        logger.assert_debug(
            subElement < getNumberOfSubElements(),
            "asked for subElement %, but the % has only % subElements",
            subElement, getName(), getNumberOfSubElements());
        return &Geometry::ReferenceTriangularPrism::Instance();
    }

    std::vector<Geometry::PointReference<3>> getNewNodeLocations(
        const PointReference<3>&) const override final {
        return {};
    }

    std::vector<const RefinementMapping*> getCodim1RefinementMaps()
        const override final {
        return std::vector<const RefinementMapping*>{
            RefinementMapForTriangle0::instance(),
            RefinementMapForTriangle0::instance(),
            RefinementMapForSquare0::instance(),
            RefinementMapForSquare0::instance(),
            RefinementMapForSquare0::instance()};
    }

    std::vector<std::size_t> getCodim1LocalNodeIndices(
        std::size_t localFaceNumber) const override final {
        switch (localFaceNumber) {
            case 0:
                return {0, 2, 1};
            case 1:
                return {3, 4, 5};
            case 2:
                return {2, 0, 5, 3};
            case 3:
                return {0, 1, 3, 4};
            case 4:
                return {1, 2, 4, 5};
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
    RefinementMapForTriangularPrism0() = default;
    RefinementMapForTriangularPrism0(const RefinementMapForTriangularPrism0&) =
        delete;
};

class RefinementMapForTriangularPrism1 : public RefinementMapping {
   public:
    static const RefinementMapForTriangularPrism1* instance() {
        static RefinementMapForTriangularPrism1 theInstance;
        return &theInstance;
    }

    std::string getName() const override final { return "split in 2"; }

    PointReference<3> refinementTransform(
        std::size_t subElementIndex,
        const PointReference<3>& p) const override final {
        logger.assert_debug(
            subElementIndex < 2,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        return {(p[0] + subElementIndex) / 2., p[1], p[2]};
    }

    LinearAlgebra::SmallMatrix<3, 3> getRefinementMappingMatrixL(
        std::size_t subElementIndex,
        const PointReference<3>& p) const override final {
        logger.assert_debug(
            subElementIndex < 2,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        // the bonus column is probably {subElementIndex / 2, 0, 0, 1}
        return LinearAlgebra::SmallMatrix<3, 3>{
            {{{0.5, 0., 0.}}, {{0., 1., 0.}}, {{0., 0., 1.}}}};
    }

    LinearAlgebra::SmallMatrix<3, 3> getRefinementMappingMatrixR(
        std::size_t subElementIndex,
        const PointReference<3>& p) const override final {
        logger.assert_debug(
            subElementIndex < 2,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        // the bonus column is most likely {-subelementIndex, 0, 0, 1}
        return LinearAlgebra::SmallMatrix<3, 3>{
            {{{2., 0., 0.}}, {{0., 1., 0.}}, {{0., 0., 1.}}}};
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
            return std::vector<std::size_t>{0, 6, 2, 3, 7, 5};
        } else {
            return std::vector<std::size_t>{6, 1, 2, 7, 4, 5};
        }
    }

    Geometry::ReferenceGeometry* getBigElementReferenceGeometry()
        const override final {
        return &Geometry::ReferenceTriangularPrism::Instance();
    }

    Geometry::ReferenceGeometry* getSubElementReferenceGeometry(
        std::size_t subElement) const override final {
        logger.assert_debug(
            subElement < getNumberOfSubElements(),
            "asked for subElement %, but the % has only % subElements",
            subElement, getName(), getNumberOfSubElements());
        return &Geometry::ReferenceTriangularPrism::Instance();
    }

    std::vector<Geometry::PointReference<3>> getNewNodeLocations(
        const PointReference<3>&) const override final {
        return {{0.5, 0., -1.}, {0.5, 0., 1.}};
    }

    std::vector<const RefinementMapping*> getCodim1RefinementMaps()
        const override final {
        return std::vector<const RefinementMapping*>{
            RefinementMapForTriangle2::instance(),
            RefinementMapForTriangle1::instance(),
            RefinementMapForSquare0::instance(),
            RefinementMapForSquare1::instance(),
            RefinementMapForSquare0::instance()};
    }

    std::vector<std::size_t> getCodim1LocalNodeIndices(
        std::size_t localFaceNumber) const override final {
        switch (localFaceNumber) {
            case 0:
                return {0, 2, 1, 6};
            case 1:
                return {3, 4, 5, 7};
            case 2:
                return {2, 0, 5, 3};
            case 3:
                return {0, 1, 3, 4, 6, 7};
            case 4:
                return {1, 2, 4, 5};
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
                return std::make_tuple(subFaceIndex, face);
            case 2:
                return std::make_tuple(0, face);
            case 3:
                return std::make_tuple(subFaceIndex, face);
            case 4:
                return std::make_tuple(1, face);
            default:
                logger(ERROR, "something broke, please check the assertions");
        }
        return {};
    }

   private:
    RefinementMapForTriangularPrism1() = default;
    RefinementMapForTriangularPrism1(const RefinementMapForTriangularPrism1&) =
        delete;
};

class RefinementMapForTriangularPrism2 : public RefinementMapping {
   public:
    static const RefinementMapForTriangularPrism2* instance() {
        static RefinementMapForTriangularPrism2 theInstance;
        return &theInstance;
    }

    std::string getName() const override final { return "split in 2"; }

    PointReference<3> refinementTransform(
        std::size_t subElementIndex,
        const PointReference<3>& p) const override final {
        logger.assert_debug(
            subElementIndex < 2,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        return {p[0], (p[1] + subElementIndex) / 2., p[2]};
    }

    LinearAlgebra::SmallMatrix<3, 3> getRefinementMappingMatrixL(
        std::size_t subElementIndex,
        const PointReference<3>& p) const override final {
        logger.assert_debug(
            subElementIndex < 2,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        // the bonus column is probably {subElementIndex / 2, 0, 0, 1}
        return LinearAlgebra::SmallMatrix<3, 3>{
            {{{1., 0., 0.}}, {{0., 0.5, 0.}}, {{0., 0., 1.}}}};
    }

    LinearAlgebra::SmallMatrix<3, 3> getRefinementMappingMatrixR(
        std::size_t subElementIndex,
        const PointReference<3>& p) const override final {
        logger.assert_debug(
            subElementIndex < 2,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        // the bonus column is most likely {-subelementIndex, 0, 0, 1}
        return LinearAlgebra::SmallMatrix<3, 3>{
            {{{1., 0., 0.}}, {{0., 2., 0.}}, {{0., 0., 1.}}}};
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
            return std::vector<std::size_t>{0, 1, 6, 3, 4, 7};
        } else {
            return std::vector<std::size_t>{6, 1, 2, 7, 4, 5};
        }
    }

    Geometry::ReferenceGeometry* getBigElementReferenceGeometry()
        const override final {
        return &Geometry::ReferenceTriangularPrism::Instance();
    }

    Geometry::ReferenceGeometry* getSubElementReferenceGeometry(
        std::size_t subElement) const override final {
        logger.assert_debug(
            subElement < getNumberOfSubElements(),
            "asked for subElement %, but the % has only % subElements",
            subElement, getName(), getNumberOfSubElements());
        return &Geometry::ReferenceTriangularPrism::Instance();
    }

    std::vector<Geometry::PointReference<3>> getNewNodeLocations(
        const PointReference<3>&) const override final {
        return {{0., 0.5, -1.}, {0., 0.5, 1.}};
    }

    std::vector<const RefinementMapping*> getCodim1RefinementMaps()
        const override final {
        return std::vector<const RefinementMapping*>{
            RefinementMapForTriangle1::instance(),
            RefinementMapForTriangle2::instance(),
            RefinementMapForSquare1::instance(),
            RefinementMapForSquare0::instance(),
            RefinementMapForSquare0::instance()};
    }

    std::vector<std::size_t> getCodim1LocalNodeIndices(
        std::size_t localFaceNumber) const override final {
        switch (localFaceNumber) {
            case 0:
                return {0, 2, 1, 6};
            case 1:
                return {3, 4, 5, 7};
            case 2:
                return {2, 0, 5, 3, 6, 7};
            case 3:
                return {0, 1, 3, 4};
            case 4:
                return {1, 2, 4, 5};
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
                return std::make_tuple(subFaceIndex, face);
            case 2:
                return std::make_tuple(subFaceIndex, face);
            case 3:
                return std::make_tuple(0, face);
            case 4:
                return std::make_tuple(1, face);
            default:
                logger(ERROR, "something broke, please check the assertions");
        }
        return {};
    }

   private:
    RefinementMapForTriangularPrism2() = default;
    RefinementMapForTriangularPrism2(const RefinementMapForTriangularPrism2&) =
        delete;
};

class RefinementMapForTriangularPrism3 : public RefinementMapping {
   public:
    static const RefinementMapForTriangularPrism3* instance() {
        static RefinementMapForTriangularPrism3 theInstance;
        return &theInstance;
    }

    std::string getName() const override final { return "split in 2"; }

    PointReference<3> refinementTransform(
        std::size_t subElementIndex,
        const PointReference<3>& p) const override final {
        logger.assert_debug(
            subElementIndex < 2,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        if (subElementIndex == 0) {
            return {0.5 * p[0], 0.5 * p[0] + p[1], p[2]};
        } else {
            return {p[0] + 0.5 * p[1], 0.5 * p[1], p[2]};
        }
    }

    LinearAlgebra::SmallMatrix<3, 3> getRefinementMappingMatrixL(
        std::size_t subElementIndex,
        const PointReference<3>& p) const override final {
        logger.assert_debug(
            subElementIndex < 2,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        if (subElementIndex == 0) {
            // the bonus column might be {0, 0, 0, 1}
            return LinearAlgebra::SmallMatrix<3, 3>{
                {{{0.5, 0.5, 0.}}, {{0., 1., 0.}}, {{0., 0., 1.}}}};
        } else {
            // the bonus column might be {0, 0, 0, 1}
            return LinearAlgebra::SmallMatrix<3, 3>{
                {{{1., 0., 0.}}, {{0.5, 0.5, 0.}}, {{0., 0., 1.}}}};
        }
    }

    LinearAlgebra::SmallMatrix<3, 3> getRefinementMappingMatrixR(
        std::size_t subElementIndex,
        const PointReference<3>& p) const override final {
        logger.assert_debug(
            subElementIndex < 2,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        if (subElementIndex == 0) {
            // the bonus column might be {0, 0, 0, 1}
            return LinearAlgebra::SmallMatrix<3, 3>{
                {{{2., -1., 0.}}, {{0., 1., 0.}}, {{0., 0., 1.}}}};
        } else {
            // the bonus column might be {0, 0, 0, 1}
            return LinearAlgebra::SmallMatrix<3, 3>{
                {{{1., 0., 0.}}, {{-1., 2., 0.}}, {{0., 0., 1.}}}};
        }
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
            return std::vector<std::size_t>{0, 6, 2, 3, 7, 5};
        } else {
            return std::vector<std::size_t>{0, 1, 6, 3, 4, 7};
        }
    }

    Geometry::ReferenceGeometry* getBigElementReferenceGeometry()
        const override final {
        return &Geometry::ReferenceTriangularPrism::Instance();
    }

    Geometry::ReferenceGeometry* getSubElementReferenceGeometry(
        std::size_t subElement) const override final {
        logger.assert_debug(
            subElement < getNumberOfSubElements(),
            "asked for subElement %, but the % has only % subElements",
            subElement, getName(), getNumberOfSubElements());
        return &Geometry::ReferenceTriangularPrism::Instance();
    }

    std::vector<Geometry::PointReference<3>> getNewNodeLocations(
        const PointReference<3>&) const override final {
        return {{0.5, 0.5, -1.}, {0.5, 0.5, 1.}};
    }

    std::vector<const RefinementMapping*> getCodim1RefinementMaps()
        const override final {
        return std::vector<const RefinementMapping*>{
            RefinementMapForTriangle3::instance(),
            RefinementMapForTriangle3::instance(),
            RefinementMapForSquare0::instance(),
            RefinementMapForSquare0::instance(),
            RefinementMapForSquare1::instance()};
    }

    std::vector<std::size_t> getCodim1LocalNodeIndices(
        std::size_t localFaceNumber) const override final {
        switch (localFaceNumber) {
            case 0:
                return {0, 2, 1, 6};
            case 1:
                return {3, 4, 5, 7};
            case 2:
                return {2, 0, 5, 3};
            case 3:
                return {0, 1, 3, 4};
            case 4:
                return {1, 2, 4, 5, 6, 7};
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
                return std::make_tuple(subFaceIndex, face);
            case 2:
                return std::make_tuple(0, face);
            case 3:
                return std::make_tuple(1, face);
            case 4:
                return std::make_tuple(subFaceIndex, face);
            default:
                logger(ERROR, "something broke, please check the assertions");
        }
        return {};
    }

   private:
    RefinementMapForTriangularPrism3() = default;
    RefinementMapForTriangularPrism3(const RefinementMapForTriangularPrism3&) =
        delete;
};

class RefinementMapForTriangularPrism4 : public RefinementMapping {
   public:
    static const RefinementMapForTriangularPrism4* instance() {
        static RefinementMapForTriangularPrism4 theInstance;
        return &theInstance;
    }

    std::string getName() const override final { return "split in 4"; }

    PointReference<3> refinementTransform(
        std::size_t subElementIndex,
        const PointReference<3>& p) const override final {
        switch (subElementIndex) {
            case 0:
                return {0.5 * p[0], 0.5 * p[1], p[2]};
            case 1:
                return {0.5 + 0.5 * p[0], 0.5 * p[1], p[2]};
            case 2:
                return {0.5 * p[0], 0.5 + 0.5 * p[1], p[2]};
            case 3:
                return {0.5 - 0.5 * p[0], 0.5 - 0.5 * p[1], p[2]};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return p;
        }
    }

    LinearAlgebra::SmallMatrix<3, 3> getRefinementMappingMatrixL(
        std::size_t subElementIndex,
        const PointReference<3>& p) const override final {
        switch (subElementIndex) {
            case 0:
                // the bonus column might be {0, 0, 0, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{0.5, 0., 0.}}, {{0., 0.5, 0.}}, {{0., 0., 1.}}}};
            case 1:
                // the bonus column might be {0.5, 0, 0, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{0.5, 0., 0.}}, {{0., 0.5, 0.}}, {{0., 0., 1.}}}};
            case 2:
                // the bonus column might be {0, 0.5, 0, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{0.5, 0., 0.}}, {{0., 0.5, 0.}}, {{0., 0., 1.}}}};
            case 3:
                // the bonus column might be {0.5, 0.5, 0, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{-0.5, 0., 0.}}, {{0., -0.5, 0.}}, {{0., 0., 1.}}}};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return LinearAlgebra::SmallMatrix<3, 3>();
        }
    }

    LinearAlgebra::SmallMatrix<3, 3> getRefinementMappingMatrixR(
        std::size_t subElementIndex,
        const PointReference<3>& p) const override final {
        switch (subElementIndex) {
            case 0:
                // the bonus column might be {0, 0, 0, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{2., 0., 0.}}, {{0., 2., 0.}}, {{0., 0., 1.}}}};
            case 1:
                // the bonus column might be {-1, 0, 0, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{2., 0., 0.}}, {{0., 2., 0.}}, {{0., 0., 1.}}}};
            case 2:
                // the bonus column might be {0, -1, 0, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{2., 0., 0.}}, {{0., 2., 0.}}, {{0., 0., 1.}}}};
            case 3:
                // the bonus column might be {-1, -1, 0, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{-2., 0., 0.}}, {{0., -2., 0.}}, {{0., 0., 1.}}}};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return LinearAlgebra::SmallMatrix<3, 3>();
        }
    }

    std::size_t getNumberOfNewNodes() const override final { return 6; }

    std::size_t getNumberOfSubElements() const override final { return 4; }

    std::vector<std::size_t> getSubElementLocalNodeIndices(
        std::size_t subElementIndex) const override final {
        switch (subElementIndex) {
            case 0:
                return std::vector<std::size_t>{0, 6, 7, 3, 9, 10};
            case 1:
                return std::vector<std::size_t>{6, 1, 8, 9, 4, 11};
            case 2:
                return std::vector<std::size_t>{7, 8, 2, 10, 11, 5};
            case 3:
                return std::vector<std::size_t>{8, 7, 6, 11, 10, 9};
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
        return &Geometry::ReferenceTriangularPrism::Instance();
    }

    Geometry::ReferenceGeometry* getSubElementReferenceGeometry(
        std::size_t subElement) const override final {
        logger.assert_debug(
            subElement < getNumberOfSubElements(),
            "asked for subElement %, but the % has only % subElements",
            subElement, getName(), getNumberOfSubElements());
        return &Geometry::ReferenceTriangularPrism::Instance();
    }

    std::vector<Geometry::PointReference<3>> getNewNodeLocations(
        const PointReference<3>&) const override final {
        return {{0.5, 0., -1.}, {0., 0.5, -1.}, {0.5, 0.5, -1.},
                {0.5, 0., 1.},  {0., 0.5, 1.},  {0.5, 0.5, 1.}};
    }

    std::vector<const RefinementMapping*> getCodim1RefinementMaps()
        const override final {
        return std::vector<const RefinementMapping*>{
            RefinementMapForTriangle4::instance(),
            RefinementMapForTriangle4::instance(),
            RefinementMapForSquare1::instance(),
            RefinementMapForSquare1::instance(),
            RefinementMapForSquare1::instance()};
    }

    std::vector<std::size_t> getCodim1LocalNodeIndices(
        std::size_t localFaceNumber) const override final {
        switch (localFaceNumber) {
            case 0:
                return {0, 2, 1, 7, 6, 8};
            case 1:
                return {3, 4, 5, 9, 10, 11};
            case 2:
                return {2, 0, 5, 3, 7, 10};
            case 3:
                return {0, 1, 3, 4, 6, 9};
            case 4:
                return {1, 2, 4, 5, 8, 11};
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
                return std::make_tuple(subFaceIndex, face);
            case 2: {
                std::size_t elementIndex;
                std::size_t faceIndex;
                std::tie(elementIndex, faceIndex) =
                    getCodim1RefinementMaps()[0]
                        ->getSubElementAndLocalFaceIndex(1, subFaceIndex);
                return std::make_tuple(elementIndex, faceIndex + 2);
            }
            case 3: {
                std::size_t elementIndex;
                std::size_t faceIndex;
                std::tie(elementIndex, faceIndex) =
                    getCodim1RefinementMaps()[0]
                        ->getSubElementAndLocalFaceIndex(0, subFaceIndex);
                return std::make_tuple(elementIndex, faceIndex + 2);
            }
            case 4: {
                std::size_t elementIndex;
                std::size_t faceIndex;
                std::tie(elementIndex, faceIndex) =
                    getCodim1RefinementMaps()[0]
                        ->getSubElementAndLocalFaceIndex(2, subFaceIndex);
                return std::make_tuple(elementIndex, faceIndex + 2);
            }
            default:
                logger(ERROR, "something broke, please check the assertions");
        }
        return {};
    }

   private:
    RefinementMapForTriangularPrism4() = default;
    RefinementMapForTriangularPrism4(const RefinementMapForTriangularPrism4&) =
        delete;
};

class RefinementMapForTriangularPrism5 : public RefinementMapping {
   public:
    static const RefinementMapForTriangularPrism5* instance() {
        static RefinementMapForTriangularPrism5 theInstance;
        return &theInstance;
    }

    std::string getName() const override final { return "split in 2"; }

    PointReference<3> refinementTransform(
        std::size_t subElementIndex,
        const PointReference<3>& p) const override final {
        switch (subElementIndex) {
            case 0:
                return {p[0], p[1], 0.5 * p[2] - 0.5};
            case 1:
                return {p[0], p[1], 0.5 * p[2] + 0.5};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return p;
        }
    }

    LinearAlgebra::SmallMatrix<3, 3> getRefinementMappingMatrixL(
        std::size_t subElementIndex,
        const PointReference<3>& p) const override final {
        switch (subElementIndex) {
            case 0:
                // the bonus column might be {0, 0, -0.5, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{1., 0., 0.}}, {{0., 1., 0.}}, {{0., 0., 0.5}}}};
            case 1:
                // the bonus column might be {0, 0, 0.5, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{1., 0., 0.}}, {{0., 1., 0.}}, {{0., 0., 0.5}}}};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return LinearAlgebra::SmallMatrix<3, 3>();
        }
    }

    LinearAlgebra::SmallMatrix<3, 3> getRefinementMappingMatrixR(
        std::size_t subElementIndex,
        const PointReference<3>& p) const override final {
        switch (subElementIndex) {
            case 0:
                // the bonus column might be {0, 0, 1, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{1., 0., 0.}}, {{0., 1., 0.}}, {{0., 0., 2.}}}};
            case 1:
                // the bonus column might be {0, 0, -1, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{1., 0., 0.}}, {{0., 1., 0.}}, {{0., 0., 2.}}}};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return LinearAlgebra::SmallMatrix<3, 3>();
        }
    }

    std::size_t getNumberOfNewNodes() const override final { return 3; }

    std::size_t getNumberOfSubElements() const override final { return 2; }

    std::vector<std::size_t> getSubElementLocalNodeIndices(
        std::size_t subElementIndex) const override final {
        switch (subElementIndex) {
            case 0:
                return std::vector<std::size_t>{0, 1, 2, 6, 7, 8};
            case 1:
                return std::vector<std::size_t>{6, 7, 8, 3, 4, 5};
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
        return &Geometry::ReferenceTriangularPrism::Instance();
    }

    Geometry::ReferenceGeometry* getSubElementReferenceGeometry(
        std::size_t subElement) const override final {
        logger.assert_debug(
            subElement < getNumberOfSubElements(),
            "asked for subElement %, but the % has only % subElements",
            subElement, getName(), getNumberOfSubElements());
        return &Geometry::ReferenceTriangularPrism::Instance();
    }

    std::vector<Geometry::PointReference<3>> getNewNodeLocations(
        const PointReference<3>&) const override final {
        return {{0., 0., 0.}, {0., 1., 0.}, {1., 0., 0.}};
    }

    std::vector<const RefinementMapping*> getCodim1RefinementMaps()
        const override final {
        return std::vector<const RefinementMapping*>{
            RefinementMapForTriangle0::instance(),
            RefinementMapForTriangle0::instance(),
            RefinementMapForSquare2::instance(),
            RefinementMapForSquare2::instance(),
            RefinementMapForSquare2::instance()};
    }

    std::vector<std::size_t> getCodim1LocalNodeIndices(
        std::size_t localFaceNumber) const override final {
        switch (localFaceNumber) {
            case 0:
                return {0, 2, 1};
            case 1:
                return {3, 4, 5};
            case 2:
                return {2, 0, 5, 3, 8, 6};
            case 3:
                return {0, 1, 3, 4, 6, 7};
            case 4:
                return {1, 2, 4, 5, 7, 8};
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
                return std::make_tuple(0, face);
            case 1:
                return std::make_tuple(1, face);
            case 2:
                return std::make_tuple(subFaceIndex, face);
            case 3:
                return std::make_tuple(subFaceIndex, face);
            case 4:
                return std::make_tuple(subFaceIndex, face);
            default:
                logger(ERROR, "something broke, please check the assertions");
        }
        return {};
    }

   private:
    RefinementMapForTriangularPrism5() = default;
    RefinementMapForTriangularPrism5(const RefinementMapForTriangularPrism5&) =
        delete;
};

class RefinementMapForTriangularPrism6 : public RefinementMapping {
   public:
    static const RefinementMapForTriangularPrism6* instance() {
        static RefinementMapForTriangularPrism6 theInstance;
        return &theInstance;
    }

    std::string getName() const override final { return "split in 3"; }

    PointReference<3> refinementTransform(
        std::size_t subElementIndex,
        const PointReference<3>& p) const override final {
        switch (subElementIndex) {
            case 0:
                return {(p[1] + 1.) * (5. + p[0]) / 24.,
                        (1. - p[0]) * (5. - p[1]) / 24., p[2]};
            case 1:
                return {((p[0] + 2.) * (2. - p[1]) + 3.) / 12.,
                        (p[1] + 1.) * (5. + p[0]) / 24., p[2]};
            case 2:
                return {(1. - p[0]) * (5. - p[1]) / 24.,
                        ((p[0] + 2.) * (2. - p[1]) + 3.) / 12., p[2]};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return p;
        }
    }

    LinearAlgebra::SmallMatrix<3, 3> getRefinementMappingMatrixL(
        std::size_t subElementIndex,
        const PointReference<3>& p) const override final {
        switch (subElementIndex) {
            case 0:
                // the bonus column should be {5./24., 5./24, 0, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{(p[1] + 1.) / 24., (p[1] - 5.) / 24., 0.}},
                     {{(p[0] + 5.) / 24., (p[0] - 1.) / 24., 0.}},
                     {{0., 0., 1.}}}};
            case 1:
                // the bonus column should be {14./24., 5./24., 0, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{(2. - p[1]) / 12., (p[1] + 1.) / 24., 0.}},
                     {{-(p[0] + 2.) / 12., (p[0] + 5.) / 24., 0.}},
                     {{0., 0., 1.}}}};
            case 2:
                // the bonus column should be {5./24., 14./24., 0, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{(p[1] - 5.) / 24., (2. - p[1]) / 12., 0.}},
                     {{(p[0] - 1.) / 24., -(p[0] + 2.) / 12., 0.}},
                     {{0., 0., 1.}}}};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return LinearAlgebra::SmallMatrix<3, 3>();
        }
    }

    LinearAlgebra::SmallMatrix<3, 3> getRefinementMappingMatrixR(
        std::size_t subElementIndex,
        const PointReference<3>& p) const override final {
        switch (subElementIndex) {
            case 0:
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{(p[0] - 1.) * 4. / (p[0] - p[1] + 4.),
                       (5. - p[1]) * 4. / (p[0] - p[1] + 4.), 0.}},
                     {{-(p[0] + 5.) * 4. / (p[0] - p[1] + 4.),
                       (p[1] + 1.) * 4. / (p[0] - p[1] + 4.), 0.}},
                     {{0., 0., 1.}}}};
            case 1:
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{(p[0] + 5.) * 4. / (p[0] - p[1] + 4.),
                       -(p[1] + 1.) * 4. / (p[0] - p[1] + 4.), 0.}},
                     {{(2. * p[0] + 4.) * 4. / (p[0] - p[1] + 4.),
                       (4. - 2. * p[1]) * 4. / (p[0] - p[1] + 4.), 0.}},
                     {{0., 0., 1.}}}};
            case 2:
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{-(2. * p[0] + 4.) * 4. / (p[0] - p[1] + 4.),
                       (2. * p[1] - 4.) * 4. / (p[0] - p[1] + 4.), 0.}},
                     {{(1. - p[0]) * 4. / (p[0] - p[1] + 4.),
                       (p[1] - 5.) * 4. / (p[0] - p[1] + 4.), 0.}},
                     {{0., 0., 1.}}}};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return LinearAlgebra::SmallMatrix<3, 3>();
        }
    }

    std::size_t getNumberOfNewNodes() const override final { return 8; }

    std::size_t getNumberOfSubElements() const override final { return 3; }

    std::vector<std::size_t> getSubElementLocalNodeIndices(
        std::size_t subElementIndex) const override final {
        switch (subElementIndex) {
            case 0:
                return std::vector<std::size_t>{7, 0, 9, 6, 11, 3, 13, 10};
            case 1:
                return std::vector<std::size_t>{6, 1, 9, 8, 10, 4, 13, 12};
            case 2:
                return std::vector<std::size_t>{8, 2, 9, 7, 12, 5, 13, 11};
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
        return &Geometry::ReferenceTriangularPrism::Instance();
    }

    Geometry::ReferenceGeometry* getSubElementReferenceGeometry(
        std::size_t subElement) const override final {
        logger.assert_debug(
            subElement < getNumberOfSubElements(),
            "asked for subElement %, but the % has only % subElements",
            subElement, getName(), getNumberOfSubElements());
        return &Geometry::ReferenceCube::Instance();
    }

    std::vector<Geometry::PointReference<3>> getNewNodeLocations(
        const PointReference<3>&) const override final {
        return {{0.5, 0., -1.},          {0., 0.5, -1.},        {0.5, 0.5, -1.},
                {1. / 3., 1. / 3., -1.}, {0.5, 0., 1.},         {0., 0.5, 1.},
                {0.5, 0.5, 1.},          {1. / 3., 1. / 3., 1.}};
    }

    std::vector<const RefinementMapping*> getCodim1RefinementMaps()
        const override final {
        return std::vector<const RefinementMapping*>{
            RefinementMapForTriangle5::instance(),
            RefinementMapForTriangle5::instance(),
            RefinementMapForSquare1::instance(),
            RefinementMapForSquare1::instance(),
            RefinementMapForSquare1::instance()};
    }

    std::vector<std::size_t> getCodim1LocalNodeIndices(
        std::size_t localFaceNumber) const override final {
        switch (localFaceNumber) {
            case 0:
                return {0, 2, 1, 7, 6, 8, 9};
            case 1:
                return {3, 4, 5, 10, 11, 12, 13};
            case 2:
                return {2, 0, 5, 3, 7, 11};
            case 3:
                return {0, 1, 3, 4, 6, 10};
            case 4:
                return {1, 2, 4, 5, 8, 12};
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
                return std::make_tuple(subFaceIndex, face);
            case 2: {
                std::size_t elementIndex;
                std::size_t faceIndex;
                std::tie(elementIndex, faceIndex) =
                    getCodim1RefinementMaps()[0]
                        ->getSubElementAndLocalFaceIndex(1, subFaceIndex);
                return std::make_tuple(elementIndex, faceIndex + 2);
            }
            case 3: {
                std::size_t elementIndex;
                std::size_t faceIndex;
                std::tie(elementIndex, faceIndex) =
                    getCodim1RefinementMaps()[0]
                        ->getSubElementAndLocalFaceIndex(0, subFaceIndex);
                return std::make_tuple(elementIndex, faceIndex + 2);
            }
            case 4: {
                std::size_t elementIndex;
                std::size_t faceIndex;
                std::tie(elementIndex, faceIndex) =
                    getCodim1RefinementMaps()[0]
                        ->getSubElementAndLocalFaceIndex(2, subFaceIndex);
                return std::make_tuple(elementIndex, faceIndex + 2);
            }
            default:
                logger(ERROR, "something broke, please check the assertions");
        }
        return {};
    }

   private:
    RefinementMapForTriangularPrism6() = default;
    RefinementMapForTriangularPrism6(const RefinementMapForTriangularPrism6&) =
        delete;
};

class RefinementMapForTriangularPrism7 : public RefinementMapping {
   public:
    static const RefinementMapForTriangularPrism7* instance() {
        static RefinementMapForTriangularPrism7 theInstance;
        return &theInstance;
    }

    std::string getName() const override final {
        return "split in a haxahedron and a triangular prism";
    }

    PointReference<3> refinementTransform(
        std::size_t subElementIndex,
        const PointReference<3>& p) const override final {
        switch (subElementIndex) {
            case 0:
                return {(1. + p[0]) / 4., (p[1] + 1.) * (3. - p[0]) / 8., p[2]};
            case 1:
                return {(1. + p[0]) / 2., p[1] / 2., p[2]};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return p;
        }
    }

    LinearAlgebra::SmallMatrix<3, 3> getRefinementMappingMatrixL(
        std::size_t subElementIndex,
        const PointReference<3>& p) const override final {
        switch (subElementIndex) {
            case 0:
                // the bonus column should be {.25, .375, 0, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{.25, -(p[1] + 1.) / 8., 0.}},
                     {{0., (3. - p[0]) / 8., 0.}},
                     {{0., 0., 1.}}}};
            case 1:
                // the bonus column should be {.5, 0, 0, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{.5, 0., 0.}}, {{0., .5, 0.}}, {{0., 0., 1.}}}};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return LinearAlgebra::SmallMatrix<3, 3>();
        }
    }

    LinearAlgebra::SmallMatrix<3, 3> getRefinementMappingMatrixR(
        std::size_t subElementIndex,
        const PointReference<3>& p) const override final {
        switch (subElementIndex) {
            case 0:
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{4., (p[1] + 1.) * 4. / (3. - p[0]), 0.}},
                     {{0., 8. / (3. - p[0]), 0.}},
                     {{0., 0., 1.}}}};
            case 1:
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{2., 0., 0.}}, {{0., 2., 0.}}, {{0., 0., 1.}}}};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return LinearAlgebra::SmallMatrix<3, 3>();
        }
    }

    std::size_t getNumberOfNewNodes() const override final { return 4; }

    std::size_t getNumberOfSubElements() const override final { return 2; }

    std::vector<std::size_t> getSubElementLocalNodeIndices(
        std::size_t subElementIndex) const override final {
        switch (subElementIndex) {
            case 0:
                return std::vector<std::size_t>{0, 6, 2, 7, 3, 8, 5, 9};
            case 1:
                return std::vector<std::size_t>{6, 1, 7, 8, 4, 9};
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
        return &Geometry::ReferenceTriangularPrism::Instance();
    }

    Geometry::ReferenceGeometry* getSubElementReferenceGeometry(
        std::size_t subElement) const override final {
        switch (subElement) {
            case 0:
                return &Geometry::ReferenceCube::Instance();
            case 1:
                return &Geometry::ReferenceTriangularPrism::Instance();
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElement, getName(), getNumberOfSubElements());
                return nullptr;
        }
    }

    std::vector<Geometry::PointReference<3>> getNewNodeLocations(
        const PointReference<3>&) const override final {
        return {{0.5, 0., -1.}, {0.5, 0.5, -1.}, {0.5, 0., 1.}, {0.5, 0.5, 1.}};
    }

    std::vector<const RefinementMapping*> getCodim1RefinementMaps()
        const override final {
        return std::vector<const RefinementMapping*>{
            RefinementMapForTriangle7::instance(),
            RefinementMapForTriangle6::instance(),
            RefinementMapForSquare0::instance(),
            RefinementMapForSquare1::instance(),
            RefinementMapForSquare1::instance()};
    }

    std::vector<std::size_t> getCodim1LocalNodeIndices(
        std::size_t localFaceNumber) const override final {
        switch (localFaceNumber) {
            case 0:
                return {0, 2, 1, 6, 7};
            case 1:
                return {3, 4, 5, 8, 9};
            case 2:
                return {2, 0, 5, 3};
            case 3:
                return {0, 1, 3, 4, 6, 8};
            case 4:
                return {1, 2, 4, 5, 7, 9};
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
                return std::make_tuple(subFaceIndex, face);
            case 2: {
                std::size_t elementIndex;
                std::size_t faceIndex;
                std::tie(elementIndex, faceIndex) =
                    getCodim1RefinementMaps()[0]
                        ->getSubElementAndLocalFaceIndex(1, subFaceIndex);
                return std::make_tuple(elementIndex, faceIndex + 2);
            }
            case 3: {
                std::size_t elementIndex;
                std::size_t faceIndex;
                std::tie(elementIndex, faceIndex) =
                    getCodim1RefinementMaps()[0]
                        ->getSubElementAndLocalFaceIndex(0, subFaceIndex);
                return std::make_tuple(elementIndex, faceIndex + 2);
            }
            case 4: {
                std::size_t elementIndex;
                std::size_t faceIndex;
                std::tie(elementIndex, faceIndex) =
                    getCodim1RefinementMaps()[0]
                        ->getSubElementAndLocalFaceIndex(2, subFaceIndex);
                return std::make_tuple(elementIndex, faceIndex + 2);
            }
            default:
                logger(ERROR, "something broke, please check the assertions");
        }
        return {};
    }

   private:
    RefinementMapForTriangularPrism7() = default;
    RefinementMapForTriangularPrism7(const RefinementMapForTriangularPrism7&) =
        delete;
};

class RefinementMapForTriangularPrism8 : public RefinementMapping {
   public:
    static const RefinementMapForTriangularPrism8* instance() {
        static RefinementMapForTriangularPrism8 theInstance;
        return &theInstance;
    }

    std::string getName() const override final {
        return "split in a haxahedron and a triangular prism";
    }

    PointReference<3> refinementTransform(
        std::size_t subElementIndex,
        const PointReference<3>& p) const override final {
        switch (subElementIndex) {
            case 0:
                return {(p[0] + 1.) * (3. - p[1]) / 8., (1. + p[1]) / 4., p[2]};
            case 1:
                return {p[0] / 2., (1 + p[1]) / 2., p[2]};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return p;
        }
    }

    LinearAlgebra::SmallMatrix<3, 3> getRefinementMappingMatrixL(
        std::size_t subElementIndex,
        const PointReference<3>& p) const override final {
        switch (subElementIndex) {
            case 0:
                // the bonus column should be {.375, .25, 0, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{(3. - p[1]) / 8., 0., 0.}},
                     {{-(p[0] + 1.) / 8., .25, 0.}},
                     {{0., 0., 1.}}}};
            case 1:
                // the bonus column should be {0, .5, 0, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{.5, 0., 0.}}, {{0., .5, 0.}}, {{0., 0., 1.}}}};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return LinearAlgebra::SmallMatrix<3, 3>();
        }
    }

    LinearAlgebra::SmallMatrix<3, 3> getRefinementMappingMatrixR(
        std::size_t subElementIndex,
        const PointReference<3>& p) const override final {
        switch (subElementIndex) {
            case 0:
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{8. / (3. - p[1]), 0., 0.}},
                     {{(p[0] + 1.) * 4. / (3. - p[1]), 4., 0.}},
                     {{0., 0., 1.}}}};
            case 1:
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{2., 0., 0.}}, {{0., 2., 0.}}, {{0., 0., 1.}}}};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return LinearAlgebra::SmallMatrix<3, 3>();
        }
    }

    std::size_t getNumberOfNewNodes() const override final { return 4; }

    std::size_t getNumberOfSubElements() const override final { return 2; }

    std::vector<std::size_t> getSubElementLocalNodeIndices(
        std::size_t subElementIndex) const override final {
        switch (subElementIndex) {
            case 0:
                return std::vector<std::size_t>{0, 1, 6, 7, 3, 4, 8, 9};
            case 1:
                return std::vector<std::size_t>{6, 7, 2, 8, 9, 5};
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
        return &Geometry::ReferenceTriangularPrism::Instance();
    }

    Geometry::ReferenceGeometry* getSubElementReferenceGeometry(
        std::size_t subElement) const override final {
        switch (subElement) {
            case 0:
                return &Geometry::ReferenceCube::Instance();
            case 1:
                return &Geometry::ReferenceTriangularPrism::Instance();
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElement, getName(), getNumberOfSubElements());
                return nullptr;
        }
    }

    std::vector<Geometry::PointReference<3>> getNewNodeLocations(
        const PointReference<3>&) const override final {
        return {{0., 0.5, -1.}, {0.5, 0.5, -1.}, {0., 0.5, 1.}, {0.5, 0.5, 1.}};
    }

    std::vector<const RefinementMapping*> getCodim1RefinementMaps()
        const override final {
        return std::vector<const RefinementMapping*>{
            RefinementMapForTriangle6::instance(),
            RefinementMapForTriangle7::instance(),
            RefinementMapForSquare1::instance(),
            RefinementMapForSquare0::instance(),
            RefinementMapForSquare1::instance()};
    }

    std::vector<std::size_t> getCodim1LocalNodeIndices(
        std::size_t localFaceNumber) const override final {
        switch (localFaceNumber) {
            case 0:
                return {0, 2, 1, 6, 7};
            case 1:
                return {3, 4, 5, 8, 9};
            case 2:
                return {2, 0, 5, 3, 6, 8};
            case 3:
                return {0, 1, 3, 4};
            case 4:
                return {1, 2, 4, 5, 7, 9};
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
                return std::make_tuple(subFaceIndex, face);
            case 2: {
                std::size_t elementIndex;
                std::size_t faceIndex;
                std::tie(elementIndex, faceIndex) =
                    getCodim1RefinementMaps()[0]
                        ->getSubElementAndLocalFaceIndex(1, subFaceIndex);
                return std::make_tuple(elementIndex, faceIndex + 2);
            }
            case 3: {
                std::size_t elementIndex;
                std::size_t faceIndex;
                std::tie(elementIndex, faceIndex) =
                    getCodim1RefinementMaps()[0]
                        ->getSubElementAndLocalFaceIndex(0, subFaceIndex);
                return std::make_tuple(elementIndex, faceIndex + 2);
            }
            case 4: {
                std::size_t elementIndex;
                std::size_t faceIndex;
                std::tie(elementIndex, faceIndex) =
                    getCodim1RefinementMaps()[0]
                        ->getSubElementAndLocalFaceIndex(2, subFaceIndex);
                return std::make_tuple(elementIndex, faceIndex + 2);
            }
            default:
                logger(ERROR, "something broke, please check the assertions");
        }
        return {};
    }

   private:
    RefinementMapForTriangularPrism8() = default;
    RefinementMapForTriangularPrism8(const RefinementMapForTriangularPrism8&) =
        delete;
};

class RefinementMapForTriangularPrism9 : public RefinementMapping {
   public:
    static const RefinementMapForTriangularPrism9* instance() {
        static RefinementMapForTriangularPrism9 theInstance;
        return &theInstance;
    }

    std::string getName() const override final {
        return "split in a haxahedron and a triangular prism";
    }

    PointReference<3> refinementTransform(
        std::size_t subElementIndex,
        const PointReference<3>& p) const override final {
        switch (subElementIndex) {
            case 0:
                return {p[0] / 2., p[1] / 2., p[2]};
            case 1:
                return {(1. - p[1]) * (3. + p[0]) / 8.,
                        (1. + p[1]) * (3 + p[0]) / 8., p[2]};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return p;
        }
    }

    LinearAlgebra::SmallMatrix<3, 3> getRefinementMappingMatrixL(
        std::size_t subElementIndex,
        const PointReference<3>& p) const override final {
        switch (subElementIndex) {
            case 0:
                // the bonus column should be {0, 0, 0, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{.5, 0., 0.}}, {{0., .5, 0.}}, {{0., 0., 1.}}}};
            case 1:
                // the bonus column should be {.375, .375, 0, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{(1. - p[1]) / 8., (1. + p[1]) / 8., 0.}},
                     {{-(3. + p[0]) / 8., (3. + p[0]) / 8., 0.}},
                     {{0., 0., 1.}}}};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return LinearAlgebra::SmallMatrix<3, 3>();
        }
    }

    LinearAlgebra::SmallMatrix<3, 3> getRefinementMappingMatrixR(
        std::size_t subElementIndex,
        const PointReference<3>& p) const override final {
        switch (subElementIndex) {
            case 0:
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{2., 0., 0.}}, {{0., 2., 0.}}, {{0., 0., 1.}}}};
            case 1:
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{4., -(1. + p[1]) * 4. / (p[0] + 3.), 0.}},
                     {{4., (1. - p[1]) * 4. / (p[0] + 3), 0.}},
                     {{0., 0., 1.}}}};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return LinearAlgebra::SmallMatrix<3, 3>();
        }
    }

    std::size_t getNumberOfNewNodes() const override final { return 4; }

    std::size_t getNumberOfSubElements() const override final { return 2; }

    std::vector<std::size_t> getSubElementLocalNodeIndices(
        std::size_t subElementIndex) const override final {
        switch (subElementIndex) {
            case 0:
                return std::vector<std::size_t>{0, 6, 7, 3, 8, 9};
            case 1:
                return std::vector<std::size_t>{6, 1, 7, 2, 8, 4, 9, 5};
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
        return &Geometry::ReferenceTriangularPrism::Instance();
    }

    Geometry::ReferenceGeometry* getSubElementReferenceGeometry(
        std::size_t subElement) const override final {
        switch (subElement) {
            case 0:
                return &Geometry::ReferenceTriangularPrism::Instance();
            case 1:
                return &Geometry::ReferenceCube::Instance();
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElement, getName(), getNumberOfSubElements());
                return nullptr;
        }
    }

    std::vector<Geometry::PointReference<3>> getNewNodeLocations(
        const PointReference<3>&) const override final {
        return {{0.5, 0., -1.}, {0., 0.5, -1.}, {0.5, 0., 1.}, {0., 0.5, 1.}};
    }

    std::vector<const RefinementMapping*> getCodim1RefinementMaps()
        const override final {
        return std::vector<const RefinementMapping*>{
            RefinementMapForTriangle8::instance(),
            RefinementMapForTriangle8::instance(),
            RefinementMapForSquare1::instance(),
            RefinementMapForSquare1::instance(),
            RefinementMapForSquare0::instance()};
    }

    std::vector<std::size_t> getCodim1LocalNodeIndices(
        std::size_t localFaceNumber) const override final {
        switch (localFaceNumber) {
            case 0:
                return {0, 2, 1, 7, 6};
            case 1:
                return {3, 4, 5, 8, 9};
            case 2:
                return {2, 0, 5, 3, 7, 9};
            case 3:
                return {0, 1, 3, 4, 6, 8};
            case 4:
                return {1, 2, 4, 5};
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
                return std::make_tuple(subFaceIndex, face);
            case 2: {
                std::size_t elementIndex;
                std::size_t faceIndex;
                std::tie(elementIndex, faceIndex) =
                    getCodim1RefinementMaps()[0]
                        ->getSubElementAndLocalFaceIndex(1, subFaceIndex);
                return std::make_tuple(elementIndex, faceIndex + 2);
            }
            case 3: {
                std::size_t elementIndex;
                std::size_t faceIndex;
                std::tie(elementIndex, faceIndex) =
                    getCodim1RefinementMaps()[0]
                        ->getSubElementAndLocalFaceIndex(0, subFaceIndex);
                return std::make_tuple(elementIndex, faceIndex + 2);
            }
            case 4: {
                std::size_t elementIndex;
                std::size_t faceIndex;
                std::tie(elementIndex, faceIndex) =
                    getCodim1RefinementMaps()[0]
                        ->getSubElementAndLocalFaceIndex(2, subFaceIndex);
                return std::make_tuple(elementIndex, faceIndex + 2);
            }
            default:
                logger(ERROR, "something broke, please check the assertions");
        }
        return {};
    }

   private:
    RefinementMapForTriangularPrism9() = default;
    RefinementMapForTriangularPrism9(const RefinementMapForTriangularPrism9&) =
        delete;
};
}  // namespace Geometry

#endif /* KERNEL_GEOMETRY_MAPPINGS_REFINEMENTMAPSFORTRIANGULARPRISM_H_ */
