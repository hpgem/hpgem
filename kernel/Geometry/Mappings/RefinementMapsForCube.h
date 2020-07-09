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

#ifndef HPGEM_KERNEL_REFINEMENTMAPSFORCUBE_H
#define HPGEM_KERNEL_REFINEMENTMAPSFORCUBE_H

#include "RefinementMapping.h"
#include "RefinementMapsForSquare.h"
#include "Geometry/PointReference.h"
#include "Geometry/ReferenceCube.h"

namespace Geometry {
/**
 * Stores all the refinement mappings for a square. Indices for nodes, edges and
 * faces are omitted to prevent excessive cluttering
 *
 *
 *           o---------------o               o-------o-------o
 *          /|              /|              /|      /|      /|
 *        /  |            /  |            /  |    /  |    /  |
 *       /   |           /   |           /   |   /   |   /   |
 *      /    |          /    |          /    |  /    |  /    |
 *    /      |        /      |        /      |/      |/      |
 *   o---------------o       |       o-------o-------o       |
 *   |       |       |       |       |       |       |       |
 *   |       o-------|-------o       |       |-------|-------o
 *   |      /        |      /        |      /|      /|      /
 *   |    /          |    /          |    /  |    /  |    /
 *   |   /       0   |   /           |   /   |   /   |   /
 *   |  /            |  /            |  /  0 |  / 1  |  /
 *   |/              |/              |/      |/      |/
 *   o---------------o               o-------o-------o
 *
 *      (identity)                   (vertical split)
 *         new                        legacy index 0
 *       index 0                         index 1
 *
 *
 *           o---------------o               o---------------o
 *          /|              /|              /|       1      /|
 *        /  |       1    /  |            /  |            /  |
 *       /   |           /   |           o---------------o   |
 *      /    o----------/----o          /|   |          /|   |
 *    /     /|        /     /|        /  |   |        /  |   |
 *   o---------------o    /  |       o---------------o   |   |
 *   |   /   |       |   /   |       |   |   |       |   |   |
 *   |  /    o-------|--/----o       |   |   o-------|---|---o
 *   |/     /        |/     /        |   |  /        |   |  /
 *   o---------------o    /          |   |/          |   |/
 *   |   /           |   /           |   o-----------|---o
 *   |  /       0    |  /            |  /   0        |  /
 *   |/              |/              |/              |/
 *   o---------------o               o---------------o
 *
 *  (horizontal split)                (vertical split)
 *    legacy index 1                   legacy index 2
 *       index 2                          index 3
 *
 *
 *           o-------o-------o               o-------o-------o
 *          /|      /|      /|              /|  2   /|   3  /|
 *        /  |  2 /  | 3  /  |            /  |    /  |    /  |
 *       /   |   /   |   /   |           o-------o-------o   |
 *      /    o--/----|--/----o          /|   |  /|   |  /|   |
 *    /     /|/     /|/     /|        /  |   |/  |   |/  |   |
 *   o-------o-------o    /  |       o-------o-------o   |   |
 *   |   /   |   /   |   /   |       |   |   |   |   |   |   |
 *   |  /    |--/----|--/----o       |   |   |---|---|---|---o
 *   |/     /|/     /|/     /        |   |  /|   |  /|   |  /
 *   o-------o-------o    /          |   |/  |   |/  |   |/
 *   |   /   |   /   |   /           |   o---|---o---|---o
 *   |  /    |  /    |  /            |  /    |  /    |  /
 *   |/   0  |/  1   |/              |/   0  |/   1  |/
 *   o-------o-------o               o-------o-------o
 *
 *     (4way split)                  (2x vertical split)
 *    legacy index 3                   legacy index 4
 *       index 4                          index 5
 *
 *
 *           o---------------o               o-------o-------o
 *          /|              /|              /|  6   /|  7   /|
 *        /  |            /  |            /  |    /  |    /  |
 *       o---------------o 3 |           o-------o-------o   |
 *      /|   o----------/|---o          /|   o--/|---o--/|---o
 *    /  |  /|        /  |  /|        /  |  /|/  |  /|/  |  /|
 *   o---------------o   |/  |       o-------o-------o 3 |/  |
 *   | 1 o-----------|---o 2 |       | 2 o---|---o---|---o 5 |
 *   |  /|   o-------|--/|---o       |  /|   |--/|---|--/|---o
 *   |/  |  /        |/  |  /        |/  |  /|/ 4|  /|/  |  /
 *   o---------------o   |/          o-------o-------o   |/
 *   | 0 o-----------|---o           |   o---|---o---|---o
 *   |  /            |  /            |  /    |  /    |  /
 *   |/              |/              |/  0   |/   1  |/
 *   o---------------o               o-------o-------o
 *
 *     (4way split)                     (8way split)
 *    legacy index 5                   legacy index 6
 *       index 6                          index 7
 *
 *
 *
 * The legacy implementation used to have an additional row and column in its
 * transformation matrix The final row contains a 1 on the diagonal and 0's
 * elsewhere, the final column is mapping-depended and will be provided in
 * comments with the appropriate function. This column seems to be the constant
 * component of the transformation
 *
 *
 * The legacy implementation of refinementTransform also provides the indices
 * 22-25, moreover the refinements 20-25 of the triangular prism appear to be
 * intended for cubes and matching indices appear to belong to matching
 * transformations. The p[2] direction is never refined, so pictures show the
 * front face only. The lines in the pictures are intended to be straight (old
 * comments included) index 20:
 *     // We use this for data transfer between a coarse element, which is
 *     // a hexahedron, and its subelements, which are three triangular-prisms.
 * o-----------o
 *     // The hexahedron itself is a subelement of a coarser triangular-prism |
 * 1  _ -   |
 *     // which element is refined in x_0 direction. | _ -       | element 0: o
 * _     2   | pMap = {2 * p[0] - 1, p[1] - 1, p[2]} |   - _     | element 1: |
 * 0    - _ | pMap = {2 * p[0] - 1, p[0] + p[1], p[2]} o-----------o element 2:
 *         pMap = {1 - 2 * p[0], 1 - p[0] - 2 * p[1], p[2]}
 * index 21:
 *     // We use this for data transfer between a coarse element, which is
 *     // a hexahedron, and its subelements, which are three triangular-prisms.
 * o-----------o
 *     // The hexahedron itself is a subelement of a coarser triangular-prism |\
 * /|
 *     // element which is refined in x_1 direction. | \   2   / | element 0: |
 * \     /  | pMap = {p[0] - 1, 2 * p[1] - 1, p[2]} |   \   /   | element 1: | 0
 * \ / 1  | pMap = {p[0] + p[1], 2 * p[1] - 1, p[2]} o-----o-----o element 2:
 *         pMap = {1 - 2 * p[0] - p[1], 1 - 2 * p[1], p[2]}
 * index 22:
 *     // We use this for data transfer between a coarse element, which is
 *     // a hexahedron, and its triangular/hexahedron subelement. The coarse
 * element o-----------o
 *     // is split into two hexahedrons and one triangular-prism. |     1     |
 *     // The parent hexahedron itself is a subsubelement of a coarser o--- ___
 * |
 *     // triangular-prism element which is refined in x_0 direction. | 2  ___
 * ---o element 0: o---        | pMap = {p[0], (-7 - p[0] + p[1] * (5. - p[0]))
 * / 12, p[2]}                |     0     | element 1: o-----------o pMap =
 * {p[0], (7 + p[0] + p[1] * (5. - p[0])) / 12, p[2]} element 2: pMap = {1 - 2 *
 * p[0], (1 - p[0] - 2 * p[1]) / 3, p[2]} index 23:
 *     // We use this for data transfer between a coarse element, which is
 *     // a hexahedron, and its triangular/hexahedron subelement. The coarse
 * element o---o---o---o
 *     // is split into two hexahedrons and one triangular-prism. |   |   |   |
 *     // The parent hexahedron itself is a subsubelement of a coarser |   | 2 |
 * |
 *     // triangular-prism element which is refined in x_1 direction. | 0  | |
 * 1 | element 0: |    | |    | pMap = {(-7 - p[1] + p[0] * (5. - p[1])) / 12,
 * p[1], p[2]}                |     |     | element 1: o-----o-----o pMap = {(7
 * + p[0] + p[1] * (5. - p[0])) / 12, p[1], p[2]} element 2: pMap = {(1 - 2 *
 * p[0]- p[1]) / 3, 1 - 2 * p[1], p[2]} index 24:
 *     // We use this for data transfer between a coarse element, which is
 *     // a hexahedron, and its triangular/hexahedron subelement. The coarse
 * element o-----------o
 *     // is split into one hexahedron and two triangular-prisms. | 1 __ --   |
 *     // The parent hexahedron itself is a subsubelement of a coarser o - |
 *     // triangular-prism element which is refined in x_0 direction. |      2 |
 *     element 0: o - __      | pMap = {2 * p[0] - 1, (2 * p[1] - 3) / 3, p[2]}
 * | 0    -- _ | element 1: o-----------o pMap = {2 * p[0] - 1, (1 + 2 * p[0] +
 * 2 * p[1]) / 3, p[2]} element 2: (originally 0) pMap = {p[0], p[1] * (2 +
 * p[0]) / 3, p[2]} index 25:
 *     // We use this for data transfer between a coarse element, which is
 *     // a hexahedron, and its triangular/hexahedron subelement. The coarse
 * element o-----------o
 *     // is split into one hexahedron and two triangular-prisms. |\         /|
 *     // The parent hexahedron itself is a subsubelement of a coarser | | | |
 *     // triangular-prism element which is refined in x_1 direction. | |   2 |
 * | element 0: |  \     /  | pMap = {(2 * p[0] - 3) / 3, 2 * p[1] - 1, p[2]} |
 * 0 |   | 1 | element 1: o---o---o---o pMap = {(1 + 2 * p[0] + 2 * p[1]) / 3, 2
 * * p[1] - 1, p[2]} element 2: (originally 0) pMap = {p[0] * (2 + p[1]) / 3,
 * p[1], p[2]}
 */

class RefinementMapForCube0 : public RefinementMapping {
   public:
    static const RefinementMapForCube0* instance() {
        static RefinementMapForCube0 theInstance;
        return &theInstance;
    }

    std::string getName() const final { return "Identity map"; }

    PointReference<3> refinementTransform(
        std::size_t subElementIndex, const PointReference<3>& p) const final {
        logger.assert_debug(
            subElementIndex == 0,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        return p;
    }

    LinearAlgebra::SmallMatrix<3, 3> getRefinementMappingMatrixL(
        std::size_t subElementIndex, const PointReference<3>& p) const final {
        logger.assert_debug(
            subElementIndex == 0,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        // the bonus column is most likely {0, 0, 0, 1}
        return LinearAlgebra::SmallMatrix<3, 3>{
            {{{1., 0., 0.}}, {{0., 1., 0.}}, {{0., 0., 1.}}}};
    }

    LinearAlgebra::SmallMatrix<3, 3> getRefinementMappingMatrixR(
        std::size_t subElementIndex, const PointReference<3>& p) const final {
        logger.assert_debug(
            subElementIndex == 0,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        // the bonus column is most likely {0, 0, 0, 1}
        return LinearAlgebra::SmallMatrix<3, 3>{
            {{{1., 0., 0.}}, {{0., 1., 0.}}, {{0., 0., 1.}}}};
    }

    std::size_t getNumberOfNewNodes() const final { return 0; }

    std::size_t getNumberOfSubElements() const final { return 1; }

    std::vector<std::size_t> getSubElementLocalNodeIndices(
        std::size_t subElementIndex) const final {
        logger.assert_debug(
            subElementIndex == 0,
            "asked for subElement %, but the % has only % subElements",
            subElementIndex, getName(), getNumberOfSubElements());
        return std::vector<std::size_t>{0, 1, 2, 3, 4, 5, 6, 7};
    }

    Geometry::ReferenceGeometry* getBigElementReferenceGeometry() const final {
        return &Geometry::ReferenceCube::Instance();
    }

    Geometry::ReferenceGeometry* getSubElementReferenceGeometry(
        std::size_t subElement) const final {
        logger.assert_debug(
            subElement == 0,
            "asked for subElement %, but the % has only % subElements",
            subElement, getName(), getNumberOfSubElements());
        return &Geometry::ReferenceCube::Instance();
    }

    std::vector<PointReference<3>> getNewNodeLocations(
        const PointReference<3>&) const final {
        return {};
    }

    std::vector<const RefinementMapping*> getCodim1RefinementMaps()
        const final {
        return std::vector<const RefinementMapping*>{
            RefinementMapForSquare0::instance(),
            RefinementMapForSquare0::instance(),
            RefinementMapForSquare0::instance(),
            RefinementMapForSquare0::instance(),
            RefinementMapForSquare0::instance(),
            RefinementMapForSquare0::instance()};
    }

    std::vector<std::size_t> getCodim1LocalNodeIndices(
        std::size_t localFaceNumber) const final {
        switch (localFaceNumber) {
            case 0:
                return {0, 1, 2, 3};
            case 1:
                return {0, 1, 4, 5};
            case 2:
                return {0, 2, 4, 6};
            case 3:
                return {1, 3, 5, 7};
            case 4:
                return {2, 3, 6, 7};
            case 5:
                return {4, 5, 6, 7};
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
    RefinementMapForCube0() = default;
    RefinementMapForCube0(const RefinementMapForCube0&) = delete;
};

class RefinementMapForCube1 : public RefinementMapping {
   public:
    static const RefinementMapForCube1* instance() {
        static RefinementMapForCube1 theInstance;
        return &theInstance;
    }

    std::string getName() const final { return "vertical split"; }

    PointReference<3> refinementTransform(
        std::size_t subElementIndex, const PointReference<3>& p) const final {
        switch (subElementIndex) {
            case 0:
                return {(p[0] - 1.) / 2., p[1], p[2]};
            case 1:
                return {(p[0] + 1.) / 2., p[1], p[2]};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return p;
        }
    }

    LinearAlgebra::SmallMatrix<3, 3> getRefinementMappingMatrixL(
        std::size_t subElementIndex, const PointReference<3>& p) const final {
        switch (subElementIndex) {
            case 0:
                // the bonus column should be {-0.5, 0, 0, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{.5, 0., 0.}}, {{0., 1., 0.}}, {{0., 0., 1.}}}};
            case 1:
                // the bonus column should be {0.5, 0, 0, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{.5, 0., 0.}}, {{0., 1., 0.}}, {{0., 0., 1.}}}};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return LinearAlgebra::SmallMatrix<3, 3>();
        }
    }

    LinearAlgebra::SmallMatrix<3, 3> getRefinementMappingMatrixR(
        std::size_t subElementIndex, const PointReference<3>& p) const final {
        switch (subElementIndex) {
            case 0:
                // the bonus column should be {1, 0, 0, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{2., 0., 0.}}, {{0., 1., 0.}}, {{0., 0., 1.}}}};
            case 1:
                // the bonus column should be {-1, 0, 0, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{2., 0., 0.}}, {{0., 1., 0.}}, {{0., 0., 1.}}}};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return LinearAlgebra::SmallMatrix<3, 3>();
        }
    }

    std::size_t getNumberOfNewNodes() const final { return 4; }

    std::size_t getNumberOfSubElements() const final { return 2; }

    std::vector<std::size_t> getSubElementLocalNodeIndices(
        std::size_t subElementIndex) const final {
        switch (subElementIndex) {
            case 0:
                return std::vector<std::size_t>{0, 8, 2, 9, 4, 10, 6, 11};
            case 1:
                return std::vector<std::size_t>{8, 1, 9, 3, 10, 5, 11, 7};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return std::vector<std::size_t>();
        }
    }

    Geometry::ReferenceGeometry* getBigElementReferenceGeometry() const final {
        return &Geometry::ReferenceCube::Instance();
    }

    Geometry::ReferenceGeometry* getSubElementReferenceGeometry(
        std::size_t subElement) const final {
        logger.assert_debug(
            subElement < 2,
            "asked for subElement %, but the % has only % subElements",
            subElement, getName(), getNumberOfSubElements());
        return &Geometry::ReferenceCube::Instance();
    }

    std::vector<PointReference<3>> getNewNodeLocations(
        const PointReference<3>&) const final {
        return {{0., -1., -1.}, {0., 1., -1.}, {0., -1., 1.}, {0., 1., 1.}};
    }

    std::vector<const RefinementMapping*> getCodim1RefinementMaps()
        const final {
        return std::vector<const RefinementMapping*>{
            RefinementMapForSquare1::instance(),
            RefinementMapForSquare1::instance(),
            RefinementMapForSquare0::instance(),
            RefinementMapForSquare0::instance(),
            RefinementMapForSquare1::instance(),
            RefinementMapForSquare1::instance()};
    }

    std::vector<std::size_t> getCodim1LocalNodeIndices(
        std::size_t localFaceNumber) const final {
        switch (localFaceNumber) {
            case 0:
                return {0, 1, 2, 3, 8, 9};
            case 1:
                return {0, 1, 4, 5, 8, 10};
            case 2:
                return {0, 2, 4, 6};
            case 3:
                return {1, 3, 5, 7};
            case 4:
                return {2, 3, 6, 7, 9, 11};
            case 5:
                return {4, 5, 6, 7, 10, 11};
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
                return std::make_tuple(subFaceIndex, face);
            case 2:
                return std::make_tuple(0, face);
            case 3:
                return std::make_tuple(1, face);
            case 4:
                return std::make_tuple(subFaceIndex, face);
            case 5:
                return std::make_tuple(subFaceIndex, face);
            default:
                logger(ERROR, "something broke, please check the assertions");
        }
        return {};
    }

   private:
    RefinementMapForCube1() = default;
    RefinementMapForCube1(const RefinementMapForCube1&) = delete;
};

class RefinementMapForCube2 : public RefinementMapping {
   public:
    static const RefinementMapForCube2* instance() {
        static RefinementMapForCube2 theInstance;
        return &theInstance;
    }

    std::string getName() const final { return "horizontal split"; }

    PointReference<3> refinementTransform(
        std::size_t subElementIndex, const PointReference<3>& p) const final {
        switch (subElementIndex) {
            case 0:
                return {p[0], (p[1] - 1.) / 2., p[2]};
            case 1:
                return {p[0], (p[1] + 1.) / 2., p[2]};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return p;
        }
    }

    LinearAlgebra::SmallMatrix<3, 3> getRefinementMappingMatrixL(
        std::size_t subElementIndex, const PointReference<3>& p) const final {
        switch (subElementIndex) {
            case 0:
                // the bonus column should be {0, -0.5, 0, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{1., 0., 0.}}, {{0., .5, 0.}}, {{0., 0., 1.}}}};
            case 1:
                // the bonus column should be {0, 0.5, 0, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{1., 0., 0.}}, {{0., .5, 0.}}, {{0., 0., 1.}}}};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return LinearAlgebra::SmallMatrix<3, 3>();
        }
    }

    LinearAlgebra::SmallMatrix<3, 3> getRefinementMappingMatrixR(
        std::size_t subElementIndex, const PointReference<3>& p) const final {
        switch (subElementIndex) {
            case 0:
                // the bonus column should be {0, 1, 0, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{1., 0., 0.}}, {{0., 2., 0.}}, {{0., 0., 1.}}}};
            case 1:
                // the bonus column should be {0, -1, 0, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{1., 0., 0.}}, {{0., 2., 0.}}, {{0., 0., 1.}}}};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return LinearAlgebra::SmallMatrix<3, 3>();
        }
    }

    std::size_t getNumberOfNewNodes() const final { return 4; }

    std::size_t getNumberOfSubElements() const final { return 2; }

    std::vector<std::size_t> getSubElementLocalNodeIndices(
        std::size_t subElementIndex) const final {
        switch (subElementIndex) {
            case 0:
                return std::vector<std::size_t>{0, 1, 8, 9, 4, 5, 10, 11};
            case 1:
                return std::vector<std::size_t>{8, 9, 2, 3, 10, 11, 6, 7};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return std::vector<std::size_t>();
        }
    }

    Geometry::ReferenceGeometry* getBigElementReferenceGeometry() const final {
        return &Geometry::ReferenceCube::Instance();
    }

    Geometry::ReferenceGeometry* getSubElementReferenceGeometry(
        std::size_t subElement) const final {
        logger.assert_debug(
            subElement < 2,
            "asked for subElement %, but the % has only % subElements",
            subElement, getName(), getNumberOfSubElements());
        return &Geometry::ReferenceCube::Instance();
    }

    std::vector<PointReference<3>> getNewNodeLocations(
        const PointReference<3>&) const final {
        return {{-1., 0., -1.}, {1., 0., -1.}, {-1., 0., 1.}, {1., 0., 1.}};
    }

    std::vector<const RefinementMapping*> getCodim1RefinementMaps()
        const final {
        return std::vector<const RefinementMapping*>{
            RefinementMapForSquare2::instance(),
            RefinementMapForSquare0::instance(),
            RefinementMapForSquare1::instance(),
            RefinementMapForSquare1::instance(),
            RefinementMapForSquare0::instance(),
            RefinementMapForSquare2::instance()};
    }

    std::vector<std::size_t> getCodim1LocalNodeIndices(
        std::size_t localFaceNumber) const final {
        switch (localFaceNumber) {
            case 0:
                return {0, 1, 2, 3, 8, 9};
            case 1:
                return {0, 1, 4, 5};
            case 2:
                return {0, 2, 4, 6, 8, 10};
            case 3:
                return {1, 3, 5, 7, 9, 11};
            case 4:
                return {2, 3, 6, 7};
            case 5:
                return {4, 5, 6, 7, 10, 11};
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
                return std::make_tuple(subFaceIndex, face);
            case 3:
                return std::make_tuple(subFaceIndex, face);
            case 4:
                return std::make_tuple(1, face);
            case 5:
                return std::make_tuple(subFaceIndex, face);
            default:
                logger(ERROR, "something broke, please check the assertions");
        }
        return {};
    }

   private:
    RefinementMapForCube2() = default;
    RefinementMapForCube2(const RefinementMapForCube2&) = delete;
};

class RefinementMapForCube3 : public RefinementMapping {
   public:
    static const RefinementMapForCube3* instance() {
        static RefinementMapForCube3 theInstance;
        return &theInstance;
    }

    std::string getName() const final { return "2way split"; }

    PointReference<3> refinementTransform(
        std::size_t subElementIndex, const PointReference<3>& p) const final {
        switch (subElementIndex) {
            case 0:
                return {p[0], p[1], (p[2] - 1.) / 2.};
            case 1:
                return {p[0], p[1], (p[2] + 1.) / 2.};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return p;
        }
    }

    LinearAlgebra::SmallMatrix<3, 3> getRefinementMappingMatrixL(
        std::size_t subElementIndex, const PointReference<3>& p) const final {
        switch (subElementIndex) {
            case 0:
                // the bonus column should be {0, 0, -0.5, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{1., 0., 0.}}, {{0., 1., 0.}}, {{0., 0., 0.5}}}};
            case 1:
                // the bonus column should be {0, 0, 0.5, 1}
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
        std::size_t subElementIndex, const PointReference<3>& p) const final {
        switch (subElementIndex) {
            case 0:
                // the bonus column should be {0, 0, 1, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{1., 0., 0.}}, {{0., 1., 0.}}, {{0., 0., 2.}}}};
            case 1:
                // the bonus column should be {0, 0, -1, 1}
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

    std::size_t getNumberOfNewNodes() const final { return 4; }

    std::size_t getNumberOfSubElements() const final { return 2; }

    std::vector<std::size_t> getSubElementLocalNodeIndices(
        std::size_t subElementIndex) const final {
        switch (subElementIndex) {
            case 0:
                return std::vector<std::size_t>{0, 1, 2, 3, 8, 9, 10, 11};
            case 1:
                return std::vector<std::size_t>{8, 9, 10, 11, 4, 5, 6, 7};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return std::vector<std::size_t>();
        }
    }

    Geometry::ReferenceGeometry* getBigElementReferenceGeometry() const final {
        return &Geometry::ReferenceCube::Instance();
    }

    Geometry::ReferenceGeometry* getSubElementReferenceGeometry(
        std::size_t subElement) const final {
        logger.assert_debug(
            subElement < 2,
            "asked for subElement %, but the % has only % subElements",
            subElement, getName(), getNumberOfSubElements());
        return &Geometry::ReferenceCube::Instance();
    }

    std::vector<PointReference<3>> getNewNodeLocations(
        const PointReference<3>&) const final {
        return {{-1., -1., 0.}, {1., -1., 0.}, {-1., 1., 0.}, {1., 1., 0.}};
    }

    std::vector<const RefinementMapping*> getCodim1RefinementMaps()
        const final {
        return std::vector<const RefinementMapping*>{
            RefinementMapForSquare0::instance(),
            RefinementMapForSquare2::instance(),
            RefinementMapForSquare2::instance(),
            RefinementMapForSquare2::instance(),
            RefinementMapForSquare2::instance(),
            RefinementMapForSquare0::instance()};
    }

    std::vector<std::size_t> getCodim1LocalNodeIndices(
        std::size_t localFaceNumber) const final {
        switch (localFaceNumber) {
            case 0:
                return {0, 1, 2, 3};
            case 1:
                return {0, 1, 4, 5, 8, 9};
            case 2:
                return {0, 2, 4, 6, 8, 10};
            case 3:
                return {1, 3, 5, 7, 9, 11};
            case 4:
                return {2, 3, 6, 7, 10, 11};
            case 5:
                return {4, 5, 6, 7};
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
            case 3:
                return std::make_tuple(subFaceIndex, face);
            case 4:
                return std::make_tuple(subFaceIndex, face);
            case 5:
                return std::make_tuple(1, face);
            default:
                logger(ERROR, "something broke, please check the assertions");
        }
        return {};
    }

   private:
    RefinementMapForCube3() = default;
    RefinementMapForCube3(const RefinementMapForCube3&) = delete;
};

class RefinementMapForCube4 : public RefinementMapping {
   public:
    static const RefinementMapForCube4* instance() {
        static RefinementMapForCube4 theInstance;
        return &theInstance;
    }

    std::string getName() const final { return "4way split"; }

    PointReference<3> refinementTransform(
        std::size_t subElementIndex, const PointReference<3>& p) const final {
        switch (subElementIndex) {
            case 0:
                return {(p[0] - 1.) / 2., (p[1] - 1.) / 2., p[2]};
            case 1:
                return {(p[0] + 1.) / 2., (p[1] - 1.) / 2., p[2]};
            case 2:
                return {(p[0] - 1.) / 2., (p[1] + 1.) / 2., p[2]};
            case 3:
                return {(p[0] + 1.) / 2., (p[1] + 1.) / 2., p[2]};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return p;
        }
    }

    LinearAlgebra::SmallMatrix<3, 3> getRefinementMappingMatrixL(
        std::size_t subElementIndex, const PointReference<3>& p) const final {
        switch (subElementIndex) {
            case 0:
                // the bonus column should be {-0.5, -0.5, 0, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{.5, 0., 0.}}, {{0., .5, 0.}}, {{0., 0., 1.}}}};
            case 1:
                // the bonus column should be {0.5, -0.5, 0, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{.5, 0., 0.}}, {{0., .5, 0.}}, {{0., 0., 1.}}}};
            case 2:
                // the bonus column should be {-0.5, 0.5, 0, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{.5, 0., 0.}}, {{0., .5, 0.}}, {{0., 0., 1.}}}};
            case 3:
                // the bonus column should be {0.5, 0.5, 0, 1}
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
        std::size_t subElementIndex, const PointReference<3>& p) const final {
        switch (subElementIndex) {
            case 0:
                // the bonus column should be {1, 1, 0, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{2., 0., 0.}}, {{0., 2., 0.}}, {{0., 0., 1.}}}};
            case 1:
                // the bonus column should be {-1, 1, 0, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{2., 0., 0.}}, {{0., 2., 0.}}, {{0., 0., 1.}}}};
            case 2:
                // the bonus column should be {1, -1, 0, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{2., 0., 0.}}, {{0., 2., 0.}}, {{0., 0., 1.}}}};
            case 3:
                // the bonus column should be {-1, -1, 0, 1}
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

    std::size_t getNumberOfNewNodes() const final { return 10; }

    std::size_t getNumberOfSubElements() const final { return 4; }

    Geometry::ReferenceGeometry* getBigElementReferenceGeometry() const final {
        return &Geometry::ReferenceCube::Instance();
    }

    Geometry::ReferenceGeometry* getSubElementReferenceGeometry(
        std::size_t subElement) const final {
        logger.assert_debug(
            subElement < getNumberOfSubElements(),
            "asked for subElement %, but the % has only % subElements",
            subElement, getName(), getNumberOfSubElements());
        return &Geometry::ReferenceCube::Instance();
    }

    std::vector<std::size_t> getSubElementLocalNodeIndices(
        std::size_t subElementIndex) const final {
        switch (subElementIndex) {
            case 0:
                return std::vector<std::size_t>{0, 8, 10, 12, 4, 13, 15, 17};
            case 1:
                return std::vector<std::size_t>{8, 1, 12, 11, 13, 5, 17, 16};
            case 2:
                return std::vector<std::size_t>{10, 12, 2, 9, 15, 17, 6, 14};
            case 3:
                return std::vector<std::size_t>{12, 11, 9, 3, 17, 16, 14, 7};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return std::vector<std::size_t>();
        }
    }

    std::vector<PointReference<3>> getNewNodeLocations(
        const PointReference<3>&) const final {
        return {{0., -1., -1.}, {0., 1., -1.}, {-1., 0., -1.}, {1., 0., -1.},
                {0., 0., -1.},  {0., -1., 1.}, {0., 1., 1.},   {-1., 0., 1.},
                {1., 0., 1.},   {0., 0., 1.}};
    }

    std::vector<const RefinementMapping*> getCodim1RefinementMaps()
        const final {
        return std::vector<const RefinementMapping*>{
            RefinementMapForSquare3::instance(),
            RefinementMapForSquare1::instance(),
            RefinementMapForSquare1::instance(),
            RefinementMapForSquare1::instance(),
            RefinementMapForSquare1::instance(),
            RefinementMapForSquare3::instance()};
    }

    std::vector<std::size_t> getCodim1LocalNodeIndices(
        std::size_t localFaceNumber) const final {
        switch (localFaceNumber) {
            case 0:
                return {0, 1, 2, 3, 8, 10, 11, 9, 12};
            case 1:
                return {0, 1, 4, 5, 8, 13};
            case 2:
                return {0, 2, 4, 6, 10, 15};
            case 3:
                return {1, 3, 5, 7, 11, 16};
            case 4:
                return {2, 3, 6, 7, 9, 14};
            case 5:
                return {4, 5, 6, 7, 13, 15, 16, 14, 17};
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
                return std::make_tuple(subFaceIndex, face);
            case 2:
                return std::make_tuple(subFaceIndex * 2, face);
            case 3:
                return std::make_tuple(subFaceIndex * 2 + 1, face);
            case 4:
                return std::make_tuple(subFaceIndex + 2, face);
            case 5:
                return std::make_tuple(subFaceIndex, face);
            default:
                logger(ERROR, "something broke, please check the assertions");
        }
        return {};
    }

   private:
    RefinementMapForCube4() = default;
    RefinementMapForCube4(const RefinementMapForCube4&) = delete;
};

class RefinementMapForCube5 : public RefinementMapping {
   public:
    static const RefinementMapForCube5* instance() {
        static RefinementMapForCube5 theInstance;
        return &theInstance;
    }

    std::string getName() const final { return "4way split"; }

    PointReference<3> refinementTransform(
        std::size_t subElementIndex, const PointReference<3>& p) const final {
        switch (subElementIndex) {
            case 0:
                return {(p[0] - 1.) / 2., p[1], (p[2] - 1.) / 2.};
            case 1:
                return {(p[0] + 1.) / 2., p[1], (p[2] - 1.) / 2.};
            case 2:
                return {(p[0] - 1.) / 2., p[1], (p[2] + 1.) / 2.};
            case 3:
                return {(p[0] + 1.) / 2., p[1], (p[2] + 1.) / 2.};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return p;
        }
    }

    LinearAlgebra::SmallMatrix<3, 3> getRefinementMappingMatrixL(
        std::size_t subElementIndex, const PointReference<3>& p) const final {
        switch (subElementIndex) {
            case 0:
                // the bonus column should be {-0.5, 0, -0.5, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{.5, 0., 0.}}, {{0., 1., 0.}}, {{0., 0., .5}}}};
            case 1:
                // the bonus column should be {0.5, 0, -0.5, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{.5, 0., 0.}}, {{0., 1., 0.}}, {{0., 0., .5}}}};
            case 2:
                // the bonus column should be {-0.5, 0, 0.5, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{.5, 0., 0.}}, {{0., 1., 0.}}, {{0., 0., .5}}}};
            case 3:
                // the bonus column should be {0.5, 0, 0.5, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{.5, 0., 0.}}, {{0., 1., 0.}}, {{0., 0., .5}}}};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return LinearAlgebra::SmallMatrix<3, 3>();
        }
    }

    LinearAlgebra::SmallMatrix<3, 3> getRefinementMappingMatrixR(
        std::size_t subElementIndex, const PointReference<3>& p) const final {
        switch (subElementIndex) {
            case 0:
                // the bonus column should be {1, 0, 1, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{2., 0., 0.}}, {{0., 1., 0.}}, {{0., 0., 2.}}}};
            case 1:
                // the bonus column should be {-1, 0, 1, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{2., 0., 0.}}, {{0., 1., 0.}}, {{0., 0., 2.}}}};
            case 2:
                // the bonus column should be {1, 0, -1, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{2., 0., 0.}}, {{0., 1., 0.}}, {{0., 0., 2.}}}};
            case 3:
                // the bonus column should be {-1, 0, -1, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{2., 0., 0.}}, {{0., 1., 0.}}, {{0., 0., 2.}}}};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return LinearAlgebra::SmallMatrix<3, 3>();
        }
    }

    std::size_t getNumberOfNewNodes() const final { return 10; }

    std::size_t getNumberOfSubElements() const final { return 4; }

    std::vector<std::size_t> getSubElementLocalNodeIndices(
        std::size_t subElementIndex) const final {
        switch (subElementIndex) {
            case 0:
                return std::vector<std::size_t>{0, 8, 2, 9, 12, 16, 14, 17};
            case 1:
                return std::vector<std::size_t>{8, 1, 9, 3, 16, 13, 17, 15};
            case 2:
                return std::vector<std::size_t>{12, 16, 14, 17, 4, 10, 6, 11};
            case 3:
                return std::vector<std::size_t>{16, 13, 17, 15, 10, 5, 11, 7};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return std::vector<std::size_t>();
        }
    }

    Geometry::ReferenceGeometry* getBigElementReferenceGeometry() const final {
        return &Geometry::ReferenceCube::Instance();
    }

    Geometry::ReferenceGeometry* getSubElementReferenceGeometry(
        std::size_t subElement) const final {
        logger.assert_debug(
            subElement < getNumberOfSubElements(),
            "asked for subElement %, but the % has only % subElements",
            subElement, getName(), getNumberOfSubElements());
        return &Geometry::ReferenceCube::Instance();
    }

    std::vector<PointReference<3>> getNewNodeLocations(
        const PointReference<3>&) const final {
        return {{0., -1., -1.}, {0., -1., 1.}, {-1., -1., 0.}, {1., -1., 0.},
                {0., -1., 0.},  {0., 1., -1.}, {0., 1., 1.},   {-1., 1., 0.},
                {1., 1., 0.},   {0., 1., 0.}};
    }

    std::vector<const RefinementMapping*> getCodim1RefinementMaps()
        const final {
        return std::vector<const RefinementMapping*>{
            RefinementMapForSquare1::instance(),
            RefinementMapForSquare3::instance(),
            RefinementMapForSquare2::instance(),
            RefinementMapForSquare2::instance(),
            RefinementMapForSquare3::instance(),
            RefinementMapForSquare1::instance()};
    }

    std::vector<std::size_t> getCodim1LocalNodeIndices(
        std::size_t localFaceNumber) const final {
        switch (localFaceNumber) {
            case 0:
                return {0, 1, 2, 3, 8, 13};
            case 1:
                return {0, 1, 4, 5, 8, 10, 11, 9, 12};
            case 2:
                return {0, 2, 4, 6, 10, 15};
            case 3:
                return {1, 3, 5, 7, 11, 16};
            case 4:
                return {2, 3, 6, 7, 9, 14};
            case 5:
                return {4, 5, 6, 7, 13, 15, 16, 14, 17};
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
                return std::make_tuple(subFaceIndex, face);
            case 2:
                return std::make_tuple(subFaceIndex * 2, face);
            case 3:
                return std::make_tuple(subFaceIndex * 2 + 1, face);
            case 4:
                return std::make_tuple(subFaceIndex, face);
            case 5:
                return std::make_tuple(subFaceIndex + 2, face);
            default:
                logger(ERROR, "something broke, please check the assertions");
        }
        return {};
    }

   private:
    RefinementMapForCube5() = default;
    RefinementMapForCube5(const RefinementMapForCube5&) = delete;
};

class RefinementMapForCube6 : public RefinementMapping {
   public:
    static const RefinementMapForCube6* instance() {
        static RefinementMapForCube6 theInstance;
        return &theInstance;
    }

    std::string getName() const final { return "4way split"; }

    PointReference<3> refinementTransform(
        std::size_t subElementIndex, const PointReference<3>& p) const final {
        switch (subElementIndex) {
            case 0:
                return {p[0], (p[1] - 1.) / 2., (p[2] - 1.) / 2.};
            case 1:
                return {p[0], (p[1] + 1.) / 2., (p[2] - 1.) / 2.};
            case 2:
                return {p[0], (p[1] - 1.) / 2., (p[2] + 1.) / 2.};
            case 3:
                return {p[0], (p[1] + 1.) / 2., (p[2] + 1.) / 2.};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return p;
        }
    }

    LinearAlgebra::SmallMatrix<3, 3> getRefinementMappingMatrixL(
        std::size_t subElementIndex, const PointReference<3>& p) const final {
        switch (subElementIndex) {
            case 0:
                // the bonus column should be {0, -0.5, -0.5, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{1., 0., 0.}}, {{0., .5, 0.}}, {{0., 0., .5}}}};
            case 1:
                // the bonus column should be {0, 0.5, -0.5, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{1., 0., 0.}}, {{0., .5, 0.}}, {{0., 0., .5}}}};
            case 2:
                // the bonus column should be {0, -0.5, 0.5, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{1., 0., 0.}}, {{0., .5, 0.}}, {{0., 0., .5}}}};
            case 3:
                // the bonus column should be {0, 0.5, 0.5, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{1., 0., 0.}}, {{0., .5, 0.}}, {{0., 0., .5}}}};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return LinearAlgebra::SmallMatrix<3, 3>();
        }
    }

    LinearAlgebra::SmallMatrix<3, 3> getRefinementMappingMatrixR(
        std::size_t subElementIndex, const PointReference<3>& p) const final {
        switch (subElementIndex) {
            case 0:
                // the bonus column should be {0, 1, 1, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{1., 0., 0.}}, {{0., 2., 0.}}, {{0., 0., 2.}}}};
            case 1:
                // the bonus column should be {0, -1, 1, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{1., 0., 0.}}, {{0., 2., 0.}}, {{0., 0., 2.}}}};
            case 2:
                // the bonus column should be {0, 1, -1, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{1., 0., 0.}}, {{0., 2., 0.}}, {{0., 0., 2.}}}};
            case 3:
                // the bonus column should be {0, -1, -1, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{1., 0., 0.}}, {{0., 2., 0.}}, {{0., 0., 2.}}}};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return LinearAlgebra::SmallMatrix<3, 3>();
        }
    }

    std::size_t getNumberOfNewNodes() const final { return 10; }

    std::size_t getNumberOfSubElements() const final { return 4; }

    std::vector<std::size_t> getSubElementLocalNodeIndices(
        std::size_t subElementIndex) const final {
        switch (subElementIndex) {
            case 0:
                return std::vector<std::size_t>{0, 1, 8, 9, 12, 13, 16, 17};
            case 1:
                return std::vector<std::size_t>{8, 9, 2, 3, 16, 17, 14, 15};
            case 2:
                return std::vector<std::size_t>{12, 13, 16, 17, 4, 5, 10, 11};
            case 3:
                return std::vector<std::size_t>{16, 17, 14, 15, 10, 11, 6, 7};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return std::vector<std::size_t>();
        }
    }

    Geometry::ReferenceGeometry* getBigElementReferenceGeometry() const final {
        return &Geometry::ReferenceCube::Instance();
    }

    Geometry::ReferenceGeometry* getSubElementReferenceGeometry(
        std::size_t subElement) const final {
        logger.assert_debug(
            subElement == 0,
            "asked for subElement %, but the % has only % subElements",
            subElement, getName(), getNumberOfSubElements());
        return &Geometry::ReferenceCube::Instance();
    }

    std::vector<PointReference<3>> getNewNodeLocations(
        const PointReference<3>&) const final {
        return {{-1., 0., -1.}, {-1., 0., 1.}, {-1., -1., 0.}, {-1., 1., 0.},
                {-1., 0., 0.},  {1., 0., -1.}, {1., 0., 1.},   {1., -1., 0.},
                {1., 1., 0.},   {1., 0., 0.}};
    }

    std::vector<const RefinementMapping*> getCodim1RefinementMaps()
        const final {
        return std::vector<const RefinementMapping*>{
            RefinementMapForSquare2::instance(),
            RefinementMapForSquare2::instance(),
            RefinementMapForSquare3::instance(),
            RefinementMapForSquare3::instance(),
            RefinementMapForSquare2::instance(),
            RefinementMapForSquare2::instance()};
    }

    std::vector<std::size_t> getCodim1LocalNodeIndices(
        std::size_t localFaceNumber) const final {
        switch (localFaceNumber) {
            case 0:
                return {0, 1, 2, 3, 8, 13};
            case 1:
                return {0, 1, 4, 5, 10, 15};
            case 2:
                return {0, 2, 4, 6, 8, 10, 11, 9, 12};
            case 3:
                return {1, 3, 5, 7, 13, 15, 16, 14, 17};
            case 4:
                return {2, 3, 6, 7, 11, 16};
            case 5:
                return {4, 5, 6, 7, 9, 14};
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
                return std::make_tuple(subFaceIndex, face);
            case 3:
                return std::make_tuple(subFaceIndex, face);
            case 4:
                return std::make_tuple(subFaceIndex * 2 + 1, face);
            case 5:
                return std::make_tuple(subFaceIndex + 2, face);
            default:
                logger(ERROR, "something broke, please check the assertions");
        }
        return {};
    }

   private:
    RefinementMapForCube6() = default;
    RefinementMapForCube6(const RefinementMapForCube6&) = delete;
};

class RefinementMapForCube7 : public RefinementMapping {
   public:
    static const RefinementMapForCube7* instance() {
        static RefinementMapForCube7 theInstance;
        return &theInstance;
    }

    std::string getName() const final { return "8way split"; }

    PointReference<3> refinementTransform(
        std::size_t subElementIndex, const PointReference<3>& p) const final {
        switch (subElementIndex) {
            case 0:
                return {(p[0] - 1.) / 2., (p[1] - 1.) / 2., (p[2] - 1.) / 2.};
            case 1:
                return {(p[0] + 1.) / 2., (p[1] - 1.) / 2., (p[2] - 1.) / 2.};
            case 2:
                return {(p[0] - 1.) / 2., (p[1] + 1.) / 2., (p[2] - 1.) / 2.};
            case 3:
                return {(p[0] + 1.) / 2., (p[1] + 1.) / 2., (p[2] - 1.) / 2.};
            case 4:
                return {(p[0] - 1.) / 2., (p[1] - 1.) / 2., (p[2] + 1.) / 2.};
            case 5:
                return {(p[0] + 1.) / 2., (p[1] - 1.) / 2., (p[2] + 1.) / 2.};
            case 6:
                return {(p[0] - 1.) / 2., (p[1] + 1.) / 2., (p[2] + 1.) / 2.};
            case 7:
                return {(p[0] + 1.) / 2., (p[1] + 1.) / 2., (p[2] + 1.) / 2.};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return p;
        }
    }

    LinearAlgebra::SmallMatrix<3, 3> getRefinementMappingMatrixL(
        std::size_t subElementIndex, const PointReference<3>& p) const final {
        switch (subElementIndex) {
            case 0:
                // the bonus column should be {-0.5, -0.5, -0.5, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{.5, 0., 0.}}, {{0., .5, 0.}}, {{0., 0., .5}}}};
            case 1:
                // the bonus column should be {0.5, -0.5, -0.5, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{.5, 0., 0.}}, {{0., .5, 0.}}, {{0., 0., .5}}}};
            case 2:
                // the bonus column should be {-0.5, 0.5, -0.5, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{.5, 0., 0.}}, {{0., .5, 0.}}, {{0., 0., .5}}}};
            case 3:
                // the bonus column should be {0.5, 0.5, -0.5, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{.5, 0., 0.}}, {{0., .5, 0.}}, {{0., 0., .5}}}};
            case 4:
                // the bonus column should be {-0.5, -0.5, 0.5, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{.5, 0., 0.}}, {{0., .5, 0.}}, {{0., 0., .5}}}};
            case 5:
                // the bonus column should be {0.5, -0.5, 0.5, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{.5, 0., 0.}}, {{0., .5, 0.}}, {{0., 0., .5}}}};
            case 6:
                // the bonus column should be {-0.5, 0.5, 0.5, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{.5, 0., 0.}}, {{0., .5, 0.}}, {{0., 0., .5}}}};
            case 7:
                // the bonus column should be {0.5, 0.5, 0.5, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{.5, 0., 0.}}, {{0., .5, 0.}}, {{0., 0., .5}}}};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return LinearAlgebra::SmallMatrix<3, 3>();
        }
    }

    LinearAlgebra::SmallMatrix<3, 3> getRefinementMappingMatrixR(
        std::size_t subElementIndex, const PointReference<3>& p) const final {
        switch (subElementIndex) {
            case 0:
                // the bonus column should be {1, 1, 1, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{2., 0., 0.}}, {{0., 2., 0.}}, {{0., 0., 2.}}}};
            case 1:
                // the bonus column should be {-1, 1, 1, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{2., 0., 0.}}, {{0., 2., 0.}}, {{0., 0., 2.}}}};
            case 2:
                // the bonus column should be {1, -1, 1, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{2., 0., 0.}}, {{0., 2., 0.}}, {{0., 0., 2.}}}};
            case 3:
                // the bonus column should be {-1, -1, 1, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{2., 0., 0.}}, {{0., 2., 0.}}, {{0., 0., 2.}}}};
            case 4:
                // the bonus column should be {1, 1, -1, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{2., 0., 0.}}, {{0., 2., 0.}}, {{0., 0., 2.}}}};
            case 5:
                // the bonus column should be {-1, 1, -1, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{2., 0., 0.}}, {{0., 2., 0.}}, {{0., 0., 2.}}}};
            case 6:
                // the bonus column should be {1, -1, -1, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{2., 0., 0.}}, {{0., 2., 0.}}, {{0., 0., 2.}}}};
            case 7:
                // the bonus column should be {-1, -1, -1, 1}
                return LinearAlgebra::SmallMatrix<3, 3>{
                    {{{2., 0., 0.}}, {{0., 2., 0.}}, {{0., 0., 2.}}}};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return LinearAlgebra::SmallMatrix<3, 3>();
        }
    }

    std::size_t getNumberOfNewNodes() const final { return 19; }

    std::size_t getNumberOfSubElements() const final { return 8; }

    std::vector<std::size_t> getSubElementLocalNodeIndices(
        std::size_t subElementIndex) const final {
        switch (subElementIndex) {
            case 0:
                return std::vector<std::size_t>{0, 8, 12, 20, 16, 22, 24, 26};
            case 1:
                return std::vector<std::size_t>{8, 1, 20, 13, 22, 17, 26, 25};
            case 2:
                return std::vector<std::size_t>{12, 20, 2, 9, 24, 26, 18, 23};
            case 3:
                return std::vector<std::size_t>{20, 13, 9, 3, 26, 25, 23, 19};
            case 4:
                return std::vector<std::size_t>{16, 22, 24, 26, 4, 10, 14, 21};
            case 5:
                return std::vector<std::size_t>{22, 17, 26, 25, 10, 5, 21, 15};
            case 6:
                return std::vector<std::size_t>{24, 26, 18, 23, 14, 21, 6, 11};
            case 7:
                return std::vector<std::size_t>{26, 25, 23, 19, 21, 15, 11, 7};
            default:
                logger(
                    ERROR,
                    "asked for subElement %, but the % has only % subElements",
                    subElementIndex, getName(), getNumberOfSubElements());
                return std::vector<std::size_t>();
        }
    }

    Geometry::ReferenceGeometry* getBigElementReferenceGeometry() const final {
        return &Geometry::ReferenceCube::Instance();
    }

    Geometry::ReferenceGeometry* getSubElementReferenceGeometry(
        std::size_t subElement) const final {
        logger.assert_debug(
            subElement < getNumberOfSubElements(),
            "asked for subElement %, but the % has only % subElements",
            subElement, getName(), getNumberOfSubElements());
        return &Geometry::ReferenceCube::Instance();
    }

    std::vector<PointReference<3>> getNewNodeLocations(
        const PointReference<3>&) const final {
        return {{0., -1., -1.}, {0., 1., -1.}, {0., -1., 1.}, {0., 1., 1.},
                {-1., 0., -1.}, {1., 0., -1.}, {-1., 0., 1.}, {1., 0., 1.},
                {-1., -1., 0.}, {1., -1., 0.}, {-1., 1., 0.}, {1., 1., 0.},
                {-1., 0., 0.},  {1., 0., 0.},  {0., -1., 0.}, {0., 1., 0.},
                {0., 0., -1.},  {0., 0., 1.},  {0., 0., 0.}};
    }

    std::vector<const RefinementMapping*> getCodim1RefinementMaps()
        const final {
        return std::vector<const RefinementMapping*>{
            RefinementMapForSquare3::instance(),
            RefinementMapForSquare3::instance(),
            RefinementMapForSquare3::instance(),
            RefinementMapForSquare3::instance(),
            RefinementMapForSquare3::instance(),
            RefinementMapForSquare3::instance()};
    }

    std::vector<std::size_t> getCodim1LocalNodeIndices(
        std::size_t localFaceNumber) const final {
        switch (localFaceNumber) {
            case 0:
                return {0, 1, 2, 3, 8, 12, 13, 9, 24};
            case 1:
                return {0, 1, 4, 5, 8, 16, 17, 10, 22};
            case 2:
                return {0, 2, 4, 6, 16, 12, 14, 18, 20};
            case 3:
                return {1, 3, 5, 7, 17, 13, 15, 19, 21};
            case 4:
                return {2, 3, 6, 7, 18, 9, 11, 19, 23};
            case 5:
                return {4, 5, 6, 7, 10, 14, 15, 11, 25};
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
                return std::make_tuple(subFaceIndex + 2 * (subFaceIndex / 2),
                                       face);
            case 2:
                return std::make_tuple(subFaceIndex * 2, face);
            case 3:
                return std::make_tuple(subFaceIndex * 2 + 1, face);
            case 4:
                return std::make_tuple(
                    subFaceIndex + 2 * (subFaceIndex / 2) + 2, face);
            case 5:
                return std::make_tuple(subFaceIndex + 4, face);
            default:
                logger(ERROR, "something broke, please check the assertions");
        }
        return {};
    }

   private:
    RefinementMapForCube7() = default;
    RefinementMapForCube7(const RefinementMapForCube7&) = delete;
};
}  // namespace Geometry

#endif  // HPGEM_KERNEL_REFINEMENTMAPSFORCUBE_H
