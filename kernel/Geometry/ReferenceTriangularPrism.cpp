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
#include "ReferenceTriangularPrism.h"
#include "ReferenceTriangle.h"
#include "ReferenceSquare.h"
#include "ReferenceLine.h"
#include "Geometry/PointReference.h"
#include "Mappings/MappingToRefFaceToTriangularPrism.h"
#include "LinearAlgebra/MiddleSizeMatrix.h"

namespace hpgem {

namespace Geometry {
std::size_t ReferenceTriangularPrism::localNodeIndexes_[5][4] = {
    {0, 2, 1}, {3, 4, 5}, {2, 0, 5, 3}, {0, 1, 3, 4}, {1, 2, 4, 5}};

std::size_t ReferenceTriangularPrism::localNodesOnEdge_[9][2] = {
    {0, 1}, {0, 2}, {1, 2}, {3, 4}, {3, 5}, {4, 5}, {0, 3}, {1, 4}, {2, 5}};

ReferenceTriangularPrism::ReferenceTriangularPrism()
    : ReferenceGeometry(ReferenceGeometryType::TRIANGULARPRISM,
                        "ReferenceTriangularPrism"),
      referenceGeometryCodim1TrianglePtr_(&ReferenceTriangle::Instance()),
      referenceGeometryCodim1SquarePtr_(&ReferenceSquare::Instance()),
      referenceGeometryCodim2Ptr_(&ReferenceLine::Instance()),
      points_(6) {
    points_[0] = {0., 0., -1.};
    points_[1] = {1., 0., -1.};
    points_[2] = {0., 1., -1.};
    points_[3] = {0., 0., 1.};
    points_[4] = {1., 0., 1.};
    points_[5] = {0., 1., 1.};
    center_ = {1. / 3., 1. / 3., 0.};

    mappingsFaceToTriangularPrism_[0] =
        &MappingToRefFaceToTriangularPrism0::Instance();
    mappingsFaceToTriangularPrism_[1] =
        &MappingToRefFaceToTriangularPrism1::Instance();
    mappingsFaceToTriangularPrism_[2] =
        &MappingToRefFaceToTriangularPrism2::Instance();
    mappingsFaceToTriangularPrism_[3] =
        &MappingToRefFaceToTriangularPrism3::Instance();
    mappingsFaceToTriangularPrism_[4] =
        &MappingToRefFaceToTriangularPrism4::Instance();
}

bool ReferenceTriangularPrism::isInternalPoint(
    const PointReference<3>& p) const {
    logger.assert_debug(p.size() == 3,
                        "The dimension of the reference point is incorrect");
    return ((-1. <= p[2]) && (1. >= p[2]) && (p[0] >= 0.) && (p[0] <= 1.) &&
            (p[1] >= 0.) && (p[1] <= 1. - p[0]));
}

std::ostream& operator<<(std::ostream& os,
                         const ReferenceTriangularPrism& prism) {
    os << prism.getName() << " = ( ";
    auto it = prism.points_.begin();
    auto end = prism.points_.end();

    for (; it != end; ++it) {
        os << (*it) << '\t';
    }
    os << ')' << std::endl;

    return os;
}

// ================================== Codimension 0
// ============================================

std::size_t ReferenceTriangularPrism::getCodim0MappingIndex(
    const std::vector<std::size_t>& list1,
    const std::vector<std::size_t>& list2) const {
    logger(FATAL,
           "ReferenceTriangularPrism::getCodim0MappingIndex: T.p to t.p "
           "mappings do not exist.\n");
    return 0;
}

const MappingReferenceToReference<0>*
    ReferenceTriangularPrism::getCodim0MappingPtr(const std::size_t i) const {
    logger(FATAL,
           "ReferenceTetrahedron::getCodim0MappingPtr: T.p to T.p mappings do "
           "not exist.\n");
    return nullptr;
}

// ================================== Codimension 1
// ============================================

std::vector<std::size_t> ReferenceTriangularPrism::getCodim1EntityLocalIndices(
    const std::size_t faceIndex) const {
    if (faceIndex < 2) {
        return std::vector<std::size_t>(localNodeIndexes_[faceIndex],
                                        localNodeIndexes_[faceIndex] + 3);
    }
    if (faceIndex < 5) {
        return std::vector<std::size_t>(localNodeIndexes_[faceIndex],
                                        localNodeIndexes_[faceIndex] + 4);
    } else {
        logger(ERROR,
               "ReferenceTriangularPrism::getCodim1EntityLocalIndices: Index "
               "out of range. T.p has 5 faces.\n");
    }
    std::vector<std::size_t> dummy(1);
    return dummy;
}

const ReferenceGeometry* ReferenceTriangularPrism::getCodim1ReferenceGeometry(
    const std::size_t faceIndex) const {
    if (faceIndex < 2) {
        return referenceGeometryCodim1TrianglePtr_;
    }
    if (faceIndex < 5) {
        return referenceGeometryCodim1SquarePtr_;
    } else {
        logger(ERROR,
               "ReferenceTriangularPrism::getCodim1ReferenceGeometry: Index "
               "out of range. T.p has 5 faces.\n");
    }
    return nullptr;
}

const BoundaryFaceMapping* ReferenceTriangularPrism::getCodim1MappingPtr(
    const std::size_t faceIndex) const {
    logger.assert_debug((faceIndex < 5),
                        "Asked for a square point index larger than 3. There "
                        "are only 4 nodes in a square!.\n");
    return mappingsFaceToTriangularPrism_[faceIndex];
}

// ================================== Codimension 2
// ============================================

std::vector<std::size_t> ReferenceTriangularPrism::getCodim2EntityLocalIndices(
    const std::size_t edgeIndex) const {
    logger.assert_debug((edgeIndex < 9),
                        "ReferenceTriangularPrism::getCodim2EntityLocalIndices "
                        "Index out of range. T.p has only 9 edges.\n");
    return std::vector<std::size_t>(localNodesOnEdge_[edgeIndex],
                                    localNodesOnEdge_[edgeIndex] + 2);
}

const ReferenceGeometry* ReferenceTriangularPrism::getCodim2ReferenceGeometry(
    const std::size_t edgeIndex) const {
    logger.assert_debug((edgeIndex < 9),
                        "ReferenceTriangularPrism::getCodim2ReferenceGeometry "
                        "Index out of range. T.p has only 9 edges.\n");
    return referenceGeometryCodim2Ptr_;
}

const MappingReferenceToReference<2>*
    ReferenceTriangularPrism::getCodim2MappingPtr(
        const std::size_t faceIndex) const {
    logger(FATAL,
           "ReferenceTriangularPrism::getCodim2MappingPtr: Line to TP mappings "
           "do not exist.\n");
    return nullptr;
}

// ================================== Codimension 3
// ============================================

std::vector<std::size_t> ReferenceTriangularPrism::getCodim3EntityLocalIndices(
    const std::size_t nodeIndex) const {
    logger.assert_debug(
        (nodeIndex < 6),
        "ReferenceTriangularPrism::Index out of range. TP has only 6 nodes.\n");
    return std::vector<std::size_t>(1, nodeIndex);
}
}  // namespace Geometry

}  // namespace hpgem
