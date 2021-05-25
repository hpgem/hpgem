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

#include "ReferenceLine.h"
#include "ReferenceSquare.h"
#include "ReferenceCube.h"
#include "Mappings/MappingReferenceToReference.h"
#include "Mappings/MappingToRefSquareToCube.h"
#include "Mappings/MappingToRefCubeToCube.h"
#include "LinearAlgebra/MiddleSizeMatrix.h"
#include "Geometry/PointReference.h"
#include "ReferenceGeometry.h"

namespace hpgem {

namespace Geometry {

std::size_t ReferenceCube::localNodeIndexes_[6][4] = {
    {0, 1, 2, 3}, {0, 1, 4, 5}, {0, 2, 4, 6},
    {1, 3, 5, 7}, {2, 3, 6, 7}, {4, 5, 6, 7}};

std::size_t ReferenceCube::localNodesOnEdge_[12][2] = {
    {0, 1}, {2, 3}, {4, 5}, {6, 7}, {0, 2}, {1, 3},
    {4, 6}, {5, 7}, {0, 4}, {1, 5}, {2, 6}, {3, 7},
};

ReferenceCube::ReferenceCube()
    : ReferenceGeometry(ReferenceGeometryType::CUBE, "ReferenceCube"),
      referenceGeometryCodim1Ptr_(&ReferenceSquare::Instance()),
      referenceGeometryCodim2Ptr_(&ReferenceLine::Instance()),
      points_(8) {
    points_[0] = {-1., -1., -1.};
    points_[1] = {1., -1., -1.};
    points_[2] = {-1., 1., -1.};
    points_[3] = {1., 1., -1.};
    points_[4] = {-1., -1., 1.};
    points_[5] = {1., -1., 1.};
    points_[6] = {-1., 1., 1.};
    points_[7] = {1., 1., 1.};
    center_ = {0., 0., 0.};

    mappingsSquareToCube_[0] = &MappingToRefSquareToCube0::Instance();
    mappingsSquareToCube_[1] = &MappingToRefSquareToCube1::Instance();
    mappingsSquareToCube_[2] = &MappingToRefSquareToCube2::Instance();
    mappingsSquareToCube_[3] = &MappingToRefSquareToCube3::Instance();
    mappingsSquareToCube_[4] = &MappingToRefSquareToCube4::Instance();
    mappingsSquareToCube_[5] = &MappingToRefSquareToCube5::Instance();

    mappingsCubeToCube_[0] = &MappingToRefCubeToCube0::Instance();
    mappingsCubeToCube_[1] = &MappingToRefCubeToCube1::Instance();
    mappingsCubeToCube_[2] = &MappingToRefCubeToCube2::Instance();
    mappingsCubeToCube_[3] = &MappingToRefCubeToCube3::Instance();
    mappingsCubeToCube_[4] = &MappingToRefCubeToCube4::Instance();
    mappingsCubeToCube_[5] = &MappingToRefCubeToCube5::Instance();
    mappingsCubeToCube_[6] = &MappingToRefCubeToCube6::Instance();
    mappingsCubeToCube_[7] = &MappingToRefCubeToCube7::Instance();
}

bool ReferenceCube::isInternalPoint(const PointReference<3>& p) const {
    logger.assert_debug(p.size() == 3,
                        "Passed a point with the wrong dimension");
    return ((p[0] >= -1.) && (p[0] <= 1.) && (p[1] >= -1.) && (p[1] <= 1.) &&
            (p[2] >= -1.) && (p[2] <= 1.));
}

std::ostream& operator<<(std::ostream& os, const ReferenceCube& cube) {
    os << cube.getName() << " =( ";
    auto it = cube.points_.begin();
    auto end = cube.points_.end();

    for (; it != end; ++it) {
        os << (*it) << '\t';
    }
    os << ')' << std::endl;

    return os;
}

// ================================== Codimension 0
// ============================================

std::size_t ReferenceCube::getCodim0MappingIndex(
    const std::vector<std::size_t>& list1,
    const std::vector<std::size_t>& list2) const {
    if (list1.size() == 8 && list2.size() == 8) {
        if ((list1[0] == list2[0]) && (list1[4] == list2[4])) {
            if ((list1[1] == list2[1])) return 0;

            return 7;
        } else if ((list1[0] == list2[1]) && (list1[4] == list2[5])) {
            if (list1[1] == list2[0]) return 5;

            return 3;
        } else if ((list1[0] == list2[2]) && (list1[4] == list2[6])) {
            if (list1[2] == list2[0]) return 4;

            return 1;
        } else if ((list1[0] == list2[3]) && (list1[4] == list2[7])) {
            if ((list1[1] == list2[1])) return 6;

            return 2;
        }
    } else {
        logger(ERROR,
               "Number of node indexes was different than 8, so this is not a "
               "cube.\n");
    }
    logger(ERROR,
           "in ReferenceCube, we should not get to the end without returning!. "
           "\n");
    return -1UL;
}

const MappingReferenceToReference<0>* ReferenceCube::getCodim0MappingPtr(
    const std::size_t i) const {
    logger.assert_debug((i < 8), "ERROR: Cube50.\n");
    return mappingsCubeToCube_[i];
}

// ================================== Codimension 1
// ============================================

const MappingReferenceToReference<1>* ReferenceCube::getCodim1MappingPtr(
    const std::size_t faceIndex) const {
    logger.assert_debug((faceIndex < 6), "Cube100.\n");
    return mappingsSquareToCube_[faceIndex];
}

const ReferenceGeometry* ReferenceCube::getCodim1ReferenceGeometry(
    const std::size_t e) const {
    logger.assert_debug((e < 8), "Cube150.\n");
    return referenceGeometryCodim1Ptr_;
}

std::vector<std::size_t> ReferenceCube::getCodim1EntityLocalIndices(
    const std::size_t i) const {
    logger.assert_debug((i < 6), "Cube75.\n");
    return std::vector<std::size_t>(localNodeIndexes_[i],
                                    localNodeIndexes_[i] + 4);
}

// ================================== Codimension 2
// ============================================

const MappingReferenceToReference<2>* ReferenceCube::getCodim2MappingPtr(
    const std::size_t lineIndex) const {
    logger.assert_debug(lineIndex < getNumberOfCodim2Entities(),
                        "Asked for line %, but a cube only has % lines",
                        lineIndex, getNumberOfCodim2Entities());
    return nullptr;
}

const ReferenceGeometry* ReferenceCube::getCodim2ReferenceGeometry(
    const std::size_t e) const {
    logger.assert_debug((e < 12), "Cube150.\n");
    return referenceGeometryCodim2Ptr_;
}

std::vector<std::size_t> ReferenceCube::getCodim2EntityLocalIndices(
    const std::size_t i) const {
    logger.assert_debug((i < 12), "Cube200.\n");
    return std::vector<std::size_t>(localNodesOnEdge_[i],
                                    localNodesOnEdge_[i] + 2);
}
}  // namespace Geometry

}  // namespace hpgem
