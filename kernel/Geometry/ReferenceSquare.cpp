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

#include "ReferenceSquare.h"
#include "ReferenceLine.h"
#include "Mappings/MappingToRefLineToSquare.h"
#include "Mappings/MappingToRefSquareToSquare.h"
#include "Geometry/PointReference.h"
#include "Geometry/ReferencePoint.h"
#include "LinearAlgebra/MiddleSizeMatrix.h"
#include "Logger.h"

namespace hpgem {

namespace Geometry {
/* Behold the reference square:
 *
 * (-1,+1) 2---3---3 (+1,+1)
 *         |       |
 *         1       2
 *         |       |
 * (-1,-1) 0---0---1 (+1,-1)
 *
 */
std::size_t ReferenceSquare::localNodeIndexes_[4][2] = {
    {0, 1},
    {0, 2},
    {1, 3},
    {2, 3},
};

ReferenceSquare::ReferenceSquare()
    : ReferenceGeometry(ReferenceGeometryType::SQUARE, "ReferenceSquare"),
      referenceGeometryCodim1Ptr_(&ReferenceLine::Instance()),
      points_(4) {
    // See MappingLineToSquare.h for further info.                 Ref.Line
    // Ref.Sqr.Side
    mappingsLineToSquare_[0] =
        &MappingToRefLineToSquare0::Instance();  // x         -> 0:( x,-1)
    mappingsLineToSquare_[1] =
        &MappingToRefLineToSquare1::Instance();  // x         -> 1:(-1, x)
    mappingsLineToSquare_[2] =
        &MappingToRefLineToSquare2::Instance();  // x         -> 2:( 1, x)
    mappingsLineToSquare_[3] =
        &MappingToRefLineToSquare3::Instance();  // x         -> 3:( x, 1)

    // See Mapping SquareToSquare for further info.                      Ref.
    // Sq.    Ref. Sq.
    mappingsSquareToSquare_[0] =
        &MappingToRefSquareToSquare0::Instance();  // (x,y)    -> (x,y)
    mappingsSquareToSquare_[1] =
        &MappingToRefSquareToSquare1::Instance();  // (x,y)    -> (-y,x)
    mappingsSquareToSquare_[2] =
        &MappingToRefSquareToSquare2::Instance();  // (x,y)    -> (-x,-y)
    mappingsSquareToSquare_[3] =
        &MappingToRefSquareToSquare3::Instance();  // (x,y)    -> (y,-x)
    mappingsSquareToSquare_[4] =
        &MappingToRefSquareToSquare4::Instance();  // (x,y)    -> (x,-y)
    mappingsSquareToSquare_[5] =
        &MappingToRefSquareToSquare5::Instance();  // (x,y)    -> (-x,y)
    mappingsSquareToSquare_[6] =
        &MappingToRefSquareToSquare6::Instance();  // (x,y)    -> (-y,-x)
    mappingsSquareToSquare_[7] =
        &MappingToRefSquareToSquare7::Instance();  // (x,y)    -> (y,x)

    // We set the actual coordinates (see top comment for drawing).

    points_[0] = {-1., -1.};
    points_[1] = {1., -1.};
    points_[2] = {-1., 1.};
    points_[3] = {1., 1.};
    center_ = {0., 0.};
}

bool ReferenceSquare::isInternalPoint(const PointReference<2>& p) const {
    logger.assert_debug(p.size() == 2,
                        "The passed reference point has the wrong dimension");
    return ((p[0] >= -1.) && (p[0] <= 1.) && (p[1] >= -1.) && (p[1] <= 1.));
}

std::ostream& operator<<(std::ostream& os, const ReferenceSquare& square) {
    os << square.getName() << "={";
    auto it = square.points_.begin();
    auto end = square.points_.end();

    for (; it != end; ++it) {
        os << (*it) << '\t';
    }
    os << '}' << std::endl;

    return os;
}

// ================================== Codimension 0
// ============================================
std::size_t ReferenceSquare::getCodim0MappingIndex(
    const std::vector<std::size_t>& list1,
    const std::vector<std::size_t>& list2) const {
    if (list1.size() == 4 && list2.size() == 4) {
        if (list1[0] == list2[0]) {
            if (list1[1] == list2[1]) return 0;
            if (list1[3] == list2[3])
                return 7;
            else
                logger(FATAL,
                       "reference square indexes were given in impossible "
                       "order.\n");
        } else if (list1[0] == list2[1]) {
            if (list1[1] == list2[0]) return 5;
            if (list1[1] == list2[3])
                return 3;
            else
                logger(FATAL,
                       "reference square indexes were given in impossible "
                       "order.\n");
        } else if (list1[0] == list2[2]) {
            if (list1[2] == list2[0]) return 4;
            if (list1[2] == list2[3])
                return 1;
            else
                logger(FATAL,
                       "reference square indexes were given in impossible "
                       "order.\n");
        } else {
            if (list1[1] == list2[1]) return 6;
            if (list1[1] == list2[2])
                return 2;
            else
                logger(FATAL,
                       "reference square indexes were given in impossible "
                       "order.\n");
        }
    } else {
        logger(FATAL,
               "Number of nodes of reference square was larger than 4. \n");
    }
    return 0;
}

const MappingReferenceToReference<0>* ReferenceSquare::getCodim0MappingPtr(
    const std::size_t i) const {
    logger.assert_debug((i < 8),
                        "ERROR: Asked for a mappingSquareToSquare larger than "
                        "7. There are only 8.\n");
    return mappingsSquareToSquare_[i];
}
// ================================== Codimension 1
// ============================================
std::vector<std::size_t> ReferenceSquare::getCodim1EntityLocalIndices(
    const std::size_t faceIndex) const {
    logger.assert_debug(
        faceIndex < 4,
        "ERROR: A square has only 4 edges, while edge % is requested",
        faceIndex);
    return std::vector<std::size_t>(localNodeIndexes_[faceIndex],
                                    localNodeIndexes_[faceIndex] + 2);
}

const ReferenceGeometry* ReferenceSquare::getCodim1ReferenceGeometry(
    const std::size_t faceIndex) const {
    logger.assert_debug((faceIndex < 4),
                        "ERROR: Asked for a square face index larger than 3. "
                        "There are only 4 faces in a square.\n");
    return referenceGeometryCodim1Ptr_;
}
const MappingReferenceToReference<1>* ReferenceSquare::getCodim1MappingPtr(
    const std::size_t faceIndex) const {
    logger.assert_debug((faceIndex < 4),
                        "ERROR: Asked for a square point index larger than 3. "
                        "There are only 4 nodes in a square.\n");
    return mappingsLineToSquare_[faceIndex];
}

const ReferenceGeometry* ReferenceSquare::getCodim2ReferenceGeometry(
    const std::size_t) const {
    return &Geometry::ReferencePoint::Instance();
}
}  // namespace Geometry

}  // namespace hpgem
