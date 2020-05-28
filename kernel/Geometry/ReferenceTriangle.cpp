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

#include "ReferenceTriangle.h"
#include "ReferenceLine.h"
#include "Mappings/MappingToRefLineToTriangle.h"
#include "Mappings/MappingToRefTriangleToTriangle.h"
#include "Geometry/PointReference.h"
#include "Geometry/ReferencePoint.h"
#include "Logger.h"

namespace Geometry {
/* The ordering of the vertex and faces in a triangle:
 *
 *   (0,1) 2
 *         | \
 *         1   2
 *         |     \
 *   (0,0) 0---0---1 (1,0)
 *
 */
std::size_t ReferenceTriangle::localNodeIndexes_[3][2] = {
    {0, 1}, {0, 2}, {1, 2}};

ReferenceTriangle::ReferenceTriangle()
    : ReferenceSimplex(ReferenceGeometryType::TRIANGLE, "ReferenceTriangle"),
      referenceGeometryCodim1Ptr_(&ReferenceLine::Instance()) {
    // See MappingLineToTriangle.h for further info. Ref.Line->Ref.Tr.Side
    mappingsLineToTriangle_[0] =
        &MappingToRefLineToTriangle0::Instance();  // x -> 0:((1+x)/2,0)
    mappingsLineToTriangle_[1] =
        &MappingToRefLineToTriangle1::Instance();  // x -> 1:(0,(1+x)/2)
    mappingsLineToTriangle_[2] =
        &MappingToRefLineToTriangle2::Instance();  // x -> 2:((1-x)/2,(1+x)/2)

    // See Mapping TriangleToTriangle for further info. Ref.Tr. Ref. Tr.
    mappingsTriangleToTriangle_[0] =
        &MappingToRefTriangleToTriangle0::Instance();  // (x,y) -> (x,y)
    mappingsTriangleToTriangle_[1] =
        &MappingToRefTriangleToTriangle1::Instance();  // (x,y) -> (-y,x)
    mappingsTriangleToTriangle_[2] =
        &MappingToRefTriangleToTriangle2::Instance();  // (x,y  -> (-x,-y)
    mappingsTriangleToTriangle_[3] =
        &MappingToRefTriangleToTriangle3::Instance();  // (x,y) -> (y,x)
    mappingsTriangleToTriangle_[4] =
        &MappingToRefTriangleToTriangle4::Instance();  // (x,y) -> (x,-y)
    mappingsTriangleToTriangle_[5] =
        &MappingToRefTriangleToTriangle5::Instance();  // (x,y) -> (-x,y)
}

bool ReferenceTriangle::isInternalPoint(const PointReference<2>& p) const {
    logger.assert_debug(p.size() == 2,
                        "The dimension of the reference point is incorrect");
    return ((p[0] >= 0.) && (p[0] <= 1.) && (p[1] >= 0.) &&
            (p[1] <= 1. - p[0]));
}

// ================================== Codimension 0
// ============================================
std::size_t ReferenceTriangle::getCodim0MappingIndex(
    const std::vector<std::size_t>& list1,
    const std::vector<std::size_t>& list2) const {
    if (list1.size() == 3 && list2.size() == 3) {
        if (list1[0] == list2[0]) {
            if (list1[1] == list2[1])
                return 0;  // 0.1.2.
            
                return 1;  // 0.2.1.
        } else {
            if (list1[0] == list2[1]) {
                if (list1[1] == list2[2])
                    return 5;  // 1.2.0.
                
                    return 3;  // 1.0.2.
            } else {
                if (list1[1] == list2[1])
                    return 4;  // 2.1.0.
                
                    return 2;  // 2.0.1.
            }
        }
    } else {
        logger(ERROR,
               "number of nodes of reference triangle was larger than 3.\n");
    }
    return 0;
}

const MappingReferenceToReference<0>* ReferenceTriangle::getCodim0MappingPtr(
    const std::size_t i) const {
    logger.assert_debug((i < 6),
                        "ERROR: Asked for a mappingTriangleToTriangle larger "
                        "than 5. There are only 6.\n");
    return mappingsTriangleToTriangle_[i];
}
// ================================== Codimension 1
// ============================================
std::vector<std::size_t> ReferenceTriangle::getCodim1EntityLocalIndices(
    const std::size_t faceIndex) const {
    logger.assert_debug(
        faceIndex < 3, "A triangle has only 3 edges, while edge % is requested",
        faceIndex);
    return std::vector<std::size_t>(localNodeIndexes_[faceIndex],
                                    localNodeIndexes_[faceIndex] + 2);
}

const ReferenceGeometry* ReferenceTriangle::getCodim1ReferenceGeometry(
    const std::size_t faceIndex) const {
    logger.assert_debug((faceIndex < 3),
                        "ERROR: Asked for a triangle face index larger than 2. "
                        "There are only 3 faces in a triangle.\n");
    return referenceGeometryCodim1Ptr_;
}
const MappingReferenceToReference<1>* ReferenceTriangle::getCodim1MappingPtr(
    const std::size_t faceIndex) const {
    logger.assert_debug((faceIndex < 3),
                        "ERROR: Asked for a triangle point index larger than "
                        "3. There are only 3 nodes in a triangle.\n");
    return mappingsLineToTriangle_[faceIndex];
}

const ReferenceGeometry* ReferenceTriangle::getCodim2ReferenceGeometry(
    const std::size_t) const {
    return &Geometry::ReferencePoint::Instance();
}

}  // namespace Geometry
