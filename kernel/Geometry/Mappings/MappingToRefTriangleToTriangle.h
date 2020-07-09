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

#ifndef HPGEM_KERNEL_MAPPINGTOREFTRIANGLETOTRIANGLE_H
#define HPGEM_KERNEL_MAPPINGTOREFTRIANGLETOTRIANGLE_H

#include "MappingReferenceToReference.h"

namespace Geometry {
/*
 * The ordering of the vertex and faces in a triangle:
 *
 *   (0,1) 2
 *         | \
 *         1   2
 *         |     \
 *   (0,0) 0---0---1 (1,0)
 *
 *
 * This implements the linear mappings of a ReferenceTriangle [0,1]^2 onto
 * itself. There are 6 possible mappings.
 *
 *
 *      index 0: (x,y)->(x,y)     ((0)(1)(2) Identity)
 *      index 1: (x,y)->(y,x)     ((0)(12) Rotation with respect to x = y)
 *      index 2: (x,y)->(1-x-y,x) ((012) Counter-clockwise rotation)
 *      index 3: (x,y)->(1-x-y,y) ((01)(2) Something weird)
 *      index 4: (x,y)->(x,1-x-y) ((02)(1) Something weird)
 *      index 5: (x,y)->(y,1-x-y) ((021) Clockwise rotation)
 */

// ~~~ index 0
// ~~~==============================================================================
class MappingToRefTriangleToTriangle0 : public MappingReferenceToReference<0> {
   public:
    static const MappingToRefTriangleToTriangle0& Instance();
    PointReference<2> transform(
        const Geometry::PointReference<2>& p1) const final;
    Jacobian<2, 2> calcJacobian(const Geometry::PointReference<2>&) const final;
    std::size_t getTargetDimension() const final { return 2; }
    MappingToRefTriangleToTriangle0(const MappingToRefTriangleToTriangle0&) =
        delete;
    MappingToRefTriangleToTriangle0& operator=(
        const MappingToRefTriangleToTriangle0&) = delete;

   private:
    MappingToRefTriangleToTriangle0();
    std::map<const PointReference<2>*, const PointReference<2>*>
        transformedCoordinates;
};
// ~~~ index 1
// ~~~==============================================================================
class MappingToRefTriangleToTriangle1 : public MappingReferenceToReference<0> {
   public:
    static const MappingToRefTriangleToTriangle1& Instance();
    PointReference<2> transform(
        const Geometry::PointReference<2>& p1) const final;
    Jacobian<2, 2> calcJacobian(const Geometry::PointReference<2>&) const final;
    std::size_t getTargetDimension() const final { return 2; }
    MappingToRefTriangleToTriangle1(const MappingToRefTriangleToTriangle1&) =
        delete;
    MappingToRefTriangleToTriangle1& operator=(
        const MappingToRefTriangleToTriangle1&) = delete;

   private:
    MappingToRefTriangleToTriangle1();
    std::map<const PointReference<2>*, const PointReference<2>*>
        transformedCoordinates;
};
// ~~~ index 2
// ~~~==============================================================================
class MappingToRefTriangleToTriangle2 : public MappingReferenceToReference<0> {
   public:
    static const MappingToRefTriangleToTriangle2& Instance();
    PointReference<2> transform(
        const Geometry::PointReference<2>& p1) const final;
    Jacobian<2, 2> calcJacobian(const Geometry::PointReference<2>&) const final;
    std::size_t getTargetDimension() const final { return 2; }
    MappingToRefTriangleToTriangle2(const MappingToRefTriangleToTriangle2&) =
        delete;
    MappingToRefTriangleToTriangle2& operator=(
        const MappingToRefTriangleToTriangle2&) = delete;

   private:
    MappingToRefTriangleToTriangle2();
    std::map<const PointReference<2>*, const PointReference<2>*>
        transformedCoordinates;
};
// ~~~ index 3
// ~~~==============================================================================
class MappingToRefTriangleToTriangle3 : public MappingReferenceToReference<0> {
   public:
    static const MappingToRefTriangleToTriangle3& Instance();
    PointReference<2> transform(
        const Geometry::PointReference<2>& p1) const final;
    Jacobian<2, 2> calcJacobian(const Geometry::PointReference<2>&) const final;
    std::size_t getTargetDimension() const final { return 2; }
    MappingToRefTriangleToTriangle3(const MappingToRefTriangleToTriangle3&) =
        delete;
    MappingToRefTriangleToTriangle3& operator=(
        const MappingToRefTriangleToTriangle3&) = delete;

   private:
    MappingToRefTriangleToTriangle3();
    std::map<const PointReference<2>*, const PointReference<2>*>
        transformedCoordinates;
};
// ~~~ index 4
// ~~~==============================================================================
class MappingToRefTriangleToTriangle4 : public MappingReferenceToReference<0> {
   public:
    static const MappingToRefTriangleToTriangle4& Instance();
    PointReference<2> transform(
        const Geometry::PointReference<2>& p1) const final;
    Jacobian<2, 2> calcJacobian(const Geometry::PointReference<2>&) const final;
    std::size_t getTargetDimension() const final { return 2; }
    MappingToRefTriangleToTriangle4(const MappingToRefTriangleToTriangle4&) =
        delete;
    MappingToRefTriangleToTriangle4& operator=(
        const MappingToRefTriangleToTriangle4&) = delete;

   private:
    MappingToRefTriangleToTriangle4();
    std::map<const PointReference<2>*, const PointReference<2>*>
        transformedCoordinates;
};
// ~~~ index 5
// ~~~==============================================================================
class MappingToRefTriangleToTriangle5 : public MappingReferenceToReference<0> {
   public:
    static const MappingToRefTriangleToTriangle5& Instance();
    PointReference<2> transform(
        const Geometry::PointReference<2>& p1) const final;
    Jacobian<2, 2> calcJacobian(const Geometry::PointReference<2>&) const final;
    std::size_t getTargetDimension() const final { return 2; }
    MappingToRefTriangleToTriangle5(const MappingToRefTriangleToTriangle5&) =
        delete;
    MappingToRefTriangleToTriangle5& operator=(
        const MappingToRefTriangleToTriangle5&) = delete;

   private:
    MappingToRefTriangleToTriangle5();
    std::map<const PointReference<2>*, const PointReference<2>*>
        transformedCoordinates;
};
}  // namespace Geometry
#endif // HPGEM_KERNEL_MAPPINGTOREFTRIANGLETOTRIANGLE_H
