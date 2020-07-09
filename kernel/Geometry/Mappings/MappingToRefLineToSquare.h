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

#ifndef HPGEM_KERNEL_MAPPINGTOREFLINETOSQUARE_H
#define HPGEM_KERNEL_MAPPINGTOREFLINETOSQUARE_H

#include "MappingReferenceToReference.h"

namespace hpgem {

namespace Geometry {
/*
 * The ordering of the vertex and faces in a square:
 *
 * (-1,+1) 2---3---3 (+1,+1)
 *         |       |
 *         1       2
 *         |       |
 * (-1,-1) 0---0---1 (1,-1)
 *
 * This maps the reference line [-1,1] to the square shown above. The mappings
 * are defined as follows:
 *
 *      faceindex 0: x -> (x,-1)
 *      faceindex 1: x -> (-1,x)
 *      faceindex 2: x -> (1,x)
 *      faceindex 3: x -> (x,1)
 *
 * The mapping can in principle be defined freely, but to simplify finding the
 * right RefFace2RefFaceMapping they are chosen to preserve the ordering of the
 * vertices (Ordering by coordinate). This makes it possible to use the
 * FaceNodeList from the Reference elements to find the RF2RFMapping. This will
 * also eliminate the need for mirror type of RF2RFMappings, since these will
 * not occur when the ordering is preserved (I think).
 */

// ~~~ index 0
// ~~~==============================================================================
class MappingToRefLineToSquare0 : public MappingReferenceToReference<1> {
   public:
    static const MappingToRefLineToSquare0& Instance();
    PointReference<2> transform(
        const Geometry::PointReference<1>& p1) const final;
    Jacobian<1, 2> calcJacobian(const Geometry::PointReference<1>&) const final;
    std::size_t getTargetDimension() const final { return 2; }
    MappingToRefLineToSquare0(const MappingToRefLineToSquare0&);
    MappingToRefLineToSquare0& operator=(const MappingToRefLineToSquare0&);

   private:
    MappingToRefLineToSquare0();
    std::map<const PointReference<1>*, const PointReference<2>*>
        transformedCoordinates;
};
// ~~~ index 1
// ~~~==============================================================================
class MappingToRefLineToSquare1 : public MappingReferenceToReference<1> {
   public:
    static const MappingToRefLineToSquare1& Instance();
    PointReference<2> transform(
        const Geometry::PointReference<1>& p1) const final;
    Jacobian<1, 2> calcJacobian(const Geometry::PointReference<1>&) const final;
    std::size_t getTargetDimension() const final { return 2; }
    MappingToRefLineToSquare1(const MappingToRefLineToSquare1&) = delete;
    MappingToRefLineToSquare1& operator=(const MappingToRefLineToSquare1&) =
        delete;

   private:
    MappingToRefLineToSquare1();
    std::map<const PointReference<1>*, const PointReference<2>*>
        transformedCoordinates;
};
// ~~~ index 2
// ~~~==============================================================================
class MappingToRefLineToSquare2 : public MappingReferenceToReference<1> {
   public:
    static const MappingToRefLineToSquare2& Instance();
    PointReference<2> transform(
        const Geometry::PointReference<1>& p1) const final;
    Jacobian<1, 2> calcJacobian(const Geometry::PointReference<1>&) const final;
    std::size_t getTargetDimension() const final { return 2; }
    MappingToRefLineToSquare2(const MappingToRefLineToSquare2&) = delete;
    MappingToRefLineToSquare1& operator=(const MappingToRefLineToSquare2&) =
        delete;

   private:
    MappingToRefLineToSquare2();
    std::map<const PointReference<1>*, const PointReference<2>*>
        transformedCoordinates;
};
// ~~~ index 3
// ~~~==============================================================================
class MappingToRefLineToSquare3 : public MappingReferenceToReference<1> {
   public:
    static const MappingToRefLineToSquare3& Instance();
    PointReference<2> transform(
        const Geometry::PointReference<1>& p1) const final;
    Jacobian<1, 2> calcJacobian(const Geometry::PointReference<1>&) const final;
    std::size_t getTargetDimension() const final { return 2; }
    MappingToRefLineToSquare3(const MappingToRefLineToSquare3&) = delete;
    MappingToRefLineToSquare3& operator=(const MappingToRefLineToSquare3&) =
        delete;

   private:
    MappingToRefLineToSquare3();
    std::map<const PointReference<1>*, const PointReference<2>*>
        transformedCoordinates;
};
}  // namespace Geometry
}  // namespace hpgem

#endif  // HPGEM_KERNEL_MAPPINGTOREFLINETOSQUARE_H
