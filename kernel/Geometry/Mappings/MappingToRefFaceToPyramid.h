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

#ifndef HPGEM_KERNEL_MAPPINGTOREFFACETOPYRAMID_H
#define HPGEM_KERNEL_MAPPINGTOREFFACETOPYRAMID_H

#include "MappingReferenceToReference.h"

namespace Geometry {
// ~~~ index 0
// ~~~==============================================================================
class MappingToRefFaceToPyramid0 : public MappingReferenceToReference<1> {
   public:
    static const MappingToRefFaceToPyramid0& Instance();
    PointReference<3> transform(
        const Geometry::PointReference<2>& p1) const final;
    Jacobian<2, 3> calcJacobian(const Geometry::PointReference<2>&) const final;
    std::size_t getTargetDimension() const final { return 3; }
    MappingToRefFaceToPyramid0(const MappingToRefFaceToPyramid0&) = delete;
    MappingToRefFaceToPyramid0& operator=(const MappingToRefFaceToPyramid0&) =
        delete;

   private:
    MappingToRefFaceToPyramid0();
    std::map<const PointReference<2>*, const PointReference<3>*>
        transformedCoordinates;
};

// ~~~ index 1
// ~~~==============================================================================

class MappingToRefFaceToPyramid1 : public MappingReferenceToReference<1> {
   public:
    static const MappingToRefFaceToPyramid1& Instance();
    PointReference<3> transform(
        const Geometry::PointReference<2>& p1) const final;
    Jacobian<2, 3> calcJacobian(const Geometry::PointReference<2>&) const final;
    std::size_t getTargetDimension() const final { return 3; }
    MappingToRefFaceToPyramid1(const MappingToRefFaceToPyramid1&) = delete;
    MappingToRefFaceToPyramid1& operator=(const MappingToRefFaceToPyramid1&) =
        delete;

   private:
    MappingToRefFaceToPyramid1();
    std::map<const PointReference<2>*, const PointReference<3>*>
        transformedCoordinates;
};

// ~~~ index 2
// ~~~==============================================================================

class MappingToRefFaceToPyramid2 : public MappingReferenceToReference<1> {
   public:
    static const MappingToRefFaceToPyramid2& Instance();
    PointReference<3> transform(
        const Geometry::PointReference<2>& p1) const final;
    Jacobian<2, 3> calcJacobian(const Geometry::PointReference<2>&) const final;
    std::size_t getTargetDimension() const final { return 3; }
    MappingToRefFaceToPyramid2(const MappingToRefFaceToPyramid2&) = delete;
    MappingToRefFaceToPyramid1& operator=(const MappingToRefFaceToPyramid2&) =
        delete;

   private:
    MappingToRefFaceToPyramid2();
    std::map<const PointReference<2>*, const PointReference<3>*>
        transformedCoordinates;
};

// ~~~ index 3
// ~~~==============================================================================

class MappingToRefFaceToPyramid3 : public MappingReferenceToReference<1> {
   public:
    static const MappingToRefFaceToPyramid3& Instance();
    PointReference<3> transform(
        const Geometry::PointReference<2>& p1) const final;
    Jacobian<2, 3> calcJacobian(const Geometry::PointReference<2>&) const final;
    std::size_t getTargetDimension() const final { return 3; }
    MappingToRefFaceToPyramid3(const MappingToRefFaceToPyramid3&) = delete;
    MappingToRefFaceToPyramid3& operator=(const MappingToRefFaceToPyramid3&) =
        delete;

   private:
    MappingToRefFaceToPyramid3();
    std::map<const PointReference<2>*, const PointReference<3>*>
        transformedCoordinates;
};

// ~~~ index 4
// ~~~==============================================================================

class MappingToRefFaceToPyramid4 : public MappingReferenceToReference<1> {
   public:
    static const MappingToRefFaceToPyramid4& Instance();
    PointReference<3> transform(
        const Geometry::PointReference<2>& p1) const final;
    Jacobian<2, 3> calcJacobian(const Geometry::PointReference<2>&) const final;
    std::size_t getTargetDimension() const final { return 3; }
    MappingToRefFaceToPyramid4(const MappingToRefFaceToPyramid4&) = delete;
    MappingToRefFaceToPyramid4& operator=(const MappingToRefFaceToPyramid4&) =
        delete;

   private:
    MappingToRefFaceToPyramid4();
    std::map<const PointReference<2>*, const PointReference<3>*>
        transformedCoordinates;
};

}  // namespace Geometry
#endif // HPGEM_KERNEL_MAPPINGTOREFFACETOPYRAMID_H
