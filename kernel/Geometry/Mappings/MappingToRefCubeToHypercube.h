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

#ifndef HPGEM_KERNEL_MAPPINGTOREFCUBETOHYPERCUBE_H
#define HPGEM_KERNEL_MAPPINGTOREFCUBETOHYPERCUBE_H

#include "MappingReferenceToReference.h"

namespace Geometry {

// ~~~ index 0
// ~~~==============================================================================
class MappingToRefCubeToHypercube0 : public MappingReferenceToReference<1> {
   public:
    static const MappingToRefCubeToHypercube0& Instance();
    PointReference<4> transform(
        const Geometry::PointReference<3>& p1) const final;
    Jacobian<3, 4> calcJacobian(const Geometry::PointReference<3>&) const final;
    std::size_t getTargetDimension() const final { return 4; }
    MappingToRefCubeToHypercube0(const MappingToRefCubeToHypercube0&) = delete;
    MappingToRefCubeToHypercube0& operator=(
        const MappingToRefCubeToHypercube0&) = delete;

   private:
    MappingToRefCubeToHypercube0();
    std::map<const PointReference<3>*, const PointReference<4>*>
        transformedCoordinates;
};

// ~~~ index 1
// ~~~==============================================================================

class MappingToRefCubeToHypercube1 : public MappingReferenceToReference<1> {
   public:
    static const MappingToRefCubeToHypercube1& Instance();
    PointReference<4> transform(
        const Geometry::PointReference<3>& p1) const final;
    Jacobian<3, 4> calcJacobian(const Geometry::PointReference<3>&) const final;
    std::size_t getTargetDimension() const final { return 4; }
    MappingToRefCubeToHypercube1(const MappingToRefCubeToHypercube1&) = delete;
    MappingToRefCubeToHypercube1& operator=(
        const MappingToRefCubeToHypercube1&) = delete;

   private:
    MappingToRefCubeToHypercube1();
    std::map<const PointReference<3>*, const PointReference<4>*>
        transformedCoordinates;
};

// ~~~ index 2
// ~~~==============================================================================

class MappingToRefCubeToHypercube2 : public MappingReferenceToReference<1> {
   public:
    static const MappingToRefCubeToHypercube2& Instance();
    PointReference<4> transform(
        const Geometry::PointReference<3>& p1) const final;
    Jacobian<3, 4> calcJacobian(const Geometry::PointReference<3>&) const final;
    std::size_t getTargetDimension() const final { return 4; }
    MappingToRefCubeToHypercube2(const MappingToRefCubeToHypercube2&) = delete;
    MappingToRefCubeToHypercube1& operator=(
        const MappingToRefCubeToHypercube2&) = delete;

   private:
    MappingToRefCubeToHypercube2();
    std::map<const PointReference<3>*, const PointReference<4>*>
        transformedCoordinates;
};

// ~~~ index 3
// ~~~==============================================================================

class MappingToRefCubeToHypercube3 : public MappingReferenceToReference<1> {
   public:
    static const MappingToRefCubeToHypercube3& Instance();
    PointReference<4> transform(
        const Geometry::PointReference<3>& p1) const final;
    Jacobian<3, 4> calcJacobian(const Geometry::PointReference<3>&) const final;
    std::size_t getTargetDimension() const final { return 4; }
    MappingToRefCubeToHypercube3(const MappingToRefCubeToHypercube3&) = delete;
    MappingToRefCubeToHypercube3& operator=(
        const MappingToRefCubeToHypercube3&) = delete;

   private:
    MappingToRefCubeToHypercube3();
    std::map<const PointReference<3>*, const PointReference<4>*>
        transformedCoordinates;
};

// ~~~ index 4
// ~~~==============================================================================

class MappingToRefCubeToHypercube4 : public MappingReferenceToReference<1> {
   public:
    static const MappingToRefCubeToHypercube4& Instance();
    PointReference<4> transform(
        const Geometry::PointReference<3>& p1) const final;
    Jacobian<3, 4> calcJacobian(const Geometry::PointReference<3>&) const final;
    std::size_t getTargetDimension() const final { return 4; }
    MappingToRefCubeToHypercube4(const MappingToRefCubeToHypercube4&) = delete;
    MappingToRefCubeToHypercube4& operator=(
        const MappingToRefCubeToHypercube4&) = delete;

   private:
    MappingToRefCubeToHypercube4();
    std::map<const PointReference<3>*, const PointReference<4>*>
        transformedCoordinates;
};

// ~~~ index 5
// ~~~==============================================================================

class MappingToRefCubeToHypercube5 : public MappingReferenceToReference<1> {
   public:
    static const MappingToRefCubeToHypercube5& Instance();
    PointReference<4> transform(
        const Geometry::PointReference<3>& p1) const final;
    Jacobian<3, 4> calcJacobian(const Geometry::PointReference<3>&) const final;
    std::size_t getTargetDimension() const final { return 4; }
    MappingToRefCubeToHypercube5(const MappingToRefCubeToHypercube5&) = delete;
    MappingToRefCubeToHypercube5& operator=(
        const MappingToRefCubeToHypercube5&) = delete;

   private:
    MappingToRefCubeToHypercube5();
    std::map<const PointReference<3>*, const PointReference<4>*>
        transformedCoordinates;
};

// ~~~ index 6
// ~~~==============================================================================

class MappingToRefCubeToHypercube6 : public MappingReferenceToReference<1> {
   public:
    static const MappingToRefCubeToHypercube6& Instance();
    PointReference<4> transform(
        const Geometry::PointReference<3>& p1) const final;
    Jacobian<3, 4> calcJacobian(const Geometry::PointReference<3>&) const final;
    std::size_t getTargetDimension() const final { return 4; }
    MappingToRefCubeToHypercube6(const MappingToRefCubeToHypercube6&) = delete;
    MappingToRefCubeToHypercube6& operator=(
        const MappingToRefCubeToHypercube6&) = delete;

   private:
    MappingToRefCubeToHypercube6();
    std::map<const PointReference<3>*, const PointReference<4>*>
        transformedCoordinates;
};

// ~~~ index 7
// ~~~==============================================================================

class MappingToRefCubeToHypercube7 : public MappingReferenceToReference<1> {
   public:
    static const MappingToRefCubeToHypercube7& Instance();
    PointReference<4> transform(
        const Geometry::PointReference<3>& p1) const final;
    Jacobian<3, 4> calcJacobian(const Geometry::PointReference<3>&) const final;
    std::size_t getTargetDimension() const final { return 4; }
    MappingToRefCubeToHypercube7(const MappingToRefCubeToHypercube7&) = delete;
    MappingToRefCubeToHypercube7& operator=(
        const MappingToRefCubeToHypercube7&) = delete;

   private:
    MappingToRefCubeToHypercube7();
    std::map<const PointReference<3>*, const PointReference<4>*>
        transformedCoordinates;
};
}  // namespace Geometry
#endif // HPGEM_KERNEL_MAPPINGTOREFCUBETOHYPERCUBE_H
