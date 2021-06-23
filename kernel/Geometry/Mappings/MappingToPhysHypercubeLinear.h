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

#ifndef HPGEM_KERNEL_MAPPINGTOPHYSHYPERCUBELINEAR_H
#define HPGEM_KERNEL_MAPPINGTOPHYSHYPERCUBELINEAR_H

#include "MappingReferenceToPhysical.h"
#include "Geometry/PointPhysical.h"

namespace hpgem {

namespace Geometry {
/*!
 * "In geometry, a hypercube is an n-dimensional analogue of a square (n = 2)
 * and a cube (n = 3)" -Wikipedia.
 *
 * This templated class defines the linear mappings between hypercubes, in the
 * corresponding dimension. See the comments in the Physical<Hypercube>.cpp
 * files to know the order of the vertex of each hypercube, an order which is
 * kept by the mappings.
 */

template <std::size_t DIM>
class MappingToPhysHypercubeLinear;

// ~~~ Dimension 1
// ~~~==========================================================================
template <>
class MappingToPhysHypercubeLinear<1>
    : public MappingReferenceToPhysicalDim<1> {

   public:
    MappingToPhysHypercubeLinear(const PhysicalGeometry<1> *const &);
    MappingToPhysHypercubeLinear(const MappingToPhysHypercubeLinear<1> &other) =
        default;
    PointPhysical<1> transform(const PointReference<1> &) const final;
    PointReference<1> inverseTransform(const PointPhysical<1> &) const final;
    Jacobian<1, 1> calcJacobian(const PointReference<1> &) const final;
    void reinit() final;
    std::size_t getTargetDimension() const final { return 1; }

   private:
    bool isValidPoint(const PointReference<1> &) const;
    double mid, slope;
};

// ~~~ Dimension 2
// ~~~==========================================================================
template <>
class MappingToPhysHypercubeLinear<2>
    : public MappingReferenceToPhysicalDim<2> {
   public:
    // Constructor.
    MappingToPhysHypercubeLinear(const PhysicalGeometry<2> *const &);
    MappingToPhysHypercubeLinear(const MappingToPhysHypercubeLinear<2> &other) =
        default;
    PointPhysical<2> transform(const PointReference<2> &) const final;
    PointReference<2> inverseTransform(const PointPhysical<2> &) const final;
    Jacobian<2, 2> calcJacobian(const PointReference<2> &) const final;
    void reinit() final;
    std::size_t getTargetDimension() const final { return 2; }

   private:
    bool isValidPoint(const PointReference<2> &) const;
    PointPhysical<2> a0, a1, a2, a12;
};

// ~~~ Dimension 3
// ~~~==========================================================================
template <>
class MappingToPhysHypercubeLinear<3>
    : public MappingReferenceToPhysicalDim<3> {
   public:
    // Constructor.
    MappingToPhysHypercubeLinear(const PhysicalGeometry<3> *const &);
    MappingToPhysHypercubeLinear(const MappingToPhysHypercubeLinear<3> &other) =
        default;
    PointPhysical<3> transform(const PointReference<3> &) const final;
    PointReference<3> inverseTransform(const PointPhysical<3> &) const final;
    Jacobian<3, 3> calcJacobian(const PointReference<3> &) const final;
    void reinit() final;
    std::size_t getTargetDimension() const final { return 3; }

   private:
    bool isValidPoint(const PointReference<3> &) const;
    PointPhysical<3> a0, a1, a2, a3, a12, a23, a13, a123;
};

// ~~~ Dimension 4
// ~~~==========================================================================
template <>
class MappingToPhysHypercubeLinear<4>
    : public MappingReferenceToPhysicalDim<4> {
   public:
    // Constructor.
    MappingToPhysHypercubeLinear(const PhysicalGeometry<4> *const &);
    MappingToPhysHypercubeLinear(const MappingToPhysHypercubeLinear<4> &other) =
        default;
    PointPhysical<4> transform(const PointReference<4> &) const final;
    PointReference<4> inverseTransform(const PointPhysical<4> &) const final;
    Jacobian<4, 4> calcJacobian(const PointReference<4> &) const final;
    void reinit() final;
    std::size_t getTargetDimension() const final { return 4; }

   private:
    bool isValidPoint(const PointReference<4> &) const;
    PointPhysical<4> abar, a0, a1, a2, a3, a01, a02, a03, a12, a13, a23, a012,
        a013, a123, a230, a0123;
};
}  // namespace Geometry
}  // namespace hpgem

#endif  // HPGEM_KERNEL_MAPPINGTOPHYSHYPERCUBELINEAR_H
