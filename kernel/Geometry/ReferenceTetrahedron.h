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
//
#ifndef ____ReferenceTetrahedron__
#define ____ReferenceTetrahedron__

#include <iostream>

#include "ReferenceSimplex_Impl.h"
#include <vector>

namespace Geometry {
/* The ordering of the vertex and faces in a tetrahedron:
 *
 * 3 o
 *   |\
 *   |  \o 2
 *   |  / \
 *   |/     \
 * 0 o--------o 1
 *
 */
class ReferenceTetrahedron : public ReferenceSimplex<3> {
   public:
    static ReferenceTetrahedron& Instance() {
        static ReferenceTetrahedron theInstance;
        return theInstance;
    }

    ReferenceTetrahedron(const ReferenceTetrahedron& copy) = delete;

    //! (see ReferenceGeometry.h)
    bool isInternalPoint(const PointReference<3>& point) const final;

    // ================================== Codimension 0
    // ========================================

    //! (see MappingCodimensions.h)
    std::size_t getCodim0MappingIndex(
        const std::vector<std::size_t>&,
        const std::vector<std::size_t>&) const final;

    //! (see MappingCodimensions.h)
    const MappingReferenceToReference<0>* getCodim0MappingPtr(
        const std::size_t) const final;

    using MappingCodimensions::getCodim0MappingPtr;

    // ================================== Codimension 1
    // ========================================

    //! (see MappingCodimensions.h)
    std::size_t getNumberOfCodim1Entities() const final { return 4; }

    //! (see MappingCodimensions.h)
    std::vector<std::size_t> getCodim1EntityLocalIndices(
        const std::size_t) const final;

    //! (see MappingCodimensions.h)
    const MappingReferenceToReference<1>* getCodim1MappingPtr(
        const std::size_t) const final;

    //! (see MappingCodimensions.h)
    const ReferenceGeometry* getCodim1ReferenceGeometry(
        const std::size_t) const final;

    // ================================== Codimension 2
    // ========================================

    //! (see MappingCodimensions.h)
    std::size_t getNumberOfCodim2Entities() const final { return 6; }

    //! (see MappingCodimensions.h)
    std::vector<std::size_t> getCodim2EntityLocalIndices(
        const std::size_t) const final;

    //! (see MappingCodimensions.h)
    const MappingReferenceToReference<2>* getCodim2MappingPtr(
        const std::size_t) const final;

    //! (see MappingCodimensions.h)
    const ReferenceGeometry* getCodim2ReferenceGeometry(
        const std::size_t) const final;

    // ================================== Codimension 3
    // ========================================

    //! (see MappingCodimensions.h)
    std::size_t getNumberOfCodim3Entities() const final { return 4; }

    //! (see MappingCodimensions.h)
    std::vector<std::size_t> getCodim3EntityLocalIndices(
        const std::size_t) const final;

   private:
    ReferenceTetrahedron();

    //! Local node indexes contains the numbering of the vertex of the shape,
    //! ordered by faces. See top comment for the corresponding numbering.
    static std::size_t localNodeIndexes_[4][3];
    static std::size_t localNodesOnEdge_[6][2];

    //! Codimension 1 mappings, from a square to a tetrahedron face. (used to
    //! map a coordinate from a face to an element)
    const MappingReferenceToReference<1>* mappingsTriangleToTetrahedron_[4];
    // const MappingReferenceToReference<0>*
    // mappingsTetrahedronToTetrahedron_[1];

    //! Pointer to the Codimension 1 reference geometry.
    ReferenceGeometry* const referenceGeometryCodim1Ptr_;
    ReferenceGeometry* const referenceGeometryCodim2Ptr_;

    //! List of valid quadrature rules for this reference geometry
    std::vector<QuadratureRules::GaussQuadratureRule*> lstGaussQuadratureRules_;
};
}  // namespace Geometry

#endif /* defined(____ReferenceTetrahedron__) */
