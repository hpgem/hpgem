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

#ifndef HPGEM_KERNEL_REFERENCETRIANGULARPRISM_H
#define HPGEM_KERNEL_REFERENCETRIANGULARPRISM_H

#include <iostream>

#include "ReferenceGeometry.h"
#include <vector>

namespace hpgem {

namespace Geometry {

/* The ordering of the vertex and faces in a triangularPrism:
 *
 *     5 o
 *      /|\
 *    /  |  \
 * 2 o   |    \
 *   |\  |3     \ 4
 *   |  \o--------o
 *   |  / \      /
 *   |/     \  /
 * 0 o--------o 1
 *
 * Note that the axes are oriented differently than for the other reference
 * geometries
 */
class ReferenceTriangularPrism : public ReferenceGeometry {
   public:
    static ReferenceTriangularPrism& Instance() {
        static ReferenceTriangularPrism theInstance;
        return theInstance;
    }

    ReferenceTriangularPrism(const ReferenceTriangularPrism& copy) = delete;

    //! (see ReferenceGeometry.h)
    bool isInternalPoint(const PointReference<3>& point) const final;

    /// Output routine.
    friend std::ostream& operator<<(std::ostream& os,
                                    const ReferenceTriangularPrism& point);

    const PointReferenceBase& getCenter() const final { return center_; }

    std::size_t getNumberOfNodes() const final { return 6; }

    const PointReferenceBase& getReferenceNodeCoordinate(
        const std::size_t& i) const final {
        logger.assert_debug(i < getNumberOfNodes(),
                            "Asked for node %, but there are only % nodes", i,
                            getNumberOfNodes());
        return points_[i];
    }

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
    std::size_t getNumberOfCodim1Entities() const final { return 5; }

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
    std::size_t getNumberOfCodim2Entities() const final { return 9; }

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
    std::size_t getNumberOfCodim3Entities() const final { return 6; }

    //! (see MappingCodimensions.h)
    std::vector<std::size_t> getCodim3EntityLocalIndices(
        const std::size_t) const final;

   private:
    ReferenceTriangularPrism();

    //! Local node indexes contains the numbering of the vertex of the shape,
    //! ordered by faces. See top comment for the corresponding numbering.
    static std::size_t localNodeIndexes_[5][4];
    static std::size_t localNodesOnEdge_[9][2];

    //! Codimension 0 mappings, from triangular prisms. Only used when the face
    //! of an element is a triangular prism
    // const MappingReferenceToReference<0>*
    // mappingsTriangularPrismToTriangularPrism_[1];

    //! Codimension 1 mappings, from a square or triangle to a triangular prism
    //! face.
    const MappingReferenceToReference<1>* mappingsFaceToTriangularPrism_[5];

    //! Pointer to the Codimension 1 reference geometry.
    ReferenceGeometry* const referenceGeometryCodim1TrianglePtr_;
    ReferenceGeometry* const referenceGeometryCodim1SquarePtr_;
    ReferenceGeometry* const referenceGeometryCodim2Ptr_;

    //! List of valid quadrature rules for this reference geometry
    std::vector<QuadratureRules::GaussQuadratureRule*> lstGaussQuadratureRules_;

    std::vector<PointReference<3> > points_;

    PointReference<3> center_;
};
}  // namespace Geometry

}  // namespace hpgem

#endif  // HPGEM_KERNEL_REFERENCETRIANGULARPRISM_H
