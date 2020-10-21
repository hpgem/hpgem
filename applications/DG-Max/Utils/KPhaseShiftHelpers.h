/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2020, University of Twente
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
#ifndef HPGEM_KPHASESHIFTHELPERS_H
#define HPGEM_KPHASESHIFTHELPERS_H

#include "Geometry/PointPhysical.h"
#include "Base/Face.h"

namespace DGMax {

using namespace hpgem;


// Several helper functions needed for implemting KPhase shifts

// Given a point x on a (periodic boundary) face, which has coordinates x_l and
// x_r as seen from the left and right face respectively. Compute the difference
// x_l - x_r,

/// \brief Compute the coordinate jump of a face.
///
/// For a point P on a (periodic boundary) the coordinate may be different when
/// computed from the left or right face. This function computes the jump in this
/// coordinate: x_l - x_r, where x_l/x_r is the coordinate of point P as seen
/// from the left/right face.
/// \tparam DIM
/// \param face
/// \return
template <std::size_t DIM>
LinearAlgebra::SmallVector<DIM> boundaryFaceShift(const Base::Face* face) {
    logger.assert_always(face->isInternal(), "Internal face boundary");
    const Geometry::PointReference<DIM - 1>& p =
        face->getReferenceGeometry()->getCenter();
    const Geometry::PointPhysical<DIM> pLeftPhys =
        face->getPtrElementLeft()->referenceToPhysical(
            face->mapRefFaceToRefElemL(p));
    const Geometry::PointPhysical<DIM> pRightPhys =
        face->getPtrElementRight()->referenceToPhysical(
            face->mapRefFaceToRefElemR(p));
    return pLeftPhys.getCoordinates() - pRightPhys.getCoordinates();
}

///  Compute the physical coordinate of a node as seen from a specific element.
template <std::size_t DIM>
Geometry::PointPhysical<DIM> getCoordinate(const Base::Element* element,
                                           const Base::Node* node) {
    return element->getPhysicalGeometry()->getLocalNodeCoordinates(
        element->getLocalId(node));
}

/// Compute the physical coordinate of the middle of an edge as seen from a
/// specific element.
template <std::size_t DIM>
Geometry::PointPhysical<DIM> getCoordinate(const Base::Element* element,
                                           const Base::Edge* edge) {
    std::size_t localEdgeId = element->getLocalId(edge);
    std::vector<std::size_t> nodeIds =
        element->getReferenceGeometry()->getCodim2EntityLocalIndices(
            localEdgeId);
    const Geometry::PointPhysical<DIM> n1 =
        element->getPhysicalGeometry()->getLocalNodeCoordinates(nodeIds[0]);
    const Geometry::PointPhysical<DIM> n2 =
        element->getPhysicalGeometry()->getLocalNodeCoordinates(nodeIds[1]);
    return 0.5 * (n1 + n2);
}

/// Compute the physical coordinate of the center of a face as seen from a
/// specific element.
template <std::size_t DIM>
Geometry::PointPhysical<DIM> getCoordinate(const Base::Element* element,
                                           const Base::Face* face) {
    bool isLeft = face->getPtrElementLeft() == element;
    logger.assert_debug(isLeft || face->getPtrElementRight() == element,
                        "Not a neighbouring element");
    const Geometry::PointReference<DIM - 1> faceCenter =
        face->getReferenceGeometry()->getCenter();
    const Geometry::PointReference<DIM> elementFaceCenter =
        isLeft ? face->mapRefFaceToRefElemL(faceCenter)
               : face->mapRefFaceToRefElemR(faceCenter);
    return element->referenceToPhysical(elementFaceCenter);
}


}

#endif  // HPGEM_KPHASESHIFTHELPERS_H
