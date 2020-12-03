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
#include "CGDGMatrixKPhaseShiftBuilder.h"

#include "Utilities/GlobalIndexing.h"
#include "KPhaseShiftHelpers.h"

namespace DGMax {

template <std::size_t DIM>
KPhaseShifts<DIM> CGDGMatrixKPhaseShiftBuilder<DIM>::build(
    const Utilities::GlobalIndexing& cgIndexing,
    const Utilities::GlobalIndexing& dgIndexing) const {

    std::vector<DGMax::KPhaseShiftBlock<DIM>> result;
    std::set<const Base::Edge*> boundaryEdges;
    std::set<const Base::Node*> boundaryNodes;
    Base::MeshManipulatorBase* mesh = cgIndexing.getMesh();

    auto end = mesh->faceColEnd();
    for (Base::TreeIterator<Base::Face*> it = mesh->faceColBegin(); it != end;
         ++it) {
        Geometry::FaceType type = (*it)->getFaceType();
        if (type == Geometry::FaceType::PERIODIC_SUBDOMAIN_BC ||
            type == Geometry::FaceType::PERIODIC_BC) {

            if ((*it)->isOwnedByCurrentProcessor()) {
                this->addFacePhaseShifts(*it, cgIndexing, dgIndexing, result);
            }

            // This is a periodic boundary face, so all the nodes are also
            // placed on the periodic boundary.
            for (const Base::Node* node : (*it)->getNodesList()) {
                if (node->isOwnedByCurrentProcessor()) {
                    boundaryNodes.emplace(node);
                }
            }

            // Additionally, all the edges are also on the periodic boundary.
            const Base::Element* element = (*it)->getPtrElementLeft();
            const Geometry::ReferenceGeometry* referenceGeometry =
                element->getReferenceGeometry();
            std::vector<std::size_t> nodeIds(
                referenceGeometry->getCodim1EntityLocalIndices(
                    element->getLocalId(*it)));
            std::vector<std::size_t> edgeNodeIds;
            for (std::size_t i = 0; i < element->getNumberOfEdges(); ++i) {
                edgeNodeIds = referenceGeometry->getCodim2EntityLocalIndices(i);
                bool found1(false);
                bool found2(false);
                for (std::size_t j = 0; j < nodeIds.size(); ++j) {
                    found1 |= nodeIds[j] == edgeNodeIds[0];
                    found2 |= nodeIds[j] == edgeNodeIds[1];
                }
                if (found1 && found2 &&
                    element->getEdge(i)->isOwnedByCurrentProcessor()) {
                    boundaryEdges.emplace(element->getEdge(i));
                }
            }
        }

        // Now that all boundary nodes and edges have been de-duplicated, add
        // their shifts.
        for (const Base::Node* node : boundaryNodes) {
            this->addNodePhaseShifts(node, cgIndexing, dgIndexing, result);
        }
        for (const Base::Edge* edge : boundaryEdges) {
            this->addEdgePhaseShifts(edge, cgIndexing, dgIndexing, result);
        }
    }
    return KPhaseShifts<DIM>(result);
}
template <std::size_t DIM>
void CGDGMatrixKPhaseShiftBuilder<DIM>::addFacePhaseShifts(
    const Base::Face* face, const Utilities::GlobalIndexing& projectorIndex,
    const Utilities::GlobalIndexing& indexing,
    std::vector<DGMax::KPhaseShiftBlock<DIM>>& out) const {
    logger.assert_debug(face->isOwnedByCurrentProcessor(),
                        "Not owned by current processor");

    // Associate the node with an element that owns it
    const Base::Element* owningElement = face->getPtrElementLeft();
    const Geometry::PointPhysical<DIM> owningCoord =
        getCoordinate<DIM>(owningElement, face);

    addElementPhaseShift(face, owningCoord, face->getPtrElementRight(),
                         projectorIndex, indexing, out);
}

template <std::size_t DIM>
void CGDGMatrixKPhaseShiftBuilder<DIM>::addEdgePhaseShifts(
    const Base::Edge* edge, const Utilities::GlobalIndexing& projectorIndex,
    const Utilities::GlobalIndexing& indexing,
    std::vector<DGMax::KPhaseShiftBlock<DIM>>& out) const {
    logger.assert_debug(edge->isOwnedByCurrentProcessor(),
                        "Not owned by current processor");

    // Associate the node with an element that owns it
    const Base::Element* owningElement = edge->getOwningElement();
    const Geometry::PointPhysical<DIM> owningCoord =
        getCoordinate<DIM>(owningElement, edge);

    for (const Base::Element* element : edge->getElements()) {
        addElementPhaseShift(edge, owningCoord, element, projectorIndex,
                             indexing, out);
    }
}

template <std::size_t DIM>
void CGDGMatrixKPhaseShiftBuilder<DIM>::addNodePhaseShifts(
    const Base::Node* node, const Utilities::GlobalIndexing& projectorIndex,
    const Utilities::GlobalIndexing& indexing,
    std::vector<DGMax::KPhaseShiftBlock<DIM>>& out) const {
    logger.assert_debug(node->isOwnedByCurrentProcessor(),
                        "Not owned by current processor");

    // Associate the node with an element that owns it
    const Base::Element* owningElement = node->getOwningElement();
    const Geometry::PointPhysical<DIM> owningCoord =
        getCoordinate<DIM>(owningElement, node);

    for (const Base::Element* element : node->getElements()) {
        addElementPhaseShift(node, owningCoord, element, projectorIndex,
                             indexing, out);
    }
}

template <std::size_t DIM>
template <typename GEOM>
void CGDGMatrixKPhaseShiftBuilder<DIM>::addElementPhaseShift(
    const GEOM* geom, const Geometry::PointPhysical<DIM>& owningCoord,
    const Base::Element* element,
    const Utilities::GlobalIndexing& projectorIndex,
    const Utilities::GlobalIndexing& indexing,
    std::vector<DGMax::KPhaseShiftBlock<DIM>>& out) const {

    // Compute the vector a_T
    const Geometry::PointPhysical<DIM> elementCoord =
        getCoordinate<DIM>(element, geom);
    LinearAlgebra::SmallVector<DIM> dx =
        owningCoord.getCoordinates() - elementCoord.getCoordinates();
    // If the owning element and the current element are not on different sides
    // of the periodic boundary, then dx = a_T = 0 and there is no phase shift.
    if (dx.l2Norm() > 1e-12) {
        const Geometry::PointReference<DIM>& center =
            element->getReferenceGeometry()->getCenter();

        if (extraShift_) {
            // Rescaling of the columns, as the shifts for the field basis
            // functions are e^{-ikx}. Where the convention is used that the
            // trial basis functions use the - sign.
            dx -= extraShift_(element);
        }

        std::size_t numProjectorDoF = geom->getLocalNumberOfBasisFunctions(1);

        // Difference in coordinates for the node => shift is needed
        std::size_t numSolutionDoF = element->getLocalNumberOfBasisFunctions(0);
        LinearAlgebra::MiddleSizeMatrix matrix(numProjectorDoF, numSolutionDoF);
        const LinearAlgebra::MiddleSizeMatrix originalMatrix =
            matrixExtractor_(element);

        // Extract the submatrix from the element matrix corresponding to the cg
        // basis functions from geom and the dg basis functions on the element.
        std::size_t localOffset = element->getBasisFunctionOffset(geom, 1);
        for (std::size_t i = 0; i < numProjectorDoF; ++i) {
            for (std::size_t j = 0; j < numSolutionDoF; ++j) {
                matrix(i, j) = originalMatrix(i + localOffset, j);
            }
        }
        // Compute their global indices
        PetscInt globalProjectorIndex = projectorIndex.getGlobalIndex(geom, 1);
        PetscInt globalElementIndex = indexing.getGlobalIndex(element, 0);

        std::vector<PetscInt> rowIndices(numProjectorDoF);
        std::vector<PetscInt> colIndices(numSolutionDoF);
        std::iota(rowIndices.begin(), rowIndices.end(), globalProjectorIndex);
        std::iota(colIndices.begin(), colIndices.end(), globalElementIndex);

        checkMatrixSize(matrix, rowIndices.size(), colIndices.size());

        out.emplace_back(DGMax::MatrixBlocks(rowIndices, colIndices, matrix),
                         dx);
    }
}

template class CGDGMatrixKPhaseShiftBuilder<2>;
template class CGDGMatrixKPhaseShiftBuilder<3>;

}  // namespace DGMax