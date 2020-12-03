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
KPhaseShifts<DIM> CGDGMatrixKPhaseShiftBuilder<DIM>::build() const {

    logger.assert_always(cgIndexing_ != nullptr, "Null CG index at build");
    logger.assert_always(dgIndexing_ != nullptr, "Null DG index at build");

    if (hermitian_) {
        logger.assert_always(
            cgIndexing_ == dgIndexing_,
            "Different index for CG and DG in Hermitian mode.");
        // Not verified that cgUnknowns_ and dgUnknowns_ are disjoint.
    }

    std::vector<DGMax::KPhaseShiftBlock<DIM>> result;
    std::set<const Base::Edge*> boundaryEdges;
    std::set<const Base::Node*> boundaryNodes;
    Base::MeshManipulatorBase* mesh = cgIndexing_->getMesh();

    auto end = mesh->faceColEnd();
    for (Base::TreeIterator<Base::Face*> it = mesh->faceColBegin(); it != end;
         ++it) {
        Geometry::FaceType type = (*it)->getFaceType();
        if (type == Geometry::FaceType::PERIODIC_SUBDOMAIN_BC ||
            type == Geometry::FaceType::PERIODIC_BC) {

            if ((*it)->getPtrElementLeft()->isOwnedByCurrentProcessor() ||
                (*it)->getPtrElementRight()->isOwnedByCurrentProcessor()) {
                this->addFacePhaseShifts(*it, result);
            }

            // This is a periodic boundary face, so all the nodes are also
            // placed on the periodic boundary.
            for (const Base::Node* node : (*it)->getNodesList()) {
                bool include = false;
                if (!hermitian_) {
                    // For the non Hermitian case the node forms the rows and
                    // should thus be owned.
                    include = node->isOwnedByCurrentProcessor();
                } else {
                    // Include not only if we own the node, but also any of its
                    // adjacent elements.
                    for (const Base::Element* element : node->getElements()) {
                        if (element->isOwnedByCurrentProcessor()) {
                            include = true;
                            break;
                        }
                    }
                }
                if (include) {
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
                if (found1 && found2) {
                    const Base::Edge* edge = element->getEdge(i);
                    // Same as for node, check if it needs to be included.
                    bool include = false;
                    if (!hermitian_) {
                        include = edge->isOwnedByCurrentProcessor();
                    } else {

                        for (const Base::Element* element1 :
                             edge->getElements()) {
                            if (element1->isOwnedByCurrentProcessor()) {
                                include = true;
                                break;
                            }
                        }
                    }
                    if (include) {
                        boundaryEdges.emplace(element->getEdge(i));
                    }
                }
            }
        }
    }

    // Now that all boundary nodes and edges have been de-duplicated, add
    // their shifts.
    for (const Base::Node* node : boundaryNodes) {
        this->addNodePhaseShifts(node, result);
    }
    for (const Base::Edge* edge : boundaryEdges) {
        this->addEdgePhaseShifts(edge, result);
    }
    return KPhaseShifts<DIM>(result);
}
template <std::size_t DIM>
void CGDGMatrixKPhaseShiftBuilder<DIM>::addFacePhaseShifts(
    const Base::Face* face,
    std::vector<DGMax::KPhaseShiftBlock<DIM>>& out) const {

    // Associate the node with an element that owns it
    const Base::Element* owningElement = face->getPtrElementLeft();
    const Geometry::PointPhysical<DIM> owningCoord =
        getCoordinate<DIM>(owningElement, face);

    // Dirty shortcut: The face DoFs are associated with the Left face.
    // Therefore, the left face does not need any shifts and can be skipped.

    if (face->isOwnedByCurrentProcessor() ||
        face->getPtrElementRight()->isOwnedByCurrentProcessor()) {
        addElementPhaseShift(face, owningCoord, face->getPtrElementRight(),
                             out);
    }
}

template <std::size_t DIM>
void CGDGMatrixKPhaseShiftBuilder<DIM>::addEdgePhaseShifts(
    const Base::Edge* edge,
    std::vector<DGMax::KPhaseShiftBlock<DIM>>& out) const {

    // Associate the node with an element that owns it
    const Base::Element* owningElement = edge->getOwningElement();
    const Geometry::PointPhysical<DIM> owningCoord =
        getCoordinate<DIM>(owningElement, edge);

    for (const Base::Element* element : edge->getElements()) {
        if (edge->isOwnedByCurrentProcessor() ||
            element->isOwnedByCurrentProcessor()) {
            addElementPhaseShift(edge, owningCoord, element, out);
        }
    }
}

template <std::size_t DIM>
void CGDGMatrixKPhaseShiftBuilder<DIM>::addNodePhaseShifts(
    const Base::Node* node,
    std::vector<DGMax::KPhaseShiftBlock<DIM>>& out) const {

    // Associate the node with an element that owns it
    const Base::Element* owningElement = node->getOwningElement();
    const Geometry::PointPhysical<DIM> owningCoord =
        getCoordinate<DIM>(owningElement, node);

    for (const Base::Element* element : node->getElements()) {
        if (node->isOwnedByCurrentProcessor() ||
            element->isOwnedByCurrentProcessor()) {
            addElementPhaseShift(node, owningCoord, element, out);
        }
    }
}

template <std::size_t DIM>
template <typename GEOM>
void CGDGMatrixKPhaseShiftBuilder<DIM>::addElementPhaseShift(
    const GEOM* geom, const Geometry::PointPhysical<DIM>& owningCoord,
    const Base::Element* element,
    std::vector<DGMax::KPhaseShiftBlock<DIM>>& out) const {

    if (hermitian_) {
        logger.assert_always(geom->isOwnedByCurrentProcessor() ||
                                 element->isOwnedByCurrentProcessor(),
                             "Hermitian case and owning neither side");
    } else {
        logger.assert_debug(
            geom->isOwnedByCurrentProcessor(),
            "Non hermitian case and not owning the geometry part");
    }

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

        const LinearAlgebra::MiddleSizeMatrix originalMatrix =
            matrixExtractor_(element);

        // The descriptors of each side
        SideDescriptor<GEOM> rows;
        rows.unknowns = cgIndexing_->getIncludedUnknowns();
        rows.desiredUnknowns = cgUnknowns_;
        rows.geom = geom;
        rows.element = element;

        SideDescriptor<Base::Element> cols;
        cols.unknowns = dgIndexing_->getIncludedUnknowns();
        cols.desiredUnknowns = dgUnknowns_;
        cols.geom = element;
        cols.element = element;

        // Compute their global indices
        std::vector<PetscInt> rowIndices = rows.globalIndices(cgIndexing_);
        std::vector<PetscInt> colIndices = cols.globalIndices(dgIndexing_);

        LinearAlgebra::MiddleSizeMatrix matrix1 =
            extractSubMatrix<GEOM, Base::Element>(rows, cols, originalMatrix);
        checkMatrixSize(matrix1, rowIndices.size(), colIndices.size());

        LinearAlgebra::MiddleSizeMatrix matrix2;
        if (hermitian_) {
            matrix2 = extractSubMatrix<Base::Element, GEOM>(cols, rows,
                                                            originalMatrix);
            checkMatrixSize(matrix2, colIndices.size(), rowIndices.size());
        }
        if (hermitian_ && element->isOwnedByCurrentProcessor() &&
            geom->isOwnedByCurrentProcessor()) {
            // Owning both blocks -> insert both
            out.emplace_back(
                DGMax::MatrixBlocks(rowIndices, colIndices, matrix1, matrix2),
                dx);
        } else {
            if (hermitian_ && !geom->isOwnedByCurrentProcessor()) {
                std::swap(rowIndices, colIndices);
                std::swap(matrix1, matrix2);
                dx = -dx;
            }
            out.emplace_back(
                DGMax::MatrixBlocks(rowIndices, colIndices, matrix1), dx);
        }
    }
}

template <std::size_t DIM>
void CGDGMatrixKPhaseShiftBuilder<DIM>::setCGUnknowns(
    const std::vector<std::size_t>& cgUnknowns) {
    checkUnknowns(cgUnknowns, cgIndexing_);
    cgUnknowns_ = cgUnknowns;
    std::sort(cgUnknowns_.begin(), cgUnknowns_.end());
}

template <std::size_t DIM>
void CGDGMatrixKPhaseShiftBuilder<DIM>::setDGUnknowns(
    const std::vector<std::size_t>& dgUnknowns) {
    checkUnknowns(dgUnknowns, dgIndexing_);
    dgUnknowns_ = dgUnknowns;
    std::sort(dgUnknowns_.begin(), dgUnknowns_.end());
}

template class CGDGMatrixKPhaseShiftBuilder<2>;
template class CGDGMatrixKPhaseShiftBuilder<3>;

}  // namespace DGMax