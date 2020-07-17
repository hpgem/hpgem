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
#include "KPhaseShift.h"

#include "LinearAlgebra/MiddleSizeMatrix.h"
#include "Utilities/GlobalIndexing.h"

namespace DGMax {

template <std::size_t DIM>
void KPhaseShiftBlock<DIM>::apply(LinearAlgebra::SmallVector<DIM> k,
                        std::vector<PetscScalar>& storage, Mat mat) const {
    if (shiftNeeded(k)) {
        // Load values in storage
        blocks_.loadBlocks(storage);
        // Apply phase shift
        std::size_t blockSize = blocks_.getBlockSize();
        const std::complex<double> phase =
            std::exp(std::complex<double>(0, k * dx_));

        for (std::size_t i = 0; i < blockSize; ++i) {
            storage[i] *= phase;
        }
        if (blocks_.isPair()) {
            const std::complex<double> antiPhase =
                std::exp(std::complex<double>(0, -(k * dx_)));
            for (std::size_t i = 0; i < blockSize; ++i) {
                storage[i + blockSize] *= antiPhase;
            }
        }
        // Insert the now phase shifted blocks into the matrix
        blocks_.insertBlocks(storage, mat);
    }
}

/// Phase Shifts
template <std::size_t DIM>
void KPhaseShifts<DIM>::apply(LinearAlgebra::SmallVector<DIM> k,
                              Mat mat) const {
    PetscErrorCode error = MatSetOption(mat, MAT_ROW_ORIENTED, PETSC_FALSE);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    std::vector<PetscScalar> storage;
    for (const KPhaseShiftBlock<DIM>& block : blocks_) {
        block.apply(k, storage, mat);
    }

    // Go back to the default row orientation
    error = MatSetOption(mat, MAT_ROW_ORIENTED, PETSC_TRUE);
    CHKERRABORT(PETSC_COMM_WORLD, error);

    error = MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    error = MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);
    CHKERRABORT(PETSC_COMM_WORLD, error);
}

/// Helper functions for the builders

// Given a point x on a (periodic boundary) face, which has coordinates x_l and
// x_r as seen from the left and right face respectively. Compute the difference
// x_l - x_r,
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

/// Builders

template <std::size_t DIM>
KPhaseShifts<DIM> FaceMatrixKPhaseShiftBuilder<DIM>::build(
    const Utilities::GlobalIndexing& indexing) const {
    std::vector<KPhaseShiftBlock<DIM>> result;
    const Base::MeshManipulatorBase* mesh = indexing.getMesh();

    auto end = mesh->faceColEnd(Base::IteratorType::GLOBAL);
    for (auto it = mesh->faceColBegin(Base::IteratorType::GLOBAL); it != end;
         ++it) {
        // With ghost faces some faces may not have two elements associated with
        // them. Calculating the boundary shift for these faces is not needed
        // nor possible, so skip them.
        if (!(*it)->isInternal()) {
            continue;
        }
        Geometry::FaceType type = (*it)->getFaceType();
        if (type == Geometry::FaceType::PERIODIC_SUBDOMAIN_BC ||
            type == Geometry::FaceType::PERIODIC_BC) {
            // Phase shift needed for the periodic BC.
            if (!(*it)->getPtrElementLeft()->isOwnedByCurrentProcessor() &&
                !(*it)->getPtrElementRight()->isOwnedByCurrentProcessor()) {
                // For each face two blocks need to be 'shifted' in the
                // stiffness matrix (that are Hermitian counterparts). Each
                // processor only shifts the rows in the matrix that it owns.
                // Shifting for a face is thus only necessary when owning an
                // element to either side.
                continue;
            }
            result.emplace_back(facePhaseShift(*it, indexing));
        } else if ((boundaryFaceShift<DIM>(*it)).l2Norm() > 1e-3) {
            // This seems to be a periodic boundary that was not marked
            // correctly. If we don't own either side this is not a problem for
            // this processor.
            logger.assert_always(
                !(*it)->getPtrElementLeft()->isOwnedByCurrentProcessor() &&
                    !(*it)->getPtrElementRight()->isOwnedByCurrentProcessor(),
                "Incorrectly marked periodic boundary");
        }
    }
    return KPhaseShifts<DIM>(result);
}

template <std::size_t DIM>
KPhaseShiftBlock<DIM> FaceMatrixKPhaseShiftBuilder<DIM>::facePhaseShift(
    const Base::Face* face, const Utilities::GlobalIndexing& indexing) const {

    // Extra data
    LinearAlgebra::SmallVector<DIM> dx = boundaryFaceShift<DIM>(face);
    if (extraShift_) {
        dx += extraShift_(face);
    }
    LinearAlgebra::MiddleSizeMatrix block1, block2;
    std::tie(block1, block2) = matrixExtractor_(face);

    // Determine element ownership //
    /////////////////////////////////

    // The face is a periodic boundary, but could also be a subdomain boundary.
    // For a non-subdomain boundary the current processor is responsible for
    // both off diagonal blocks in the face matrix. On a subdomain boundary it
    // is only responsible of the off diagonal block that has locally owned rows
    // (i.e. test functions).
    bool isSubdomainBoundary;
    // Element that this processor owns
    const Base::Element* ownedElement = face->getPtrElementLeft();
    // Element that this processor might own if this is not a subdomain boundary
    const Base::Element* otherElement = face->getPtrElementRight();
    if (!ownedElement->isOwnedByCurrentProcessor()) {
        isSubdomainBoundary = true;
        std::swap(ownedElement, otherElement);
        std::swap(block1, block2);
        // dx is the jump in coordinate to go from other->owned, swapping owned
        // and other means that we get the inverse jump.
        dx = -dx;

    } else {
        isSubdomainBoundary = !otherElement->isOwnedByCurrentProcessor();
    }
    logger.assert_debug(ownedElement->isOwnedByCurrentProcessor(),
                        "Not owning either side of the face!");

    // Compute vectors with the global indices //
    /////////////////////////////////////////////
    std::vector<PetscInt> rowIndices(
        ownedElement->getNumberOfBasisFunctions(0));
    std::iota(rowIndices.begin(), rowIndices.end(),
              indexing.getGlobalIndex(ownedElement, 0));

    std::vector<PetscInt> colIndices(
        otherElement->getNumberOfBasisFunctions(0));
    std::iota(colIndices.begin(), colIndices.end(),
              indexing.getGlobalIndex(otherElement, 0));

    // Construct the result //
    //////////////////////////

    if (isSubdomainBoundary) {
        return KPhaseShiftBlock<DIM>(DGMax::MatrixBlocks(rowIndices, colIndices, block1),
                           dx);
    } else {
        return KPhaseShiftBlock<DIM>(
            DGMax::MatrixBlocks(rowIndices, colIndices, block1, block2), dx);
    }
}

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

        out.emplace_back(DGMax::MatrixBlocks(rowIndices, colIndices, matrix),
                         dx);
    }
}

template class KPhaseShifts<2>;
template class KPhaseShifts<3>;

template class FaceMatrixKPhaseShiftBuilder<2>;
template class FaceMatrixKPhaseShiftBuilder<3>;

template class CGDGMatrixKPhaseShiftBuilder<2>;
template class CGDGMatrixKPhaseShiftBuilder<3>;

}  // namespace DGMax