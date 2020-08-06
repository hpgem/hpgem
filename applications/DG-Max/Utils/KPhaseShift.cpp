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

/// Helper class for extractSubMatrix() describing either the rows or columns.
///
/// \tparam GEOM
template <typename GEOM>
struct SideDescriptor {

    SideDescriptor()
        : unknowns(), desiredUnknowns(), geom(nullptr), element(nullptr){};

    /// All the unknowns that form all the row/column of the local matrix
    std::vector<std::size_t> unknowns;
    /// The subset of the unknowns that needs to be extracted
    std::vector<std::size_t> desiredUnknowns;
    /// The geometry part for which to extract the entries.
    const GEOM* geom;
    /// The element to which the side of the matrix corresponds, geom should be
    /// adjacent to this element (not checked).
    const Base::Element* element;

    /// Compute the indices in the local (input) matrix that form the resulting
    /// submatrix. In matlab syntax:
    /// result = input(rows.localIndices(), cols.localIndices())
    std::vector<std::size_t> localIndices() const {
        // Index in the local matrix where, for each unknown, the basis
        // functions for the geometrical part start.
        std::vector<std::size_t> unknownOffsets(element->getNumberOfUnknowns() +
                                                1);
        // First compute the offsets of the unknowns
        for (const std::size_t& unknown : unknowns) {
            unknownOffsets[unknown + 1] =
                element->getNumberOfBasisFunctions(unknown);
        }
        std::partial_sum(unknownOffsets.begin(), unknownOffsets.end(),
                         unknownOffsets.begin());
        // Add the local offset for the unknown
        for (const std::size_t& unknown : unknowns) {
            unknownOffsets[unknown] +=
                element->getBasisFunctionOffset(geom, unknown);
        }

        // Compute the result mapping
        std::vector<std::size_t> result(numberOfEntries());
        std::size_t resultOffset = 0;
        for (const std::size_t& unknown : desiredUnknowns) {
            std::size_t numDoFs = geom->getLocalNumberOfBasisFunctions(unknown);
            for (std::size_t i = 0; i < numDoFs; ++i) {
                result[i + resultOffset] = unknownOffsets[unknown] + i;
            }
            resultOffset += numDoFs;
        }
        return result;
    }

    /// Global indices corresponding to the extracted block.
    ///
    /// \param indexing The GlobalIndex of this side, assumes that
    /// indexing->getIncludedUnknowns() == unknowns (not checked).
    ///
    /// \return The global indices for the extracted block.
    std::vector<PetscInt> globalIndices(
        const Utilities::GlobalIndexing* indexing) {
        const std::size_t size = numberOfEntries();
        std::vector<PetscInt> result(size);
        auto start = result.begin();
        for (const std::size_t& unknown : unknowns) {
            std::size_t numDoFs = geom->getLocalNumberOfBasisFunctions(unknown);
            std::size_t offset = indexing->getGlobalIndex(geom, unknown);
            std::iota(start, start + numDoFs, offset);
            start += numDoFs;
        }
        logger.assert_debug(result.size() == size, "Incorrect size given");
        return result;
    }

   private:
    std::size_t numberOfEntries() const {
        std::size_t result = 0;
        for (const std::size_t& unknown : unknowns) {
            result += geom->getLocalNumberOfBasisFunctions(unknown);
        }
        return result;
    }
};

/// Given a local matrix, extract the submatrix that corresponds to the basis
/// functions that correspond to a certain geometrical part and from these basis
/// functions only those for a certain unknown.
template <typename RGeom, typename CGeom>
LinearAlgebra::MiddleSizeMatrix extractSubMatrix(
    const SideDescriptor<RGeom>& rows, const SideDescriptor<CGeom>& cols,
    const LinearAlgebra::MiddleSizeMatrix& localMatrix) {

    std::vector<std::size_t> rowIndices = rows.localIndices();
    std::vector<std::size_t> colIndices = cols.localIndices();
    std::size_t numRows = rowIndices.size();
    std::size_t numCols = colIndices.size();

    LinearAlgebra::MiddleSizeMatrix result(numRows, numCols);

    for (std::size_t i = 0; i < numRows; ++i) {
        for (std::size_t j = 0; j < numCols; ++j) {
            result(i, j) = localMatrix(rowIndices[i], colIndices[j]);
        }
    }
    return result;
}

template <std::size_t DIM>
void KPhaseShiftBlock<DIM>::apply(LinearAlgebra::SmallVector<DIM> k,
                                  std::vector<PetscScalar>& storage,
                                  Mat mat) const {
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

//////////////
// BUILDERS //
//////////////

void checkUnknowns(const std::vector<std::size_t>& unknowns,
                   const Utilities::GlobalIndexing* indexing) {
    if (indexing == nullptr) {
        logger(WARN, "Setting unknowns without indexing, no verification done");
        return;
    }
    // Check whether the unknowns are included in the index;
    for (const std::size_t unknown : unknowns) {
        const std::vector<std::size_t>& includedUnknowns =
            indexing->getIncludedUnknowns();
        logger.assert_always(
            std::find(includedUnknowns.begin(), includedUnknowns.end(),
                      unknown) != includedUnknowns.end(),
            "Unknown % is not available in the index", unknown);
    }
}

template <std::size_t DIM>
void FaceMatrixKPhaseShiftBuilder<DIM>::setUnknowns(
    const std::vector<std::size_t>& unknowns) {
    checkUnknowns(unknowns, indexing_);
    unknowns_ = unknowns;
    std::sort(unknowns_.begin(), unknowns_.end());
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

//

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

void checkMatrixSize(const LinearAlgebra::MiddleSizeMatrix& mat,
                     std::size_t rows, std::size_t cols) {
    logger.assert_always(mat.getNumberOfRows() == rows,
                         "Incorrect number of rows, expected % got %", rows,
                         mat.getNumberOfRows());
    logger.assert_always(mat.getNumberOfColumns() == cols,
                         "Incorrect number of columns, expected % got %", cols,
                         mat.getNumberOfColumns());
}

/// Builders

template <std::size_t DIM>
KPhaseShifts<DIM> FaceMatrixKPhaseShiftBuilder<DIM>::build() const {
    std::vector<KPhaseShiftBlock<DIM>> result;
    logger.assert_always(indexing_ != nullptr, "No indexing set");
    const Base::MeshManipulatorBase* mesh = indexing_->getMesh();

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
            result.emplace_back(facePhaseShift(*it));
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
    const Base::Face* face) const {

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

    // Compute the actual content of the shift //
    /////////////////////////////////////////////

    SideDescriptor<Base::Element> rows, cols;
    rows.unknowns = indexing_->getIncludedUnknowns();
    rows.desiredUnknowns = unknowns_;
    rows.element = ownedElement;
    rows.geom = ownedElement;

    cols.unknowns = indexing_->getIncludedUnknowns();
    cols.desiredUnknowns = unknowns_;
    cols.element = otherElement;
    cols.geom = otherElement;

    std::vector<PetscInt> rowIndices = rows.globalIndices(indexing_);
    std::vector<PetscInt> colIndices = cols.globalIndices(indexing_);
    block1 = extractSubMatrix(rows, cols, block1);
    checkMatrixSize(block1, rowIndices.size(), colIndices.size());
    if (!isSubdomainBoundary) {
        block2 = extractSubMatrix(cols, rows, block2);
        checkMatrixSize(block2, colIndices.size(), rowIndices.size());
    }

    // Construct the result //
    //////////////////////////

    if (isSubdomainBoundary) {
        return KPhaseShiftBlock<DIM>(
            DGMax::MatrixBlocks(rowIndices, colIndices, block1), dx);
    } else {
        return KPhaseShiftBlock<DIM>(
            DGMax::MatrixBlocks(rowIndices, colIndices, block1, block2), dx);
    }
}

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
                for (unsigned long nodeId : nodeIds) {
                    found1 |= nodeId == edgeNodeIds[0];
                    found2 |= nodeId == edgeNodeIds[1];
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

template class KPhaseShifts<2>;
template class KPhaseShifts<3>;

template class FaceMatrixKPhaseShiftBuilder<2>;
template class FaceMatrixKPhaseShiftBuilder<3>;

template class CGDGMatrixKPhaseShiftBuilder<2>;
template class CGDGMatrixKPhaseShiftBuilder<3>;

}  // namespace DGMax