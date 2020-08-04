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

inline void checkMatrixSize(const LinearAlgebra::MiddleSizeMatrix mat,
                            std::size_t rows, std::size_t cols) {
    logger.assert_always(mat.getNumberOfRows() == rows,
                         "Incorrect number of rows, expected % got %", rows,
                         mat.getNumberOfRows());
    logger.assert_always(mat.getNumberOfColumns() == cols,
                         "Incorrect number of columns, expected % got %", cols,
                         mat.getNumberOfColumns());
}

inline void checkUnknowns(const std::vector<std::size_t>& unknowns,
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

// Given a point x on a (periodic boundary) face, which has coordinates x_l and
// x_r as seen from the left and right face respectively. Compute the difference
// x_l - x_r,

/// \brief Compute the coordinate jump of a face.
///
/// For a point P on a (periodic boundary) the coordinate may be different when
/// computed from the left or right face. This function computes the jump in
/// this coordinate: x_l - x_r, where x_l/x_r is the coordinate of point P as
/// seen from the left/right face. \tparam DIM \param face \return
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

}  // namespace DGMax

#endif  // HPGEM_KPHASESHIFTHELPERS_H
