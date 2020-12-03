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
#include "FaceKPhaseShiftBuilder.h"

#include "Utilities/GlobalIndexing.h"
#include "KPhaseShiftHelpers.h"

namespace DGMax {

template class FaceMatrixKPhaseShiftBuilder<2>;
template class FaceMatrixKPhaseShiftBuilder<3>;

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

    // Using the functionality from GlobalIndexing is only safe for DG type
    // basis functions, as they are only located on the elements of the mesh.
    // This is not checked and left up to the user.
    std::vector<PetscInt> rowIndices = indexing.getGlobalIndices(ownedElement);
    std::vector<PetscInt> colIndices = indexing.getGlobalIndices(otherElement);

    checkMatrixSize(block1, rowIndices.size(), colIndices.size());
    checkMatrixSize(block2, colIndices.size(), rowIndices.size());

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

}  // namespace DGMax
