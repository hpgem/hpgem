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

#include "Base/MpiContainer.h"
#include "GlobalMatrix.h"
#include "SparsityEstimator.h"
#include <vector>
#include "Base/MeshManipulatorBase.h"
#include "Base/Edge.h"
#include "Base/Face.h"
#include "Base/Element.h"
#include "LinearAlgebra/MiddleSizeVector.h"
#include "Geometry/PointPhysical.h"
#include "Geometry/PhysicalGeometry.h"
#include "Geometry/ReferenceGeometry.h"
#include "Geometry/PointReference.h"
#include "Base/Mesh.h"
#include "Logger.h"
#include <numeric>

namespace Utilities {

GlobalMatrix::GlobalMatrix(const GlobalIndexing& rowIndexing,
                           const GlobalIndexing& columnIndexing,
                           int elementMatrixID, int faceMatrixID)
    : meshLevel_(-2),
      rowIndexing_(rowIndexing),
      columnIndexing_(columnIndexing),
      elementMatrixID_(elementMatrixID)
      // Use pointer comparison to see if the row & colum indices are the same
      ,
      symmetricIndexing_(&rowIndexing == &columnIndexing),
      faceMatrixID_(faceMatrixID) {
    logger.assert_always(columnIndexing.getMesh() == nullptr ||
                             rowIndexing.getMesh() == nullptr ||
                             columnIndexing.getMesh() == rowIndexing.getMesh(),
                         "Row and column indexing based on different meshes.");
}

void GlobalMatrix::getMatrixBCEntries(const Base::Face* face,
                                      std::size_t& numberOfEntries,
                                      std::vector<int>& entries) {
    logger.assert_debug(face != nullptr, "Invalid face passed");
    if (rowIndexing_.getMesh() == nullptr) {
        entries.clear();
        numberOfEntries = 0;
        return;
    }
    // Face basis functions
    const std::vector<std::size_t>& unknowns =
        rowIndexing_.getIncludedUnknowns();
    for (std::size_t unknown : unknowns) {
        std::size_t nBasis = face->getLocalNumberOfBasisFunctions(unknown);
        int elementBasis0 = rowIndexing_.getGlobalIndex(face, unknown);
        for (std::size_t basisId = 0; basisId < nBasis; ++basisId) {
            entries.push_back(elementBasis0 + basisId);
        }
        numberOfEntries += nBasis;
    }
    // Edges around the face

    std::vector<std::size_t> nodeEntries =
        face->getPtrElementLeft()
            ->getPhysicalGeometry()
            ->getGlobalFaceNodeIndices(face->localFaceNumberLeft());
    std::vector<std::size_t> edgeIndex(2);
    for (std::size_t i = 0; i < face->getPtrElementLeft()->getNumberOfEdges();
         ++i) {
        edgeIndex = face->getPtrElementLeft()
                        ->getReferenceGeometry()
                        ->getCodim2EntityLocalIndices(i);
        edgeIndex[0] =
            face->getPtrElementLeft()->getPhysicalGeometry()->getNodeIndex(
                edgeIndex[0]);
        edgeIndex[1] =
            face->getPtrElementLeft()->getPhysicalGeometry()->getNodeIndex(
                edgeIndex[1]);
        bool firstFound(false), secondFound(false);
        for (unsigned long nodeEntrie : nodeEntries) {
            if (nodeEntrie == edgeIndex[0]) firstFound = true;
            if (nodeEntrie == edgeIndex[1]) secondFound = true;
        }
        if (firstFound && secondFound) {
            const Base::Edge* edge = face->getPtrElementLeft()->getEdge(i);
            for (std::size_t unknown : unknowns) {
                std::size_t nEdgeBasis =
                    edge->getLocalNumberOfBasisFunctions(unknown);
                int edgeBasis0 = rowIndexing_.getGlobalIndex(edge, unknown);
                for (std::size_t basisId = 0; basisId < nEdgeBasis; ++basisId) {
                    entries.push_back(edgeBasis0 + basisId);
                }
                numberOfEntries += nEdgeBasis;
            }
        }
    }
    // Nodes around the face
    if (rowIndexing_.getMesh()->dimension() > 1) {
        nodeEntries =
            face->getPtrElementLeft()
                ->getPhysicalGeometry()
                ->getLocalFaceNodeIndices(face->localFaceNumberLeft());
        for (std::size_t i : nodeEntries) {
            const Base::Node* node = face->getPtrElementLeft()->getNode(i);
            for (std::size_t unknown : unknowns) {
                std::size_t nNodeBasis =
                    node->getLocalNumberOfBasisFunctions(unknown);
                int nodeBasis0 = rowIndexing_.getGlobalIndex(node, unknown);
                for (std::size_t basisId = 0; basisId < nNodeBasis; ++basisId) {
                    entries.push_back(nodeBasis0 + basisId);
                }
                numberOfEntries += nNodeBasis;
            }
        }
    }
}

#if defined(HPGEM_USE_ANY_PETSC)

GlobalPetscMatrix::GlobalPetscMatrix(const GlobalIndexing& rowIndexing,
                                     const GlobalIndexing& columnIndexing,
                                     int elementMatrixID, int faceMatrixID)
    : GlobalMatrix(rowIndexing, columnIndexing, elementMatrixID, faceMatrixID) {
    PetscBool petscRuns;
    PetscInitialized(&petscRuns);
    logger.assert_debug(
        petscRuns == PETSC_TRUE,
        "Early call, firstly the command line arguments should be parsed");

    // Dummy call to always have an valid matrix in A_.
    MatCreateSeqAIJ(PETSC_COMM_SELF, 0, 0, 0, PETSC_NULL, &A_);

    reinit();
}

GlobalPetscMatrix::~GlobalPetscMatrix() {
    int ierr = MatDestroy(&A_);
    // giving error about Petsc has generated inconsistent data and likely
    // memory corruption in heap
    CHKERRV(ierr);
}

GlobalPetscMatrix::operator Mat() {
    if (HPGEM_LOGLEVEL >= Log::DEBUG) {
        MatChop(A_, 1e-13);
        MatScale(A_, 24.);
        MatView(A_, PETSC_VIEWER_STDOUT_WORLD);
        MatScale(A_, 1. / 24.);
    }
    return A_;
}

void GlobalPetscMatrix::reinit() {
    if (rowIndexing_.getMesh() == nullptr ||
        columnIndexing_.getMesh() == nullptr) {
        // Uninitialized index for at least the rows or columns. So replace it
        // with an empty matrix.
        PetscErrorCode error;
        error = MatDestroy(&A_);
        CHKERRABORT(PETSC_COMM_WORLD, error);
        error = MatCreateSeqAIJ(PETSC_COMM_SELF, 0, 0, 0, PETSC_NULL, &A_);
        CHKERRABORT(PETSC_COMM_WORLD, error);
        return;
    } else {
        logger.assert_always(
            rowIndexing_.getMesh() == columnIndexing_.getMesh(),
            "Row and column indexing from different meshes");
        createMat();
        assemble();
    }
}

void GlobalPetscMatrix::createMat() {
    PetscErrorCode ierr = MatDestroy(&A_);
    CHKERRV(ierr);

    // Create matrix
    ierr = MatCreate(PETSC_COMM_WORLD, &A_);
    CHKERRV(ierr);
    ierr = MatSetType(A_, MATMPIAIJ);
    CHKERRV(ierr);

    // Set sizes based on indexing
    const std::size_t numberOfLocalRows =
        rowIndexing_.getNumberOfLocalBasisFunctions();
    const std::size_t numberOfLocalColumns =
        columnIndexing_.getNumberOfLocalBasisFunctions();
    MatSetSizes(A_, numberOfLocalRows, numberOfLocalColumns, PETSC_DETERMINE,
                PETSC_DETERMINE);

    if (rowIndexing_.getMesh() != nullptr) {
        std::vector<PetscInt> numberOfPositionsPerRow;
        std::vector<PetscInt> offDiagonalPositionsPerRow;
        logger.assert_always(rowIndexing_.getMesh() != nullptr, "Null mesh");
        SparsityEstimator estimator(rowIndexing_, columnIndexing_);
        estimator.computeSparsityEstimate(numberOfPositionsPerRow,
                                          offDiagonalPositionsPerRow,
                                          faceMatrixID_ >= 0);
        // Zeros are passed for the ignored arguments.
        ierr = MatMPIAIJSetPreallocation(A_, 0, numberOfPositionsPerRow.data(),
                                         0, offDiagonalPositionsPerRow.data());
        CHKERRV(ierr);
    } else {
        logger.assert_always(
            numberOfLocalRows == 0 && numberOfLocalColumns == 0,
            "Degrees of freedom without a mesh");
        ierr = MatMPIAIJSetPreallocation(A_, 0, nullptr, 0, nullptr);
        CHKERRV(ierr);
    }

    // Most options can only be set after the preallocation is done
    MatSetOption(A_, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);  // performance
    // While the estimates should be correct, the zeroing of a row without
    // anything on the diagonal will cause an error with the previous option but
    // without the following.
    MatSetOption(A_, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);
    // Enable the following to test that no new allocations were needed
    // MatSetOption(A_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
}

void GlobalPetscMatrix::assemble() {
    int ierr = MatZeroEntries(A_);
    CHKERRV(ierr);

    const Base::MeshManipulatorBase* mesh = rowIndexing_.getMesh();
    if (mesh == nullptr) {
        // Empty matrix, nothing to assemble
        return;
    }
    // MediumSizeMatrix uses column oriented storage
    ierr = MatSetOption(A_, MAT_ROW_ORIENTED, PETSC_FALSE);
    CHKERRV(ierr);
    if (elementMatrixID_ >= 0) {
        std::vector<PetscInt> localToGlobalRow;
        std::vector<PetscInt> localToGlobalColumn;
        const std::vector<PetscInt>& localToGlobalColumnRef =
            symmetricIndexing_ ? localToGlobalRow : localToGlobalColumn;
        for (Base::Element* element : mesh->getElementsList()) {
            rowIndexing_.getGlobalIndices(element, localToGlobalRow);
            if (!symmetricIndexing_) {
                columnIndexing_.getGlobalIndices(element, localToGlobalColumn);
            }
            const LinearAlgebra::MiddleSizeMatrix& localMatrix =
                element->getElementMatrix(elementMatrixID_);
            logger.assert_debug(
                localMatrix.getNumberOfRows() == localToGlobalRow.size() &&
                    localMatrix.getNumberOfColumns() ==
                        localToGlobalColumnRef.size(),
                "Incorrect element matrix size");
            ierr = MatSetValues(
                A_, localToGlobalRow.size(), localToGlobalRow.data(),
                localToGlobalColumnRef.size(), localToGlobalColumnRef.data(),
                localMatrix.data(), ADD_VALUES);
            CHKERRV(ierr);
        }
    }

    LinearAlgebra::MiddleSizeMatrix faceMatrix;
    if (faceMatrixID_ >= 0) {
        std::vector<PetscInt> localToGlobalRow;
        std::vector<PetscInt> localToGlobalColumn;
        const std::vector<PetscInt>& localToGlobalColumnRef =
            symmetricIndexing_ ? localToGlobalRow : localToGlobalColumn;
        for (Base::Face* face : mesh->getFacesList()) {
            if (!face->isOwnedByCurrentProcessor()) continue;

            rowIndexing_.getGlobalIndices(face, localToGlobalRow);
            if (!symmetricIndexing_) {
                columnIndexing_.getGlobalIndices(face, localToGlobalColumn);
            }
            faceMatrix = face->getFaceMatrixMatrix(faceMatrixID_);
            logger.assert_debug(
                faceMatrix.getNumberOfRows() == localToGlobalRow.size() &&
                    faceMatrix.getNumberOfColumns() ==
                        localToGlobalColumnRef.size(),
                "Incorrect face matrix size");
            ierr = MatSetValues(
                A_, localToGlobalRow.size(), localToGlobalRow.data(),
                localToGlobalColumnRef.size(), localToGlobalColumnRef.data(),
                faceMatrix.data(), ADD_VALUES);
            CHKERRV(ierr);
        }
    }
    // Reset the matrix to row oriented storage
    ierr = MatSetOption(A_, MAT_ROW_ORIENTED, PETSC_TRUE);
    CHKERRV(ierr);

    // Matrix assembly
    ierr = MatAssemblyBegin(A_, MAT_FINAL_ASSEMBLY);
    CHKERRV(ierr);
    ierr = MatAssemblyEnd(A_, MAT_FINAL_ASSEMBLY);
    CHKERRV(ierr);
}

void GlobalPetscMatrix::printMatInfo(MatInfoType type, std::ostream& stream) {
    MatInfo info;
    PetscErrorCode error = MatGetInfo(A_, type, &info);
    CHKERRABORT(PETSC_COMM_WORLD, error);
    stream << "Blocksize " << info.block_size << ", Nonzero " << info.nz_used
           << " used " << info.nz_unneeded << " unused."
           << " Assembled " << info.assemblies << " mallocs " << info.mallocs
           << std::endl;
}

void GlobalPetscMatrix::writeMatlab(const std::string& fileName) {
    PetscViewer viewer;
    PetscErrorCode err;

    err = PetscViewerASCIIOpen(PETSC_COMM_WORLD, fileName.c_str(), &viewer);
    CHKERRABORT(PETSC_COMM_WORLD, err);
    err = PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
    CHKERRABORT(PETSC_COMM_WORLD, err);
    err = MatView(A_, viewer);
    CHKERRABORT(PETSC_COMM_WORLD, err);
    err = PetscViewerDestroy(&viewer);
    CHKERRABORT(PETSC_COMM_WORLD, err);
}
#endif
}  // namespace Utilities
