/*
 This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
 http://hpgem.org/about-the-code/team
 
 This code is distributed using BSD 3-Clause License. A copy of which can found below.
 
 
 Copyright (c) 2014, University of Twente
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "Base/MpiContainer.h"
#include "GlobalMatrix.h"
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

namespace Utilities
{
    
    GlobalMatrix::GlobalMatrix(Base::MeshManipulatorBase* theMesh, int elementMatrixID, int faceMatrixID)
            : meshLevel_(-2), elementMatrixID_(elementMatrixID), faceMatrixID_(faceMatrixID), theMesh_(theMesh)
    {
        logger.assert_debug(theMesh != nullptr, "Invalid mesh passed");
    }
    
    void GlobalMatrix::getMatrixBCEntries(const Base::Face* face, std::size_t& numberOfEntries, std::vector<int>& entries)
    {
        logger.assert_debug(face != nullptr, "Invalid face passed");
        // Face basis functions
        std::size_t nFaceBasisLocal = face->getTotalLocalNumberOfBasisFunctions();
        std::size_t nUnknowns = face->getPtrElementLeft()->getNumberOfUnknowns();
        for (std::size_t unknown = 0; unknown < nUnknowns; ++unknown)
        {
            std::size_t nBasis = face->getLocalNumberOfBasisFunctions(unknown);
            int elementBasis0 = indexing_.getGlobalIndex(face, unknown);
            for (std::size_t basisId = 0; basisId < nBasis; ++basisId)
            {
                entries.push_back(elementBasis0 + basisId);
            }
            numberOfEntries += nBasis;
        }
        // Edges around the face

        std::vector<std::size_t> nodeEntries = face->getPtrElementLeft()->getPhysicalGeometry()->getGlobalFaceNodeIndices(face->localFaceNumberLeft());
        std::vector<std::size_t> edgeIndex(2);
        for (std::size_t i = 0; i < face->getPtrElementLeft()->getNumberOfEdges(); ++i)
        {
            edgeIndex = face->getPtrElementLeft()->getReferenceGeometry()->getCodim2EntityLocalIndices(i);
            edgeIndex[0] = face->getPtrElementLeft()->getPhysicalGeometry()->getNodeIndex(edgeIndex[0]);
            edgeIndex[1] = face->getPtrElementLeft()->getPhysicalGeometry()->getNodeIndex(edgeIndex[1]);
            bool firstFound(false), secondFound(false);
            for (std::size_t j = 0; j < nodeEntries.size(); ++j)
            {
                if (nodeEntries[j] == edgeIndex[0])
                    firstFound = true;
                if (nodeEntries[j] == edgeIndex[1])
                    secondFound = true;
            }
            if (firstFound && secondFound)
            {
                const Base::Edge* edge = face->getPtrElementLeft()->getEdge(i);
                for (std::size_t unknown = 0; unknown < nUnknowns; ++unknown)
                {
                    std::size_t nEdgeBasis = edge->getLocalNumberOfBasisFunctions(unknown);
                    int edgeBasis0 = indexing_.getGlobalIndex(edge, unknown);
                    for (std::size_t basisId = 0; basisId < nEdgeBasis; ++basisId)
                    {
                        entries.push_back(edgeBasis0 + basisId);
                    }
                    numberOfEntries += nEdgeBasis;
                }
            }
        }
        // Nodes around the face
        if (theMesh_->dimension() > 1)
        {
            nodeEntries = face->getPtrElementLeft()->getPhysicalGeometry()->getLocalFaceNodeIndices(
                    face->localFaceNumberLeft());
            for (std::size_t i : nodeEntries)
            {
                const Base::Node *node = face->getPtrElementLeft()->getNode(i);
                for (std::size_t unknown = 0; unknown < nUnknowns; ++unknown)
                {
                    std::size_t nNodeBasis = node->getLocalNumberOfBasisFunctions(unknown);
                    int nodeBasis0 = indexing_.getGlobalIndex(node, unknown);
                    for (std::size_t basisId = 0; basisId < nNodeBasis; ++basisId)
                    {
                        entries.push_back(nodeBasis0 + basisId);
                    }
                    numberOfEntries += nNodeBasis;
                }
            }
        }
    }

#if defined(HPGEM_USE_ANY_PETSC)
    
    GlobalPetscMatrix::GlobalPetscMatrix(Base::MeshManipulatorBase* theMesh, int elementMatrixID, int faceMatrixID)
            : GlobalMatrix(theMesh, elementMatrixID, faceMatrixID)
    {
        logger.assert_debug(theMesh != nullptr, "Invalid mesh passed");
        PetscBool petscRuns;
        PetscInitialized(&petscRuns);
        logger.assert_debug(petscRuns == PETSC_TRUE, "Early call, firstly the command line arguments should be parsed");
        //temporary
        MatCreateSeqAIJ(PETSC_COMM_SELF, 1, 1, 1, PETSC_NULL, &A_);
        
        reAssemble();
    }
    
    GlobalPetscMatrix::~GlobalPetscMatrix()
    {
        int ierr = MatDestroy(&A_);
        //giving error about Petsc has generated inconsistent data and likely memory corruption in heap
        CHKERRV(ierr);
    }
    
    GlobalPetscMatrix::operator Mat()
    {
        if(HPGEM_LOGLEVEL>=Log::DEBUG)
        {
            MatChop(A_, 1e-13);
            MatScale(A_, 24.);
            MatView(A_, PETSC_VIEWER_STDOUT_WORLD);
            MatScale(A_, 1. / 24.);
        }
        return A_;
    }
    
    void GlobalPetscMatrix::reset()
    {
        int ierr = MatZeroEntries(A_);
        CHKERRV(ierr);
        ierr = MatSetOption(A_, MAT_ROW_ORIENTED, PETSC_FALSE);
        CHKERRV(ierr);
        
        LinearAlgebra::MiddleSizeMatrix elementMatrix;
        
        if (elementMatrixID_ >= 0)
        {
            std::vector<PetscInt> localToGlobal;
            for (Base::Element* element : theMesh_->getElementsList())
            {
                indexing_.getGlobalIndices(element, localToGlobal);
                elementMatrix = element->getElementMatrix(elementMatrixID_);
                logger.assert_debug(elementMatrix.getNumberOfRows() == elementMatrix.getNumberOfColumns()
                    && elementMatrix.getNumberOfRows() == localToGlobal.size(),
                    "Incorrect element matrix size");
                logger(DEBUG, "%", elementMatrix * 24.);
                ierr = MatSetValues(A_, localToGlobal.size(), localToGlobal.data(), localToGlobal.size(), localToGlobal.data(), elementMatrix.data(), ADD_VALUES);
                CHKERRV(ierr);
            }
        }
        
        LinearAlgebra::MiddleSizeMatrix faceMatrix;
        
        if (faceMatrixID_ >= 0)
        {
            std::vector<PetscInt> localToGlobal;
            for (Base::Face* face : theMesh_->getFacesList())
            {
                if (!face->isOwnedByCurrentProcessor())
                    continue;

                indexing_.getGlobalIndices(face, localToGlobal);
                faceMatrix = face->getFaceMatrixMatrix(faceMatrixID_);
                logger.assert_debug(faceMatrix.getNumberOfRows() == faceMatrix.getNumberOfColumns()
                                    && faceMatrix.getNumberOfRows() == localToGlobal.size(),
                                    "Incorrect face matrix size");
                logger(DEBUG, "%", faceMatrix * 24.);
                ierr = MatSetValues(A_, localToGlobal.size(), localToGlobal.data(), localToGlobal.size(), localToGlobal.data(), faceMatrix.data(), ADD_VALUES);
                CHKERRV(ierr);
            }
        }

        ierr = MatSetOption(A_, MAT_ROW_ORIENTED, PETSC_TRUE);
        CHKERRV(ierr);

        ierr = MatAssemblyBegin(A_, MAT_FINAL_ASSEMBLY);
        ierr = MatAssemblyEnd(A_, MAT_FINAL_ASSEMBLY);
        
        CHKERRV(ierr);
    }

    /// Compute the global DoFs that have support on an element like. Unlike
    /// getGlobalIndices from GlobalIndexing the dofs are here split into those owned by
    /// the current processor and those that are not.
    void addLocalBasisFunctions(const Base::Element* element, std::size_t dim, const GlobalIndexing& indexing,
            std::set<std::size_t>& owned, std::set<std::size_t>& notOwned)
    {
        logger.assert_debug(element->isOwnedByCurrentProcessor(), "Non owned element");
        const std::size_t unknowns = element->getNumberOfUnknowns();
        {
            // For the element
            std::set<std::size_t> &toChange = element->isOwnedByCurrentProcessor() ? owned : notOwned;
            for (std::size_t unknown = 0; unknown < unknowns; ++unknown)
            {
                const std::size_t offset = indexing.getGlobalIndex(element, unknown);
                const std::size_t localBasisFunctions = element->getLocalNumberOfBasisFunctions(unknown);
                for (std::size_t i = 0; i < localBasisFunctions; ++i)
                {
                    toChange.insert(offset + i);
                }
            }
        }
        // For the faces
        for (const Base::Face* face : element->getFacesList())
        {
            std::set<std::size_t>& toChange = face->isOwnedByCurrentProcessor() ? owned : notOwned;
            for (std::size_t unknown = 0; unknown < unknowns; ++unknown)
            {
                const std::size_t offset = indexing.getGlobalIndex(face, unknown);
                const std::size_t localBasisFunctions = face->getLocalNumberOfBasisFunctions(unknown);
                for (std::size_t i = 0; i < localBasisFunctions; ++i)
                {
                    toChange.insert(offset + i);
                }
            }
        }
        // For the edges
        for (const Base::Edge* edge : element->getEdgesList())
        {
            std::set<std::size_t>& toChange = edge->isOwnedByCurrentProcessor() ? owned : notOwned;
            for (std::size_t unknown = 0; unknown < unknowns; ++unknown)
            {
                const std::size_t offset = indexing.getGlobalIndex(edge, unknown);
                const std::size_t localBasisFunctions = edge->getLocalNumberOfBasisFunctions(unknown);
                for (std::size_t i = 0; i < localBasisFunctions; ++i)
                {
                    toChange.insert(offset + i);
                }
            }
        }
        // For the nodes
        if (dim > 1)
        {
            // For the edges
            for (const Base::Node* node : element->getNodesList())
            {
                std::set<std::size_t>& toChange = node->isOwnedByCurrentProcessor() ? owned : notOwned;
                for (std::size_t unknown = 0; unknown < unknowns; ++unknown)
                {
                    const std::size_t offset = indexing.getGlobalIndex(node, unknown);
                    const std::size_t localBasisFunctions = node->getLocalNumberOfBasisFunctions(unknown);
                    for (std::size_t i = 0; i < localBasisFunctions; ++i)
                    {
                        toChange.insert(offset + i);
                    }
                }
            }
        }
    }

    void computeConnectivity(Base::MeshManipulatorBase& mesh, GlobalIndexing& indexing,
            std::vector<PetscInt>& nonZerosPerRowDiag, std::vector<PetscInt>& nonZerosPerRowOffDiag)
    {
        const std::size_t totalNumberOfDoF = indexing.getNumberOfLocalBasisFunctions();
        const std::size_t nUknowns = indexing.getNumberOfUnknowns();
        // Starting position, all zeros
        nonZerosPerRowDiag.assign(totalNumberOfDoF, 0);
        nonZerosPerRowOffDiag.assign(totalNumberOfDoF, 0);
        // Sets of global indices of owned and unowned DoFs. By using sets these are automatically deduplicated.
        std::set<std::size_t> localDoFs;
        std::set<std::size_t> nonLocalDoFs;

        for (const Base::Element* element : mesh.getElementsList())
        {
            localDoFs.clear();
            nonLocalDoFs.clear();
            // Add local basis functions
            addLocalBasisFunctions(element, mesh.dimension(), indexing, localDoFs, nonLocalDoFs);

            // Face-face coupling
            for(const Base::Face* face : element->getFacesList())
            {
                if (!face->isInternal())
                    continue;
                const Base::Element* other = face->getPtrElementLeft() == element
                        ? face->getPtrElementRight() : face->getPtrElementLeft();
                addLocalBasisFunctions(other, mesh.dimension(), indexing, localDoFs, nonLocalDoFs);
            }

            // With all the global indices counted, allocate them.
            for (std::size_t unknown = 0; unknown < nUknowns; ++unknown)
            {
                std::size_t offset = indexing.getGlobalIndex(element, unknown);
                std::size_t nbasis = element->getLocalNumberOfBasisFunctions(unknown);
                for (std::size_t i = 0; i < nbasis; ++i)
                {
                    nonZerosPerRowDiag[i + offset] = localDoFs.size();
                    nonZerosPerRowOffDiag[i + offset] = nonLocalDoFs.size();
                }
            }
        }
        // Faces
        for (const Base::Face* face : mesh.getFacesList())
        {
            logger.assert_debug(face->isOwnedByCurrentProcessor(), "Non owned face");
            localDoFs.clear();
            nonLocalDoFs.clear();
            std::vector<const Base::Element*> faceElements = {face->getPtrElementLeft()};
            if (face->isInternal())
                faceElements.push_back(face->getPtrElementRight());
            // Compute all DoFs such that either have support on the elements adjacent to
            // the face, or have support on one of the faces of the adjacent elements.
            for (const Base::Element* element : faceElements)
            {
                // Add basis functions from the left/right element
                addLocalBasisFunctions(element, mesh.dimension(), indexing, localDoFs, nonLocalDoFs);
                // Loop over all adjacent faces
                for (const Base::Face* otherFace : element->getFacesList())
                {
                    if (otherFace == face)
                        continue; // Both left & right element of the original face are already visited
                    if (!face->isInternal())
                        continue;
                    const Base::Element* otherElement = face->getPtrElementLeft() == element
                            ? face->getPtrElementRight() : face->getPtrElementLeft();
                    // The basis function from the face may couple through any other face of the elements adjacent to
                    // the face
                    addLocalBasisFunctions(otherElement, mesh.dimension(), indexing, localDoFs, nonLocalDoFs);
                }
            }
            // Store the information
            for (std::size_t unknown = 0; unknown < nUknowns; ++unknown)
            {
                const std::size_t offset = indexing.getGlobalIndex(face, unknown);
                const std::size_t nbasis = face->getLocalNumberOfBasisFunctions(unknown);
                for (std::size_t i = 0; i < nbasis; ++i)
                {
                    nonZerosPerRowDiag[offset + i] = localDoFs.size();
                    nonZerosPerRowOffDiag[offset + i] = nonLocalDoFs.size();
                }
            }
        }
        // Edges
        for (const Base::Edge* edge : mesh.getEdgesList())
        {
            logger.assert_debug(edge->isOwnedByCurrentProcessor(), "Non owned edge");
            localDoFs.clear();
            nonLocalDoFs.clear();
            // Compute all DoFs such that either have support on the elements adjacent to
            // the edge, or have support on one of the faces of the adjacent elements.
            for (const Base::Element* element : edge->getElements())
            {
                // Add basis functions that have support on the adjacent element
                addLocalBasisFunctions(element, mesh.dimension(), indexing, localDoFs, nonLocalDoFs);
                // Loop over the adjacent elements
                for (const Base::Face* face : element->getFacesList())
                {
                    if (!face->isInternal())
                        continue;
                    const Base::Element* otherElement = face->getPtrElementLeft() == element
                            ? face->getPtrElementRight() : face->getPtrElementLeft();
                    addLocalBasisFunctions(otherElement, mesh.dimension(), indexing, localDoFs, nonLocalDoFs);
                }
            }
            // Store the information
            for (std::size_t unknown = 0; unknown < nUknowns; ++unknown)
            {
                const std::size_t offset = indexing.getGlobalIndex(edge, unknown);
                const std::size_t nbasis = edge->getLocalNumberOfBasisFunctions(unknown);
                for (std::size_t i = 0; i < nbasis; ++i)
                {
                    nonZerosPerRowDiag[offset + i] = localDoFs.size();
                    nonZerosPerRowOffDiag[offset + i] = nonLocalDoFs.size();
                }
            }
        }
        // Nodes
        if (mesh.dimension() > 1)
        {
            for (const Base::Node* node : mesh.getNodesList())
            {
                logger.assert_debug(node->isOwnedByCurrentProcessor(), "Non owned node");
                localDoFs.clear();
                nonLocalDoFs.clear();
                // Compute all DoFs such that either have support on the elements adjacent to
                // the face, or have support on one of the faces of the adjacent elements.
                for (const Base::Element* element : node->getElements())
                {
                    addLocalBasisFunctions(element, mesh.dimension(), indexing, localDoFs, nonLocalDoFs);
                    // Loop over the adjacent faces
                    for (const Base::Face* face : element->getFacesList())
                    {
                        if (!face->isInternal())
                            continue;
                        const Base::Element* otherElement = face->getPtrElementLeft() == element
                                                            ? face->getPtrElementRight() : face->getPtrElementLeft();
                        addLocalBasisFunctions(otherElement, mesh.dimension(), indexing, localDoFs, nonLocalDoFs);
                    }
                }
                // Store the information
                for (std::size_t unknown = 0; unknown < nUknowns; ++unknown)
                {
                    const std::size_t offset = indexing.getGlobalIndex(node, unknown);
                    const std::size_t nbasis = node->getLocalNumberOfBasisFunctions(unknown);
                    for (std::size_t i = 0; i < nbasis; ++i)
                    {
                        nonZerosPerRowDiag[offset + i] = localDoFs.size();
                        nonZerosPerRowOffDiag[offset + i] = nonLocalDoFs.size();
                    }
                }
            }
        }
    }

    ///\todo figure out a nice way to keep local data local
    //debug note: GlobalPetscVector 'independently' chooses an ordering for the degrees of freedom, but hpGEM assumes both orderings to be the same
    void GlobalPetscMatrix::reAssemble()
    {
        PetscErrorCode  ierr = MatDestroy(&A_);
        CHKERRV(ierr);

        indexing_.reset(theMesh_, GlobalIndexing::BLOCKED_PROCESSOR);

        const std::size_t totalNumberOfDOF = indexing_.getNumberOfLocalBasisFunctions();
        const std::size_t nUnknowns = indexing_.getNumberOfUnknowns();

        //now construct the only bit of data where PETSc expects a local numbering...
        std::vector<PetscInt> numberOfPositionsPerRow(totalNumberOfDOF, 0);
        std::vector<PetscInt> offDiagonalPositionsPerRow(totalNumberOfDOF, 0);
        computeConnectivity(*theMesh_, indexing_, numberOfPositionsPerRow, offDiagonalPositionsPerRow);
        
        ierr = MatCreateAIJ(PETSC_COMM_WORLD, totalNumberOfDOF, totalNumberOfDOF, PETSC_DETERMINE, PETSC_DETERMINE, -1, numberOfPositionsPerRow.data(), 0, offDiagonalPositionsPerRow.data(), &A_);
        CHKERRV(ierr);
        MatSetOption(A_, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE); //performance
        // While the estimates should be correct, the zeroing of a row without anything on the diagonal will cause an
        // error with the previous option but without the following.
        MatSetOption(A_, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);
        // Enable the following to test that no new allocations were needed
        // MatSetOption(A_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
        ierr = MatSetUp(A_);
        CHKERRV(ierr);
        reset();
    }

    void GlobalPetscMatrix::printMatInfo(MatInfoType type, std::ostream &stream)
    {
        MatInfo info;
        PetscErrorCode error = MatGetInfo(A_, type, &info);
        CHKERRABORT(PETSC_COMM_WORLD, error);
        stream
            << "Blocksize " << info.block_size
            << ", Nonzero " << info.nz_used << " used " << info.nz_unneeded << " unused."
            << " Assembled " << info.assemblies << " mallocs " << info.mallocs<<std::endl;
    }

    void GlobalPetscMatrix::writeMatlab(const std::string &fileName)
    {
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
}

