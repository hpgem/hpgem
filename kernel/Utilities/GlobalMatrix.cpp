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
#include "Base/ElementCacheData.h"
#include "LinearAlgebra/MiddleSizeVector.h"
#include "Base/FaceCacheData.h"
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
        
        for (Base::Element* element : theMesh_->getElementsList())
        {
            std::size_t nElementBasisTotal = element->getTotalNumberOfBasisFunctions();

            for (std::size_t unknown = 0; unknown < nUnknowns; ++unknown)
            {
                // Note, we assume here that the basis functions for a single
                // unknown are laid out consecutively.
                std::size_t localIndex0 = indexing_.getProcessorLocalIndex(element, unknown);
                for (std::size_t basisId = 0; basisId < element->getLocalNumberOfBasisFunctions(unknown); ++basisId)
                {
                    numberOfPositionsPerRow[localIndex0 + basisId] += nElementBasisTotal;
                }
            }

            for (int i = 0; i < element->getReferenceGeometry()->getNumberOfCodim1Entities(); ++i)
            {
                //conforming contributions
                for (std::size_t unknown = 0; unknown < nUnknowns; ++unknown)
                {
                    // Note, we assume here that the basis functions for a single
                    // unknown are laid out consecutively.
                    const Base::Face* face = element->getFace(i);
                    if (!face->isOwnedByCurrentProcessor())
                        continue;
                    std::size_t localIndex0 = indexing_.getProcessorLocalIndex(face, unknown);
                    for (std::size_t basisId = 0; basisId < face->getLocalNumberOfBasisFunctions(unknown); ++basisId)
                    {
                        numberOfPositionsPerRow[localIndex0 + basisId] += nElementBasisTotal;
                    }
                }
            }
            for (int i = 0; i < element->getNumberOfEdges(); ++i)
            {
                for (std::size_t unknown = 0; unknown < nUnknowns; ++unknown)
                {
                    // Note, we assume here that the basis functions for a single
                    // unknown are laid out consecutively.
                    const Base::Edge* edge = element->getEdge(i);
                    if (!edge->isOwnedByCurrentProcessor())
                        continue;
                    int localIndex0 = indexing_.getProcessorLocalIndex(edge, unknown);
                    for (std::size_t basisId = 0; basisId < edge->getLocalNumberOfBasisFunctions(unknown); ++basisId)
                    {
                        numberOfPositionsPerRow[localIndex0 + basisId] += nElementBasisTotal;
                    }
                }
            }
            if (theMesh_->dimension() > 1)
            {
                for (int i = 0; i < element->getNumberOfNodes(); ++i)
                {
                    for (std::size_t unknown = 0; unknown < nUnknowns; ++unknown)
                    {
                        // Note, we assume here that the basis functions for a single
                        // unknown are laid out consecutively.
                        const Base::Node *node = element->getNode(i);
                        if (!node->isOwnedByCurrentProcessor())
                            continue;
                        int localIndex0 = indexing_.getProcessorLocalIndex(node, unknown);
                        for (std::size_t basisId = 0;
                             basisId < node->getLocalNumberOfBasisFunctions(unknown); ++basisId)
                        {
                            numberOfPositionsPerRow[localIndex0 + basisId] += nElementBasisTotal;
                        }
                    }
                }
            }
        }
        
        for(Base::Face* face : theMesh_->getFacesList())
        {
            if(face->isInternal()
                    && (face->getPtrElementLeft()->isOwnedByCurrentProcessor()
                        || face->getPtrElementRight()->isOwnedByCurrentProcessor()
                    ))
            {
                // Choose the correct vector to modify.
                std::vector<int>& changeVec = (face->getFaceType() == Geometry::FaceType::SUBDOMAIN_BOUNDARY || face->getFaceType() == Geometry::FaceType::PERIODIC_SUBDOMAIN_BC)
                        ? offDiagonalPositionsPerRow : numberOfPositionsPerRow;
                std::size_t nDuplicates = 0;
                std::vector<int> duplicates;
                getMatrixBCEntries(face, nDuplicates, duplicates);

                for(Base::Element* element : {face->getPtrElementLeft(), face->getPtrElementRight()})
                {
                    std::size_t nElementBasisTotal = element->getTotalNumberOfBasisFunctions();
                    if (element->isOwnedByCurrentProcessor())
                    {
                        for (std::size_t unknown = 0; unknown < nUnknowns; ++unknown)
                        {
                            // Note, we assume here that the basis functions for a single
                            // unknown are laid out consecutively.
                            std::size_t localIndex0 = indexing_.getProcessorLocalIndex(element, unknown);
                            for (std::size_t basisId = 0;
                                 basisId < element->getLocalNumberOfBasisFunctions(unknown); ++basisId)
                            {
                                changeVec[localIndex0 + basisId] += nElementBasisTotal - nDuplicates;
                            }
                        }
                    }
                    for (int i = 0; i < element->getReferenceGeometry()->getNumberOfCodim1Entities(); ++i)
                    {
                        const Base::Face* face = element->getFace(i);
                        if (face->isOwnedByCurrentProcessor())
                        {
                            //conforming contributions
                            for (std::size_t unknown = 0; unknown < nUnknowns; ++unknown)
                            {
                                // Note, we assume here that the basis functions for a single
                                // unknown are laid out consecutively.
                                int localIndex0 = indexing_.getProcessorLocalIndex(face, unknown);
                                for (std::size_t basisId = 0;
                                     basisId < face->getLocalNumberOfBasisFunctions(unknown); ++basisId)
                                {
                                    changeVec[localIndex0 + basisId] += nElementBasisTotal - nDuplicates;
                                }
                            }
                        }
                    }
                    for (int i = 0; i < element->getNumberOfEdges(); ++i)
                    {
                        const Base::Edge* edge = element->getEdge(i);
                        if (edge->isOwnedByCurrentProcessor())
                        {
                            for (std::size_t unknown = 0; unknown < nUnknowns; ++unknown)
                            {
                                // Note, we assume here that the basis functions for a single
                                // unknown are laid out consecutively.
                                int localIndex0 = indexing_.getProcessorLocalIndex(edge, unknown);
                                for (std::size_t basisId = 0;
                                     basisId < edge->getLocalNumberOfBasisFunctions(unknown); ++basisId)
                                {
                                    changeVec[localIndex0 + basisId] += nElementBasisTotal - nDuplicates;
                                }
                            }
                        }
                    }
                    if (theMesh_->dimension() > 1)
                    {
                        for (int i = 0; i < element->getNumberOfNodes(); ++i)
                        {
                            const Base::Node *node = element->getNode(i);
                            if (node->isOwnedByCurrentProcessor())
                            {
                                for (std::size_t unknown = 0; unknown < nUnknowns; ++unknown)
                                {
                                    // Note, we assume here that the basis functions for a single
                                    // unknown are laid out consecutively.
                                    int localIndex0 = indexing_.getProcessorLocalIndex(node, unknown);
                                    for (std::size_t basisId = 0;
                                         basisId < node->getLocalNumberOfBasisFunctions(unknown); ++basisId)
                                    {
                                        changeVec[localIndex0 + basisId] += nElementBasisTotal - nDuplicates;
                                    }
                                }
                            }
                        }
                    }
                }
                for(int globalIndex : duplicates)
                {
                    int localIndex = indexing_.globalToProcessorLocalIndex(globalIndex);
                    if (localIndex != -1)
                    {
                        changeVec[localIndex] -= face->getPtrElementLeft()->getTotalNumberOfBasisFunctions()
                                - nDuplicates;
                    }
                }
            }
        }

        for (int i = 0; i < totalNumberOfDOF; ++i)
        {
            if (numberOfPositionsPerRow[i] > totalNumberOfDOF)
            {
                numberOfPositionsPerRow[i] = totalNumberOfDOF; //a row cant have more nonzero entries than the number of columns
            }
            if (numberOfPositionsPerRow[i] + offDiagonalPositionsPerRow[i] > totalNumberOfDOF)
            {
                offDiagonalPositionsPerRow[i] = totalNumberOfDOF - numberOfPositionsPerRow[i]; //a row cant have more nonzero entries than the number of columns
            }
        }
        
        ierr = MatCreateAIJ(PETSC_COMM_WORLD, totalNumberOfDOF, totalNumberOfDOF, PETSC_DETERMINE, PETSC_DETERMINE, -1, numberOfPositionsPerRow.data(), 0, offDiagonalPositionsPerRow.data(), &A_);
        CHKERRV(ierr);
        MatSetOption(A_, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE); //performance
        MatSetOption(A_, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE); //the estimate is known to be wrong for mixed element cases and conforming parallel cases
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

