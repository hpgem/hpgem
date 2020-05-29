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
#include "GlobalVector.h"
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
#include <Logger.h>
#include <numeric>
#if defined(HPGEM_USE_ANY_PETSC)
#include "petscis.h"
#endif

namespace Utilities {

GlobalVector::GlobalVector(const GlobalIndexing& indexing, int elementVectorID,
                           int faceVectorID)
    : meshLevel_(-2),
      indexing_(indexing),
      elementVectorID_(elementVectorID),
      faceVectorID_(faceVectorID),
      startPositionsOfElementsInTheVector_() {}

#if defined(HPGEM_USE_ANY_PETSC)

GlobalPetscVector::GlobalPetscVector(const GlobalIndexing& indexing,
                                     int elementVectorID, int faceVectorID)
    : GlobalVector(indexing, elementVectorID, faceVectorID) {
    PetscBool petscRuns;
    PetscInitialized(&petscRuns);
    logger.assert_debug(
        petscRuns == PETSC_TRUE,
        "Early call, firstly the command line arguments should be parsed");
    VecCreateSeq(PETSC_COMM_SELF, 0, &b_);

    reinit();
}

GlobalPetscVector::~GlobalPetscVector() {
    int ierr = VecDestroy(&b_);
    CHKERRV(ierr);
}

GlobalPetscVector::operator Vec() {
    if (HPGEM_LOGLEVEL >= Log::DEBUG) {
        VecChop(b_, 1e-13);
        VecScale(b_, 9.);
        VecView(b_, PETSC_VIEWER_STDOUT_WORLD);
        VecScale(b_, 1. / 9.);
    }
    return b_;
}

void GlobalPetscVector::createVec() {
    int ierr = VecDestroy(&b_);
    CHKERRV(ierr);
    ierr = VecCreateMPI(PETSC_COMM_WORLD,
                        indexing_.getNumberOfLocalBasisFunctions(),
                        PETSC_DETERMINE, &b_);
    CHKERRV(ierr);
}

void GlobalPetscVector::zeroVector() {
    int ierr = VecZeroEntries(b_);
    CHKERRV(ierr);
}

void GlobalPetscVector::reinit() {
    createVec();
    assemble();
}

void GlobalPetscVector::assemble() {
    zeroVector();

    const Base::MeshManipulatorBase* mesh = indexing_.getMesh();
    if (mesh == nullptr) {
        // Without any mesh there is nothing to assemble
        return;
    }

    std::vector<PetscInt> elementToGlobal(0);

    if (elementVectorID_ >= 0) {
        for (Base::Element* element : mesh->getElementsList()) {
            indexing_.getGlobalIndices(element, elementToGlobal);
            const LinearAlgebra::MiddleSizeVector& elementVector =
                element->getElementVector(elementVectorID_);
            logger.assert_debug(elementVector.size() == elementToGlobal.size(),
                                "Incorrectly sized element vector.");
            int ierr =
                VecSetValues(b_, elementToGlobal.size(), elementToGlobal.data(),
                             elementVector.data(), ADD_VALUES);
            CHKERRV(ierr);
        }
    }

    LinearAlgebra::MiddleSizeVector faceVector;
    std::vector<PetscInt> faceToGlobal(0);
    if (faceVectorID_ >= 0) {
        for (Base::Face* face : mesh->getFacesList()) {
            if (!face->isOwnedByCurrentProcessor()) continue;

            faceToGlobal.clear();
            indexing_.getGlobalIndices(face, faceToGlobal);
            faceVector = face->getFaceVector(faceVectorID_);
            logger.assert_debug(faceVector.size() == faceToGlobal.size(),
                                "Incorrectly sized face vector");
            int ierr =
                VecSetValues(b_, faceToGlobal.size(), faceToGlobal.data(),
                             faceVector.data(), ADD_VALUES);
            CHKERRV(ierr);
        }
    }

    int ierr = VecAssemblyBegin(b_);
    ierr = VecAssemblyEnd(b_);
    CHKERRV(ierr);
}

void GlobalPetscVector::constructFromTimeIntegrationVector(
    std::size_t timeIntegrationVectorId, std::size_t solutionVar) {
    zeroVector();

    const Base::MeshManipulatorBase* mesh = indexing_.getMesh();
    logger.assert_always(mesh != nullptr, "No mesh to for the GlobalVector");

    LinearAlgebra::MiddleSizeVector elementData;
    std::vector<PetscInt> localToGlobal;
    for (Base::Element* element : mesh->getElementsList()) {
        indexing_.getGlobalIndices(element, localToGlobal);
        elementData.resize(element->getTotalNumberOfBasisFunctions());
        for (std::size_t j = 0; j < element->getNumberOfUnknowns(); ++j) {
            for (std::size_t i = 0; i < element->getNumberOfBasisFunctions(j);
                 ++i) {

                if (j == solutionVar) {
                    elementData[element->convertToSingleIndex(i, solutionVar)] =
                        element->getTimeIntegrationData(timeIntegrationVectorId,
                                                        solutionVar, i);
                } else {
                    localToGlobal[element->convertToSingleIndex(i, j)] = -1;
                }
            }
        }
        int ierr = VecSetValues(b_, localToGlobal.size(), localToGlobal.data(),
                                elementData.data(), INSERT_VALUES);
        CHKERRV(ierr);
    }

    int ierr = VecAssemblyBegin(b_);
    ierr = VecAssemblyEnd(b_);
    CHKERRV(ierr);
}

void GlobalPetscVector::constructFromTimeIntegrationVector(
    std::size_t timeIntegrationVectorId) {
    zeroVector();
    const Base::MeshManipulatorBase* mesh = indexing_.getMesh();
    if (mesh == nullptr) {
        logger(WARN,
               "Construction from time integration vector without a mesh");
        return;
    }

    LinearAlgebra::MiddleSizeVector elementData;
    std::vector<PetscInt> localToGlobal;
    for (Base::Element* element : mesh->getElementsList()) {
        std::size_t numberOfBasisFunctions =
            element->getTotalNumberOfBasisFunctions();
        indexing_.getGlobalIndices(element, localToGlobal);
        int ierr = VecSetValues(
            b_, numberOfBasisFunctions, localToGlobal.data(),
            element->getTimeIntegrationVector(timeIntegrationVectorId).data(),
            INSERT_VALUES);
        CHKERRV(ierr);
    }

    int ierr = VecAssemblyBegin(b_);
    ierr = VecAssemblyEnd(b_);
    CHKERRV(ierr);
}

void GlobalPetscVector::writeTimeIntegrationVector(
    std::size_t timeIntegrationVectorId) {
    PetscScalar* data;

    VecScatter scatter;
    Vec localB;
    // create a local vector...
    VecScatterCreateToAll(b_, &scatter, &localB);
    // we dont need the scatter context, just the vector
    VecScatterDestroy(&scatter);

    const Base::MeshManipulatorBase* mesh = indexing_.getMesh();
    logger.assert_always(mesh != nullptr, "No mesh to for the GlobalVector");

    std::vector<PetscInt> positions;
    std::size_t totalPositions = 0;
    for (Base::Element* element : mesh->getElementsList()) {
        totalPositions += element->getTotalNumberOfBasisFunctions();
    }
    positions.reserve(totalPositions);
    std::vector<PetscInt> newPositions;
    for (Base::Element* element : mesh->getElementsList()) {
        indexing_.getGlobalIndices(element, newPositions);
        for (auto& a : newPositions) {
            positions.push_back(a);
        }
    }

    IS scatterIS;
    ISCreateGeneral(PETSC_COMM_SELF, positions.size(), positions.data(),
                    PETSC_COPY_VALUES, &scatterIS);
    ISSortRemoveDups(scatterIS);
    VecScatterCreate(b_, scatterIS, localB, scatterIS, &scatter);
    VecScatterBegin(scatter, b_, localB, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(scatter, b_, localB, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterDestroy(&scatter);
    ISDestroy(&scatterIS);

    int ierr = VecGetArray(localB, &data);
    CHKERRV(ierr);
    // Question: wat doet deze lijn
    for (Base::MeshManipulatorBase::ConstElementIterator it =
             mesh->elementColBegin();
         it != mesh->elementColEnd(); ++it) {
        // Create vector for the local data that has to be written
        LinearAlgebra::MiddleSizeVector localData(
            (*it)->getTotalNumberOfBasisFunctions());
        std::size_t runningTotal = 0;
        for (std::size_t index = 0; index < (*it)->getNumberOfUnknowns();
             ++index)  // for every variable iV
        {
            std::size_t nElementBasis =
                (*it)->getLocalNumberOfBasisFunctions(index);
            int elementBasis0 = indexing_.getGlobalIndex((*it), index);
            for (std::size_t i = 0; i < nElementBasis;
                 ++i)  // Get the local basis functions of the element
            {
                // Copy the values from data to the localData and update the
                // running number
                localData[runningTotal] = std::real(data[elementBasis0 + i]);
                ++runningTotal;
            }

            for (std::size_t i = 0;
                 i < (*it)->getPhysicalGeometry()->getNumberOfFaces();
                 ++i)  // for all faces of the element
            {
                std::size_t nFaceBasis =
                    (*it)->getFace(i)->getLocalNumberOfBasisFunctions(index);
                int faceBasis0 =
                    indexing_.getGlobalIndex((*it)->getFace(i), index);
                for (std::size_t j = 0; j < nFaceBasis;
                     ++j)  // get local basis functions of a face
                {
                    // Copy the values from data to the localData and update the
                    // running number
                    localData[runningTotal] = std::real(data[faceBasis0 + j]);
                    ++runningTotal;
                }
            }
            for (std::size_t i = 0; i < (*it)->getNumberOfEdges();
                 ++i)  // For all edges of the element
            {
                std::size_t nEdgeBasis =
                    (*it)->getEdge(i)->getLocalNumberOfBasisFunctions(index);
                int edgeBasis0 =
                    indexing_.getGlobalIndex((*it)->getEdge(i), index);
                for (std::size_t j = 0; j < nEdgeBasis;
                     ++j)  // Get the local basis function of an edge
                {
                    // Copy the values from data to the localData and update the
                    // running number

                    localData[runningTotal] = std::real(data[edgeBasis0 + j]);
                    ++runningTotal;
                }
            }
            if (mesh->dimension() > 1)  // There are no nodes in a 1D problem
            {
                for (std::size_t i = 0; i < (*it)->getNumberOfNodes();
                     ++i)  // For all nodes
                {
                    std::size_t nNodeBasis =
                        (*it)->getNode(i)->getLocalNumberOfBasisFunctions(
                            index);
                    int nodeBasis0 =
                        indexing_.getGlobalIndex((*it)->getNode(i), index);
                    for (std::size_t j = 0; j < nNodeBasis;
                         ++j)  // Get the local number of basis function of a
                               // node
                    {
                        // Copy the values from data to the localData and update
                        // the running number
                        localData[runningTotal] =
                            std::real(data[nodeBasis0 + j]);
                        ++runningTotal;
                    }
                }
            }
        }
        // Put the localData in the element
        logger.assert_debug(localData.size() == runningTotal,
                            "not enough info to fill the vector");
        (*it)->setTimeIntegrationVector(timeIntegrationVectorId, localData);
    }
    ierr = VecRestoreArray(localB, &data);
    VecDestroy(&localB);
    CHKERRV(ierr);
}

void GlobalPetscVector::writeTimeIntegrationVector(
    std::size_t timeIntegrationVectorId, std::size_t unknown) {
    PetscScalar* data;

    VecScatter scatter;
    Vec localB;
    // create a local vector...
    VecScatterCreateToAll(b_, &scatter, &localB);
    // we dont need the scatter context, just the vector
    VecScatterDestroy(&scatter);

    const Base::MeshManipulatorBase* mesh = indexing_.getMesh();
    logger.assert_always(mesh != nullptr, "No mesh to for the GlobalVector");

    std::vector<PetscInt> positions;
    std::size_t totalPositions = 0;
    for (Base::Element* element : mesh->getElementsList()) {
        totalPositions += element->getTotalNumberOfBasisFunctions();
    }
    positions.reserve(totalPositions);
    std::vector<PetscInt> localPositions;
    for (Base::Element* element : mesh->getElementsList()) {
        indexing_.getGlobalIndices(element, localPositions);
        for (auto& a : localPositions) {
            positions.push_back(a);
        }
    }

    IS scatterIS;
    ISCreateGeneral(PETSC_COMM_SELF, positions.size(), positions.data(),
                    PETSC_COPY_VALUES, &scatterIS);
    ISSortRemoveDups(scatterIS);
    VecScatterCreate(b_, scatterIS, localB, scatterIS, &scatter);
    VecScatterBegin(scatter, b_, localB, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(scatter, b_, localB, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterDestroy(&scatter);
    ISDestroy(&scatterIS);

    int ierr = VecGetArray(localB, &data);
    CHKERRV(ierr);
    for (Base::MeshManipulatorBase::ConstElementIterator it =
             mesh->elementColBegin();
         it != mesh->elementColEnd(); ++it) {
        LinearAlgebra::MiddleSizeVector localData(
            (*it)->getTotalNumberOfBasisFunctions());
        std::size_t runningTotal = 0;
        for (std::size_t index = 0; index < (*it)->getNumberOfUnknowns();
             ++index) {
            std::size_t nElementBasis =
                (*it)->getLocalNumberOfBasisFunctions(index);
            int elementBasis0 = indexing_.getGlobalIndex((*it), index);
            for (std::size_t i = 0; i < nElementBasis; ++i) {
                localData[runningTotal] = std::real(data[elementBasis0 + i]);
                ++runningTotal;
            }

            for (std::size_t i = 0;
                 i < (*it)->getPhysicalGeometry()->getNumberOfFaces(); ++i) {
                std::size_t nFaceBasis =
                    (*it)->getFace(i)->getLocalNumberOfBasisFunctions(index);
                int faceBasis0 =
                    indexing_.getGlobalIndex((*it)->getFace(i), index);
                for (std::size_t j = 0; j < nFaceBasis; ++j) {
                    localData[runningTotal] = std::real(data[faceBasis0 + j]);
                    ++runningTotal;
                }
            }
            for (std::size_t i = 0; i < (*it)->getNumberOfEdges(); ++i) {
                std::size_t nEdgeBasis =
                    (*it)->getEdge(i)->getLocalNumberOfBasisFunctions(index);
                int edgeBasis0 =
                    indexing_.getGlobalIndex((*it)->getEdge(i), index);
                for (std::size_t j = 0; j < nEdgeBasis; ++j) {
                    localData[runningTotal] = std::real(data[edgeBasis0 + j]);
                    ++runningTotal;
                }
            }
            if (mesh->dimension() > 1) {
                for (std::size_t i = 0; i < (*it)->getNumberOfNodes(); ++i) {
                    std::size_t nNodeBasis =
                        (*it)->getNode(i)->getLocalNumberOfBasisFunctions(
                            index);
                    int nodeBasis0 =
                        indexing_.getGlobalIndex((*it)->getNode(i), index);
                    for (std::size_t j = 0; j < nNodeBasis; ++j) {
                        localData[runningTotal] =
                            std::real(data[nodeBasis0 + j]);
                        ++runningTotal;
                    }
                }
            }
        }
        logger.assert_debug(localData.size() == runningTotal,
                            "not enough info to fill the vector");
        LinearAlgebra::MiddleSizeVector singleUnknownData(
            (*it)->getNumberOfBasisFunctions());
        for (std::size_t i = 0; i < (*it)->getNumberOfBasisFunctions(); ++i) {
            singleUnknownData[i] =
                localData[(*it)->convertToSingleIndex(i, unknown)];
        }
        (*it)->setTimeIntegrationSubvector(timeIntegrationVectorId, unknown,
                                           singleUnknownData);
    }
    ierr = VecRestoreArray(localB, &data);
    VecDestroy(&localB);
    CHKERRV(ierr);
}
#endif

}  // namespace Utilities
