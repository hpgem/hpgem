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

#include "HpgemAPIBase.h"

#include "Base/Element.h"

#include "Geometry/PointPhysical.h"
#include "LinearAlgebra/MiddleSizeVector.h"
#include "ConfigurationData.h"
namespace hpgem {
namespace Base {
template <std::size_t DIM>
HpgemAPIBase<DIM>::HpgemAPIBase(GlobalData* const global,
                                const ConfigurationData* config)
    : meshes_(), globalData_(global), configData_(config) {
    if (!parse_isDone()) {
        logger(WARN,
               "Warning: Command line arguments have not been parsed.\n"
               "  Please call Base::parse_options(argc, argv); first.\n"
               "  This application may not behave as intended.");
    }
}

// Destructor, destructs the meshes, configData_ and globalData_
template <std::size_t DIM>
HpgemAPIBase<DIM>::~HpgemAPIBase() {
    for (std::size_t i = 0; i < meshes_.size(); ++i) delete meshes_[i];
    delete configData_;
    delete globalData_;
}

template <std::size_t DIM>
bool HpgemAPIBase<DIM>::initialiseMeshMover(
    const MeshMoverBase<DIM>* meshMoverBase, std::size_t meshID) {
    meshes_[meshID]->setMeshMover(meshMoverBase);
    return true;
}

template <std::size_t DIM>
std::size_t HpgemAPIBase<DIM>::addMesh(const std::string& fileName,
                                       std::size_t numberOfElementMatrixes,
                                       std::size_t numberOfElementVectors,
                                       std::size_t numberOfFaceMatrixes,
                                       std::size_t numberOfFaceVectors) {
    std::size_t numberOfMeshes = meshes_.size();
    MeshManipulator<DIM>* mesh = new MeshManipulator<DIM>(
        configData_, numberOfElementMatrixes, numberOfElementVectors,
        numberOfFaceMatrixes, numberOfFaceVectors);
    mesh->readMesh(fileName);
    mesh->getElementsList();
    meshes_.push_back(mesh);
    logger(INFO, "HpgemAPIBase::addMesh read a mesh.");
    return numberOfMeshes;
}

template <std::size_t DIM>
void HpgemAPIBase<DIM>::synchronize(const std::size_t timeIntegrationVectorId) {
#ifdef HPGEM_USE_MPI

    // Now, set it up.
    Base::MeshManipulator<DIM>* meshManipulator = this->meshes_[0];
    Base::Submesh& mesh = meshManipulator->getMesh().getSubmesh();

    const auto& pushes = mesh.getPushElements();
    const auto& pulls = mesh.getPullElements();

    // receive first for lower overhead
    for (const auto& it : pulls) {
        for (Base::Element* ptrElement : it.second) {
            logger.assert_debug(
                ptrElement->getTimeIntegrationVector(timeIntegrationVectorId)
                        .size() == ptrElement->getTotalNumberOfBasisFunctions(),
                "Size of time integration vector % is wrong: % instead of %.",
                timeIntegrationVectorId,
                ptrElement->getTimeIntegrationVector(timeIntegrationVectorId)
                    .size(),
                ptrElement->getTotalNumberOfBasisFunctions());

            Base::MPIContainer::Instance().receive(
                ptrElement->getTimeIntegrationVector(timeIntegrationVectorId),
                it.first, ptrElement->getID());
        }
    }
    for (const auto& it : pushes) {
        for (Base::Element* ptrElement : it.second) {
            logger.assert_debug(
                ptrElement->getTimeIntegrationVector(timeIntegrationVectorId)
                        .size() == ptrElement->getTotalNumberOfBasisFunctions(),
                "Size of time integration vector % is wrong: % instead of %.",
                timeIntegrationVectorId,
                ptrElement->getTimeIntegrationVector(timeIntegrationVectorId)
                    .size(),
                ptrElement->getTotalNumberOfBasisFunctions());

            Base::MPIContainer::Instance().send(
                ptrElement->getTimeIntegrationVector(timeIntegrationVectorId),
                it.first, ptrElement->getID());
        }
    }
    Base::MPIContainer::Instance().sync();

#endif
}

template <std::size_t DIM>
void HpgemAPIBase<DIM>::setNumberOfTimeIntegrationVectorsGlobally(
    std::size_t numberOfTimeIntegrationVectors) {
    if (HpgemAPIBase<DIM>::meshes_.size() == 0) {
        logger(ERROR,
               "Error no mesh created : You need to create at least one mesh "
               "to set the number of time integration vectors");
    }
    for (auto mesh : this->meshes_) {
        for (Base::Element* ptrElement :
             mesh->getElementsList(IteratorType::GLOBAL)) {
            ptrElement->setNumberOfTimeIntegrationVectors(
                numberOfTimeIntegrationVectors);
        }
    }
}

template <std::size_t DIM>
void HpgemAPIBase<DIM>::copyTimeIntegrationToTimeLevelData(
    std::size_t timeIntegrationVectorId, std::size_t timeLevel,
    std::size_t meshId) {
    for (Base::Element* ptrElement : this->meshes_[meshId]->getElementsList()) {
        ptrElement->getTimeLevelDataVector(timeLevel) =
            ptrElement->getTimeIntegrationVector(timeIntegrationVectorId);
    }
}

template <std::size_t DIM>
void HpgemAPIBase<DIM>::copyTimeLevelToTimeIntegrationData(
    std::size_t timeLevel, std::size_t timeIntegrationVectorId,
    std::size_t meshId) {
    for (Base::Element* ptrElement : this->meshes_[meshId]->getElementsList()) {
        ptrElement->getTimeIntegrationVector(timeIntegrationVectorId) =
            ptrElement->getTimeLevelDataVector(timeLevel);
    }
}

template <std::size_t DIM>
typename HpgemAPIBase<DIM>::ConstElementIterator
    HpgemAPIBase<DIM>::elementColBegin(std::size_t mId) const {
    return meshes_[mId]->elementColBegin();
}

template <std::size_t DIM>
typename HpgemAPIBase<DIM>::ConstElementIterator
    HpgemAPIBase<DIM>::elementColEnd(std::size_t mId) const {
    return meshes_[mId]->elementColEnd();
}

template <std::size_t DIM>
typename HpgemAPIBase<DIM>::ElementIterator HpgemAPIBase<DIM>::elementColBegin(
    std::size_t mId) {
    return meshes_[mId]->elementColBegin();
}

template <std::size_t DIM>
typename HpgemAPIBase<DIM>::ElementIterator HpgemAPIBase<DIM>::elementColEnd(
    std::size_t mId) {
    return meshes_[mId]->elementColEnd();
}

template <std::size_t DIM>
typename HpgemAPIBase<DIM>::ConstFaceIterator HpgemAPIBase<DIM>::faceColBegin(
    std::size_t mId) const {
    return meshes_[mId]->faceColBegin();
}

template <std::size_t DIM>
typename HpgemAPIBase<DIM>::ConstFaceIterator HpgemAPIBase<DIM>::faceColEnd(
    std::size_t mId) const {
    return meshes_[mId]->faceColEnd();
}

template <std::size_t DIM>
typename HpgemAPIBase<DIM>::FaceIterator HpgemAPIBase<DIM>::faceColBegin(
    std::size_t mId) {
    return meshes_[mId]->faceColBegin();
}

template <std::size_t DIM>
typename HpgemAPIBase<DIM>::FaceIterator HpgemAPIBase<DIM>::faceColEnd(
    std::size_t mId) {
    return meshes_[mId]->faceColEnd();
}

}  // namespace Base
}  // namespace hpgem