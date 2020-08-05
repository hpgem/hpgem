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

#include "PhysicalElement.h"
#include "Face.h"
#include "Element.h"

#include "Logger.h"
#include <iostream>
#include "Integration/QuadratureRules/GaussQuadratureRule.h"
#include "Geometry/ReferenceGeometry.h"
#include "Geometry/PointReference.h"
#include "LinearAlgebra/MiddleSizeVector.h"
#include "L2Norm.h"
#include "Node.h"
#include "Geometry/PhysicalGeometry.h"
#include "Geometry/PointPhysical.h"
#include "TreeIterator.h"

namespace hpgem {

namespace Base {

class Face;

/// \details The user does not need to worry about the construction of faces.
/// This is done by mesh-generators. For example the interface HpgemAPIBase can
/// be used to create meshes.
Face::Face(Element* ptrElemL, const std::size_t& localFaceNumberL,
           Element* ptrElemR, const std::size_t& localFaceNumberR,
           std::size_t faceID, std::size_t numberOfFaceMatrixes,
           std::size_t numberOfFaceVectors)
    : FaceGeometry(ptrElemL, localFaceNumberL, ptrElemR, localFaceNumberR),
      FaceData(ptrElemL->getTotalNumberOfBasisFunctions() +
                   ptrElemR->getTotalNumberOfBasisFunctions(),
               numberOfFaceMatrixes, numberOfFaceVectors),
      elementLeft_(ptrElemL),
      elementRight_(ptrElemR),
      numberOfConformingDOFOnTheFace_(std::vector<std::size_t>(1, 0)),
      faceID_(faceID) {
    logger.assert_debug(ptrElemL != nullptr, "Invalid element passed");
    logger.assert_debug(ptrElemR != nullptr,
                        "Error: passing a boundary face to the constructor for "
                        "internal faces!");
    createQuadratureRules();
    ptrElemL->setFace(localFaceNumberL, this);
    ptrElemR->setFace(localFaceNumberR, this);
    numberOfConformingDOFOnTheFace_.resize(ptrElemL->getNumberOfUnknowns(), 0);

    std::vector<std::size_t> leftNodes, rightNodes;
    std::vector<std::size_t> localLeftNodes =
        ptrElemL->getPhysicalGeometry()->getLocalFaceNodeIndices(
            localFaceNumberL);
    std::vector<std::size_t> localRightNodes =
        ptrElemR->getPhysicalGeometry()->getLocalFaceNodeIndices(
            localFaceNumberR);

    {
        // Determine if the face is on a Periodic boundary by looking at the
        // global node indices.
        std::vector<std::size_t> globalLeftNodes =
            ptrElemL->getPhysicalGeometry()->getGlobalFaceNodeIndices(
                localFaceNumberL);
        std::vector<std::size_t> globalRightNodes =
            ptrElemR->getPhysicalGeometry()->getGlobalFaceNodeIndices(
                localFaceNumberR);
        std::sort(globalLeftNodes.begin(), globalLeftNodes.end());
        std::sort(globalRightNodes.begin(), globalRightNodes.end());
        if (globalLeftNodes != globalRightNodes) {
            setFaceType(Geometry::FaceType::PERIODIC_BC);
        }
    }

    for (std::size_t i = 0; i < getReferenceGeometry()->getNumberOfNodes();
         ++i) {
        leftNodes.push_back(ptrElemL->getNode(localLeftNodes[i])->getID());
        rightNodes.push_back(ptrElemR->getNode(localRightNodes[i])->getID());
    }
    initialiseFaceToFaceMapIndex(leftNodes, rightNodes);
}

Face::Face(Element* ptrElemL, const std::size_t& localFaceNumberL,
           const Geometry::FaceType& faceType, std::size_t faceID,
           std::size_t numberOfFaceMatrixes, std::size_t numberOfFaceVectors)
    : FaceGeometry(ptrElemL, localFaceNumberL, faceType),
      FaceData(ptrElemL->getTotalNumberOfBasisFunctions(), numberOfFaceMatrixes,
               numberOfFaceVectors),
      elementLeft_(ptrElemL),
      elementRight_(nullptr),
      numberOfConformingDOFOnTheFace_(std::vector<std::size_t>(1, 0)),
      faceID_(faceID) {
    logger.assert_debug(ptrElemL != nullptr, "Invalid element passed");
    createQuadratureRules();
    ptrElemL->setFace(localFaceNumberL, this);
    numberOfConformingDOFOnTheFace_.resize(ptrElemL->getNumberOfUnknowns(), 0);
}

Face::Face(const Face& other, Element* elementL, const std::size_t localFaceL,
           Element* elementR, const std::size_t localFaceR)
    : FaceGeometry(other, elementL, localFaceL, elementR, localFaceR),
      FaceData(other),
      elementLeft_(elementL),
      elementRight_(elementR),
      quadratureRule_(other.quadratureRule_),
      numberOfConformingDOFOnTheFace_(other.numberOfConformingDOFOnTheFace_),
      faceID_(other.faceID_) {
    logger.assert_debug(elementL != nullptr, "Invalid element passed");
    logger(DEBUG, "Coupling (left) face % to element %", faceID_,
           elementL->getID());
    elementL->setFace(localFaceL, this);
    logger.assert_debug(elementL->getNumberOfFaces() > 0,
                        "Element does not contain any face!");
    if (elementR != nullptr) {
        elementR->setFace(localFaceR, this);
        logger(DEBUG, "Coupling (right) face % to element %", faceID_,
               elementR->getID());
    }
}

void Face::createQuadratureRules() {
    // order of quadrature rules:
    std::size_t rightOrder =
        (elementRight_ == nullptr
             ? 0
             : elementRight_->getGaussQuadratureRule()->order());
    std::size_t leftOrder = elementLeft_->getGaussQuadratureRule()->order();
    if (leftOrder >= rightOrder) {
        quadratureRule_ = elementLeft_->getReferenceGeometry()
                              ->getCodim1ReferenceGeometry(localFaceNumberLeft_)
                              ->getGaussQuadratureRule(leftOrder);
    } else {
        logger(DEBUG, "again..... Face<DIM>::createQuadratureRules(): % %.",
               leftOrder, rightOrder);
        quadratureRule_ =
            elementRight_->getReferenceGeometry()
                ->getCodim1ReferenceGeometry(localFaceNumberRight_)
                ->getGaussQuadratureRule(rightOrder);
    }
}

std::size_t Face::getNumberOfBasisFunctions() const {
    if (isInternal()) {
        return getPtrElementLeft()->getNumberOfBasisFunctions() +
               getPtrElementRight()->getNumberOfBasisFunctions();
    }
    return getPtrElementLeft()->getNumberOfBasisFunctions();
}

std::size_t Face::getNumberOfBasisFunctions(std::size_t unknown) const {
    // Check on unknown delegated to the Elements
    if (isInternal()) {
        return getPtrElementLeft()->getNumberOfBasisFunctions(unknown) +
               getPtrElementRight()->getNumberOfBasisFunctions(unknown);
    }
    return getPtrElementLeft()->getNumberOfBasisFunctions(unknown);
}

/// Get the time integration vectors from both elements and concatenate them.
/// Note that we assume that the data is stored as column "vectors".
LinearAlgebra::MiddleSizeVector Face::getTimeIntegrationVector(
    std::size_t timeIntegrationVectorId, std::size_t unknown) const {
    LinearAlgebra::MiddleSizeVector resLeft =
        getPtrElementLeft()->getTimeIntegrationSubvector(
            timeIntegrationVectorId, unknown);
    if (isInternal()) {
        std::size_t numberOfBasisFunctions = getNumberOfBasisFunctions(unknown);
        std::size_t numberOfBasisFunctionsLeft =
            getPtrElementLeft()->getNumberOfBasisFunctions(unknown);
        resLeft.resize(numberOfBasisFunctions);
        LinearAlgebra::MiddleSizeVector resRight =
            getPtrElementRight()->getTimeIntegrationSubvector(
                timeIntegrationVectorId, unknown);
        for (std::size_t i = numberOfBasisFunctionsLeft;
             i < numberOfBasisFunctions; ++i) {
            resLeft[i] = resRight[i - numberOfBasisFunctionsLeft];
        }
    }
    return resLeft;
}

/// \param[in] side The side of the face.
/// \param[in] varId The index corresponding to the variable.
/// \param[in] scalarBasisFunctionId The index corresponding to the
/// basisfunction.
std::size_t Face::convertToSingleIndex(Side side,
                                       std::size_t scalarBasisFunctionId,
                                       std::size_t varId) const {
    logger.assert_debug(varId < getPtrElementLeft()->getNumberOfUnknowns(),
                        "Asked for unknown %, but there are only % unknowns",
                        varId, getPtrElementLeft()->getNumberOfUnknowns());
    if (side == Side::LEFT) {
        std::size_t number = 0;
        for (std::size_t i = 0; i < varId; ++i) {
            number += getPtrElementLeft()->getNumberOfBasisFunctions(i);
        }
        logger.assert_debug(
            scalarBasisFunctionId <
                getPtrElementLeft()->getNumberOfBasisFunctions(varId),
            "Asked for basis function %, but there are only % basis functions",
            scalarBasisFunctionId,
            getPtrElementLeft()->getNumberOfBasisFunctions(varId));
        return (number + scalarBasisFunctionId);
    }
    logger.assert_debug(isInternal(),
                        "boundary faces only have a \"left\" element");
    std::size_t number = 0;
    for (std::size_t i = 0; i < varId; ++i) {
        number += getPtrElementRight()->getNumberOfBasisFunctions(i);
    }
    logger.assert_debug(
        scalarBasisFunctionId <
            getPtrElementRight()->getNumberOfBasisFunctions(varId),
        "Asked for basis function %, but there are only % basis functions",
        scalarBasisFunctionId,
        getPtrElementRight()->getNumberOfBasisFunctions(varId));
    std::size_t nDOFLeft =
        getPtrElementLeft()->getTotalLocalNumberOfBasisFunctions();
    return nDOFLeft + number + scalarBasisFunctionId;
}

Side Face::getSide(std::size_t faceBasisFunctionId) const {
    std::size_t nDOFLeft =
        getPtrElementLeft()->getTotalNumberOfBasisFunctions();
    if (faceBasisFunctionId < nDOFLeft) {
        return Side::LEFT;
    }
    logger.assert_debug(
        faceBasisFunctionId <
            nDOFLeft +
                (isInternal()
                     ? getPtrElementRight()->getTotalNumberOfBasisFunctions()
                     : 0),
        "The index for the face basis (vector)function (%) is larger than "
        "the number of basis (vector)functions at the adjacent elements "
        "(%)",
        faceBasisFunctionId,
        nDOFLeft + (isInternal()
                        ? getPtrElementRight()->getTotalNumberOfBasisFunctions()
                        : 0));
    return Side::RIGHT;
}

std::size_t Face::getElementBasisFunctionId(
    std::size_t faceBasisFunctionId) const {
    std::size_t nDOFLeft =
        getPtrElementLeft()->getTotalNumberOfBasisFunctions();
    if (faceBasisFunctionId < nDOFLeft) {
        return faceBasisFunctionId;
    }
    logger.assert_debug(
        faceBasisFunctionId <
            nDOFLeft +
                (isInternal()
                     ? getPtrElementRight()->getTotalNumberOfBasisFunctions()
                     : 0),
        "The index for the face basis (vector)function (%) is larger than "
        "the number of basis (vector)functions at the adjacent elements "
        "(%)",
        faceBasisFunctionId,
        nDOFLeft + (isInternal()
                        ? getPtrElementRight()->getTotalNumberOfBasisFunctions()
                        : 0));
    return faceBasisFunctionId - nDOFLeft;
}

void Face::addElement(Element* ptrElementR, std::size_t localFaceNumberR) {
    logger.assert_debug(ptrElementR != nullptr,
                        "Error: passing a boundary face to the constructor for "
                        "internal faces!");
    logger.assert_debug(!isInternal(),
                        "can only add an extra element if there is only one");
    elementRight_ = ptrElementR;
    rightElementGeom_ = ptrElementR;
    localFaceNumberRight_ = localFaceNumberR;
    // deliberatly bypass the check that boundary faces should remain boundary
    // faces since this routine in intended to turn a subdomain boundary into a
    // internal face
    faceType_ = Geometry::FaceType::INTERNAL;
    ptrElementR->setFace(localFaceNumberR, this);

    std::vector<std::size_t> leftNodes, rightNodes;
    std::vector<std::size_t> localLeftNodes =
        elementLeft_->getPhysicalGeometry()->getLocalFaceNodeIndices(
            localFaceNumberLeft());
    std::vector<std::size_t> localRightNodes =
        ptrElementR->getPhysicalGeometry()->getLocalFaceNodeIndices(
            localFaceNumberR);

    {
        // Determine if the face is on a Periodic boundary by looking at the
        // global node indices.
        std::vector<std::size_t> globalLeftNodes =
            elementLeft_->getPhysicalGeometry()->getGlobalFaceNodeIndices(
                localFaceNumberLeft_);
        std::vector<std::size_t> globalRightNodes =
            ptrElementR->getPhysicalGeometry()->getGlobalFaceNodeIndices(
                localFaceNumberR);
        std::sort(globalLeftNodes.begin(), globalLeftNodes.end());
        std::sort(globalRightNodes.begin(), globalRightNodes.end());
        if (globalLeftNodes != globalRightNodes) {
            setFaceType(Geometry::FaceType::PERIODIC_BC);
        }
    }

    for (std::size_t i = 0; i < getReferenceGeometry()->getNumberOfNodes();
         ++i) {
        leftNodes.push_back(elementLeft_->getNode(localLeftNodes[i])->getID());
        rightNodes.push_back(ptrElementR->getNode(localRightNodes[i])->getID());
    }
    initialiseFaceToFaceMapIndex(leftNodes, rightNodes);
}

Element* Face::getRootElement() {
    auto root = elementLeft_->getPositionInTree();
    while (!root->isRoot()) {
        root = root->getParent();
    }
    return *root->getIterator(TreeTraversalMethod::ALLLEVEL);
}

const std::vector<Base::Node*> Face::getNodesList() const {
    std::vector<std::size_t> localFaceNodeIDs =
        elementLeft_->getPhysicalGeometry()->getLocalFaceNodeIndices(
            localFaceNumberLeft_);
    std::size_t nNodesAtFace = localFaceNodeIDs.size();

    std::vector<Base::Node*> ptrNodesAtFace(nNodesAtFace);
    for (std::size_t j = 0; j < nNodesAtFace; j++) {
        std::size_t localNodeIndex = localFaceNodeIDs[j];
        ptrNodesAtFace[j] = elementLeft_->getNode(localNodeIndex);
    }

    return ptrNodesAtFace;
}
}  // namespace Base

}  // namespace hpgem
