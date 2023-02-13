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

#include "Node.h"
#include "Element.h"
#include "Face.h"
#include "LinearAlgebra/MiddleSizeVector.h"
#include <algorithm>

namespace hpgem {

void Base::Node::addElement(Element *element, std::size_t localNodeNumber) {
    logger.assert_debug(std::find(elements_.begin(), elements_.end(),
                                  element) == elements_.end(),
                        "Trying to add the same Element (%) twice", *element);
    elements_.push_back(element);
    localNodeNumbers_.push_back(localNodeNumber);
    element->setNode(localNodeNumber, this);
    if (numberOfConformingDOFOnTheNode_.size() == 0) {
        numberOfConformingDOFOnTheNode_.resize(element->getNumberOfUnknowns(),
                                               0);
    } else {
        logger.assert_debug(numberOfConformingDOFOnTheNode_.size() ==
                                element->getNumberOfUnknowns(),
                            "The element thinks there are % unknowns, but the "
                            "node thinks there are % unknowns",
                            element->getNumberOfUnknowns(),
                            numberOfConformingDOFOnTheNode_.size());
    }
}

Base::Element *Base::Node::getElement(std::size_t i) {
    logger.assert_debug(i < getNumberOfElements(),
                        "asked for element %, but there are only % elements", i,
                        getNumberOfElements());
    return elements_[i];
}

const Base::Element *Base::Node::getElement(std::size_t i) const {
    logger.assert_debug(i < getNumberOfElements(),
                        "asked for element %, but there are only % elements", i,
                        getNumberOfElements());
    return elements_[i];
}

std::size_t Base::Node::getNumberOfElements() const { return elements_.size(); }

const Base::Element *Base::Node::getRootElement() const {
    logger.assert_debug(
        getNumberOfElements() > 0,
        "Add at least one element before queriyng about neighbouring elements");
    auto root = elements_[0]->getPositionInTree();
    while (!root->isRoot()) {
        root = root->getParent();
    }
    return root->getData();
}

const std::vector<Base::Face *> Base::Node::getFaces() const {
    const std::vector<Base::Element *> &ptrElementsAtNode = elements_;
    std::vector<Base::Face *> ptrFacesAtNode;

    std::size_t nElementsAtNode = ptrElementsAtNode.size();
    for (std::size_t i = 0; i < nElementsAtNode; i++) {
        std::size_t localNodeIndex = localNodeNumbers_[i];

        Geometry::PhysicalGeometryBase *ptrPhysicalGeometry =
            ptrElementsAtNode[i]->getPhysicalGeometry();
        std::vector<Base::Face *> ptrFacesAtElement =
            ptrElementsAtNode[i]->getFacesList();
        std::size_t nFacesAtElement = ptrFacesAtElement.size();

        for (std::size_t j = 0; j < nFacesAtElement; j++) {
            bool faceIsInVector = false;
            bool faceIsAtNode = false;
            for (auto &k : ptrFacesAtNode) {
                if (ptrFacesAtElement[j]->getID() == k->getID()) {
                    faceIsInVector = true;
                    break;
                }
            }

            if (!faceIsInVector) {
                std::vector<std::size_t> nodesAtFaceLocalIDs =
                    ptrPhysicalGeometry->getLocalFaceNodeIndices(j);
                std::size_t nNodesAtFace = nodesAtFaceLocalIDs.size();
                for (std::size_t k = 0; k < nNodesAtFace; k++) {
                    if (nodesAtFaceLocalIDs[k] == localNodeIndex) {
                        faceIsAtNode = true;
                    }
                }
                if (faceIsAtNode) {
                    ptrFacesAtNode.push_back(ptrFacesAtElement[j]);
                }
            }
        }
    }

    return ptrFacesAtNode;
}

long long Base::Node::getTopologicalNodeNumber(std::size_t i) const {
    logger.assert_debug(i < getNumberOfElements(),
                        "Asked for element %, but there are only % elements", i,
                        getNumberOfElements());
    return elements_[i]->getReferenceGeometry()->getTopologicalLocalIndex(
        localNodeNumbers_[i]);
}

bool Base::Node::isOwnedByCurrentProcessor() const {
    return elements_.size() > 0 && elements_[0]->isOwnedByCurrentProcessor();
}

Base::Element *Base::Node::getOwningElement() const {
#if HPGEM_ASSERTS
    if (!isOwnedByCurrentProcessor()) {
        // The node is part of the boundary layer of ghost elements. On the
        // outer side of this layer there are nodes where the list of elements
        // is not complete because those elements are outside the ghost layer.
        // One of these missing Elements could be the owner. It is only
        // guaranteed that the list of Elements is correct, when at least one of
        // the elements is owned by the current processor.
        logger.assert_debug(
            std::any_of(elements_.begin(), elements_.end(),
                        [](const Base::Element *el) {
                            return el->isOwnedByCurrentProcessor();
                        }),
            "Owning element of this node may be inaccurate");
    }
#endif
    return elements_.empty() ? nullptr : elements_[0];
}

}  // namespace hpgem
