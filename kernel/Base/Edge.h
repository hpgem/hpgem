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

#ifndef HPGEM_KERNEL_EDGE_H
#define HPGEM_KERNEL_EDGE_H

#include <vector>
#include <cstdlib>

#include "Logger.h"
#include "TreeEntry.h"

namespace hpgem {

namespace Base {

class Element;

/**
 * generic class that contains entities of codimension 2 or greater that are not
 * vertices. At the moment no integration takes place on edges, so they don't
 * care about their own shape or position. They do know what elements are nearby
 * so they can connect edge-based conforming degrees of freedom to the proper
 * elements. \todo 4D support
 */
class Edge {
   public:
    explicit Edge(std::size_t ID)
        : numberOfConformingDOFOnTheEdge_(std::vector<std::size_t>(0, 0)),
          ID_(ID) {}

    // copy constructor: this is not intended for use and is therefore deleted.
    Edge(const Edge& other) = delete;
    Edge& operator=(const Edge& other) = delete;

    void addElement(Element* element, std::size_t edgeNumber);

    std::size_t getLocalNumberOfBasisFunctions() const {
        std::size_t number = numberOfConformingDOFOnTheEdge_[0];
        for (std::size_t index : numberOfConformingDOFOnTheEdge_)
            logger.assert_debug(
                index == number,
                "number of basis functions is different for different unknown");
        return numberOfConformingDOFOnTheEdge_[0];
    }

    std::size_t getLocalNumberOfBasisFunctions(std::size_t unknown) const {
        // TODO: LC, numberOfConformingDOFOnTheNode_ might be smaller than
        // the number of unknowns (as that is not known here). Thus we might
        // index beyond the number of unknowns.
        if (unknown >= numberOfConformingDOFOnTheEdge_.size()) return 0;
        return numberOfConformingDOFOnTheEdge_[unknown];
    }

    std::size_t getTotalLocalNumberOfBasisFunctions() const {
        std::size_t result = 0;
        for (auto nbasis : numberOfConformingDOFOnTheEdge_) result += nbasis;
        return result;
    }

    ///\deprecated Does not conform naming conventions, use
    /// getLocalNumberOfBasisFunctions instead
    std::size_t getLocalNrOfBasisFunctions() const {
        return getLocalNumberOfBasisFunctions();
    }

    std::size_t getID() const { return ID_; }

    ///\deprecated Does not conform naming conventions, use getNumberOfElements
    /// instead
    std::size_t getNrOfElements() { return getNumberOfElements(); }

    std::size_t getNumberOfElements() const;

    Element* getElement(std::size_t i);

    /// get the root element that is the (indirect) parent of one of the
    /// adjacent elements
    const Element* getRootElement() const;

    std::vector<Element*> getElements();
    const std::vector<Element*> getElements() const;

    ///\deprecated Does not conform naming conventions, use getEdgeNumber
    /// instead
    std::size_t getEdgeNr(std::size_t i) { return getEdgeNumber(i); }

    std::size_t getEdgeNumber(std::size_t i) {
        logger.assert_debug(
            i < getNumberOfElements(),
            "Asked for element %, but there are only % elements", i,
            getNumberOfElements());
        return localEdgeNumbers_[i];
    }

    std::size_t getOrientation(std::size_t i) {
        logger.assert_debug(
            i < getNumberOfElements(),
            "Asked for element %, but there are only % elements", i,
            getNumberOfElements());
        return orientation_[i];
    }

    ///\deprecated Does not conform naming conventions, use
    /// setLocalNumberOfBasisFunctions instead
    void setLocalNrOfBasisFunctions(std::size_t number) {
        setLocalNumberOfBasisFunctions(number);
    }

    void setLocalNumberOfBasisFunctions(std::size_t number) {
        for (std::size_t unknown = 0;
             unknown < numberOfConformingDOFOnTheEdge_.size(); ++unknown) {
            setLocalNumberOfBasisFunctions(number, unknown);
        }
    }

    void setLocalNumberOfBasisFunctions(std::size_t number,
                                        std::size_t unknown) {
        logger.assert_debug(unknown < numberOfConformingDOFOnTheEdge_.size(),
                            "Setting unknown % but there are only %", unknown,
                            numberOfConformingDOFOnTheEdge_.size());
        numberOfConformingDOFOnTheEdge_[unknown] = number;
    }

    void setPositionInTree(const TreeEntry<Edge*>* position) {
        logger.assert_debug(
            position->getData() == this,
            "Trying to set the position of another edge as this edge");
        positionInTheTree_ = position;
    }

    const TreeEntry<Edge*>* getPositionInTree() const {
        return positionInTheTree_;
    }

    bool isOwnedByCurrentProcessor() const;

    /// The element owning this edge, only valid if the edge is owned by the
    /// current processor
    Element* getOwningElement() const;

   private:
    std::vector<Element*> elements_;
    std::vector<std::size_t> localEdgeNumbers_;
    std::vector<std::size_t> orientation_;

    const TreeEntry<Edge*>* positionInTheTree_;

    std::vector<std::size_t> numberOfConformingDOFOnTheEdge_;
    std::size_t ID_;
};

} /* namespace Base */

}  // namespace hpgem

#endif  // HPGEM_KERNEL_EDGE_H
