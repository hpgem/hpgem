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
#include "Edge.h"

#include "Element.h"
#include "Geometry/ReferenceGeometry.h"
#include "Geometry/PhysicalGeometry.h"
#include "LinearAlgebra/MiddleSizeMatrix.h"
#include "LinearAlgebra/MiddleSizeVector.h"
#include "ElementCacheData.h"
#include "Geometry/PointReference.h"
#include "Geometry/PointPhysical.h"
#include <vector>

namespace Base
{

    void Edge::addElement(Element* element, std::size_t edgeNumber)
    {
        logger.assert_debug(element != nullptr, "Invalid element detected");
        elements_.push_back(element);
        localEdgeNumbers_.push_back(edgeNumber);
        element->setEdge(edgeNumber, this);
        std::vector<std::size_t> indices = element->getReferenceGeometry()->getCodim2EntityLocalIndices(edgeNumber);
        indices[0] = element->getPhysicalGeometry()->getNodeIndex(indices[0]);
        indices[1] = element->getPhysicalGeometry()->getNodeIndex(indices[1]);
        orientation_.push_back((indices[0] < indices[1] ? 0 : 1));
        
        if(numberOfConformingDOFOnTheEdge_.size() == 0)
        {
            numberOfConformingDOFOnTheEdge_.resize(element->getNumberOfUnknowns(), 0);
        }
        else
        {
            logger.assert_debug(numberOfConformingDOFOnTheEdge_.size() == element->getNumberOfUnknowns(),
                                "The element thinks there are % unknowns, but the edge thinks there are % unknowns", element->getNumberOfUnknowns(),
                                numberOfConformingDOFOnTheEdge_.size());
        }
    }
    
    std::size_t Edge::getNumberOfElements() const
    {
        return elements_.size();
    }
    
    Element* Edge::getElement(std::size_t i)
    {
        logger.assert_debug(i < getNumberOfElements(), "Asked for element %, but there are only % elements", i, getNumberOfElements());
        return elements_[i];
    }
    
    std::vector<Element*> Edge::getElements()
    {
        return elements_;
    }

    const Element *Edge::getRootElement() const
    {
        logger.assert_debug(getNumberOfElements() > 0, "Add at least one element before queriyng about neighbouring elements");
        auto root = elements_[0]->getPositionInTree();
        while(!root->isRoot())
        {
            root = root->getParent();
        }
        return root->getData();
    }

}
