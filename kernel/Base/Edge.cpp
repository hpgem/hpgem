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
#include "Edge.hpp"

#include "Element.hpp"
#include "Geometry/ReferenceGeometry.hpp"
#include "Geometry/PhysicalGeometry.hpp"
#include "LinearAlgebra/Matrix.hpp"
#include "LinearAlgebra/NumericalVector.hpp"
#include "ElementCacheData.hpp"
#include "Geometry/PointReference.hpp"
#include "Geometry/PointPhysical.hpp"
#include <vector>

namespace Base {

/*Edge::Edge(std::vector<Element*>& elements,std::vector<std::size_t> localEdgeNrs, std::size_t ID) :
		ID_(ID), elements_(elements), localEdgeNrs_(localEdgeNrs), nrOfConformingDOFOnTheEdge_(0), orientation_(elements_.size()) {
	std::vector<std::size_t> indices(2);
	for (int i = 0; i < elements_.size(); ++i) {
		elements_[i]->setEdge(localEdgeNrs_[i], this);
		elements_[i]->getReferenceGeometry()->getCodim2EntityLocalIndices(localEdgeNrs_[i], indices);
		indices[0] = elements_[i]->getPhysicalGeometry()->getNodeIndex(indices[0]);
		indices[1] = elements_[i]->getPhysicalGeometry()->getNodeIndex(indices[1]);
		orientation_[i] = (indices[0] < indices[1]) ? 0 : 1;
	}
}*/

    void Edge::addElement(Element* element, std::size_t edgeNr)
    {
        elements_.push_back(element);
        localEdgeNrs_.push_back(edgeNr);
        element->setEdge(edgeNr,this);
        std::vector<std::size_t> indices(2);
        element->getReferenceGeometry()->getCodim2EntityLocalIndices(edgeNr, indices);
        indices[0] = element->getPhysicalGeometry()->getNodeIndex(indices[0]);
        indices[1] = element->getPhysicalGeometry()->getNodeIndex(indices[1]);
        orientation_.push_back((indices[0] < indices [1] ? 0 : 1));
    }

	int Edge::getNrOfElements() {
		return elements_.size();
	}

	Element* Edge::getElement(int i) {
		return elements_[i];
	}

}
