/*
 * Edge.cpp
 *
 *  Created on: Feb 20, 2014
 *      Author: brinkf
 */

#include <Edge.hpp>

namespace Base {

Edge::Edge(std::vector<Element*>& elements,std::vector<unsigned int> localEdgeNrs, unsigned int ID) :
		ID_(ID), elements_(elements), localEdgeNrs_(localEdgeNrs), nrOfConformingDOFOnTheEdge_(0), orientation_(elements_.size()) {
	std::vector<unsigned int> indices(2);
	for (int i = 0; i < elements_.size(); ++i) {
		elements_[i]->setEdge(localEdgeNrs_[i], this);
		elements_[i]->getReferenceGeometry()->getCodim2EntityLocalIndices(localEdgeNrs_[i], indices);
		indices[0] = elements_[i]->getPhysicalGeometry()->getNodeIndex(indices[0]);
		indices[1] = elements_[i]->getPhysicalGeometry()->getNodeIndex(indices[1]);
		orientation_[i] = (indices[0] < indices[1]) ? 0 : 1;
	}
}

}
