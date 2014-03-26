/*
 * Edge.hpp
 *
 *  Created on: Feb 20, 2014
 *      Author: brinkf
 */

#include "Element.hpp"

#ifndef EDGE_HPP_
#define EDGE_HPP_

namespace Base {

/**
 * generic class that contains entities of codimension 2 or greater that are not vertexes.
 * At the moment no integration takes place on edges, so they dont care about their own shape or
 * position. They do know what elements are nearby so they can connent edge-based conforming
 * degrees of freedom to the proper elements.
 * \TODO 4D support
 */
class Edge {
public:
	Edge(std::vector< Element*>& elements, std::vector<unsigned int> localEdgeNrs,unsigned int ID):ID_(ID),elements_(elements),localEdgeNrs_(localEdgeNrs),nrOfConformingDOFOnTheEdge_(0),orientation_(elements_.size())
	{
		std::vector<unsigned int> indices(2);
		for(int i=0;i<elements_.size();++i){
			elements_[i]->setEdge(localEdgeNrs_[i],this);
			elements_[i]->getReferenceGeometry()->getCodim2EntityLocalIndices(localEdgeNrs_[i],indices);
			indices[0]=elements_[i]->getPhysicalGeometry()->getNodeIndex(indices[0]);
			indices[1]=elements_[i]->getPhysicalGeometry()->getNodeIndex(indices[1]);
			orientation_[i]=(indices[0]<indices[1])?0:1;
		}
	}
	virtual ~Edge(){}

	int                             getLocalNrOfBasisFunctions() const{return nrOfConformingDOFOnTheEdge_;}

	int getID()const{return ID_;}

	int getNrOfElements(){return elements_.size();}

	Element* getElement(int i){return elements_[i];}
	unsigned int getEdgeNr(int i){return localEdgeNrs_[i];}
	unsigned int getOrientation(int i){return orientation_[i];}

	void setLocalNrOfBasisFunctions(int number){nrOfConformingDOFOnTheEdge_=number;}

private:

	std::vector< Element*>  elements_;
	std::vector<unsigned int>    localEdgeNrs_;
	std::vector<unsigned int>    orientation_;

    unsigned int  				 nrOfConformingDOFOnTheEdge_;
    int ID_;
};

} /* namespace Base */

#endif /* EDGE_HPP_ */
