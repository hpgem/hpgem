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

#ifndef EDGE_HPP_
#define EDGE_HPP_

#include <vector>

namespace Base {

	class Element;

/**
 * generic class that contains entities of codimension 2 or greater that are not vertexes.
 * At the moment no integration takes place on edges, so they dont care about their own shape or
 * position. They do know what elements are nearby so they can connent edge-based conforming
 * degrees of freedom to the proper elements.
 * \TODO 4D support
 */
class Edge {
public:
	Edge(std::vector<Element*>& elements,std::vector<unsigned int> localEdgeNrs, unsigned int ID);
	virtual ~Edge(){}

	int getLocalNrOfBasisFunctions() const{return nrOfConformingDOFOnTheEdge_;}

	int getID()const{return ID_;}

		int getNrOfElements();

		Element* getElement(int i);
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
