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
#include <cstdlib>

#include "Logger.h"

namespace Base
{
    
    class Element;
    
    /**
     * generic class that contains entities of codimension 2 or greater that are not vertices.
     * At the moment no integration takes place on edges, so they don't care about their own shape or
     * position. They do know what elements are nearby so they can connect edge-based conforming
     * degrees of freedom to the proper elements.
     * \todo 4D support
     */
    class Edge
    {
    public:
        
        explicit Edge(std::size_t ID)
                : numberOfConformingDOFOnTheEdge_(0), ID_(ID)
        {
        }     
                
        //copy constructor: this is not intended for use and is therefore deleted.
        Edge(const Edge &other) = delete;
        Edge& operator=(const Edge &other) = delete;
        
        void addElement(Element* element, std::size_t edgeNr);

        std::size_t getLocalNrOfBasisFunctions() const
        {
            return numberOfConformingDOFOnTheEdge_;
        }
        
        std::size_t getID() const
        {
            return ID_;
        }
        
        std::size_t getNrOfElements();

        Element* getElement(std::size_t i);
        
        std::vector<Element*> getElements();

        std::size_t getEdgeNr(std::size_t i)
        {
            logger.assert(i < getNrOfElements(), "Asked for element %, but there are only % elements", i, getNrOfElements());
            return localEdgeNrs_[i];
        }
        
        std::size_t getOrientation(std::size_t i)
        {
            logger.assert(i < getNrOfElements(), "Asked for element %, but there are only % elements", i, getNrOfElements());
            return orientation_[i];
        }
        
        void setLocalNrOfBasisFunctions(std::size_t number)
        {
            numberOfConformingDOFOnTheEdge_ = number;
        }
        
    private:

        std::vector<Element*> elements_;
        std::vector<std::size_t> localEdgeNrs_;
        std::vector<std::size_t> orientation_;

        std::size_t numberOfConformingDOFOnTheEdge_;
        std::size_t ID_;
    };

} /* namespace Base */

#endif /* EDGE_HPP_ */
