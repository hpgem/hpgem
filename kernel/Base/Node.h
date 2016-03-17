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

#ifndef NODE_HPP
#define	NODE_HPP

#include <cstdlib>
#include <vector>

#include "Logger.h"

namespace Base
{
    
    class Element;
    
    /// \brief an identification token for vertices that is more likely to be the same when it should be then a PointPhysical
    /// \details Node is oblivious of its physical location, but gets its relevance from the fact that it can tell what elements
    /// are connecting in this node. This information cannot be inferred from the physical coordinates of the node. In particular
    /// this implementation can deal with arbitrary periodic connectivities along the physical 'boundary' of the domain. In that
    /// situation the location of a node is not uniquely defined and its only identifying feature is the set of elements connected
    /// to it. Elements store a PointPhysical for each node independently from this class that can be used for geometric operations.
    class Node
    {
    public:
        
        explicit Node(std::size_t ID)
            : elements_(), localNodeNumbers_(), numberOfConformingDOFOnTheNode_(0), ID_(ID) { }

        //Since individual parts of the mesh should not be copied and there is no
        //need for a copy constructor of Node while copying a whole mesh, the copy
        //constructor of Node is deleted.
        Node(const Node &other) = delete;

        void addElement(Element* element, std::size_t localNodeNumber);

        ///\deprecated Does not conform naming conventions, use getLocalNumberOfBasisFunctions instead
        std::size_t getLocalNrOfBasisFunctions() const
        {
            return getLocalNumberOfBasisFunctions();
        }
        
        std::size_t getLocalNumberOfBasisFunctions() const
        {
            return numberOfConformingDOFOnTheNode_;
        }
        
        std::size_t getID() const
        {
            return ID_;
        }
        
        ///\deprecated Does not conform naming conventions, use getNumberOfElements instead
        std::size_t getNrOfElements() const
        {
            return getNumberOfElements();
        }
        
        std::size_t getNumberOfElements() const;

        Element* getElement(std::size_t i);
        const Element* getElement(std::size_t i) const;

        ///get the root element that is the (indirect) parent of one of the adjacent elements
        const Element* getRootElement() const;

        const std::vector<Element*>& getElements() const
        {
            return elements_;
        }
        
        ///\deprecated Does not conform naming conventions, use getNodeNumber instead
        std::size_t getNodeNr(std::size_t i) const
        {
            return getNodeNumber(i);
        }
        
        ///Get the local node number for this node at the given element.
        std::size_t getNodeNumber(std::size_t i) const
        {
            logger.assert(i < getNumberOfElements(), "Asked for element %, but there are only % elements", i, getNumberOfElements());
            return localNodeNumbers_[i];
        }
        
        void setLocalNumberOfBasisFunctions(std::size_t number)
        {
            numberOfConformingDOFOnTheNode_ = number;
        }

        void setLocalNrOfBasisFunctions(std::size_t number)
        {
            setLocalNumberOfBasisFunctions(number);
        }

    private:
        
        //provide information to map back to a unique corner of the element
        std::vector<Element*> elements_;
        std::vector<std::size_t> localNodeNumbers_;

        //number of basis-functions that are associated to this node (most likely 1(conforming) or 0(DG))
        std::size_t numberOfConformingDOFOnTheNode_;
        std::size_t ID_;
    };

}
#endif	/* NODE_HPP */

