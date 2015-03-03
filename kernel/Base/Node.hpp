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

namespace Base
{

    class Element;

    ///\brief an identification token for vertices that is more likely to be the same when it should be then a PointPhysical

    class Node
    {
    public:

        explicit Node(std::size_t ID) : nrOfConformingDOFOnTheNode_(0), ID_(ID), elements_(), localNodeNrs_() { }

        virtual ~Node() { }

        void addElement(Element* element, std::size_t localNodeNr);

        std::size_t getLocalNrOfBasisFunctions() const
        {
            return nrOfConformingDOFOnTheNode_;
        }

        std::size_t getID()const
        {
            return ID_;
        }

        std::size_t getNrOfElements() const;

        Element* getElement(std::size_t i);
        const Element* getElement(std::size_t i) const;

        const std::vector<Element*> getElements() const
        {
            return elements_;
        }

        std::size_t getVertexNr(std::size_t i) const
        {
            return localNodeNrs_[i];
        }

        void setLocalNrOfBasisFunctions(std::size_t number)
        {
            nrOfConformingDOFOnTheNode_ = number;
        }
    private:

        //provide information to map back to a unique corner of the element
        std::vector<Element*> elements_;
        std::vector<std::size_t> localNodeNrs_;

        //number of basis-functions that are accosiated to this node (most likely 1(conforming) or 0(DG))
        std::size_t nrOfConformingDOFOnTheNode_;
        std::size_t ID_;
    };

}
#endif	/* NODE_HPP */

