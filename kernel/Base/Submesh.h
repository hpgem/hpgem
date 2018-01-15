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

#ifndef SUBMESH_HPP
#define	SUBMESH_HPP

#include <vector>
#include <map>
#include "LevelTree.h"

namespace Geometry
{
    template<std::size_t DIM>
    class PointPhysical;
}

namespace Base
{
    
    class Element;
    class Face;
    class Edge;
    class Node;
    template<std::size_t DIM>
    class Mesh;
    
    class Submesh
    {
    public:
        /// \brief request to receive data from an element that does not neighbor this partition
        /// \details this partition will receive data about the element on all following communication
        /// elements directly neighboring this partition are required for common flux computations
        /// have already been added by default and don't have to be added again. This method is intended
        /// to facilitate use of non-local information. Passing an element twice or passing an element
        /// that is inside your partition will not have any detrimental effects on communication times.
        /// For example, if you need information from the next-nearest neighbors of your elements
        /// you can simply pass all neighbors of all neighbors of every element in your partition.
        /// Calling this routine intermittently can be expensive, multiple calls should be clustered between
        /// synchronization steps where possible
        void addPullElement(Element* element)
        {
            logger.assert(element!=nullptr, "Invalid element passed");
            if(!std::binary_search(otherPulls_.begin(), otherPulls_.end(), element, [](Element* a, Element* b)->bool{return a->getID() < b->getID();}))
            {
                otherPulls_.push_back(element);
                std::inplace_merge(otherPulls_.begin(), otherPulls_.end() - 1, otherPulls_.end(), [](Element* a, Element* b)->bool{return a->getID() < b->getID();});
            }
        }

        //! Get const list of elements
        const LevelTree<Element*>& getElementsList() const
        {
            return elements_;
        }
        
        //! Get non-const list of elements
        LevelTree<Element*>& getElementsList()
        {
            return elements_;
        }
        
        //! Get const list of faces
        const LevelTree<Face*>& getFacesList() const
        {
            return faces_;
        }
        
        //! Get non-const list of faces
        LevelTree<Face*>& getFacesList()
        {
            return faces_;
        }
        
        const LevelTree<Edge*>& getEdgesList() const
        {
            return edges_;
        }
        
        LevelTree<Edge*>& getEdgesList()
        {
            return edges_;
        }
        
        const std::vector<Node*>& getNodesList() const
        {
            return nodes_;
        }
        
        std::vector<Node*>& getNodesList()
        {
            return nodes_;
        }
        
        const std::map<int, std::vector<Element*> > & getPullElements()
        {
            processPullRequests();
            return pullElements_;
        }
        
        const std::map<int, std::vector<Element*> > & getPushElements()
        {
            processPullRequests();
            return pushElements_;
        }

        //! adds an element to this submesh
        void add(Element* element);

        //! adds a push or pull element. Make sure to add push elements after you fill the list of element belonging to this submesh
        //! processorID is the 0 based index of the processor that will be communicated with about this element
        void addPush(Element* element, int processorID);
        void addPull(Element* element, int processorID);

        //! adds a face to this submesh
        //note that interfacial faces should appear in the submesh of both their left and right element
        void add(Face* face);

        //! adds an edge to this submesh
        //note that interfacial edges should appear in the submeshes of one their adjacent elements
        void add(Edge* edge);

        //! adds a node to this submesh
        //note that interfacial edges should appear in the submeshes of one their adjacent elements
        void add(Node* node);

    private:
        //Design note: Mesh is a friend of this class, because we want mesh to
        //access all functionality of the Submesh, not because of mesh messing
        //with private data members.
        //Form a meta-physical point of view mesh is not derived from submesh
        //(and it cannot provide the functionality of submesh), so currently this
        //is probably the best imperfect solution.
        friend Mesh<1>;
        friend Mesh<2>;
        friend Mesh<3>;
        friend Mesh<4>;
        friend std::vector<Base::Submesh>::allocator_type;

        Submesh() = default;
        Submesh(const Submesh& orig) = delete;

        //! looks up processor IDs for new pull requests and adds them to the lists of push and pull elements
        void processPullRequests();

        //! signals the submesh to prepare for a redistribution (user has to make sure non-geometric data is also redistributed properly)
        void clear();
        
        //! List of all elements.
        LevelTree<Element*> elements_;

        //! List of all faces.
        //! This contains the list of all faces connected to at least one element in this sub-domain
        LevelTree<Face*> faces_;

        //! List of all edges.
        //! This contains the list of all edges connected to at least one element in this sub-domain
        LevelTree<Edge*> edges_;

        //! List of all edges.
        //! This contains the list of all edges connected to at least one element in this sub-domain
        std::vector<Node*> nodes_;

        //! Tracks the shadow elements (that needs information form another processor each update step, instead of a computation)
        //! pullElements_[i] contains the list of all elements that need info from process i.
        //! some entries in this vector are empty lists; pullElements[get_rank()] is guaranteed to be empty
        //map because 'usually' only a few processors neighbour this processor
        std::map<int, std::vector<Element*> > pullElements_;

        //! Tracks the shadow elements of other processes (that need their elements send to another processor each update step)
        //! pushElement_[i] constains the list of all elements that send info to another process
        std::map<int, std::vector<Element*> > pushElements_;

        //! List of elements that should be a shadow element, but are not coupled to a partition yet
        std::vector<Element*> otherPulls_;
    };

}

#endif	/* SUBMESH_HPP */

