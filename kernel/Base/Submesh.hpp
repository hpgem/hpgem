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
#include <list>

namespace Geometry{
    class PointPhysical;
}

namespace Base{
    
    class Element;
    class Face;
    class Edge;
    class Mesh;

class Submesh {
private:
    //Design note: I dont want mesh messing with private data members, but I only want mesh to access the functionality of mesh
    //Form a meta-physical point of view mesh is not derived from submesh (and it cannot provede the functionality of submesh),
    //so I feel this is the best imperfect solution. If you have a better idea, feel free to change the implementation -FB
    friend Mesh;
    friend std::vector<Submesh>::allocator_type;
    
    Submesh();
    Submesh(const Submesh& orig);
    virtual ~Submesh();
    
    //! adds an element to this submesh
    void add(Element* element);
    
    //! adds a push or pull element. Make sure to add push elements after you fill this list of element belonging to this submesh
    //! processorID is the 0 based index of the processor that will be communicated with about this element
    void addPush(Element* element,int processorID);
    void addPull(Element* element,int processorID);
    
    //! adds a face to this submesh
    //note that interfacial faces should appear in the submesh of both their left and right element
    void add(Face* face);
    
    //! adds an edge to this submesh
    //note that interfacial edges should appear in the submeshes of all their adjacent elements
    void add(Edge* edge);
    
    //! adds a node to this submesh
    //note that interfacial nodes should appear in the submeshes of all their adjacent elements
    void add(Geometry::PointPhysical& node);
    
    //! Get const list of elements
    const std::list<Element*>&          getElementsList() const {return elements_; }
    //! Get non-const list of elements
    std::list<Element*>&                getElementsList() { return elements_; }

    //! Get const list of faces
    const std::list<Face*>&             getFacesList() const { return faces_; }
    //! Get non-const list of faces
    std::list<Face*>&                   getFacesList() { return faces_; }

    const std::list<Edge*>&             getEdgesList() const {return edges_;}
    std::list<Edge*>&                   getEdgesList() {return edges_;}
    
private:
    //! List of all elements. TODO: this should be replaced by the mesh-tree structure
    std::list<Element*>                  elements_;

    //! List of all faces. TODO: this should be replaced by the mesh-tree structure
    //! This contains the list of all faces connected to at least one element in this sub-domain
    std::list<Face*>                     faces_;

    //! List of all edges.
    //! This contains the list of all edges connected to at least one element in this sub-domain
    std::list< Edge*>                    edges_;
    
    //! Tracks the shadow elements (that needs information form another processor each update step, instead of a computation)
    //! pullElements_[i] contains the list of all elements that need info from process i.
    std::vector<std::list<Element*> >    pullElements_;
    
    //! Tracks the shadow elements of other processes (that need their elements send to another processor each update step)
    //! pushElement_[i] constains the list of all elements that send info to process i.
    std::vector<std::list<Element*> >    pushElements_; 
};

}

#endif	/* SUBMESH_HPP */

