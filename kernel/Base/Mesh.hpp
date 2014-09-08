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

#ifndef MESH_HPP
#define	MESH_HPP

#include <list>
#include <vector>

#include "Face.hpp"

namespace Geometry{
    class PointPhysical;
}

namespace Base{
    class Edge;
    class Face;
    class Element;

class Mesh {
public:
    Mesh();
    Mesh(const Mesh& orig);
    virtual ~Mesh();

    Element*                            addElement(const std::vector<unsigned int>& globalNodeIndexes);

    bool                                addFace(Element* leftElementPtr, unsigned int leftElementLocalFaceNo, 
                                                Element* rightElementPtr, unsigned int rightElementLocalFaceNo,
                                                const Geometry::FaceType& faceType=Geometry::WALL_BC);

    void                                addEdge(std::vector< Element*> elements, std::vector<unsigned int> localEdgeNrs);
    
    void                                addNode(Geometry::PointPhysical node);

    unsigned int                        getNumberOfElements(unsigned int meshId=0) const {return elements_.size(); }
    unsigned int                        getNumberOfFaces(unsigned int meshId=0) const {return faces_.size();}
    unsigned int                        getNumberOfEdges(unsigned int meshId=0) const {return edges_.size();}
    unsigned int                        getNumberOfNodes()const {return points_.size();}
    
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
    
    const std::vector<Geometry::PointPhysical>&  getNodes()const{return points_;}
    std::vector<Geometry::PointPhysical>&        getNodes(){return points_;}
    
    //********************************************************************************
    std::list<Element*>::const_iterator elementColBegin()const{return elements_.begin();}
    std::list<Element*>::const_iterator elementColEnd()const{return elements_.end();}

    std::list<Element*>::iterator       elementColBegin(){return elements_.begin();}
    std::list<Element*>::iterator       elementColEnd(){return elements_.end();}

    std::list<Face*>::const_iterator    faceColBegin()const{return faces_.begin();}
    std::list<Face*>::const_iterator    faceColEnd()const{return faces_.end();}

    std::list<Face*>::iterator          faceColBegin(){return faces_.begin();}
    std::list<Face*>::iterator          faceColEnd(){return faces_.end();}

    std::list< Edge*>::const_iterator   edgeColBegin()const{return edges_.begin();}
    std::list< Edge*>::const_iterator   edgeColEnd()const{return edges_.end();}

    std::list< Edge*>::iterator         edgeColBegin(){return edges_.begin();}
    std::list< Edge*>::iterator         edgeColEnd(){return edges_.end();}
    //********************************************************************************
    
private:

    unsigned int                         elementcounter_;
    unsigned int                         faceCounter_;
    unsigned int                         edgeCounter_;
    //! List of all elements. TODO: this should be replaced by the mesh-tree structure
    std::list<Element*>                  elements_;

    //! List of all faces. TODO: this should be replaced by the mesh-tree structure
    std::list<Face*>                     faces_;

    //! List of all edges.
    std::list< Edge*>                    edges_;

    //! Global vector of physical nodes.
    std::vector<Geometry::PointPhysical> points_;
};

}

#endif	/* MESH_HPP */

