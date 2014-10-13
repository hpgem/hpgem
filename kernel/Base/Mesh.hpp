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
#include "Submesh.hpp"

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

    unsigned int                        getNumberOfElements(unsigned int meshId=0) const {return getElementsList().size(); }
    unsigned int                        getNumberOfFaces(unsigned int meshId=0) const {return getFacesList().size();}
    unsigned int                        getNumberOfEdges(unsigned int meshId=0) const {return getEdgesList().size();}
    unsigned int                        getNumberOfNodes()const {return getNodes().size();}
    
    //! Get const list of elements
    const std::list<Element*>&          getElementsList() const;
    //! Get non-const list of elements
    std::list<Element*>&                getElementsList();

    //! Get const list of faces
    const std::list<Face*>&             getFacesList() const;
    //! Get non-const list of faces
    std::list<Face*>&                   getFacesList();

    const std::list<Edge*>&             getEdgesList()const ;
    std::list<Edge*>&                   getEdgesList();
    
    const std::vector<Geometry::PointPhysical>&  getNodes()const;
    std::vector<Geometry::PointPhysical>&        getNodes();
    
    //********************************************************************************
    std::list<Element*>::const_iterator elementColBegin()const{return getElementsList().begin();}
    std::list<Element*>::const_iterator elementColEnd()const{return getElementsList().end();}

    std::list<Element*>::iterator       elementColBegin(){return getElementsList().begin();}
    std::list<Element*>::iterator       elementColEnd(){return getElementsList().end();}

    std::list<Face*>::const_iterator    faceColBegin()const{return getFacesList().begin();}
    std::list<Face*>::const_iterator    faceColEnd()const{return getFacesList().end();}

    std::list<Face*>::iterator          faceColBegin(){return getFacesList().begin();}
    std::list<Face*>::iterator          faceColEnd(){return getFacesList().end();}

    std::list< Edge*>::const_iterator   edgeColBegin()const{return getEdgesList().begin();}
    std::list< Edge*>::const_iterator   edgeColEnd()const{return getEdgesList().end();}

    std::list< Edge*>::iterator         edgeColBegin(){return getEdgesList().begin();}
    std::list< Edge*>::iterator         edgeColEnd(){return getEdgesList().end();}
    //********************************************************************************
    
private:
    
    //! 'distributes' the mesh across the nodes
    //! this routine assumes all nodes generated the mesh in the same way (so no randomness or thread dependence)
    void                                split();
    
    unsigned int                        localProcessorID_;
    
    Submesh                             submeshes_;

    unsigned int                        elementcounter_;
    unsigned int                        faceCounter_;
    unsigned int                        edgeCounter_;
    //! List of all elements. TODO: this should be replaced by the mesh-tree structure
    std::list<Element*>                 elements_;

    //! List of all faces. TODO: this should be replaced by the mesh-tree structure
    std::list<Face*>                    faces_;

    //! List of all edges.
    std::list< Edge*>                   edges_;

    //! Global vector of physical nodes.
    std::vector<Geometry::PointPhysical>points_;
};

}

#endif	/* MESH_HPP */

