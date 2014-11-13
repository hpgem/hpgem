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

#include <vector>

#include "Face.hpp"
#include "Submesh.hpp"

namespace Geometry {
    class PointPhysical;
}

namespace Base {
    class Edge;
    class Face;
    class Element;

    ///\brief select if you want to iterate over the local part of the mesh or the whole mesh
    ///\bug mesh information is inherently linked to data and geometric information, but the last two are
    /// only computed/updated locally and not communicated back even when you iterate oven the whole mesh
    enum class IteratorType {
        LOCAL, GLOBAL
    };
    
    class Mesh {
    public:
        Mesh();
        Mesh(const Mesh& orig);
        virtual ~Mesh();

        Element* addElement(const std::vector<unsigned int>& globalNodeIndexes);

        bool addFace(Element* leftElementPtr, unsigned int leftElementLocalFaceNo,
                Element* rightElementPtr, unsigned int rightElementLocalFaceNo,
                const Geometry::FaceType& faceType = Geometry::WALL_BC);

        void addEdge(std::vector< Element*> elements, std::vector<unsigned int> localEdgeNrs);

        void addNode(Geometry::PointPhysical node);

        unsigned int getNumberOfElements(IteratorType part=IteratorType::LOCAL) const {
            return getElementsList(part).size();
        }

        unsigned int getNumberOfFaces(IteratorType part=IteratorType::LOCAL) const {
            return getFacesList(part).size();
        }

        unsigned int getNumberOfEdges(IteratorType part=IteratorType::LOCAL) const {
            return getEdgesList(part).size();
        }

        unsigned int getNumberOfNodes()const {
            return getNodes().size();
        }

        //! Get const list of elements
        const std::vector<Element*>& getElementsList(IteratorType part=IteratorType::LOCAL) const;
        //! Get non-const list of elements
        std::vector<Element*>& getElementsList(IteratorType part=IteratorType::LOCAL);

        //! Get const list of faces
        const std::vector<Face*>& getFacesList(IteratorType part=IteratorType::LOCAL) const;
        //! Get non-const list of faces
        std::vector<Face*>& getFacesList(IteratorType part=IteratorType::LOCAL);

        const std::vector<Edge*>& getEdgesList(IteratorType part=IteratorType::LOCAL)const;
        std::vector<Edge*>& getEdgesList(IteratorType part=IteratorType::LOCAL);

        const std::vector<Geometry::PointPhysical>& getNodes()const;
        std::vector<Geometry::PointPhysical>& getNodes();

        //********************************************************************************

        std::vector<Element*>::const_iterator elementColBegin(IteratorType part=IteratorType::LOCAL)const {
            return getElementsList(part).begin();
        }

        std::vector<Element*>::const_iterator elementColEnd(IteratorType part=IteratorType::LOCAL)const {
            return getElementsList(part).end();
        }

        std::vector<Element*>::iterator elementColBegin(IteratorType part=IteratorType::LOCAL) {
            return getElementsList(part).begin();
        }

        std::vector<Element*>::iterator elementColEnd(IteratorType part=IteratorType::LOCAL) {
            return getElementsList(part).end();
        }

        std::vector<Face*>::const_iterator faceColBegin(IteratorType part=IteratorType::LOCAL)const {
            return getFacesList(part).begin();
        }

        std::vector<Face*>::const_iterator faceColEnd(IteratorType part=IteratorType::LOCAL)const {
            return getFacesList(part).end();
        }

        std::vector<Face*>::iterator faceColBegin(IteratorType part=IteratorType::LOCAL) {
            return getFacesList(part).begin();
        }

        std::vector<Face*>::iterator faceColEnd(IteratorType part=IteratorType::LOCAL) {
            return getFacesList(part).end();
        }

        std::vector< Edge*>::const_iterator edgeColBegin(IteratorType part=IteratorType::LOCAL)const {
            return getEdgesList(part).begin();
        }

        std::vector< Edge*>::const_iterator edgeColEnd(IteratorType part=IteratorType::LOCAL)const {
            return getEdgesList(part).end();
        }

        std::vector< Edge*>::iterator edgeColBegin(IteratorType part=IteratorType::LOCAL) {
            return getEdgesList(part).begin();
        }

        std::vector< Edge*>::iterator edgeColEnd(IteratorType part=IteratorType::LOCAL) {
            return getEdgesList(part).end();
        }
        //********************************************************************************

        Submesh& getSubmesh() {
            return submeshes_;
        }

        const Submesh& getSubmesh() const {
            return submeshes_;
        }

    private:

        //! 'distributes' the mesh across the nodes
        //! this routine assumes all nodes generated the mesh in the same way (so no randomness or thread dependence)
        void split();

        bool hasToSplit_;

        unsigned int localProcessorID_;

        Submesh submeshes_;

        unsigned int elementcounter_;
        unsigned int faceCounter_;
        unsigned int edgeCounter_;
        //! List of all elements. TODO: this should be replaced by the mesh-tree structure
        std::vector<Element*> elements_;

        //! List of all faces. TODO: this should be replaced by the mesh-tree structure
        std::vector<Face*> faces_;

        //! List of all edges.
        std::vector< Edge*> edges_;

        //! Global vector of physical nodes.
        std::vector<Geometry::PointPhysical>points_;
    };

}

#endif	/* MESH_HPP */

