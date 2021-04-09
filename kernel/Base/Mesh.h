/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2014, University of Twente
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.

 3. Neither the name of the copyright holder nor the names of its contributors
 may be used to endorse or promote products derived from this software without
 specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef HPGEM_KERNEL_MESH_H
#define HPGEM_KERNEL_MESH_H
#ifdef HPGEM_USE_MPI
#include <mpi.h>
#endif
#include <vector>

#include "Face.h"
#include "Submesh.h"
#include "Geometry/PointPhysical.h"
#include "Node.h"
#include "LevelTree.h"
#include "Zone.h"

namespace hpgem {
namespace Geometry {
template <std::size_t DIM>
class PointPhysical;
}

namespace Base {
class Edge;
class Face;
class Element;

///\brief select if you want to iterate over the local part of the mesh or the
/// whole mesh \bug mesh information is inherently linked to data and geometric
/// information, but the last two are
/// only computed/updated locally and not communicated back even when you
/// iterate oven the whole mesh
enum class IteratorType { LOCAL, GLOBAL };

// class is made final so we don't have to create a v-table specifically for the
// destructor
template <std::size_t DIM>
class Mesh final {
   public:
    Mesh();

    Mesh(const Mesh& orig);

    ~Mesh();

    Element* addElement(const std::vector<std::size_t>& globalNodeIndexes,
                        std::size_t zoneId, std::size_t owner, bool owning);

    void addSubElements(Base::Element* parent,
                        const std::vector<Base::Element*> subElements);

    bool addFace(
        Element* leftElementPtr, std::size_t leftElementLocalFaceNo,
        Element* rightElementPtr, std::size_t rightElementLocalFaceNo,
        const Geometry::FaceType& faceType = Geometry::FaceType::WALL_BC);

    void addSubFaces(const Base::Face* parent,
                     const std::vector<Base::Face*> subFaces);

    Edge* addEdge();

    void addNodeCoordinate(Geometry::PointPhysical<DIM> node);

    Node* addNode();

    void clear();

    std::size_t getNumberOfElements(
        IteratorType part = IteratorType::LOCAL) const {
        return getElementsList(part).size();
    }

    std::size_t getNumberOfFaces(
        IteratorType part = IteratorType::LOCAL) const {
        return getFacesList(part).size();
    }

    std::size_t getNumberOfEdges(
        IteratorType part = IteratorType::LOCAL) const {
        return getEdgesList(part).size();
    }

    std::size_t getNumberOfNodes(
        IteratorType part = IteratorType::LOCAL) const {
        return getNodesList(part).size();
    }

    std::size_t getNumberOfNodeCoordinates() const {
        return getNodeCoordinates().size();
    }

    /// Get a vector of elements. If the IteratorType is LOCAL, get all elements
    /// on the core you're working on. If the IteratorType is GLOBAL, get all
    /// elements in the mesh. Usually an application uses only the local
    /// elements, but after for example changing a mesh, the iterator type
    /// should be global to get all elements.
    const LevelTree<Element*>& getElementsList(
        IteratorType part = IteratorType::LOCAL) const;

    /// Get a vector of elements. If the IteratorType is LOCAL, get all elements
    /// on the core you're working on. If the IteratorType is GLOBAL, get all
    /// elements in the mesh. Usually an application uses only the local
    /// elements, but after for example changing a mesh, the iterator type
    /// should be global to get all elements.
    LevelTree<Element*>& getElementsList(
        IteratorType part = IteratorType::LOCAL);

    /// Get a vector of faces. If the IteratorType is LOCAL, get all faces of
    /// the elements on the core you're working on. If the IteratorType is
    /// GLOBAL, get all faces in the mesh. Usually an application uses only the
    /// local faces.
    const LevelTree<Face*>& getFacesList(
        IteratorType part = IteratorType::LOCAL) const;

    /// Get a vector of faces. If the IteratorType is LOCAL, get all faces of
    /// the elements on the core you're working on. If the IteratorType is
    /// GLOBAL, get all faces in the mesh. Usually an application uses only the
    /// local faces.
    LevelTree<Face*>& getFacesList(IteratorType part = IteratorType::LOCAL);

    const LevelTree<Edge*>& getEdgesList(
        IteratorType part = IteratorType::LOCAL) const;
    LevelTree<Edge*>& getEdgesList(IteratorType part = IteratorType::LOCAL);

    const std::vector<Node*>& getNodesList(
        IteratorType part = IteratorType::LOCAL) const;
    std::vector<Node*>& getNodesList(IteratorType part = IteratorType::LOCAL);

    const std::vector<Geometry::PointPhysical<DIM> >& getNodeCoordinates()
        const;
    std::vector<Geometry::PointPhysical<DIM> >& getNodeCoordinates();

    //********************************************************************************

    TreeIteratorConst<Element*> elementColBegin(
        IteratorType part = IteratorType::LOCAL) const {
        return getElementsList(part).begin();
    }

    TreeIteratorConst<Element*> elementColEnd(
        IteratorType part = IteratorType::LOCAL) const {
        return getElementsList(part).end();
    }

    TreeIterator<Element*> elementColBegin(
        IteratorType part = IteratorType::LOCAL) {
        return getElementsList(part).begin();
    }

    TreeIterator<Element*> elementColEnd(
        IteratorType part = IteratorType::LOCAL) {
        return getElementsList(part).end();
    }

    TreeIteratorConst<Face*> faceColBegin(
        IteratorType part = IteratorType::LOCAL) const {
        return getFacesList(part).begin();
    }

    TreeIteratorConst<Face*> faceColEnd(
        IteratorType part = IteratorType::LOCAL) const {
        return getFacesList(part).end();
    }

    TreeIterator<Face*> faceColBegin(IteratorType part = IteratorType::LOCAL) {
        return getFacesList(part).begin();
    }

    TreeIterator<Face*> faceColEnd(IteratorType part = IteratorType::LOCAL) {
        return getFacesList(part).end();
    }

    TreeIteratorConst<Edge*> edgeColBegin(
        IteratorType part = IteratorType::LOCAL) const {
        return getEdgesList(part).begin();
    }

    TreeIteratorConst<Edge*> edgeColEnd(
        IteratorType part = IteratorType::LOCAL) const {
        return getEdgesList(part).end();
    }

    TreeIterator<Edge*> edgeColBegin(IteratorType part = IteratorType::LOCAL) {
        return getEdgesList(part).begin();
    }

    TreeIterator<Edge*> edgeColEnd(IteratorType part = IteratorType::LOCAL) {
        return getEdgesList(part).end();
    }

    std::vector<Node*>::const_iterator nodeColBegin(
        IteratorType part = IteratorType::LOCAL) const {
        return getNodesList(part).begin();
    }

    std::vector<Node*>::const_iterator nodeColEnd(
        IteratorType part = IteratorType::LOCAL) const {
        return getNodesList(part).end();
    }

    std::vector<Node*>::iterator nodeColBegin(
        IteratorType part = IteratorType::LOCAL) {
        return getNodesList(part).begin();
    }

    std::vector<Node*>::iterator nodeColEnd(
        IteratorType part = IteratorType::LOCAL) {
        return getNodesList(part).end();
    }
    //********************************************************************************

    Submesh& getSubmesh() { return submeshes_; }

    const Submesh& getSubmesh() const { return submeshes_; }

    const std::map<int, std::vector<Element*> >& getPullElements() {
        return submeshes_.getPullElements();
    }

    const std::map<int, std::vector<Element*> >& getPushElements() {
        return submeshes_.getPushElements();
    }

   private:
    Submesh submeshes_;
    //! List of all elements.
    LevelTree<Element*> elements_;

    //! List of all faces.
    LevelTree<Face*> faces_;

    //! List of all edges.
    LevelTree<Edge*> edges_;

    //! List of all nodes. (connectivity-based location of vertices)
    std::vector<Node*> nodes_;

    //! Global vector of physical nodes. (physical location of vertices)
    std::vector<Geometry::PointPhysical<DIM> > nodeCoordinates_;
};

}  // namespace Base
}  // namespace hpgem
#include "Mesh_Impl.h"

#endif  // HPGEM_KERNEL_MESH_H
