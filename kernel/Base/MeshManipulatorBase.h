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

#ifndef HPGEM_KERNEL_MESHMANIPULATORBASE_H
#define HPGEM_KERNEL_MESHMANIPULATORBASE_H

#include <vector>
#include <fstream>

#include "Geometry/FaceGeometry.h"
#include "Mesh.h"
#include "GlobalNamespaceBase.h"
#include "FE/BasisFunctionSet.h"
#include "Zone.h"

namespace hpgem {

namespace Base {
template <std::size_t DIM>
class MeshManipulator;
}

template <std::size_t DIM>
std::ostream& operator<<(std::ostream&, const Base::MeshManipulator<DIM>&);

namespace Geometry {
template <std::size_t DIM>
class PointPhysical;
template <std::size_t DIM>
class PointReference;
}  // namespace Geometry

namespace Base {
class BasisFunctionSet;
class OrientedBasisFunctionSet;
class Face;
template <std::size_t DIM>
class MeshMoverBase;
template <class V>
class LevelTree;
class Element;
struct ConfigurationData;
class Edge;

class MeshManipulatorBase {
   public:
    using CollectionOfBasisFunctionSets =
        Element::CollectionOfBasisFunctionSets;

    using ConstElementIterator = TreeIteratorConst<Element*>;
    using ElementIterator = TreeIterator<Element*>;

    using ConstFaceIterator = TreeIteratorConst<Face*>;
    using FaceIterator = TreeIterator<Face*>;

    /// idRangeBegin is the beginning of the range, from where the Element's ids
    /// should be assigned. In case of multiple meshes, one has to take care of
    /// empty intersection of those ranges!!!
    MeshManipulatorBase(const ConfigurationData* configData,
                        std::size_t dimension,
                        std::size_t numberOfElementMatrixes = 0,
                        std::size_t numberOfElementVectors = 0,
                        std::size_t numberOfFaceMatrixes = 0,
                        std::size_t numberOfFaceVectors = 0);

    // Meshes own a large number of entities like Elements. Copying them
    // requires a deep copy and is practically never required
    MeshManipulatorBase(const MeshManipulatorBase& other) = delete;

    virtual ~MeshManipulatorBase() = default;

    virtual std::size_t getNumberOfElements(
        IteratorType part = IteratorType::LOCAL) const = 0;

    virtual std::size_t getNumberOfFaces(
        IteratorType part = IteratorType::LOCAL) const = 0;

    virtual std::size_t getNumberOfEdges(
        IteratorType part = IteratorType::LOCAL) const = 0;

    virtual std::size_t getNumberOfNodes(
        IteratorType part = IteratorType::LOCAL) const = 0;

    /// *****************Iteration through the Elements*******************

    virtual ConstElementIterator elementColBegin(
        IteratorType part = IteratorType::LOCAL) const = 0;

    virtual ConstElementIterator elementColEnd(
        IteratorType part = IteratorType::LOCAL) const = 0;

    virtual ElementIterator elementColBegin(
        IteratorType part = IteratorType::LOCAL) = 0;

    virtual ElementIterator elementColEnd(
        IteratorType part = IteratorType::LOCAL) = 0;

    virtual ConstFaceIterator faceColBegin(
        IteratorType part = IteratorType::LOCAL) const = 0;

    virtual ConstFaceIterator faceColEnd(
        IteratorType part = IteratorType::LOCAL) const = 0;

    virtual FaceIterator faceColBegin(
        IteratorType part = IteratorType::LOCAL) = 0;

    virtual FaceIterator faceColEnd(
        IteratorType part = IteratorType::LOCAL) = 0;

    virtual TreeIteratorConst<Edge*> edgeColBegin(
        IteratorType part = IteratorType::LOCAL) const = 0;

    virtual TreeIteratorConst<Edge*> edgeColEnd(
        IteratorType part = IteratorType::LOCAL) const = 0;

    virtual TreeIterator<Edge*> edgeColBegin(
        IteratorType part = IteratorType::LOCAL) = 0;

    virtual TreeIterator<Edge*> edgeColEnd(
        IteratorType part = IteratorType::LOCAL) = 0;

    virtual std::vector<Node*>::const_iterator nodeColBegin(
        IteratorType part = IteratorType::LOCAL) const = 0;

    virtual std::vector<Node*>::const_iterator nodeColEnd(
        IteratorType part = IteratorType::LOCAL) const = 0;

    virtual std::vector<Node*>::iterator nodeColBegin(
        IteratorType part = IteratorType::LOCAL) = 0;

    virtual std::vector<Node*>::iterator nodeColEnd(
        IteratorType part = IteratorType::LOCAL) = 0;
    //  *****************Iteration through the Elements*******************

    //! Get const list of elements
    virtual const LevelTree<Element*>& getElementsList(
        IteratorType part = IteratorType::LOCAL) const = 0;

    //! Get non-const list of elements
    virtual LevelTree<Element*>& getElementsList(
        IteratorType part = IteratorType::LOCAL) = 0;

    //! Get const list of faces
    virtual const LevelTree<Face*>& getFacesList(
        IteratorType part = IteratorType::LOCAL) const = 0;

    //! Get non-const list of faces
    virtual LevelTree<Face*>& getFacesList(
        IteratorType part = IteratorType::LOCAL) = 0;

    virtual const LevelTree<Edge*>& getEdgesList(
        IteratorType part = IteratorType::LOCAL) const = 0;

    virtual LevelTree<Edge*>& getEdgesList(
        IteratorType part = IteratorType::LOCAL) = 0;

    virtual const std::vector<Node*>& getNodesList(
        IteratorType part = IteratorType::LOCAL) const = 0;

    virtual std::vector<Node*>& getNodesList(
        IteratorType part = IteratorType::LOCAL) = 0;

    virtual const std::map<int, std::vector<Element*>>& getPullElements() = 0;

    virtual const std::map<int, std::vector<Element*>>& getPushElements() = 0;

    std::size_t dimension() const { return dimension_; }

    const ConfigurationData* getConfigData() { return configData_; }

    /**
     * Add a zone to the mesh.
     *
     * The lifetime of the created zone (and references to it) is tied to the
     * lifetime of the mesh.
     *
     * @param name The name of the zone
     * @return A reference to the newly created zone.
     */
    Zone& addZone(std::string name);

    /**
     * @return All the zones in this mesh
     */
    const std::vector<std::unique_ptr<Zone>>& getZones() const {
        return zones_;
    }

   protected:
    const ConfigurationData* configData_;
    //! Periodicity in x-direction.
    bool periodicX_;

    //! Periodicity in y-direction.
    bool periodicY_;

    //! Periodicity in z-direction.
    bool periodicZ_;

    std::size_t numberOfElementMatrices_;
    std::size_t numberOfFaceMatrices_;
    std::size_t numberOfElementVectors_;
    std::size_t numberOfFaceVectors_;
    const std::size_t dimension_;

    /// Zones in the mesh
    // Implementation note: These are pointers so that elements can directly
    // refer to them. They are unique_ptr as they are tied to the mesh (just as
    // as for example Elements)
    std::vector<std::unique_ptr<Zone>> zones_;
};

}  // namespace Base

}  // namespace hpgem

#endif  // HPGEM_KERNEL_MESHMANIPULATORBASE_H
