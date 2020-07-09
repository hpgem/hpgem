/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2017, University of Twente
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

#ifndef HPGEM_APP_MESH_H
#define HPGEM_APP_MESH_H

#include "Base/ConstIterableWrapper.h"
#include "LinearAlgebra/SmallVector.h"
#include "elementShape.h"

// the data structures needed to represent a mesh
// they are collected in one file because they are meaningless alone and
// intimately linked they are redesigned here because the classes used in the
// rest of hpGEM need way more information than what is needed to construct and
// distribute the mesh

namespace Preprocessor {

template <std::size_t dimension, std::size_t gridDimension>
class MeshEntity;

namespace Detail {

template <std::size_t entityDimension, std::size_t dimension = entityDimension>
struct MeshEntities : public MeshEntities<entityDimension - 1, dimension> {
    std::vector<MeshEntity<entityDimension, dimension>> data;

    template <std::size_t SUB_DIM>
    std::vector<MeshEntity<SUB_DIM, dimension>>& getData() {
        static_assert(SUB_DIM <= entityDimension,
                      "Requested data is not available in this entity");
        return MeshEntities<SUB_DIM, dimension>::data;
    }

    template <std::size_t SUB_DIM>
    const std::vector<MeshEntity<SUB_DIM, dimension>>& getData() const {
        static_assert(SUB_DIM <= entityDimension,
                      "Requested data is not available in this entity");
        return MeshEntities<SUB_DIM, dimension>::data;
    }
};

template <std::size_t dimension>
struct MeshEntities<0, dimension> {
    std::vector<MeshEntity<0, dimension>> data;

    template <std::size_t SUB_DIM>
    std::vector<MeshEntity<0, dimension>>& getData() {
        static_assert(SUB_DIM == 0,
                      "Requested data is not available in this entity");
        return data;
    }

    template <std::size_t SUB_DIM>
    const std::vector<MeshEntity<0, dimension>>& getData() const {
        static_assert(SUB_DIM == 0,
                      "Requested data is not available in this entity");
        return data;
    }
};

constexpr std::size_t exp2(std::size_t power) {
    if (power == 0) {
        return 1;
    }
    if ((power / 2) * 2 == power) {
        return exp2(power / 2) * exp2(power / 2);
    } else {
        return 2 * exp2(power - 1);
    }
}
}  // namespace Detail

template <std::size_t dimension>
class Element;

template <std::size_t dimension>
class Mesh;

// Nodes, Edges and Faces behave the same so they can be described by the same
// class
/**
 * @brief While this is conceptually just an index into the mesh, it also
 * presents some functions for local inspection and manipulation of the mesh
 * @tparam gridDimension the dimension of the mesh
 */
template <std::size_t dimension, std::size_t gridDimension>
class MeshEntity {
   public:
    // the default constructor constructs an entity in a moved-from state so it
    // is only useful if you want to create a std::array or something similar
    // that is filled later
    MeshEntity() = default;

    ~MeshEntity() = default;

    MeshEntity(const MeshEntity&) = default;

    MeshEntity& operator=(const MeshEntity&) = default;

    MeshEntity(MeshEntity&& other) = default;

    MeshEntity& operator=(MeshEntity&&) = default;

    Element<gridDimension>& getElement(std::size_t i);
    const Element<gridDimension>& getElement(std::size_t i) const;

    std::size_t getElementIndex(const Element<gridDimension>& element) const;

    std::size_t getNumberOfElements() const;

    std::size_t getLocalIndex(std::size_t i) const;

    std::size_t getLocalIndex(const Element<gridDimension>& element) const;

    std::size_t getGlobalIndex() const;

    // note: this will return a vector of MeshEntities. To access individual
    // elements, call getElement(0) on the element or call getElement i for all
    // i on this entity
    std::vector<MeshEntity<gridDimension, gridDimension>> getElementsList()
        const {
        return getIncidenceList<gridDimension>();
    }

    std::vector<MeshEntity<gridDimension - 1, gridDimension>> getFacesList()
        const {
        return getIncidenceList<gridDimension - 1>();
    }

    std::size_t getNumberOfFaces() const { return getFacesList().size(); }

    std::vector<MeshEntity<1, gridDimension>> getEdgesList() const {
        return getIncidenceList<1>();
    }

    std::size_t getNumberOfEdges() const { return getEdgesList().size(); }

    std::vector<MeshEntity<0, gridDimension>> getNodesList() const {
        return getIncidenceList<0>();
    }

    std::size_t getNumberOfNodes() const { return getNodesList().size(); }

    // when d == dimension, it will just return {*this}
    template <int d>
    std::vector<MeshEntity<(d < 0 ? d + gridDimension : d), gridDimension>>
        getIncidenceList() const;

    template <int d>
    std::vector<std::size_t> getIncidenceListAsIndices() const;

    template <int d>
    std::size_t getNumberOfIncidentEntities() const {
        return getIncidenceListAsIndices<d>().size();
    }

    template <int d>
    MeshEntity<(d < 0 ? d + gridDimension : d), gridDimension>
        getIncidentEntity(std::size_t i) const {
        return mesh->template getEntity<(d < 0 ? d + gridDimension : d)>(
            getIncidenceListAsIndices<(d < 0 ? d + gridDimension : d)>()[i]);
    }

    const Mesh<gridDimension>* getMesh() const { return mesh; }

    bool operator==(const MeshEntity&) const;
    bool operator!=(const MeshEntity& other) const { return !(*this == other); }

   protected:
    friend Mesh<gridDimension>;
    MeshEntity(Mesh<gridDimension>* mesh, std::size_t entityID)
        : mesh(mesh), entityID(entityID) {}

    void addElement(std::size_t elementID, std::size_t localEntityIndex);

    Mesh<gridDimension>* mesh;
    std::size_t entityID = std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> elementIDs;
    std::vector<std::size_t> localIDs;
};

/**
 * @brief While this is conceptually just an index into the mesh, it also
 * presents some functions for local inspection and manipulation of the mesh
 * @tparam dimension the dimension of the mesh
 */
template <std::size_t dimension>
class Element : public MeshEntity<dimension, dimension> {
   public:
    // the default constructor constructs an entity in a moved-from state so it
    // is only useful if you want to create a std::array or something similar
    // that is filled later
    Element() = default;

    ~Element() = default;

    Element(const Element&) = default;

    Element(Element&&) = default;

    Element& operator=(const Element&) = default;

    Element& operator=(Element&&) = default;

    bool operator==(const Element&) const;
    bool operator!=(const Element& element) const {
        return !(*this == element);
    }

    LinearAlgebra::SmallVector<dimension> getCoordinate(
        std::size_t localIndex) const;
    std::size_t getCoordinateIndex(std::size_t localIndex) const;

    std::vector<LinearAlgebra::SmallVector<dimension>> getCoordinatesList()
        const;

    void setNodeCoordinate(std::size_t localIndex,
                           LinearAlgebra::SmallVector<dimension> newCoordinate);

    void setNode(std::size_t localIndex, std::size_t globalIndex,
                 std::size_t coordinateIndex);

    template <std::size_t d>
    std::enable_if_t<(d > 0)> setEntity(std::size_t localIndex,
                                        std::size_t globalIndex);

    using MeshEntity<dimension, dimension>::getIncidenceList;
    using MeshEntity<dimension, dimension>::getIncidenceListAsIndices;

    // get entities adjacent to both this element and the passed entity
    // e.g. adjacent faces of this element that are adjacent to the specified
    // node or edges adjacent to the element (by passing *this) the passed
    // entity must be adjacent to this element according to the zero argument
    // version
    template <int d, std::size_t entityDimension>
    std::vector<MeshEntity<(d < 0 ? d + dimension : d), dimension>>
        getIncidenceList(
            const MeshEntity<entityDimension, dimension>& entity) const;

    template <int d, std::size_t entityDimension>
    std::vector<std::size_t> getIncidenceListAsIndices(
        const MeshEntity<entityDimension, dimension>& entity) const;

    template <int d, std::size_t entityDimension>
    std::vector<std::size_t> getLocalIncidenceListAsIndices(
        const MeshEntity<entityDimension, dimension>& entity) const;

   private:
    friend Mesh<dimension>;
    Element(Mesh<dimension>* mesh, std::size_t elementID)
        : MeshEntity<dimension, dimension>(mesh, elementID) {
        this->addElement(elementID, 0);
    }

    void addNode(std::size_t globalNodeIndex, std::size_t coordinateIndex);

    template <std::size_t d>
    std::enable_if_t<(d > 0)> addEntity(std::size_t globalIndex);

    void setGeometry(const ElementShape<dimension>* shape) {
        referenceGeometry = shape;
    }

    const ElementShape<dimension>* referenceGeometry;
    // technically a map from the integer d to a list of indices of entities of
    // dimension d
    std::array<std::vector<std::size_t>, dimension> incidenceLists;
    std::vector<std::size_t> globalCoordinateIndices;
};

// the mesh class is designed for the preprocessor of hpGEM. Since a mesh will
// always have coordinates accosiated with it these are baked into the data
// structure further, there is no anticipated need for removing elements at this
// stage of the computation, so this is currently not supported
template <std::size_t dimension>
class Mesh {

    struct coordinateData {
        std::size_t nodeIndex;
        LinearAlgebra::SmallVector<dimension> coordinate;
    };

    template <std::size_t d>
    struct tag {};

   public:
    Mesh() = default;
    ~Mesh() = default;

    Mesh(const Mesh& other)
        : elementsList(other.elementsList),
          coordinates(other.coordinates),
          otherEntities(other.otherEntities) {
        for (auto& element : elementsList) {
            element.mesh = this;
        }
        copyEntities(tag<dimension>{});
    }

    Mesh(Mesh&& other)
        : elementsList(std::move(other.elementsList)),
          coordinates(std::move(other.coordinates)),
          otherEntities(std::move(other.otherEntities)) {
        for (auto& element : elementsList) {
            element.mesh = this;
        }
        copyEntities(tag<dimension>{});
    }

    Mesh& operator=(const Mesh& other) {
        elementsList = other.elementsList;
        coordinates = other.coordinates;
        otherEntities = other.otherEntities;
        for (auto& element : elementsList) {
            element.mesh = this;
        }
        copyEntities(tag<dimension>{});
    }

    Mesh& operator=(Mesh&& other) {
        elementsList = std::move(other.elementsList);
        coordinates = std::move(other.coordinates);
        otherEntities = std::move(other.otherEntities);
        for (auto& element : elementsList) {
            element.mesh = this;
        }
        copyEntities(tag<dimension>{});
    }

    std::vector<Element<dimension>>& getElements();
    const std::vector<Element<dimension>>& getElements() const;
    Element<dimension>& getElement(std::size_t i) { return getElements()[i]; }
    const Element<dimension>& getElement(std::size_t i) const {
        return getElements()[i];
    }
    std::size_t getNumberOfElements() const { return getElements().size(); }

    std::vector<MeshEntity<dimension - 1, dimension>>& getFaces();
    const std::vector<MeshEntity<dimension - 1, dimension>>& getFaces() const;
    MeshEntity<dimension - 1, dimension> getFace(std::size_t i) const {
        return getFaces()[i];
    };
    std::size_t getNumberOfFaces() const { return getFaces().size(); }

    std::vector<MeshEntity<1, dimension>>& getEdges();
    const std::vector<MeshEntity<1, dimension>>& getEdges() const;
    MeshEntity<1, dimension> getEdge(std::size_t i) const {
        return getEdges()[i];
    };
    std::size_t getNumberOfEdges() const { return getEdges().size(); }

    std::vector<MeshEntity<0, dimension>>& getNodes();
    const std::vector<MeshEntity<0, dimension>>& getNodes() const;
    MeshEntity<0, dimension> getNode(std::size_t i) const {
        return getNodes()[i];
    };
    std::size_t getNumberOfNodes() const { return getNodes().size(); }

    std::vector<coordinateData>& getNodeCoordinates();
    const std::vector<coordinateData>& getNodeCoordinates() const;

    template <int entityDimension>
    std::vector<MeshEntity<(entityDimension < 0 ? entityDimension + dimension
                                                : entityDimension),
                           dimension>>&
        getEntities();
    template <int entityDimension>
    const std::vector<MeshEntity<
        (entityDimension < 0 ? entityDimension + dimension : entityDimension),
        dimension>>&
        getEntities() const;
    template <int entityDimension>
    MeshEntity<entityDimension, dimension> getEntity(std::size_t i) const {
        return getEntities<entityDimension>()[i];
    };
    template <int entityDimension>
    std::size_t getNumberOfEntities() const {
        return entityDimension + dimension >= 0
                   ? getEntities<entityDimension>().size()
                   : 0;
    }

    void setNumberOfNodes(std::size_t number);
    void addNode();
    void addNodes(std::size_t count);

    std::size_t addNodeCoordinate(
        std::size_t nodeIndex,
        LinearAlgebra::SmallVector<dimension> coordinate);

    void addElement(std::vector<std::size_t> nodeCoordinateIDs);

    void updateCoordinate(std::size_t coordinateIndex,
                          LinearAlgebra::SmallVector<dimension> coordinate) {
        getNodeCoordinates()[coordinateIndex].coordinate = coordinate;
    }
    LinearAlgebra::SmallVector<dimension> getCoordinate(
        std::size_t coordinateIndex) const {
        return getNodeCoordinates()[coordinateIndex].coordinate;
    }

    bool isValid() const {
        logger(DEBUG, "The mesh has % elements", getNumberOfElements());
        for (auto element : getElements()) {
            if (getElement(element.getGlobalIndex()) != element) {
                logger(ERROR,
                       "The index for element % has been set incorrectly",
                       element.getGlobalIndex());
                return false;
            }
            if (!checkBoundingEntities(element, tag<dimension - 1>{})) {
                return false;
            }
        }
        return checkEntities(tag<dimension>{});
    }
    // todo: is this needed?
    void fixConnectivity();

   private:
    // do all 'faces' that bound this element know that they bound this element?
    template <std::size_t d>
    bool checkBoundingEntities(const Element<dimension>& element,
                               tag<d>) const {
        for (auto boundingEntity : element.template getIncidenceList<d>()) {
            auto candidateElements = boundingEntity.getElementsList();
            if (std::find(candidateElements.begin(), candidateElements.end(),
                          element) == candidateElements.end()) {
                logger(ERROR,
                       "The element is bounded by a %-dimensional shape that "
                       "is not near the element",
                       d);
                return false;
            }
        }
        return checkBoundingEntities(element, tag<d - 1>{});
    }
    bool checkBoundingEntities(const Element<dimension>& element,
                               tag<0>) const {
        for (auto boundingEntity : element.template getIncidenceList<0>()) {
            auto candidateElements = boundingEntity.getElementsList();
            if (std::find(candidateElements.begin(), candidateElements.end(),
                          element) == candidateElements.end()) {
                logger(ERROR,
                       "The element is bounded by a %-dimensional shape that "
                       "is not near the element",
                       0);
                return false;
            }
        }
        return true;
    }

    // do the elements bounded by this 'face' know this
    template <std::size_t d>
    bool checkEntities(tag<d>) const {
        for (auto entity : otherEntities.template getData<d>()) {
            for (std::size_t i = 0; i < entity.getNumberOfElements(); ++i) {
                if (entity !=
                    entity.getElement(i).template getIncidentEntity<d>(
                        entity.getLocalIndex(i))) {
                    logger(ERROR,
                           "This %-dimensional shape is adjacent to an element "
                           "that is not bounded by this shape",
                           d);
                    return false;
                }
            }
        }
        return checkEntities(tag<d - 1>{});
    }
    bool checkEntities(tag<0>) const {
        for (auto entity : otherEntities.template getData<0>()) {
            for (std::size_t i = 0; i < entity.getNumberOfElements(); ++i) {
                if (entity !=
                    entity.getElement(i).template getIncidentEntity<0>(
                        entity.getLocalIndex(i))) {
                    logger(ERROR,
                           "This %-dimensional shape is adjacent to an element "
                           "that is not bounded by this shape",
                           0);
                    return false;
                }
            }
        }
        return true;
    }

    template <std::size_t d>
    void fixElement(Element<dimension>& element, tag<d>);
    void fixElement(Element<dimension>& element, tag<1>);
    void fixElement(Element<dimension>& element, tag<0>) {}

    template <std::size_t d>
    void copyEntities(tag<d>) {
        for (auto& entity : otherEntities.template getData<d>()) {
            entity.mesh = this;
        }
        copyEntities(tag<d - 1>{});
    }
    void copyEntities(tag<0>) {
        for (auto& entity : otherEntities.template getData<0>()) {
            entity.mesh = this;
        }
    }

    template <std::size_t d>
    void fixEntity(Element<dimension>& element, std::size_t i);

    template <std::size_t entityDimension>
    std::size_t newEntity() {
        std::size_t newIndex =
            otherEntities.template getData<entityDimension>().size();
        otherEntities.template getData<entityDimension>().push_back(
            {this, newIndex});
        return newIndex;
    }

    const ElementShape<dimension>* findGeometry(std::size_t numberOfNodes);

    std::vector<Element<dimension>> elementsList;
    std::vector<coordinateData> coordinates;
    // re-stores slices of the elements for consistent use of getEntities<d>,
    // algorithms that need elements can un-slice by asking for the first (and
    // only) element
    Detail::MeshEntities<dimension> otherEntities;
};

// for use with file readers that can 'guess' the correct numbering of the node
// coordinates file readers that don't do this should overload this function
template <std::size_t dimension, typename file_t>
Mesh<dimension> readFile(file_t& file) {
    Mesh<dimension> result;
    for (auto nodeCoordinates : file.getNodeCoordinates()) {
        result.addNode();
        for (auto coordinate : nodeCoordinates) {
            logger.assert_debug(
                coordinate.size() == dimension,
                "The coordinates read by this reader have the wrong dimension");
            result.addNodeCoordinate(result.getNumberOfNodes() - 1,
                                     coordinate.data());
        }
    }
    for (auto element : file.getElements()) {
        result.addElement(element);
    }
    logger.assert_debug(result.isValid(), "Unspecified problem with the mesh");
    return result;
};
}  // namespace Preprocessor

#include "mesh_impl.h"

using namespace hpgem;

#endif  // HPGEM_APP_MESH_H
