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

#include "LinearAlgebra/SmallVector.h"
#include "elementShape.h"
#include "MeshSource.h"

// the data structures needed to represent a mesh
// they are collected in one file because they are meaningless alone and
// intimately linked they are redesigned here because the classes used in the
// rest of hpGEM need way more information than what is needed to construct and
// distribute the mesh

namespace Preprocessor {

template <std::size_t dimension, std::size_t gridDimension>
class MeshEntity;

template <std::size_t dimension>
class Element;

template <std::size_t dimension>
class Mesh;

/**
 * \brief Topological part of the Mesh
 *
 * A MeshEntity describes the topology (connectivity) of a single element, face,
 * edge or node of the mesh. To illustrate the connectivity let us look at a
 * Face (2D) in a 3D mesh. This face has a boundary that consist of Nodes (0D)
 * and Edges (1D) and is the boundary of one or more Elements (3D). In general,
 * a MeshEntity of dimension dimE:
 * - Is connected to MeshEntity-s with dimensions smaller than DimE which form
 *   its boundary.
 * - Is connected to MeshEntity-s with dimensions larger than DimE where it is
 *   part of the boundary of such an entity.
 * - Is connected to MeshEntity-s with dimension  equal to DimE, where it is
 *   itself.
 * The last relation may seem superfluous, but this ensures that for
 * any two entities we can determine whether they are connected or not.
 *
 * Note that a MeshEntity does not contain any geometrical information (e.g.
 * coordinate of a Node), as this is not possible in the case of periodic
 * boundaries. This information should be accessed via a connected Element,
 * which does contain geometrical information.
 *
 * @tparam entityDimension The dimension of the entity itself (e.g. 0 for a
 * Node)
 * @tparam meshDimension The dimension of the mesh.
 */
template <std::size_t entityDimension, std::size_t meshDimension>
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

    /// Get a connected Element by its local index.
    Element<meshDimension>& getElement(std::size_t i);
    const Element<meshDimension>& getElement(std::size_t i) const;

    /// Given a connected Element, what is its local index?
    std::size_t getElementIndex(const Element<meshDimension>& element) const;

    /// The number of elements that this MeshEntity is part of.
    std::size_t getNumberOfElements() const;

    /// The Local index of this MeshEntity on the i-th Element.
    std::size_t getLocalIndex(std::size_t i) const;

    /// Get the local index of this MeshEntity for an element that it is part
    /// of.
    /// \param element The element
    /// \return The local index on the element.
    std::size_t getLocalIndex(const Element<meshDimension>& element) const;

    /// The global index of this MeshEntity. For a given Mesh and
    /// entityDimension this uniquely determines the MeshEntity.
    std::size_t getGlobalIndex() const;

    std::vector<MeshEntity<meshDimension, meshDimension>> getElementsList()
        const {
        return getIncidenceList<meshDimension>();
    }

    /// List of facet-MeshEntities to which this is connected.
    /// \return
    std::vector<MeshEntity<meshDimension - 1, meshDimension>> getFacesList()
        const {
        return getIncidenceList<meshDimension - 1>();
    }

    /// Number of facet-MeshEntities to which this is connected.
    std::size_t getNumberOfFaces() const { return getFacesList().size(); }

    std::vector<MeshEntity<1, meshDimension>> getEdgesList() const {
        return getIncidenceList<1>();
    }

    std::size_t getNumberOfEdges() const { return getEdgesList().size(); }

    std::vector<MeshEntity<0, meshDimension>> getNodesList() const {
        return getIncidenceList<0>();
    }

    std::size_t getNumberOfNodes() const { return getNodesList().size(); }

    // when d == entityDimension, it will just return {*this}
    /// List of all connected entities
    template <int d>
    std::vector<MeshEntity<(d < 0 ? d + meshDimension : d), meshDimension>>
        getIncidenceList() const;

    template <int d>
    std::vector<std::size_t> getIncidenceListAsIndices() const;

    template <int d>
    std::size_t getNumberOfIncidentEntities() const {
        return getIncidenceListAsIndices<d>().size();
    }

    template <int d>
    MeshEntity<(d < 0 ? d + meshDimension : d), meshDimension>
        getIncidentEntity(std::size_t i) const {
        return mesh->template getEntity<(d < 0 ? d + meshDimension : d)>(
            getIncidenceListAsIndices<(d < 0 ? d + meshDimension : d)>()[i]);
    }

    const Mesh<meshDimension>* getMesh() const { return mesh; }

    bool operator==(const MeshEntity&) const;
    bool operator!=(const MeshEntity& other) const { return !(*this == other); }

   protected:
    friend Mesh<meshDimension>;
    MeshEntity(Mesh<meshDimension>* mesh, std::size_t entityID)
        : mesh(mesh), entityID(entityID) {}

    /// Add an element that is this MeshEntity is part of
    ///
    /// \param elementID The entityID of the element
    /// \param localEntityIndex The localIndex of this MeshEntity for the
    /// element.
    void addElement(std::size_t elementID, std::size_t localEntityIndex);

    Mesh<meshDimension>* mesh;
    /// The id of this MeshEntity
    std::size_t entityID = std::numeric_limits<std::size_t>::max();
    /// The entityIDs for the Elements that this MeshEntity is part of
    std::vector<std::size_t> elementIDs;
    /// For each of the elements that this MeshEntity is part of, the local id
    /// of this MeshEntity on the element.
    std::vector<std::size_t> localIDs;
    // NOTE: Probably the following identity should hold:
    // elements[elementIDs[i]].incidenceList[dimension][localIDs[i]] == entityID
};

/**
 * @brief While this is conceptually just an index into the mesh, it also
 * presents some functions for local inspection and manipulation of the mesh
 * @tparam dim the dimension of the mesh
 */
/**
 * \brief An Element of the Mesh.
 *
 * The MeshEntity-s of largest dimension in the Mesh are called Elements. In
 * addition to the topological information they also include geometrical
 * information. Speficially:
 *  - A polytope as shape (e.g. a triangle)
 *  - Coordinates for each of the node-type MeshEntity-s (i.e.
 *    MeshEntity<0,dim>) that form the corners of the polytopical shape.
 *
 * The shape determines not only the shape of the Element, but also of the
 * MeshEntity-s that are on the boundary of this Element. The local index of a
 * MeshEntity on this Element can be used to determine its shape.
 *
 * Note that the coordinates for each of the nodes is Element local. While most
 * connected elements to a node usually have the same coordinate, this is not
 * required and not possible with periodic boundaries.
 *
 * @tparam dim The dimension of the
 */
template <std::size_t dim>
class Element : public MeshEntity<dim, dim> {
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

    LinearAlgebra::SmallVector<dim> getCoordinate(std::size_t localIndex) const;
    std::size_t getCoordinateIndex(std::size_t localIndex) const;

    std::vector<LinearAlgebra::SmallVector<dim>> getCoordinatesList() const;

    /// Set the global coordinate of a node
    /// \param localIndex The local index of the node
    /// \param newCoordinate The new coordinate.
    void setNodeCoordinate(std::size_t localIndex,
                           LinearAlgebra::SmallVector<dim> newCoordinate);

    /// Set the node at a specific local index.
    ///
    /// \param localIndex The local index of the node
    /// \param globalIndex The global index of the node
    /// \param coordinateIndex  The global coordinate index for the point
    /// corresponding to the node.
    void setNode(std::size_t localIndex, std::size_t globalIndex,
                 std::size_t coordinateIndex);

    template <std::size_t d>
    std::enable_if_t<(d > 0)> setEntity(std::size_t localIndex,
                                        std::size_t globalIndex);

    using MeshEntity<dim, dim>::getIncidenceList;
    using MeshEntity<dim, dim>::getIncidenceListAsIndices;

    /// \brief For a MeshEntity on the boundary of this Element, compute the
    /// MeshEntity-s that are shared with this element.
    ///
    /// For a MeshEntity on the boundary (and thus is connected to this Element)
    /// computes the MeshEntity-s that are connected to this Element and are
    /// connected to the input MeshEntity. For example, for a node-MeshEntity
    /// compute the face-MeshEntity-s of this element that have that node.
    ///
    /// \tparam d The dimension of the resulting MeshEntity-s.
    /// \tparam entityDimension The dimension of the input MeshEntity
    /// \param entity The MeshEntity on the boundary.
    /// \return The shared MeshEntity-s
    template <int d, std::size_t entityDimension>
    std::vector<MeshEntity<(d < 0 ? d + dim : d), dim>> getIncidenceList(
        const MeshEntity<entityDimension, dim>& entity) const;

    /// Same as getIncidenceList(const MeshEntity<entityDimension, dim>& entity)
    /// but computing the global indices of the MeshEntity-s
    /// \tparam d The dimension of the shared MeshEntity-s
    /// \tparam entityDimension The dimension of the input MeshEntity
    /// \param entity The MeshEntity on the boundary
    /// \return The global indices of the shared MeshEntity-s
    template <int d, std::size_t entityDimension>
    std::vector<std::size_t> getIncidenceListAsIndices(
        const MeshEntity<entityDimension, dim>& entity) const;

    /// Same as getIncidenceList(const MeshEntity<entityDimension, dim>& entity)
    /// but computing the local indices of the MeshEntity-s
    /// \tparam d The dimension of the shared MeshEntity-s
    /// \tparam entityDimension The dimension of the input MeshEntity
    /// \param entity The MeshEntity on the boundary
    /// \return The local indices of the shared MeshEntity-s
    template <int d, std::size_t entityDimension>
    std::vector<std::size_t> getLocalIncidenceListAsIndices(
        const MeshEntity<entityDimension, dim>& entity) const;

    std::string getZoneName() { return this->mesh->zoneNames[zoneId]; }

   private:
    friend Mesh<dim>;
    Element(Mesh<dim>* mesh, std::size_t elementID, std::size_t zoneId)
        : MeshEntity<dim, dim>(mesh, elementID), zoneId(zoneId) {
        this->addElement(elementID, 0);
    }

    /// Add a node to this element.
    /// \param globalNodeIndex The global index of the (topological) node
    /// \param coordinateIndex The global index of the coordinate for the node.
    void addNode(std::size_t globalNodeIndex, std::size_t coordinateIndex);

    template <std::size_t d>
    std::enable_if_t<(d > 0)> addEntity(std::size_t globalIndex);

    void setGeometry(const ElementShape<dim>* shape) {
        referenceGeometry = shape;
    }

    /// The geometry of this element (e.g. a cube or tetrahedron)
    const ElementShape<dim>* referenceGeometry;
    /// The global indices of all the MeshEntities that form the boundary of
    /// this element, grouped by the dimension of the MeshEntities.
    std::array<std::vector<std::size_t>, dim> incidenceLists;
    /// The global indices of the coordinates for the corners of this element.
    std::vector<std::size_t> globalCoordinateIndices;

    std::size_t zoneId;
};

namespace Detail {

/// MeshEntities stores MeshEntity-s of all dimensions.
// Implementation note, inheritiance is used to recursively store all
// dimensions.
template <std::size_t entityDimension, std::size_t dimension = entityDimension>
struct MeshEntities : public MeshEntities<entityDimension - 1, dimension> {
    /// The entities themselves
    std::vector<MeshEntity<entityDimension, dimension>> data;

    /// Get the MeshEntity's of a specific entityDimension
    template <std::size_t SUB_DIM>
    std::vector<MeshEntity<SUB_DIM, dimension>>& getData() {
        static_assert(SUB_DIM <= entityDimension,
                      "Requested data is not available in this entity");
        return MeshEntities<SUB_DIM, dimension>::data;
    }

    /// Get the MeshEntity's of a specific entityDimension
    template <std::size_t SUB_DIM>
    const std::vector<MeshEntity<SUB_DIM, dimension>>& getData() const {
        static_assert(SUB_DIM <= entityDimension,
                      "Requested data is not available in this entity");
        return MeshEntities<SUB_DIM, dimension>::data;
    }
};

// Base case of MeshEntities, there are no MeshEntities of dimension lower than
// 0.
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

/// Simplified Mesh finite element Mesh for use with the preprocessor.
///
/// A bare minimum Mesh for use with the Preprocessor. This is conceptually the
/// same as the Mesh in hpGEM, except that it is stripped down to the bare
/// minimum topological and geometrical information needed for the tasks of the
/// preprocessor. This ensures that we can still store a large mesh for
/// partitioning.
///
/// \tparam dimension The dimension of the Mesh
template <std::size_t dimension>
class Mesh {

    /// A coordinate of a Node in the Mesh
    struct coordinateData {
        /// The node to which this belongs. Note that there maybe more than one
        /// coordinate with the same nodeIndex.
        std::size_t nodeIndex;
        /// The coordinates
        LinearAlgebra::SmallVector<dimension> coordinate;
    };

    /// Template class helper to pass a dimension d.
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
    /// Get an element by its index.
    /// \param i The global index of the element
    /// \return The element.
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

    /// Get the list of MeshEntity-s of a specific dimension
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

    /// Set the number of nodes in the Mesh. Will add or remove Nodes as
    /// necessary.
    void setNumberOfNodes(std::size_t number);
    std::size_t addNode();
    void addNodes(std::size_t count);

    std::size_t addNodeCoordinate(
        std::size_t nodeIndex,
        LinearAlgebra::SmallVector<dimension> coordinate);

    void addElement(std::vector<std::size_t> nodeCoordinateIDs,
                    const std::string& zoneName = "main");

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
    /// Recursively checks that for each MeshEntity that is adjacent to an
    /// Element, that that MeshEntity has the Element in its list of adjacent
    /// Elements.
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
    // Base case.
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
    /// Reverse of checkBoundingEntities(). For each MeshEntity, check for all
    /// pairs of (adjecent Element, localIndex) that the MeshEntity is indeed at
    /// localIndex of the Element.
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
    // Base case.
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

    std::size_t getZoneId(const std::string& zoneName);

    const ElementShape<dimension>* findGeometry(std::size_t numberOfNodes);

    std::vector<Element<dimension>> elementsList;
    std::vector<coordinateData> coordinates;
    // re-stores slices of the elements for consistent use of getEntities<d>,
    // algorithms that need elements can un-slice by asking for the first (and
    // only) element
    Detail::MeshEntities<dimension> otherEntities;
    std::vector<std::string> zoneNames;
};

// for use with file readers that can 'guess' the correct numbering of the node
// coordinates file readers that don't do this should overload this function
template <std::size_t dimension>
Mesh<dimension> readFile(MeshSource& file) {
    Mesh<dimension> result;
    logger.assert_always(dimension == file.getDimension(),
                         "Mismatching dimensions");
    for (auto nodeCoordinates : file.getNodeCoordinates()) {
        result.addNode();
        for (auto coordinate : nodeCoordinates.coordinates) {
            logger.assert_debug(
                coordinate.size() == dimension,
                "The coordinates read by this reader have the wrong dimension");
            result.addNodeCoordinate(result.getNumberOfNodes() - 1,
                                     coordinate.data());
        }
    }
    for (auto element : file.getElements()) {
        logger.assert_always(!element.zoneName.empty(),
                             "Element without a zone name");
        result.addElement(element.coordinateIds, element.zoneName);
    }
    logger.assert_debug(result.isValid(), "Unspecified problem with the mesh");
    return result;
}

template <std::size_t dimension>
Mesh<dimension> fromMeshSource(MeshSource2& file) {
    Mesh<dimension> result;
    logger.assert_always(dimension == file.getDimension(),
                         "Mismatching dimensions");
    // Mapping of the nodeId, source nodeId -> mesh nodeId
    std::map<std::size_t, std::size_t> nodeMapping;
    // Add all the nodes
    for (auto coord : file.getCoordinates()) {
        // Find the node in the mesh, or create a new one if needed
        std::size_t sourceNodeId = coord.nodeId;
        std::size_t meshNodeId;
        auto current = nodeMapping.find(sourceNodeId);
        if (current != nodeMapping.end()) {
            meshNodeId = current->second;
        } else {
            meshNodeId = result.addNode();
            nodeMapping[sourceNodeId] = meshNodeId;
        }
        // Add the actual coordinate
        logger.assert_debug(
            coord.coordinate.size() == dimension,
            "The coordinates read by this reader have the wrong dimension");
        result.addNodeCoordinate(meshNodeId, coord.coordinate);
    }
    // Add all the elements
    for (auto element : file.getElements()) {
        logger.assert_always(!element.zoneName.empty(),
                             "Element without a zone name");
        result.addElement(element.coordinateIds, element.zoneName);
    }
    logger.assert_debug(result.isValid(), "Unspecified problem with the mesh");
    return result;
}

}  // namespace Preprocessor

#include "mesh_impl.h"

#endif  // HPGEM_APP_MESH_H
