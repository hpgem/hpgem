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

// NOTE: Purely definitions (for includes and structure see ../Mesh.h)

namespace Preprocessor {

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
        EntityGId nodeIndex;
        /// The coordinates
        LinearAlgebra::SmallVector<dimension> coordinate;
    };

   public:
    Mesh() = default;
    ~Mesh() = default;

    Mesh(const Mesh& other)
        : elementsList(other.elementsList),
          coordinates(other.coordinates),
          meshEntities(other.meshEntities) {
        for (auto& element : elementsList) {
            element.mesh = this;
        }
        copyEntities(tag<dimension>{});
    }

    Mesh(Mesh&& other)
        : elementsList(std::move(other.elementsList)),
          coordinates(std::move(other.coordinates)),
          meshEntities(std::move(other.meshEntities)) {
        for (auto& element : elementsList) {
            element.mesh = this;
        }
        copyEntities(tag<dimension>{});
    }

    Mesh& operator=(const Mesh& other) {
        elementsList = other.elementsList;
        coordinates = other.coordinates;
        meshEntities = other.meshEntities;
        for (auto& element : elementsList) {
            element.mesh = this;
        }
        copyEntities(tag<dimension>{});
    }

    Mesh& operator=(Mesh&& other) {
        elementsList = std::move(other.elementsList);
        coordinates = std::move(other.coordinates);
        meshEntities = std::move(other.meshEntities);
        for (auto& element : elementsList) {
            element.mesh = this;
        }
        copyEntities(tag<dimension>{});
    }

    // Elements //
    //////////////

    std::vector<Element<dimension>>& getElements();
    const std::vector<Element<dimension>>& getElements() const;

    /// Get an element by its index.
    /// \param i The global index of the element
    /// \return The element.
    Element<dimension>& getElement(EntityGId i) { return getElements()[i.id]; }
    const Element<dimension>& getElement(EntityGId i) const {
        return getElements()[i.id];
    }
    std::size_t getNumberOfElements() const { return getElements().size(); }

    // Faces //
    ///////////

    std::vector<MeshEntity<dimension - 1, dimension>>& getFaces();
    const std::vector<MeshEntity<dimension - 1, dimension>>& getFaces() const;
    MeshEntity<dimension - 1, dimension> getFace(std::size_t i) const {
        return getFaces()[i];
    };
    std::size_t getNumberOfFaces() const { return getFaces().size(); }

    // Edges //
    ///////////

    std::vector<MeshEntity<1, dimension>>& getEdges();
    const std::vector<MeshEntity<1, dimension>>& getEdges() const;
    MeshEntity<1, dimension> getEdge(std::size_t i) const {
        return getEdges()[i];
    };
    std::size_t getNumberOfEdges() const { return getEdges().size(); }

    // Nodes //
    ///////////

    std::vector<MeshEntity<0, dimension>>& getNodes();
    const std::vector<MeshEntity<0, dimension>>& getNodes() const;
    MeshEntity<0, dimension>& getNode(EntityGId i) { return getNodes()[i.id]; };
    const MeshEntity<0, dimension>& getNode(EntityGId i) const {
        return getNodes()[i.id];
    };
    std::size_t getNumberOfNodes() const { return getNodes().size(); }

    // Coordinates //
    /////////////////

    std::vector<coordinateData>& getNodeCoordinates();
    const std::vector<coordinateData>& getNodeCoordinates() const;

    LinearAlgebra::SmallVector<dimension> getCoordinate(
        CoordId coordinateIndex) const {
        return getNodeCoordinates()[coordinateIndex.id].coordinate;
    }

    // Zones //
    ///////////

    const std::vector<std::string>& getZoneNames() const { return zoneNames; }

    // MeshEntity //
    ////////////////

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
    MeshEntity<entityDimension, dimension>& getEntity(EntityGId i) {
        return getEntities<entityDimension>()[i.id];
    };
    template <int entityDimension>
    std::size_t getNumberOfEntities() const {
        return entityDimension + dimension >= 0
                   ? getEntities<entityDimension>().size()
                   : 0;
    }

    /**
     * Get the number of entities of a certain dimension.
     *
     * This includes the unused entities (i.e. entities without an attached
     * element).
     *
     * Note the following equivalence should hold for 0 <= i <= dimension
     * getNumberOfEntities(i) == getNumberOfEntities<i>(), but the latter is
     * probably slightly faster.
     *
     * @param entityDim The dimension of the entity, must be <= dimension
     * @return The number entities
     */
    std::size_t getNumberOfEntities(std::size_t entityDim) const {
        if (entityDim == dimension) {
            return elementsList.size();
        } else {
            logger.assert_debug(
                entityDim < dimension,
                "Asking for entities with dimension higher than the mesh.");
            return getNumberOfEntities(entityDim, itag<dimension - 1>{});
        }
    }

    // Modifiers //
    ///////////////

    EntityGId addNode();
    void addNodes(std::size_t count);

    CoordId addNodeCoordinate(EntityGId nodeIndex,
                              LinearAlgebra::SmallVector<dimension> coordinate);

    void addElement(std::vector<CoordId> nodeCoordinateIDs,
                    const std::string& zoneName = "main");

    void updateCoordinate(CoordId coordinateIndex,
                          LinearAlgebra::SmallVector<dimension> coordinate) {
        getNodeCoordinates()[coordinateIndex.id].coordinate = coordinate;
    }

    /**
     * Removes unused MeshEntities, i.e. the MeshEntities without attached
     * Elements.
     */
    void removeUnusedEntities() { removeUnusedEntities(itag<dimension - 1>{}); }

    // Checks //
    ////////////

    bool isValid() const;
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
    bool checkEntities(tag<d> dimTag) const {
        EntityGId id = EntityGId(0);
        for (auto entity : meshEntities[dimTag]) {
            if (entity.getGlobalIndex() != id) {
                logger(ERROR,
                       "Incorrect id % for d-dimensional entity at position %",
                       entity.getGlobalIndex(), d, id);
            }
            id++;  // For next entity
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
    bool checkEntities(tag<0> dimTag) const {
        EntityGId id = EntityGId(0);
        for (auto entity : meshEntities[dimTag]) {
            if (entity.getGlobalIndex() != id) {
                logger(ERROR,
                       "Incorrect id % for 0-dimensional entity at position %",
                       entity.getGlobalIndex(), id);
            }
            id++;  // For next entity
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
    void fixEntity(Element<dimension>& element, std::size_t i);

    template <int d>
    void removeUnusedEntities(itag<d> dimTag);
    // Base case
    void removeUnusedEntities(itag<-1>){};

    /**
     * Helper for getNumberOfEntities(std::size_t)
     * @tparam d The recursion dimension
     * @param entityDimension The dimension to get the count of
     * @return The number of entities of entityDimension.
     */
    template <int d>
    std::size_t getNumberOfEntities(std::size_t entityDimension, itag<d>) const;
    std::size_t getNumberOfEntities(std::size_t, itag<-1>) const { return 0; }

    template <std::size_t d>
    void copyEntities(tag<d> dimTag) {
        for (auto& entity : meshEntities[dimTag]) {
            entity.mesh = this;
        }
        copyEntities(tag<d - 1>{});
    }
    void copyEntities(tag<0> dimTag) {
        for (auto& entity : meshEntities[dimTag]) {
            entity.mesh = this;
        }
    }

    template <std::size_t entityDimension>
    EntityGId newEntity() {
        EntityGId newIndex =
            EntityGId(meshEntities.template get<entityDimension>().size());
        meshEntities.template get<entityDimension>().push_back(
            {this, newIndex});
        return newIndex;
    }

    std::size_t getZoneId(const std::string& zoneName);

    const ElementShape<dimension>* findGeometry(std::size_t numberOfNodes);

    std::vector<Element<dimension>> elementsList;
    std::vector<coordinateData> coordinates;

    // Alias template to fill in the dimension of the mesh
    template <std::size_t entityDimension>
    using MeshEntityVec = std::vector<MeshEntity<entityDimension, dimension>>;
    // Note: Also stores the elements (in addition to elementsList). However,
    // here as MeshEntitiy<dim,dim>, the parentType of Element. This allows much
    // easier implementation of getEntities<dim>()
    TemplateArray<dimension + 1, MeshEntityVec> meshEntities;
    std::vector<std::string> zoneNames;
};
}  // namespace Preprocessor

#endif  // HPGEM_APP_MESH_H
