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

#include "mesh.h"

namespace Preprocessor {

template <std::size_t entityDimension, std::size_t meshDimension>
Element<meshDimension>& MeshEntity<entityDimension, meshDimension>::getElement(
    EntityLId i) {
    return mesh->getElement(elementIDs[i]);
};

template <std::size_t entityDimension, std::size_t meshDimension>
const Element<meshDimension>&
    MeshEntity<entityDimension, meshDimension>::getElement(EntityLId i) const {
    return mesh->getElement(elementIDs[i]);
}

template <std::size_t entityDimension, std::size_t meshDimension>
EntityGId MeshEntity<entityDimension, meshDimension>::getElementId(
    EntityLId i) const {
    return elementIDs[i];
}

template <std::size_t entityDimension, std::size_t meshDimension>
EntityLId MeshEntity<entityDimension, meshDimension>::getElementIndex(
    const Element<meshDimension>& element) const {
    for (std::size_t i = 0; i < elementIDs.size(); ++i) {
        if (getElement(i) == element) {
            return i;
        }
    }
    logger(ERROR, "Element % not found", element.getGlobalIndex());
    return 0;
}

template <std::size_t entityDimension, std::size_t meshDimension>
std::size_t MeshEntity<entityDimension, meshDimension>::getNumberOfElements()
    const {
    return elementIDs.size();
}

template <std::size_t entityDimension, std::size_t meshDimension>
EntityLId MeshEntity<entityDimension, meshDimension>::getLocalIndex(
    EntityLId i) const {
    return localIDs[i];
}

template <std::size_t entityDimension, std::size_t meshDimension>
EntityLId MeshEntity<entityDimension, meshDimension>::getLocalIndex(
    const Element<meshDimension>& element) const {
    for (std::size_t i = 0; i < elementIDs.size(); ++i) {
        if (getElement(i).getGlobalIndex() == element.getGlobalIndex()) {
            logger.assert_debug(getElement(i) == element,
                                "Two different elements got the same index");
            return localIDs[i];
        }
    }
    logger(ERROR, "Element % not found", element.getGlobalIndex());
    return 0;
}

template <std::size_t entityDimension, std::size_t meshDimension>
EntityGId MeshEntity<entityDimension, meshDimension>::getGlobalIndex() const {
    return entityID;
}

template <std::size_t entityDimension, std::size_t meshDimension>
template <int d>
std::vector<MeshEntity<(d < 0 ? d + meshDimension : d), meshDimension>>
    MeshEntity<entityDimension, meshDimension>::getIncidenceList() const {
    auto indices = getIncidenceListAsIndices<d>();
    std::vector<MeshEntity<(d < 0 ? d + meshDimension : d), meshDimension>>
        result(indices.size());
    for (std::size_t i = 0; i < indices.size(); ++i) {
        result[i] = mesh->template getEntity<(d < 0 ? d + meshDimension : d)>(
            indices[i]);
    }
    return result;
}

template <std::size_t entityDimension, std::size_t meshDimension>
template <int d>
std::vector<EntityGId> MeshEntity<
    entityDimension, meshDimension>::getIncidenceListAsIndices() const {
    static_assert(d + entityDimension >= 0,
                  "The codimension you are interested in is too high for the "
                  "dimension of this object");
    if (d < 0) {
        return getIncidenceListAsIndices<(d < 0 ? d + meshDimension : d)>();
    }
    if (d == entityDimension) {
        return {getGlobalIndex()};
    }
    if (d == meshDimension) {
        return std::vector<EntityGId>(elementIDs.begin(), elementIDs.end());
    } else if (d < entityDimension) {
        // easy case: all adjacent entities of this entityDimension are adjacent
        // to all adjacent elements
        return getElement(0).template getIncidenceListAsIndices<d>(*this);
    } else {
        std::vector<EntityGId> result;
        for (auto elementID : elementIDs) {
            const auto& element = mesh->getElement(elementID);
            std::vector<EntityGId> newEntries =
                element.template getIncidenceListAsIndices<d>(*this);
            result.insert(result.end(), newEntries.begin(), newEntries.end());
        }
        std::sort(result.begin(), result.end());
        result.erase(std::unique(result.begin(), result.end()), result.end());
        return result;
    }
}

template <std::size_t entityDimension, std::size_t meshDimension>
void MeshEntity<entityDimension, meshDimension>::addElement(
    EntityGId elementID, EntityLId localEntityIndex) {
    elementIDs.push_back(elementID);
    localIDs.push_back(localEntityIndex);
}

template <std::size_t entityDimension, std::size_t meshDimension>
void MeshEntity<entityDimension, meshDimension>::removeElement(
    EntityGId elementID) {
    auto entry = std::find(elementIDs.begin(), elementIDs.end(), elementID);
    logger.assert_debug(entry != elementIDs.end(),
                        "Removing a non linked element");
    std::size_t offset = entry - elementIDs.begin();

    elementIDs.erase(entry);
    localIDs.erase(localIDs.begin() + offset);
}

template <std::size_t entityDimension, std::size_t meshDimension>
bool MeshEntity<entityDimension, meshDimension>::operator==(
    const MeshEntity& other) const {
    return mesh == other.mesh && entityID == other.entityID &&
           elementIDs == other.elementIDs && localIDs == other.localIDs;
}

template <std::size_t dimension>
bool Element<dimension>::operator==(const Element& other) const {
    return static_cast<const MeshEntity<dimension, dimension>&>(*this) ==
               static_cast<const MeshEntity<dimension, dimension>&>(other) &&
           other.globalCoordinateIndices == globalCoordinateIndices &&
           other.incidenceLists == incidenceLists;
}

template <std::size_t dimension>
LinearAlgebra::SmallVector<dimension> Element<dimension>::getCoordinate(
    EntityLId localIndex) const {
    return this->mesh->getCoordinate(globalCoordinateIndices[localIndex.id]);
}

template <std::size_t dimension>
CoordId Element<dimension>::getCoordinateIndex(EntityLId localIndex) const {
    return globalCoordinateIndices[localIndex.id];
}

template <std::size_t dimension>
std::vector<LinearAlgebra::SmallVector<dimension>>
    Element<dimension>::getCoordinatesList() const {
    std::vector<LinearAlgebra::SmallVector<dimension>> result;
    result.reserve(referenceGeometry->getNumberOfNodes());
    for (std::size_t i = 0; i < referenceGeometry->getNumberOfNodes(); ++i) {
        result.push_back(getCoordinate(i));
    }
    return result;
}

template <std::size_t dimension>
void Element<dimension>::setNodeCoordinate(
    EntityLId localIndex, LinearAlgebra::SmallVector<dimension> newCoordinate) {
    CoordId globalIndex = globalCoordinateIndices[localIndex.id];
    this->mesh->updateCoordinate(globalIndex, newCoordinate);
}

template <std::size_t dimension>
void Element<dimension>::setNode(EntityLId localIndex, EntityGId globalIndex,
                                 CoordId coordinateIndex) {
    EntityGId oldNodeIndex = incidenceLists[0][localIndex.id];
    incidenceLists[0][localIndex.id] = globalIndex;
    globalCoordinateIndices[localIndex.id] = coordinateIndex;
    if (oldNodeIndex != globalIndex) {
        MeshEntity<0, dimension>& newNode = this->mesh->getNode(globalIndex);
        MeshEntity<0, dimension>& oldNode = this->mesh->getNode(oldNodeIndex);
        newNode.addElement(this->entityID, localIndex);
        oldNode.removeElement(this->entityID);
    }
}

template <std::size_t dimension>
template <std::size_t d>
std::enable_if_t<(d > 0)> Element<dimension>::setEntity(EntityLId localIndex,
                                                        EntityGId globalIndex) {
    EntityGId current = incidenceLists[d][localIndex];
    incidenceLists[d][localIndex] = globalIndex;
    if (current != globalIndex) {
        // Update happened
        MeshEntity<d, dimension>& newEntity =
            this->mesh->template getEntity<d>(globalIndex);
        newEntity.addElement(this->getGlobalIndex(), localIndex);

        MeshEntity<d, dimension>& oldEntity =
            this->mesh->template getEntity<d>(current);
        oldEntity.removeElement(this->getGlobalIndex());
    }
}

template <std::size_t dimension>
template <int d, std::size_t entityDimension>
std::vector<MeshEntity<(d < 0 ? d + dimension : d), dimension>>
    Element<dimension>::getIncidenceList(
        const MeshEntity<entityDimension, dimension>& entity) const {
    auto indices = getIncidenceListAsIndices<d>(entity);
    std::vector<MeshEntity<(d < 0 ? d + dimension : d), dimension>> result(
        indices.size());
    for (std::size_t i = 0; i < indices.size(); ++i) {
        result[i] = this->mesh->template getEntity<(d < 0 ? d + dimension : d)>(
            indices[i]);
    }
}

template <std::size_t dimension>
template <int d, std::size_t entityDimension>
std::vector<EntityGId> Element<dimension>::getIncidenceListAsIndices(
    const MeshEntity<entityDimension, dimension>& entity) const {
    auto localIndices = getLocalIncidenceListAsIndices<d>(entity);
    std::vector<EntityGId> globalIndices(localIndices.size());
    for (std::size_t i = 0; i < globalIndices.size(); ++i) {
        globalIndices[i] =
            incidenceLists[(d < 0 ? d + dimension : d)][localIndices[i]];
    }
    return globalIndices;
}

template <std::size_t dimension>
template <int d, std::size_t entityDimension>
std::vector<EntityLId> Element<dimension>::getLocalIncidenceListAsIndices(
    const MeshEntity<entityDimension, dimension>& entity) const {
    static_assert(d + dimension >= 0,
                  "The requested codimension is too high for the dimension of "
                  "this element");
    constexpr std::size_t actualDimension = (d < 0 ? d + dimension : d);
    if (incidenceLists[actualDimension].size() <
        referenceGeometry->template getNumberOfEntities<actualDimension>())
        return {};
    std::vector<std::size_t> result =
        referenceGeometry
            ->template getAdjacentEntities<entityDimension, actualDimension>(
                entity.getLocalIndex(*this));
    std::vector<EntityLId> realResult(result.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
        realResult[i] = EntityLId(result[i]);
    }
    return realResult;
}

template <std::size_t dimension>
void Element<dimension>::addNode(EntityGId globalNodeIndex,
                                 CoordId coordinateIndex) {
    incidenceLists[0].push_back(globalNodeIndex);
    globalCoordinateIndices.push_back(coordinateIndex);
}

template <std::size_t dimension>
template <std::size_t d>
std::enable_if_t<(d > 0)> Element<dimension>::addEntity(EntityGId globalIndex) {
    incidenceLists[d].push_back(globalIndex);
}

template <std::size_t dimension>
template <std::size_t d>
void Element<dimension>::remapIndices(
    const std::map<EntityGId, EntityGId>& mapping) {
    for (EntityGId& id : incidenceLists[d]) {
        auto newIndex = mapping.find(id);
        logger.assert_debug(newIndex != mapping.end(), "Non mapped index %",
                            id);
        id = newIndex->second;
    }
}

template <std::size_t dimension>
std::vector<Element<dimension>>& Mesh<dimension>::getElements() {
    return elementsList;
}

template <std::size_t dimension>
const std::vector<Element<dimension>>& Mesh<dimension>::getElements() const {
    return elementsList;
}

template <std::size_t dimension>
std::vector<MeshEntity<dimension - 1, dimension>>& Mesh<dimension>::getFaces() {
    return getEntities<dimension - 1>();
}

template <std::size_t dimension>
const std::vector<MeshEntity<dimension - 1, dimension>>&
    Mesh<dimension>::getFaces() const {
    return getEntities<dimension - 1>();
}

template <std::size_t dimension>
std::vector<MeshEntity<1, dimension>>& Mesh<dimension>::getEdges() {
    return getEntities<1>();
}

template <std::size_t dimension>
const std::vector<MeshEntity<1, dimension>>& Mesh<dimension>::getEdges() const {
    return getEntities<1>();
}

template <std::size_t dimension>
std::vector<MeshEntity<0, dimension>>& Mesh<dimension>::getNodes() {
    return getEntities<0>();
}

template <std::size_t dimension>
const std::vector<MeshEntity<0, dimension>>& Mesh<dimension>::getNodes() const {
    return getEntities<0>();
}

template <std::size_t dimension>
template <int entityDimension>
std::vector<MeshEntity<(entityDimension < 0 ? entityDimension + dimension
                                            : entityDimension),
                       dimension>>&
    Mesh<dimension>::getEntities() {
    static_assert(entityDimension + dimension >= 0,
                  "The requested codimension is too high for the dimension of "
                  "this element");
    constexpr std::size_t actualDimension =
        (entityDimension < 0 ? entityDimension + dimension : entityDimension);
    return otherEntities.template getData<actualDimension>();
};

template <std::size_t dimension>
template <int entityDimension>
const std::vector<MeshEntity<(entityDimension < 0 ? entityDimension + dimension
                                                  : entityDimension),
                             dimension>>&
    Mesh<dimension>::getEntities() const {
    static_assert(entityDimension + dimension >= 0,
                  "The requested codimension is too high for the dimension of "
                  "this element");
    constexpr std::size_t actualDimension =
        (entityDimension < 0 ? entityDimension + dimension : entityDimension);
    return otherEntities.template getData<actualDimension>();
}

template <std::size_t dimension>
std::vector<typename Mesh<dimension>::coordinateData>&
    Mesh<dimension>::getNodeCoordinates() {
    return coordinates;
}

template <std::size_t dimension>
const std::vector<typename Mesh<dimension>::coordinateData>&
    Mesh<dimension>::getNodeCoordinates() const {
    return coordinates;
}

template <std::size_t dimension>
void Mesh<dimension>::setNumberOfNodes(std::size_t number) {
    if (number < otherEntities.template getData<0>().size()) {
        otherEntities.template getData<0>().resize(number);
    } else
        addNodes(number - otherEntities.template getData<0>().size());
}

template <std::size_t dimension>
void Mesh<dimension>::addNode() {
    // TODO: Check
    EntityGId newIndex = EntityGId(otherEntities.template getData<0>().size());
    otherEntities.template getData<0>().push_back({this, newIndex});
}

template <std::size_t dimension>
void Mesh<dimension>::addNodes(std::size_t count) {
    for (std::size_t i = 0; i < count; ++i) {
        addNode();
    }
}

template <std::size_t dimension>
std::size_t Mesh<dimension>::addNodeCoordinate(
    EntityGId nodeIndex, LinearAlgebra::SmallVector<dimension> coordinate) {
    coordinates.push_back({nodeIndex, coordinate});
    return coordinates.size() - 1;
}

template <std::size_t dimension>
void Mesh<dimension>::addElement(std::vector<CoordId> nodeCoordinateIDs,
                                 const std::string& zoneName) {
    EntityGId elementID = EntityGId(elementsList.size());
    Element<dimension> newElement{this, elementID, 0};
    newElement.setGeometry(findGeometry(nodeCoordinateIDs.size()));
    for (auto coordinateID : nodeCoordinateIDs) {
        EntityGId nodeID = coordinates[coordinateID.id].nodeIndex;
        newElement.addNode(nodeID, coordinateID);
    }
    elementsList.push_back(newElement);
    otherEntities.template getData<dimension>().push_back(newElement);
    for (std::size_t i = 0; i < nodeCoordinateIDs.size(); ++i) {
        auto coordinateID = nodeCoordinateIDs[i];
        EntityGId nodeID = coordinates[coordinateID.id].nodeIndex;
        // TODO: Check
        otherEntities.template getData<0>()[nodeID.id].addElement(elementID, i);
    }
    fixElement(elementsList.back(), tag<dimension - 1>{});
}

template <std::size_t dimension>
void Mesh<dimension>::fixConnectivity() {
    for (auto& element : getElements()) {
        fixElement(element, tag<dimension - 1>{});
    }
}

template <std::size_t dimension>
template <std::size_t d>
void Mesh<dimension>::fixElement(Element<dimension>& element, tag<d>) {
    for (std::size_t i = element.incidenceLists[d].size();
         i < element.referenceGeometry->template getNumberOfEntities<d>();
         ++i) {
        fixEntity<d>(element, i);
    }
    fixElement(element, tag<d - 1>{});
}

template <std::size_t dimension>
void Mesh<dimension>::fixElement(Element<dimension>& element, tag<1>) {
    for (std::size_t i = element.incidenceLists[1].size();
         i < element.referenceGeometry->template getNumberOfEntities<1>();
         ++i) {
        fixEntity<1>(element, i);
    }
}

template <std::size_t dimension>
template <std::size_t d>
void Mesh<dimension>::fixEntity(Element<dimension>& element,
                                std::size_t index) {
    std::vector<std::size_t> localNodeIndices =
        element.referenceGeometry->template getAdjacentEntities<d, 0>(index);
    std::vector<EntityGId> globalNodeIndices(localNodeIndices.size());
    for (std::size_t i = 0; i < localNodeIndices.size(); ++i) {
        globalNodeIndices[i] =
            element.template getIncidentEntity<0>(localNodeIndices[i])
                .getGlobalIndex();
    }
    auto candidates = getNodes()[globalNodeIndices[0].id]
                          .template getIncidenceListAsIndices<d>();
    std::sort(candidates.begin(), candidates.end());
    for (std::size_t i = 1; i < globalNodeIndices.size(); ++i) {
        auto newElements = getNodes()[globalNodeIndices[i].id]
                               .template getIncidenceListAsIndices<d>();
        std::sort(newElements.begin(), newElements.end());
        std::vector<EntityGId> temp;
        temp.reserve(std::max(candidates.size(), newElements.size()));
        std::set_intersection(candidates.begin(), candidates.end(),
                              newElements.begin(), newElements.end(),
                              std::back_inserter(temp));
        candidates = std::move(temp);
    }
    logger.assert_always(candidates.size() < 2, "mesh data is not consistent");
    EntityGId entityIndex =
        (candidates.size() == 1 ? candidates[0] : newEntity<d>());
    element.template addEntity<d>(entityIndex);
    otherEntities.template getData<d>()[entityIndex.id].addElement(
        element.getGlobalIndex(), index);
}

template <std::size_t dimension>
const ElementShape<dimension>* Mesh<dimension>::findGeometry(
    std::size_t numberOfNodes) {
    for (auto shape : defaultShapes<dimension>) {
        if (shape->getNumberOfNodes() == numberOfNodes) return shape;
    }
    logger(ERROR, "There are no % dimensional default shapes with % nodes",
           dimension, numberOfNodes);
    return nullptr;
}

template <std::size_t DIM>
std::size_t Mesh<DIM>::getZoneId(const std::string& zoneName) {
    for (std::size_t i = 0; i < zoneNames.size(); ++i) {
        if (zoneNames[i] == zoneName) {
            return i;
        }
    }
    // Not found
    zoneNames.push_back(zoneName);
    return zoneNames.size() - 1;
}

template <std::size_t dimension>
template <std::size_t d>
void Mesh<dimension>::mergeEntities(
    EntityGId face1Id, EntityGId face2Id,
    const std::vector<EntityLId>& nodePermutation, tag<d>) {

    static_assert(d < dimension, "Can not merge elements");
    // Nodes should have their own overload
    static_assert(d > 0, "Wrong function");

    MeshEntity<d, dimension>& face1 = getFace(face1Id);
    MeshEntity<d, dimension>& face2 = getFace(face2Id);

    logger(INFO, "Merging faces % (elem %) and % (elem %)",
           face1.getGlobalIndex(), face1.getElementId(0),
           face2.getGlobalIndex(), face2.getElementId(0));

    // Merge nodes around the face as needed.
    std::vector<EntityGId> f1nodes =
        face1.template getIncidenceListAsIndices<0>();
    std::vector<EntityGId> f2nodes =
        face2.template getIncidenceListAsIndices<0>();
    // Mapping of Node EntityGId's from face1 -> face2.
    for (std::size_t i = 0; i < nodePermutation.size(); ++i) {
        EntityGId& node1 = f1nodes[i];
        EntityGId& node2 = f2nodes[nodePermutation[i].id];
        if (node1 != node2) {
            mergeNodes(node1, node2);
        }
    }
    mergeEntitiesInternal(face1, face2);
}
template <std::size_t dimension>
void Mesh<dimension>::mergeEntities(
    EntityGId node1, EntityGId node2,
    const std::vector<EntityLId>& nodePermutation, tag<0>) {
    logger.assert_debug(nodePermutation.size() == 1,
                        "Incorrect node permutiation size");
    mergeNodes(node1, node2);
}

template <std::size_t dimension>
template <std::size_t d>
void Mesh<dimension>::mergeEntitiesInternal(MeshEntity<d, dimension>& entity1,
                                            MeshEntity<d, dimension>& entity2) {
    static_assert(d < dimension, "Can't merge elements");
    // The case d == 0 should default to the base case.
    static_assert(d > 0, "Not the right function for merging nodes");

    if (d > 1) {
        // Merging the boundary is not straightforward, as there are multiple
        // entities on the boundary and we do not apriori know which boundary
        // entities of E1 should be merged with which boundary entities of E2.
        //
        // A necessary requirement for merging boundary entity Eb1 of E1 and Eb2
        // of E2 is that both use the same set of nodes. We assume that all
        // boundary nodes on E1 (and thus also E2) have unique EntityGId's. That
        // way we can find the pairs Eb1 and Eb2 by looking for d-1 dimensional
        // boundary entities that share the same set of EntityGIds.
        std::vector<EntityGId> boundary1 =
            entity1.template getIncidenceListAsIndices<d - 1>();
        std::vector<EntityGId> boundary2 =
            entity2.template getIncidenceListAsIndices<d - 1>();

        // Possible future optimization:
        // Remove entries that are present in both boundaries, as these are
        // already merged.

        // Mapping Nodes on boundary -> corresponding boundary entity,
        // the vectors MUST be sorted.
        std::map<std::vector<EntityGId>, EntityGId> boundary1ByNodes;
        std::map<std::vector<EntityGId>, EntityGId> boundary2ByNodes;

        for (EntityGId& bound1Id : boundary1) {
            MeshEntity<d - 1, dimension>& bound1 =
                getEntities<d - 1>()[bound1Id.id];
            std::vector<EntityGId> nodes =
                bound1.template getIncidenceListAsIndices<0>();
            std::sort(nodes.begin(), nodes.end());
            boundary1ByNodes[nodes] = bound1Id;
        }

        for (EntityGId& bound2Id : boundary2) {
            MeshEntity<d - 1, dimension>& bound2 =
                getEntities<d - 1>()[bound2Id.id];
            std::vector<EntityGId> nodes =
                bound2.template getIncidenceListAsIndices<0>();
            std::sort(nodes.begin(), nodes.end());
            boundary2ByNodes[nodes] = bound2Id;
        }
        // Merge boundary parts
        for (auto const& bound2 : boundary2ByNodes) {
            auto bound1 = boundary1ByNodes.find(bound2.first);
            logger.assert_debug(bound1 != boundary1ByNodes.end(),
                                "Could not find matching boundary part");
            // Check if the entities use a different entity on the boundary.
            if (bound1->second != bound2.second) {
                // A merge is needed
                MeshEntity<d - 1, dimension>& bound1Entity =
                    getEntities<d - 1>()[bound1->second.id];
                MeshEntity<d - 1, dimension>& bound2Entity =
                    getEntities<d - 1>()[bound2.second.id];
                mergeEntitiesInternal(bound1Entity, bound2Entity);
            }
        }
    }

    // Merge the entities themselves
    // For this we simply need to replace all references to entity2 by entity1,
    // as 0 < d < dimension there are only references in Element.

    while (entity2.getNumberOfElements() > 0) {
        EntityGId element = entity2.getElementId(0);
        EntityLId index = entity2.getLocalIndex(0);
        getElement(element).template setEntity<d>(index,
                                                  entity1.getGlobalIndex());
    }
}

template <std::size_t dimension>
void Mesh<dimension>::mergeNodes(EntityGId node1, EntityGId node2) {
    // Merge steps:
    // Coordinates: Update references from node2 -> node1
    // Elements: Replace linkage from node2 -> node1
    // Transitive, update edges -> faces

    // Update coordinates
    for (coordinateData& coordinate : coordinates) {
        if (coordinate.nodeIndex == node2) {
            coordinate.nodeIndex = node1;
        }
    }
    MeshEntity<0, dimension>& node1Ref = getNodes()[node1.id];
    MeshEntity<0, dimension>& node2Ref = getNodes()[node2.id];

    // Update elements
    for (Element<dimension>& element : elementsList) {
        std::vector<EntityGId> nodeIds =
            element.template getIncidenceListAsIndices<0>();
        for (std::size_t i = 0; i < element.getNumberOfNodes(); ++i) {
            // Relink element -> node
            if (nodeIds[i] == node2) {
                element.setNode(i, node1, element.getCoordinateIndex(i));
            }
        }
    }
}

template <std::size_t dimension>
void Mesh<dimension>::removeUnusedEntities() {
    // Recursively compact in between indices
    removeUnusedEntities(tag<dimension - 1>{});
    // entity dimension 0 is special, it requires also updating the coordinates
    {
        std::map<EntityGId, EntityGId> nodeRemapping;
        EntityGId newIndex = EntityGId(0);
        EntityGId maxId = EntityGId(getNumberOfNodes());
        std::vector<MeshEntity<0, dimension>>& nodes = getNodes();
        for (EntityGId i = EntityGId(0); i < maxId; ++i) {
            if (getNode(i).getNumberOfElements() > 0) {
                nodeRemapping[i] = newIndex;
                if (newIndex != i) {
                    // Move to new location
                    nodes[newIndex.id] = nodes[i.id];
                    // Update the id of the node
                    nodes[newIndex.id].entityID = newIndex;
                }
                newIndex++;
            }
        }
        if (newIndex != maxId) {
            // There is at least one deleted node.
            nodes.resize(newIndex.id);

            // Update elements
            EntityGId maxElemId = EntityGId(elementsList.size());
            for (EntityGId i = EntityGId(0); i < maxElemId; ++i) {
                getElement(i).template remapIndices<0>(nodeRemapping);
            }
            // Update coordinates
            for (coordinateData& coord : coordinates) {
                coord.nodeIndex = nodeRemapping[coord.nodeIndex];
            }
            logger(INFO, "Reduced the number of nodes from % to %", maxId,
                   newIndex);
        }
    }
    // Note: Coordinate remapping not implemented as:
    // 1. This is far harder (one needs to do a liveness check via elements)
    // 2. I don't think this affects the structure
}
template <std::size_t dimension>
template <std::size_t d>
void Mesh<dimension>::removeUnusedEntities(tag<d>) {
    std::map<EntityGId, EntityGId> nodeRemapping;
    EntityGId newIndex = EntityGId(0);
    EntityGId maxId = EntityGId(getNumberOfEntities<d>());
    std::vector<MeshEntity<d, dimension>>& entities = getEntities<d>();

    for (EntityGId i = EntityGId(0); i < maxId; ++i) {
        if (getEntity<d>(i).getNumberOfElements() > 0) {
            nodeRemapping[i] = newIndex;
            if (newIndex != i) {
                entities[newIndex.id] = entities[i.id];
                entities[newIndex.id].entityID = newIndex;
            }
            newIndex++;
        }
    }

    if (newIndex != maxId) {
        // There is at least one deleted entity
        entities.resize(newIndex.id);

        // Update elements
        EntityGId maxElemId = EntityGId(elementsList.size());
        for (EntityGId i = EntityGId(0); i < maxElemId; ++i) {
            getElement(i).template remapIndices<d>(nodeRemapping);
        }
        logger(INFO, "Reduced the number of dim-%-entities from % to %", d,
               maxId, newIndex);
    }
    // Recurse to lower dimensions
    removeUnusedEntities(tag<d - 1>{});
}

}  // namespace Preprocessor
