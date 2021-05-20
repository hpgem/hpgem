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

#include "ElementShapes.h"

namespace Preprocessor {

template <std::size_t entityDimension, std::size_t meshDimension>
Element<meshDimension>& MeshEntity<entityDimension, meshDimension>::getElement(
    std::size_t i) {
    return mesh->getElement(elementIDs[i]);
};

template <std::size_t entityDimension, std::size_t meshDimension>
const Element<meshDimension>&
    MeshEntity<entityDimension, meshDimension>::getElement(
        std::size_t i) const {
    return mesh->getElement(elementIDs[i]);
}

template <std::size_t entityDimension, std::size_t meshDimension>
std::size_t MeshEntity<entityDimension, meshDimension>::getElementIndex(
    const Element<meshDimension>& element) const {
    for (std::size_t i = 0; i < elementIDs.size(); ++i) {
        if (getElement(i) == element) {
            return i;
        }
    }
    logger(ERROR, "Element not found");
    return 0;
}

template <std::size_t entityDimension, std::size_t meshDimension>
std::size_t MeshEntity<entityDimension, meshDimension>::getNumberOfElements()
    const {
    return elementIDs.size();
}

template <std::size_t entityDimension, std::size_t meshDimension>
std::size_t MeshEntity<entityDimension, meshDimension>::getLocalIndex(
    std::size_t i) const {
    return localIDs[i];
}

template <std::size_t entityDimension, std::size_t meshDimension>
std::size_t MeshEntity<entityDimension, meshDimension>::getLocalIndex(
    const Element<meshDimension>& element) const {
    for (std::size_t i = 0; i < elementIDs.size(); ++i) {
        if (getElement(i).getGlobalIndex() == element.getGlobalIndex()) {
            logger.assert_debug(getElement(i) == element,
                                "Two different elements got the same index");
            return localIDs[i];
        }
    }
    logger(ERROR, "Element not found");
    return 0;
}

template <std::size_t entityDimension, std::size_t meshDimension>
std::size_t MeshEntity<entityDimension, meshDimension>::getGlobalIndex() const {
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
std::vector<std::size_t> MeshEntity<
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
        return std::vector<std::size_t>(elementIDs.begin(), elementIDs.end());
    } else if (d < entityDimension) {
        // easy case: all adjacent entities of this entityDimension are adjacent
        // to all adjacent elements
        return getElement(0).template getIncidenceListAsIndices<d>(*this);
    } else {
        std::vector<std::size_t> result;
        for (auto elementID : elementIDs) {
            const auto& element = mesh->getElement(elementID);
            std::vector<std::size_t> newEntries =
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
    std::size_t elementID, std::size_t localEntityIndex) {
    elementIDs.push_back(elementID);
    localIDs.push_back(localEntityIndex);
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
    std::size_t localIndex) const {
    return this->mesh->getCoordinate(globalCoordinateIndices[localIndex]);
}

template <std::size_t dimension>
std::size_t Element<dimension>::getCoordinateIndex(
    std::size_t localIndex) const {
    return globalCoordinateIndices[localIndex];
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
    std::size_t localIndex,
    LinearAlgebra::SmallVector<dimension> newCoordinate) {
    std::size_t globalIndex = globalCoordinateIndices[localIndex];
    this->mesh->updateCoordinate(globalIndex, newCoordinate);
}

template <std::size_t dimension>
void Element<dimension>::setNode(std::size_t localIndex,
                                 std::size_t globalIndex,
                                 std::size_t coordinateIndex) {
    incidenceLists[0][localIndex] = globalIndex;
    globalCoordinateIndices[localIndex] = coordinateIndex;
}

template <std::size_t dimension>
template <std::size_t d>
std::enable_if_t<(d > 0)> Element<dimension>::setEntity(
    std::size_t localIndex, std::size_t globalIndex) {
    incidenceLists[d][localIndex] = globalIndex;
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
std::vector<std::size_t> Element<dimension>::getIncidenceListAsIndices(
    const MeshEntity<entityDimension, dimension>& entity) const {
    auto result = getLocalIncidenceListAsIndices<d>(entity);
    for (auto& index : result) {
        index = incidenceLists[(d < 0 ? d + dimension : d)][index];
    }
    return result;
}

template <std::size_t dimension>
template <int d, std::size_t entityDimension>
std::vector<std::size_t> Element<dimension>::getLocalIncidenceListAsIndices(
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
    return result;
}

template <std::size_t dimension>
void Element<dimension>::addNode(std::size_t globalNodeIndex,
                                 std::size_t coordinateIndex) {
    incidenceLists[0].push_back(globalNodeIndex);
    globalCoordinateIndices.push_back(coordinateIndex);
}

template <std::size_t dimension>
template <std::size_t d>
std::enable_if_t<(d > 0)> Element<dimension>::addEntity(
    std::size_t globalIndex) {
    incidenceLists[d].push_back(globalIndex);
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
    return otherEntities.template get<actualDimension>();
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
    return otherEntities.template get<actualDimension>();
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
    if (number < otherEntities.template get<0>().size()) {
        otherEntities.template get<0>().resize(number);
    } else
        addNodes(number - otherEntities.template get<0>().size());
}

template <std::size_t dimension>
std::size_t Mesh<dimension>::addNode() {
    std::size_t newIndex = otherEntities.template get<0>().size();
    otherEntities.template get<0>().push_back({this, newIndex});
    return newIndex;
}

template <std::size_t dimension>
void Mesh<dimension>::addNodes(std::size_t count) {
    for (std::size_t i = 0; i < count; ++i) {
        addNode();
    }
}

template <std::size_t dimension>
std::size_t Mesh<dimension>::addNodeCoordinate(
    std::size_t nodeIndex, LinearAlgebra::SmallVector<dimension> coordinate) {
    coordinates.push_back({nodeIndex, coordinate});
    return coordinates.size() - 1;
}

template <std::size_t dimension>
void Mesh<dimension>::addElement(std::vector<std::size_t> nodeCoordinateIDs,
                                 const std::string& zoneName) {
    std::size_t elementID = elementsList.size();
    Element<dimension> newElement{this, elementID, getZoneId(zoneName)};
    newElement.setGeometry(findGeometry(nodeCoordinateIDs.size()));
    for (auto coordinateID : nodeCoordinateIDs) {
        std::size_t nodeID = coordinates[coordinateID].nodeIndex;
        newElement.addNode(nodeID, coordinateID);
    }
    elementsList.push_back(newElement);
    otherEntities.template get<dimension>().push_back(newElement);
    for (std::size_t i = 0; i < nodeCoordinateIDs.size(); ++i) {
        auto coordinateID = nodeCoordinateIDs[i];
        std::size_t nodeID = coordinates[coordinateID].nodeIndex;
        otherEntities.template get<0>()[nodeID].addElement(elementID, i);
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
    for (auto& nodeIndex : localNodeIndices) {
        nodeIndex =
            element.template getIncidentEntity<0>(nodeIndex).getGlobalIndex();
    }
    auto candidates =
        getNodes()[localNodeIndices[0]].template getIncidenceListAsIndices<d>();
    std::sort(candidates.begin(), candidates.end());
    for (std::size_t i = 1; i < localNodeIndices.size(); ++i) {
        auto newElements = getNodes()[localNodeIndices[i]]
                               .template getIncidenceListAsIndices<d>();
        std::sort(newElements.begin(), newElements.end());
        std::vector<std::size_t> temp;
        temp.reserve(std::max(candidates.size(), newElements.size()));
        std::set_intersection(candidates.begin(), candidates.end(),
                              newElements.begin(), newElements.end(),
                              std::back_inserter(temp));
        candidates = std::move(temp);
    }
    logger.assert_always(candidates.size() < 2, "mesh data is not consistent");
    std::size_t entityIndex =
        (candidates.size() == 1 ? candidates[0] : newEntity<d>());
    element.template addEntity<d>(entityIndex);
    otherEntities.template get<d>()[entityIndex].addElement(
        element.getGlobalIndex(), index);
}

template <std::size_t dimension>
const ElementShape<dimension>* Mesh<dimension>::findGeometry(
    std::size_t numberOfNodes) {
    for (auto shape : hpgemShapes.get<dimension>()) {
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
}  // namespace Preprocessor
