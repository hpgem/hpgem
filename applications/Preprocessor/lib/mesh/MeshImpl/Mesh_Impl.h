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

// NOTE: Purely implementations related to Mesh (for includes and structure
// see ../Mesh.h)

namespace Preprocessor {

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
    return meshEntities.template get<actualDimension>();
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
    return meshEntities.template get<actualDimension>();
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
EntityGId Mesh<dimension>::addNode() {
    return newEntity<0>();
}

template <std::size_t dimension>
void Mesh<dimension>::addNodes(std::size_t count) {
    for (std::size_t i = 0; i < count; ++i) {
        addNode();
    }
}

template <std::size_t dimension>
CoordId Mesh<dimension>::addNodeCoordinate(
    EntityGId nodeIndex, LinearAlgebra::SmallVector<dimension> coordinate) {
    coordinates.push_back({nodeIndex, coordinate});
    return CoordId(coordinates.size() - 1);
}

template <std::size_t dimension>
void Mesh<dimension>::addElement(std::vector<CoordId> nodeCoordinateIDs,
                                 const std::string& zoneName) {
    EntityGId elementID = EntityGId(elementsList.size());
    Element<dimension> newElement{this, elementID, getZoneId(zoneName)};
    newElement.setGeometry(findGeometry(nodeCoordinateIDs.size()));
    for (auto coordinateID : nodeCoordinateIDs) {
        EntityGId nodeID = coordinates[coordinateID.id].nodeIndex;
        newElement.addNode(nodeID, coordinateID);
    }
    elementsList.push_back(newElement);
    meshEntities.template get<dimension>().push_back(newElement);
    fixElement(elementsList.back(), tag<dimension - 1>{});
}

template <std::size_t dimension>
bool Mesh<dimension>::isValid() const {
    logger(DEBUG, "The mesh has % elements", getNumberOfElements());
    for (auto element : getElements()) {
        if (getElement(element.getGlobalIndex()) != element) {
            logger(ERROR, "The index for element % has been set incorrectly",
                   element.getGlobalIndex().id);
            return false;
        }
        if (!checkBoundingEntities(element, tag<dimension - 1>{})) {
            return false;
        }
    }
    return checkEntities(tag<dimension>{});
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
}

template <std::size_t dimension>
template <int d>
void Mesh<dimension>::removeUnusedEntities(itag<d> dimTag) {
    std::vector<EntityGId> renumbering;
    // Fill renumbering with invalid data
    std::size_t currentCount = this->meshEntities[dimTag].size();
    renumbering.resize(currentCount, EntityGId(-1));
    EntityGId newId = EntityGId(0);
    // Whether we are actually removing entities
    bool removing = false;
    for (std::size_t i = 0; i < currentCount; ++i) {
        MeshEntity<d, dimension>& entity = meshEntities[dimTag][i];
        if (entity.getNumberOfElements() > 0) {
            // Renumber
            renumbering[i] = newId;
            entity.entityID = newId;
            if (removing) {
                // Not only do we need to update the id, we also need to move
                // the entity to the corresponding location.
                meshEntities[dimTag][newId.id] = entity;
            }
            newId++;
        } else {
            removing = true;
        }
    }

    if (removing) {
        meshEntities[dimTag].resize(newId.id);

        for (Element<dimension>& element : elementsList) {
            element.renumberEntities(d, renumbering);
        }

        if (d == 0) {
            // For nodes we also need to remap the coordinate->node mapping
            for (coordinateData& coord : getNodeCoordinates()) {
                coord.nodeIndex = renumbering[coord.nodeIndex.id];
            }
        }
    }

    removeUnusedEntities(itag<d - 1>{});
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
