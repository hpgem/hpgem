/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2021, University of Twente
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
#ifndef HPGEM_MESHENTITY_IMPL_H
#define HPGEM_MESHENTITY_IMPL_H

// NOTE: Purely implementations related to MeshEntity (for includes and
// structure see ../Mesh.h)

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
EntityLId MeshEntity<entityDimension, meshDimension>::getElementIndex(
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
    logger(ERROR, "Element not found");
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
template <std::size_t d>
std::enable_if_t<d == entityDimension>
    MeshEntity<entityDimension, meshDimension>::merge(
        MeshEntity<entityDimension, meshDimension>& source, tag<d>) {
    // Merging elements is fundamentally more complicated and should not be done
    // via this method.
    static_assert(entityDimension < meshDimension,
                  "Can't use this to merge elements");

    logger.assert_debug(source.getMesh() == getMesh(),
                        "Merging from different meshes");
    logger.assert_debug(this->getGlobalIndex() != source.getGlobalIndex(),
                        "Can't merge an entity into itself");
    // Loop to move all connections
    // Internally removes the links from the source element.
    while (!source.elementIDs.empty()) {
        Element<meshDimension>& element = source.getElement(0);
        EntityLId localId = source.getLocalIndex(0);
        if (entityDimension != 0) {
            element.template setEntity<entityDimension>(localId,
                                                        getGlobalIndex());
        }
    }
}

template <std::size_t entityDimension, std::size_t dimension>
void MeshEntity<entityDimension, dimension>::merge(
    MeshEntity<entityDimension, dimension>& source, tag<0>) {
    logger.assert_debug(source.getMesh() == getMesh(),
                        "Merging from different meshes");
    logger.assert_debug(this->getGlobalIndex() != source.getGlobalIndex(),
                        "Can't merge an entity into itself");

    // Loop to move all connections
    // Internally removes the links from the source element.
    while (!source.elementIDs.empty()) {
        Element<dimension>& element = source.getElement(0);
        EntityLId localId = source.getLocalIndex(0);

        // Merging nodes requires that we also update the corresponding
        // coordinates.
        CoordId coordId = element.getCoordinateIndex(localId);
        // The coordinate should belong to either the source MeshEntity or
        // the current one. Check this for robustness.
        EntityGId& nodeIndex = mesh->getNodeCoordinates()[coordId.id].nodeIndex;
        logger.assert_always(
            nodeIndex == source.getGlobalIndex() ||
                nodeIndex == getGlobalIndex(),
            "Coordinate belongs to neither the old nor new node.");
        nodeIndex = getGlobalIndex();

        element.setNode(localId, getGlobalIndex(), coordId);
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
    EntityGId elementID, EntityLId localEntityIndex) {
    logger.assert_debug(elementIDs.size() == localIDs.size(),
                        "ElementIDs and LocalIDs have different length");
    // Simultaneous iteration
    auto elemIter = elementIDs.begin();
    auto localIter = localIDs.begin();
    for (; elemIter != elementIDs.end(); ++elemIter, ++localIter) {
        if (*elemIter == elementID && *localIter == localEntityIndex) {
            elementIDs.erase(elemIter);
            localIDs.erase(localIter);
            return;
        }
    }
    logger.assert_always(false,
                         "Element % (localID %) not found on this MeshEntity",
                         elementID, localEntityIndex);
}

template <std::size_t entityDimension, std::size_t meshDimension>
bool MeshEntity<entityDimension, meshDimension>::operator==(
    const MeshEntity& other) const {
    return mesh == other.mesh && entityID == other.entityID &&
           elementIDs == other.elementIDs && localIDs == other.localIDs;
}

}  // namespace Preprocessor

#endif  // HPGEM_MESHENTITY_IMPL_H
