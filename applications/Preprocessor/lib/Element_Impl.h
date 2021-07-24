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
#ifndef HPGEM_ELEMENT_IMPL_H
#define HPGEM_ELEMENT_IMPL_H

namespace Preprocessor {

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
    incidenceLists[0][localIndex.id] = globalIndex;
    globalCoordinateIndices[localIndex.id] = coordinateIndex;
}

template <std::size_t dimension>
template <std::size_t d>
std::enable_if_t<(d > 0)> Element<dimension>::setEntity(EntityLId localIndex,
                                                        EntityGId globalIndex) {
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
}  // namespace Preprocessor

#endif  // HPGEM_ELEMENT_IMPL_H
