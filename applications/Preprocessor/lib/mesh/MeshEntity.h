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

#ifndef HPGEM_APP_MESH_ENTITY_H
#define HPGEM_APP_MESH_ENTITY_H

#include "idtypes.h"
#include "MeshPredeclarations.h"

#include <array>
#include <vector>

namespace Preprocessor {

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
    Element<meshDimension>& getElement(EntityLId i);
    const Element<meshDimension>& getElement(EntityLId i) const;

    /// Given a connected Element, what is its local index?
    EntityLId getElementIndex(const Element<meshDimension>& element) const;

    /// The number of elements that this MeshEntity is part of.
    std::size_t getNumberOfElements() const;

    /// The Local index of this MeshEntity on the i-th Element.
    EntityLId getLocalIndex(EntityLId i) const;

    /// Get the local index of this MeshEntity for an element that it is part
    /// of.
    /// \param element The element
    /// \return The local index on the element.
    EntityLId getLocalIndex(const Element<meshDimension>& element) const;

    /// The global index of this MeshEntity. For a given Mesh and
    /// entityDimension this uniquely determines the MeshEntity.
    EntityGId getGlobalIndex() const;

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
    std::vector<EntityGId> getIncidenceListAsIndices() const;

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
    MeshEntity(Mesh<meshDimension>* mesh, EntityGId entityID)
        : mesh(mesh), entityID(entityID) {}

    /// Add an element that is this MeshEntity is part of
    ///
    /// \param elementID The entityID of the element
    /// \param localEntityIndex The localIndex of this MeshEntity for the
    /// element.
    void addElement(EntityGId elementID, EntityLId EntityLId);

    Mesh<meshDimension>* mesh;
    /// The id of this MeshEntity
    EntityGId entityID = EntityGId(std::numeric_limits<std::size_t>::max());
    /// The entityIDs for the Elements that this MeshEntity is part of
    std::vector<EntityGId> elementIDs;
    /// For each of the elements that this MeshEntity is part of, the local id
    /// of this MeshEntity on the element.
    std::vector<EntityLId> localIDs;
    // NOTE: Probably the following identity should hold:
    // elements[elementIDs[i]].incidenceList[dimension][localIDs[i]] == entityID
};

}  // namespace Preprocessor

// Include all implementations
#include "MeshEntity_Impl.h"

#endif  // HPGEM_APP_MESH_ENTITY_H