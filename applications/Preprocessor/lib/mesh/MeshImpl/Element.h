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
#ifndef HPGEM_ELEMENT_H
#define HPGEM_ELEMENT_H

// NOTE: Purely definitions (for includes and structure see ../Mesh.h)

namespace Preprocessor {

/**
 * \brief An Element of the Mesh.
 *
 * The MeshEntity-s of largest dimension in the Mesh are called Elements. In
 * addition to the topological information they also include geometrical
 * information. Specifically:
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

    LinearAlgebra::SmallVector<dim> getCoordinate(EntityLId localIndex) const;
    CoordId getCoordinateIndex(EntityLId localIndex) const;

    std::vector<LinearAlgebra::SmallVector<dim>> getCoordinatesList() const;

    /// Set the global coordinate of a node
    /// \param localIndex The local index of the node
    /// \param newCoordinate The new coordinate.
    void setNodeCoordinate(EntityLId localIndex,
                           LinearAlgebra::SmallVector<dim> newCoordinate);

    /// Set the node at a specific local index.
    ///
    /// \param localIndex The local index of the node
    /// \param globalIndex The global index of the node
    /// \param coordinateIndex  The global coordinate index for the point
    /// corresponding to the node.
    void setNode(EntityLId localIndex, EntityGId globalIndex,
                 CoordId coordinateIndex);

    template <std::size_t d>
    std::enable_if_t<(d > 0)> setEntity(EntityLId localIndex,
                                        EntityGId globalIndex);

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
    std::vector<EntityGId> getIncidenceListAsIndices(
        const MeshEntity<entityDimension, dim>& entity) const;

    /// Same as getIncidenceList(const MeshEntity<entityDimension, dim>& entity)
    /// but computing the local indices of the MeshEntity-s
    /// \tparam d The dimension of the shared MeshEntity-s
    /// \tparam entityDimension The dimension of the input MeshEntity
    /// \param entity The MeshEntity on the boundary
    /// \return The local indices of the shared MeshEntity-s
    template <int d, std::size_t entityDimension>
    std::vector<EntityLId> getLocalIncidenceListAsIndices(
        const MeshEntity<entityDimension, dim>& entity) const;

    std::string getZoneName() const { return this->mesh->zoneNames[zoneId]; }

    size_t getZoneId() const { return zoneId; }

   private:
    friend Mesh<dim>;
    Element(Mesh<dim>* mesh, EntityGId elementID, std::size_t zoneId)
        : MeshEntity<dim, dim>(mesh, elementID), zoneId(zoneId) {
        this->addElement(elementID, 0);
    }

    /// Add a node to this element.
    /// \param globalNodeIndex The global index of the (topological) node
    /// \param coordinateIndex The global index of the coordinate for the node.
    void addNode(EntityGId globalNodeIndex, CoordId coordinateIndex);

    /// Add a MeshEntity on the boundary of this Element.
    ///
    /// @tparam d The dimension of the boundary MeshEntity
    /// @param globalIndex The globalIndex of the boundary MeshEntity
    template <std::size_t d>
    std::enable_if_t<(d > 0)> addEntity(EntityGId globalIndex);

    void setGeometry(const ElementShape<dim>* shape) {
        referenceGeometry = shape;
    }

    /// The geometry of this element (e.g. a cube or tetrahedron)
    const ElementShape<dim>* referenceGeometry;
    /// The global indices of all the MeshEntities that form the boundary of
    /// this element, grouped by the dimension of the MeshEntities.
    std::array<std::vector<EntityGId>, dim> incidenceLists;
    /// The global indices of the coordinates for the corners of this element.
    std::vector<CoordId> globalCoordinateIndices;

    std::size_t zoneId;
};

}  // namespace Preprocessor

#endif  // HPGEM_ELEMENT_H
