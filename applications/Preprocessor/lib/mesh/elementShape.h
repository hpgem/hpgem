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

#ifndef HPGEM_APP_ELEMENTSHAPE_H
#define HPGEM_APP_ELEMENTSHAPE_H

#include <array>
#include <cstddef>
#include <numeric>
#include <vector>
#include "Logger.h"

#include "utils/tag.h"
#include "utils/TemplateArray.h"

namespace Preprocessor {

using namespace hpgem;

template <std::size_t dimension>
class ElementShape;

namespace Detail {

/// Part of an ElementShape
///
/// Part of an ElementShape, this could be the whole shape, or a part of the
/// boundary (i.e. a vertex, edge or face).
///
/// \tparam partDimension The dimension of the part
/// \tparam dimension The dimension of the original shape
template <std::size_t partDimension, std::size_t dimension>
class ElementShapePart {
   public:
    /// Define part of an ElementShape
    ///
    /// Each part is characterized by three things.
    ///  1. The shape of the part
    ///  2. The nodes of the parent shape used in this part
    ///  3. The mapping of node numbers of the part-shape to that of the
    ///     parent-shape.
    /// For example the first triangular face of the tetrahedron is defined by
    /// 1. shape: Triangle
    /// 2. nodes: {0,2,3}
    /// 3. mapping (triangle node -> tetrahedron node): 0 -> 0, 1 -> 2, 2 -> 3
    /// With the other three faces having a different subset of nodes, and a
    /// corresponding mapping.
    ///
    /// The mapping and subset of nodes are defined in the same variable, by
    /// listing the node numbers of the parent shape in the order of the
    /// mapping. So in the example they should be specified as {0,3,2}.
    ///
    /// WARNING: The information should match the definition in the particular
    ///   ReferenceGeometry in the kernel
    /// \param shape The geometrical shape of the part
    ///
    /// \param nodeIds The subset of nodes of the parent shape that are used in
    /// this part, in the right order.
    ElementShapePart(const ElementShape<partDimension>* shape,
                     std::vector<std::size_t> nodeIds)
        : shape(shape), adjacentShapes() {
        adjacentShapes[0] = nodeIds;
        logger.assert_always(
            shape->getNumberOfNodes() == nodeIds.size(),
            "Incorrect number of nodes provided, expected %, got %",
            shape->getNumberOfNodes(), nodeIds.size());
    };
    ElementShapePart() = default;

    const ElementShape<partDimension>* getReferenceGeometry() const { return shape; }

    std::vector<std::size_t>& getAdjacentShapes(std::size_t targetDimension) {
        logger.assert_debug(targetDimension <= dimension,
                            "Invalid target dimension");
        return adjacentShapes[targetDimension];
    }
    const std::vector<std::size_t>& getAdjacentShapes(
        std::size_t targetDimension) const {
        logger.assert_debug(targetDimension <= dimension,
                            "Invalid target dimension");
        return adjacentShapes[targetDimension];
    }

    const std::vector<std::size_t>& getNodes() const {
        return getAdjacentShapes(0);
    }

   private:
    /// The shape of this part of the boundary
    const ElementShape<partDimension>* shape;
    /// The adjacent boundary parts, as indices into the the list of boundary
    /// parts of the corresponding element. Grouped by the dimension of the
    /// specific part. This results in three cases:
    /// adjacentDim < entityDim
    ///   -> adjacentShape is part of the boundary of this boundary part.
    /// adjacentDim == entityDim
    ///    -> This boundary part itself
    /// adjacentDim > entityDim
    ///    -> This boundary part is on the boundary of the adjacentShape
    std::array<std::vector<std::size_t>, dimension + 1> adjacentShapes;
};
}  // namespace Detail

/// @brief Description of the topology of the possible elements in a mesh
///
/// The ElementShape defines the shape of the possible reference elements in
/// the mesh. Which are assumed to be polytopes (point, line, polygon,
/// polyhedra). Such a reference polytope has several properties that are used
/// in a FEM code:
///
///  1. A coordinate system, describing the geometry of the reference element
///     and used for the basis functions.
///  2. The topology of the reference element, as a polytope. That is how are
///     the parts (vertices, edges, etc.) connected to form the polytope.
///  3. A numbering of the boundary parts so that each part can be indexed and
///     their orientation.
/// In the kernel we don't need the geometry part, but we do need the topology
/// and the indices of the parts.
///
/// To define the topology of the polytope we use two properties:
///   1. It's boundary is a combination of 0, 1, ... D-1, dimensional polytopes
///   2. It has N vertices (0 dimensional polytopes).
/// The vertices are numbered 0...N-1. The intermediate polytopes on the
/// boundary are then defined by taking a subset of the vertices and providing
/// the associated polytopical shape. The ordering of the vertices then defines
/// the relative orientation.
///
/// Note: As both of these properties influence the output, they should match
/// the definition of the shapes in the kernel.
///
/// \tparam dimension The dimension of this shape
template <std::size_t dimension>
class ElementShape {

   public:
    ElementShape() = default;

    ~ElementShape() = default;

    // Disallow copying and moving for two reasons:
    // 1. It contains a self reference as part of its shapeParts, this would
    //    need a custom copying and moving routine to keep this reference
    //    correct
    // 2. Larger shapes store pointers to smaller shapes (for the boundary).
    //    By moving and copying one could easily break these pointers.

    ElementShape(const ElementShape&) = delete;

    ElementShape(ElementShape&&) noexcept = delete;

    ElementShape& operator=(const ElementShape&) = delete;

    ElementShape& operator=(ElementShape&&) noexcept = delete;

    /// Define a shape by it's boundary parts
    ///
    /// \tparam Args The types of the boundary, should be
    ///    BoundaryPart<dimension-1,dimension>, ... BoundaryPart<1,dimension>,
    ///    BoundaryPart<0,dimension>
    ///
    /// \param args The boundary parts, starting at the highest dimension (codim
    /// 1) and going down to the vertices.
    template <typename... Args>
    ElementShape(Args... args) : shapeParts({}, args...) {
        // Insert ElementShapePart corresponding to the whole element.
        auto& shapeVec = shapeParts.template get<dimension>();
        std::vector<std::size_t> vertices(getNumberOfNodes());
        std::iota(vertices.begin(), vertices.end(), 0);
        shapeVec.push_back(
            Detail::ElementShapePart<dimension, dimension>{this, vertices});
        // Compute the topological connections
        completeSubShapes();
        logger.assert_debug(checkShape(),
                            "Input generated an inconsistent shape");
    }

    std::size_t getNumberOfNodes() const { return getNumberOfEntities<0>(); }

    std::size_t getNumberOfEdges() const { return getNumberOfEntities<1>(); }

    std::size_t getNumberOfFaces() const {
        return getNumberOfEntities<dimension - 1>();
    }

    // if (entityDimension >= 0)
    //      ElementShape<entityDimension>
    // else if(entityDimension + dimension >= 0)
    //      ElementShape<entityDimension + dimension>
    // codim case else ElementShape<0>
    // make the compiler not crash while we invoke the logger
    template <int entityDimension>
    using ShapeType = const ElementShape<(
        entityDimension < 0
            ? (entityDimension + dimension < 0 ? 0
                                               : entityDimension + dimension)
            : entityDimension)>;

    template <std::size_t entityDimension>
    std::vector<std::size_t> getNodesOfEntity(std::size_t entityIndex) const {
        return getAdjacentEntities<entityDimension, 0>(entityIndex);
    }

    /**
     * @brief returns the number of entities of the specified dimension.
     */
    template <std::size_t entityDimension>
    std::enable_if_t<(entityDimension <= dimension), std::size_t>
        getNumberOfEntities() const;

    template <std::size_t entityDimension>
    std::enable_if_t<(entityDimension <= dimension),
                     const ShapeType<entityDimension>*>
        getBoundaryShape(std::size_t entityIndex) const;

    /**
     * @brief Lookup adjacency information about entities of this shape
     *
     * Given an entity on this shape find the adjacent entities
     * from a specific target dimension. For example let the current shape be a
     * triangle, and the entity be an edge of this tetrahedron. Then depending
     * on the target dimension one would get:
     *  0. The nodes connected to the edge
     *  1. The edge itself
     *  2. The triangular faces which have this edge
     *  3. The tetrahedron
     *
     * @tparam entityDimension The dimension of the entity
     * @tparam targetDimension The dimension of the entities to look to
     * @param entityIndex The index specifying the specific entity to find the
     * adjacency of.
     * @return The indices of the `targetDimension` dimensional entities that
     * are adjacent to the input entity.
     */

    template <std::size_t entityDimension, std::size_t targetDimension>
    std::enable_if_t<(entityDimension <= dimension &&
                      targetDimension <= dimension),
                     std::vector<std::size_t>>
        getAdjacentEntities(std::size_t entityIndex) const;

    bool checkShape() const { return checkBoundaryShape(itag<dimension>{}); }

   private:
    /// For a shape S on the boundary of this elementShape, do some consistency
    /// checks on the boundary of S.
    ///
    /// \tparam d The dimension of the entities to check on the boundary of S
    /// \tparam shapeDimension The dimension of S
    /// \param boundingShape The bounding shape S
    /// \param boundaryNodes The nodes for the S (subset of the nodes of the
    /// elementshape)
    ///
    /// \param currentIndex Index of S with respect to the elementShape.
    ///
    /// \return Whether correct
    template <int d, std::size_t shapeDimension>
    bool checkSingleShape(const ElementShape<shapeDimension>* boundingShape,
                          std::vector<std::size_t> boundaryNodes,
                          std::size_t currentIndex, itag<d>) const;

    template <std::size_t shapeDimension>
    bool checkSingleShape(const ElementShape<shapeDimension>* boundingShape,
                          std::vector<std::size_t> boundaryNodes,
                          std::size_t currentIndex, itag<-1>) const;

    /// Check the d-dimensional boundary of this shape
    template <int d>
    bool checkBoundaryShape(itag<d>) const;

    bool checkBoundaryShape(itag<-1>) const;

    /// Derive the topological connections from the vertices on each subshape.
    void completeSubShapes() {
        completeSubShapes(tag<dimension>{}, tag<dimension>{});
    }

    /// Derive the topological connections of the boundary parts.
    ///
    /// \tparam entityDimension The dimension of the boundary parts to derive
    /// the topological connections for.
    ///
    /// \tparam targetDimension The dimension of the boundary parts to which
    /// they are connected.
    template <std::size_t entityDimension, std::size_t targetDimension>
    void completeSubShapes(tag<entityDimension>, tag<targetDimension>);

    template <std::size_t entityDimension>
    void completeSubShapes(tag<entityDimension>, tag<0>);

    void completeSubShapes(tag<0>, tag<0>) {}

    /// Specialization to fix the dimension of the entity shape
    template <std::size_t partDimension>
    using ShapePartVec =
        std::vector<Detail::ElementShapePart<partDimension, dimension>>;

    TemplateArray<dimension + 1, ShapePartVec> shapeParts;
};

// Special case for the zero dimension shape: The point. Unlike higher
// dimensional shapes it does not have a boundary consisting of lower
// dimensional shapes.
template <>
class ElementShape<0> {
   public:
    ElementShape() = default;

    ~ElementShape() = default;

    ElementShape(const ElementShape&) = default;

    ElementShape(ElementShape&&) noexcept = default;

    ElementShape& operator=(const ElementShape&) = default;

    ElementShape& operator=(ElementShape&&) noexcept = default;

    std::size_t getNumberOfNodes() const { return 1; }

    std::size_t getNumberOfEdges() const { return 0; }

    std::size_t getNumberOfFaces() const { return 0; }

    // Templated to match the generic case.
    template <int entityDimension>
    using ShapeType = ElementShape<0>;

    template <int entityDimension>
    std::vector<std::size_t> getNodesOfEntity(std::size_t entityIndex) const {
        return getAdjacentEntities<entityDimension, 0>(entityIndex);
    }

    /**
     * @brief returns the number of entities of the specified dimension.
     * Negative dimensions are treated as codimensions
     */
    template <std::size_t entityDimension>
    std::size_t getNumberOfEntities() const {
        if (entityDimension == 0) return 1;
        return 0;
    }

    template <std::size_t entityDimension>
    ShapeType<entityDimension> getBoundaryShape(std::size_t entityIndex) const {
        logger.assert_debug(
            entityDimension == 0,
            "A point is not bounded by shapes of other dimensions");
        logger.assert_debug(
            entityIndex == 0,
            "A point consists of only 1 shape, but you asked for shape %",
            entityIndex);
        return *this;
    }

    template <std::size_t entityDimension, std::size_t targetDimension>
    std::vector<std::size_t> getAdjacentEntities(
        std::size_t entityIndex) const {
        logger.assert_debug(
            entityDimension == 0,
            "A point is not bounded by shapes of other dimensions");
        logger.assert_debug(
            entityIndex == 0,
            "A point consists of only 1 shape, but you asked for shape %",
            entityIndex);
        // The shape itself
        if (targetDimension == 0) return {0};

        return {};
    };

    // points are hardcoded to be correct
    bool checkShape() const { return true; }
};

}  // namespace Preprocessor

#include "elementShape_impl.h"

#endif  // HPGEM_APP_ELEMENTSHAPE_H
