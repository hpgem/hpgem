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

#include <cstddef>
#include <vector>
#include "Logger.h"

#include "tag.h"

namespace Preprocessor {

using namespace hpgem;

template <std::size_t dimension>
class ElementShape;

namespace Detail {

template <std::size_t entityDimension, std::size_t dimension>
struct EntityData : public EntityData<entityDimension - 1, dimension> {
    template <typename... superArgs>
    EntityData(std::vector<std::vector<std::size_t>> adjacentNodes,
               std::vector<const ElementShape<entityDimension>*> entityShapes,
               superArgs... args)
        : EntityData<entityDimension - 1, dimension>(args...),
          adjacentShapes(),
          entityShapes(entityShapes) {
        adjacentShapes[0] = adjacentNodes;
    }

    EntityData() = default;

    ~EntityData() = default;

    EntityData(const EntityData&) = default;

    EntityData(EntityData&&) noexcept = default;

    EntityData& operator=(const EntityData&) = default;

    EntityData& operator=(EntityData&&) noexcept = default;

    // - Target dimension
    // - Shape
    // -> indices
    std::array<std::vector<std::vector<std::size_t>>, dimension> adjacentShapes;
    std::vector<const ElementShape<entityDimension>*> entityShapes;

    template <std::size_t SUB_DIM>
    std::vector<const ElementShape<SUB_DIM>*>& getEntityShapes() {
        return EntityData<SUB_DIM, dimension>::entityShapes;
    }

    template <std::size_t SUB_DIM>
    const std::vector<const ElementShape<SUB_DIM>*>& getEntityShapes() const {
        return EntityData<SUB_DIM, dimension>::entityShapes;
    }

    template <std::size_t SUB_DIM>
    std::array<std::vector<std::vector<std::size_t>>, dimension>&
        getAdjacentShapes() {
        return this->EntityData<SUB_DIM, dimension>::adjacentShapes;
    }

    template <std::size_t SUB_DIM>
    const std::array<std::vector<std::vector<std::size_t>>, dimension>&
        getAdjacentShapes() const {
        return EntityData<SUB_DIM, dimension>::adjacentShapes;
    }
};

template <std::size_t dimension>
struct EntityData<0, dimension> {
    EntityData(std::vector<std::vector<std::size_t>> adjacentNodes,
               std::vector<const ElementShape<0>*> entityShapes)
        : adjacentShapes(), entityShapes(entityShapes) {
        adjacentShapes[0] = adjacentNodes;
    }

    EntityData() = default;

    ~EntityData() = default;

    EntityData(const EntityData&) = default;

    EntityData(EntityData&&) noexcept = default;

    EntityData& operator=(const EntityData&) = default;

    EntityData& operator=(EntityData&&) noexcept = default;

    std::array<std::vector<std::vector<std::size_t>>, dimension> adjacentShapes;
    std::vector<const ElementShape<0>*> entityShapes;

    template <std::size_t SUB_DIM>
    std::vector<const ElementShape<SUB_DIM>*>& getEntityShapes() {
        return EntityData<SUB_DIM, dimension>::entityShapes;
    }

    template <std::size_t SUB_DIM>
    std::array<std::vector<std::vector<std::size_t>>, dimension>&
        getAdjacentShapes() {
        return this->EntityData<SUB_DIM, dimension>::adjacentShapes;
    }

    template <std::size_t SUB_DIM>
    const std::vector<const ElementShape<SUB_DIM>*>& getEntityShapes() const {
        return EntityData<SUB_DIM, dimension>::entityShapes;
    }

    template <std::size_t SUB_DIM>
    const std::array<std::vector<std::vector<std::size_t>>, dimension>&
        getAdjacentShapes() const {
        return this->EntityData<SUB_DIM, dimension>::adjacentShapes;
    }
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

    ElementShape(const ElementShape&) = default;

    ElementShape(ElementShape&&) noexcept = default;

    ElementShape& operator=(const ElementShape&) = default;

    ElementShape& operator=(ElementShape&&) noexcept = default;

    template <typename... Args>
    ElementShape(Args... args) : shapeData(args...) {
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
    std::enable_if_t<(entityDimension == dimension), std::size_t>
        getNumberOfEntities() const;
    template <std::size_t entityDimension>
    std::enable_if_t<(entityDimension >= 0 && entityDimension < dimension),
                     std::size_t>
        getNumberOfEntities() const;

    template <std::size_t entityDimension>
    std::enable_if_t<(entityDimension == dimension),
                     const ShapeType<entityDimension>*>
        getBoundaryShape(std::size_t entityIndex) const;
    template <std::size_t entityDimension>
    std::enable_if_t<(entityDimension < dimension),
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
    std::enable_if_t<(entityDimension == dimension ||
                       targetDimension == dimension),
                     std::vector<std::size_t>>
        getAdjacentEntities(std::size_t entityIndex) const;

    template <std::size_t entityDimension, std::size_t targetDimension>
    std::enable_if_t<(entityDimension < dimension &&
                      targetDimension < dimension),
                     std::vector<std::size_t>>
        getAdjacentEntities(std::size_t entityIndex) const;

    bool checkShape() const { return checkBoundaryShape(itag<dimension - 1>{}); }

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
        completeSubShapes(tag<dimension - 1>{}, tag<dimension - 1>{});
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

    Detail::EntityData<dimension - 1, dimension> shapeData;
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

// note when debugging new shapes: invoking the logger during static
// initialisation can be a bit messy so you may need to use a debugger to see
// the error message
const ElementShape<0> point{};
const ElementShape<1> line(std::vector<std::vector<std::size_t>>{{0}, {1}},
                           std::vector<const ElementShape<0>*>{&point, &point});
const ElementShape<2> triangle(
    std::vector<std::vector<std::size_t>>{{0, 1}, {0, 2}, {1, 2}},
    std::vector<const ElementShape<1>*>{&line, &line, &line},
    std::vector<std::vector<std::size_t>>{{0}, {1}, {2}},
    std::vector<const ElementShape<0>*>{&point, &point, &point});
const ElementShape<2> square(
    std::vector<std::vector<std::size_t>>{{0, 1}, {0, 2}, {1, 3}, {2, 3}},
    std::vector<const ElementShape<1>*>{&line, &line, &line, &line},
    std::vector<std::vector<std::size_t>>{{0}, {1}, {2}, {3}},
    std::vector<const ElementShape<0>*>{&point, &point, &point, &point});
const ElementShape<3> tetrahedron(
    std::vector<std::vector<std::size_t>>{
        {0, 3, 2}, {0, 1, 3}, {0, 2, 1}, {1, 2, 3}},
    std::vector<const ElementShape<2>*>{&triangle, &triangle, &triangle,
                                        &triangle},
    std::vector<std::vector<std::size_t>>{
        {0, 1}, {0, 2}, {0, 3}, {2, 3}, {1, 3}, {1, 2}},
    std::vector<const ElementShape<1>*>{&line, &line, &line, &line, &line,
                                        &line},
    std::vector<std::vector<std::size_t>>{{0}, {1}, {2}, {3}},
    std::vector<const ElementShape<0>*>{&point, &point, &point, &point});
const ElementShape<3> cube(
    std::vector<std::vector<std::size_t>>{{0, 1, 2, 3},
                                          {0, 1, 4, 5},
                                          {0, 2, 4, 6},
                                          {1, 3, 5, 7},
                                          {2, 3, 6, 7},
                                          {4, 5, 6, 7}},
    std::vector<const ElementShape<2>*>{&square, &square, &square, &square,
                                        &square, &square},
    std::vector<std::vector<std::size_t>>{{0, 1},
                                          {2, 3},
                                          {4, 5},
                                          {6, 7},
                                          {0, 2},
                                          {1, 3},
                                          {4, 6},
                                          {5, 7},
                                          {0, 4},
                                          {1, 5},
                                          {2, 6},
                                          {3, 7}},
    std::vector<const ElementShape<1>*>{&line, &line, &line, &line, &line,
                                        &line, &line, &line, &line, &line,
                                        &line, &line},
    std::vector<std::vector<std::size_t>>{
        {0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}},
    std::vector<const ElementShape<0>*>{&point, &point, &point, &point, &point,
                                        &point, &point, &point});
const ElementShape<3> triangularPrism(
    std::vector<std::vector<std::size_t>>{
        {0, 2, 1}, {3, 4, 5}, {2, 0, 5, 3}, {0, 1, 3, 4}, {1, 2, 4, 5}},
    std::vector<const ElementShape<2>*>{&triangle, &triangle, &square, &square,
                                        &square},
    std::vector<std::vector<std::size_t>>{
        {0, 1}, {0, 2}, {1, 2}, {3, 4}, {3, 5}, {4, 5}, {0, 3}, {1, 4}, {2, 5}},
    std::vector<const ElementShape<1>*>{&line, &line, &line, &line, &line,
                                        &line, &line, &line, &line},
    std::vector<std::vector<std::size_t>>{{0}, {1}, {2}, {3}, {4}, {5}},
    std::vector<const ElementShape<0>*>{&point, &point, &point, &point, &point,
                                        &point});
const ElementShape<3> pyramid(
    std::vector<std::vector<std::size_t>>{
        {3, 4, 1, 2}, {3, 1, 0}, {2, 4, 0}, {1, 2, 0}, {4, 3, 0}},
    std::vector<const ElementShape<2>*>{&square, &triangle, &triangle,
                                        &triangle, &triangle},
    std::vector<std::vector<std::size_t>>{
        {0, 1}, {0, 2}, {0, 3}, {0, 4}, {1, 2}, {2, 4}, {4, 3}, {3, 1}},
    std::vector<const ElementShape<1>*>{&line, &line, &line, &line, &line,
                                        &line, &line, &line},
    std::vector<std::vector<std::size_t>>{{0}, {1}, {2}, {3}, {4}},
    std::vector<const ElementShape<0>*>{&point, &point, &point, &point,
                                        &point});

template <std::size_t dimension>
const std::vector<const ElementShape<dimension>*> defaultShapes;

template <>
const std::vector<const ElementShape<0>*> defaultShapes<0> = {&point};

template <>
const std::vector<const ElementShape<1>*> defaultShapes<1> = {&line};

template <>
const std::vector<const ElementShape<2>*> defaultShapes<2> = {&triangle,
                                                              &square};

template <>
const std::vector<const ElementShape<3>*> defaultShapes<3> = {
    &tetrahedron, &cube, &triangularPrism, &pyramid};

}  // namespace Preprocessor

#include <elementShape_impl.h>

#endif  // HPGEM_APP_ELEMENTSHAPE_H
