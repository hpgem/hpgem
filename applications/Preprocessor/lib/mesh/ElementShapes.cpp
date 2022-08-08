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

#include "ElementShapes.h"

namespace Preprocessor {

// Forward declare some helper functions.
// The implementation of these relies on the point and line shape, which are
// defined only afterwards.
namespace Detail {
/// Helper function to quickly generate a set of boundary vertices numbered 0 to
/// N-1
template <std::size_t dim>
std::vector<ElementShapePart<0, dim>> generateVertices(std::size_t numPoints);

/// Helper function to generate a set of boundary edges given pairs of
/// coordinates.
template <std::size_t dim>
std::vector<ElementShapePart<1, dim>> generateEdges(
    std::initializer_list<std::vector<std::size_t>> linePoints);
}  // namespace Detail

// note when debugging new shapes: invoking the logger during static
// initialisation can be a bit messy so you may need to use a debugger to see
// the error message

// NOTE ABOUT THE VALUES
// ---------------------
// These value define the hpGEM shapes, so these should match the definitions in
// the kernel. See for the actual pictures the Reference[Line,etc.] classes in
// the kernel (some additional pictures in the reference-mapping classes).
//
// For the details about the format here see the documentation of
// ElementShape(Part).

const ElementShape<0> point{};

const ElementShape<1> line{Detail::generateVertices<1>(2)};
// Quadratic line
const ElementShape<1> line2{Detail::generateVertices<1>(3)};

const ElementShape<2> triangle{
    // Boundary lines
    Detail::generateEdges<2>({{0, 1}, {0, 2}, {1, 2}}),
    // Corner vertices
    Detail::generateVertices<2>(3)};
const ElementShape<2> square{
    // Boundary lines
    Detail::generateEdges<2>({{0, 1}, {0, 2}, {1, 3}, {2, 3}}),
    // Corner vertices
    Detail::generateVertices<2>(4)};

const ElementShape<2> triangle2{
    std::vector<Detail::ElementShapePart<1, 2>>{
        Detail::ElementShapePart<1, 2>{&line2, {0, 1, 2}},
        Detail::ElementShapePart<1, 2>{&line2, {0, 3, 5}},
        Detail::ElementShapePart<1, 2>{&line2, {2, 4, 5}}},
    Detail::generateVertices<2>(6)};

const ElementShape<3> tetrahedron{
    // Faces
    std::vector<Detail::ElementShapePart<2, 3>>{
        Detail::ElementShapePart<2, 3>{&triangle, {0, 3, 2}},
        Detail::ElementShapePart<2, 3>{&triangle, {0, 1, 3}},
        Detail::ElementShapePart<2, 3>{&triangle, {0, 2, 1}},
        Detail::ElementShapePart<2, 3>{&triangle, {1, 2, 3}},
    },
    // Edges
    Detail::generateEdges<3>({{0, 1}, {0, 2}, {0, 3}, {2, 3}, {1, 3}, {1, 2}}),
    // Vertices
    Detail::generateVertices<3>(4)};

const ElementShape<3> cube{
    // Faces
    std::vector<Detail::ElementShapePart<2, 3>>{
        Detail::ElementShapePart<2, 3>{&square, {0, 1, 2, 3}},
        Detail::ElementShapePart<2, 3>{&square, {0, 1, 4, 5}},
        Detail::ElementShapePart<2, 3>{&square, {0, 2, 4, 6}},
        Detail::ElementShapePart<2, 3>{&square, {1, 3, 5, 7}},
        Detail::ElementShapePart<2, 3>{&square, {2, 3, 6, 7}},
        Detail::ElementShapePart<2, 3>{&square, {4, 5, 6, 7}},
    },
    // Edges
    Detail::generateEdges<3>({{0, 1},
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
                              {3, 7}}),
    // Vertices
    Detail::generateVertices<3>(8)};

const ElementShape<3> triangularPrism{
    // Faces
    std::vector<Detail::ElementShapePart<2, 3>>{
        Detail::ElementShapePart<2, 3>{&triangle, {0, 2, 1}},
        Detail::ElementShapePart<2, 3>{&triangle, {3, 4, 5}},
        Detail::ElementShapePart<2, 3>{&square, {2, 0, 5, 3}},
        Detail::ElementShapePart<2, 3>{&square, {0, 1, 3, 4}},
        Detail::ElementShapePart<2, 3>{&square, {1, 2, 4, 5}}},
    // Edges
    Detail::generateEdges<3>({{0, 1},
                              {0, 2},
                              {1, 2},
                              {3, 4},
                              {3, 5},
                              {4, 5},
                              {0, 3},
                              {1, 4},
                              {2, 5}}),
    // Vertices
    Detail::generateVertices<3>(6)};
const ElementShape<3> pyramid{
    // Faces
    std::vector<Detail::ElementShapePart<2, 3>>{
        Detail::ElementShapePart<2, 3>{&square, {3, 4, 1, 2}},
        Detail::ElementShapePart<2, 3>{&triangle, {3, 1, 0}},
        Detail::ElementShapePart<2, 3>{&triangle, {2, 4, 0}},
        Detail::ElementShapePart<2, 3>{&triangle, {1, 2, 0}},
        Detail::ElementShapePart<2, 3>{&triangle, {4, 3, 0}}

    },
    // Edges
    Detail::generateEdges<3>(
        {{0, 1}, {0, 2}, {0, 3}, {0, 4}, {1, 2}, {2, 4}, {4, 3}, {3, 1}}),
    // Vertices
    Detail::generateVertices<3>(5)};

const ElementShape<3> tetrahedron2{
    // Faces
    std::vector<Detail::ElementShapePart<2, 3>>{
        Detail::ElementShapePart<2, 3>{&triangle2, {0, 6, 9, 3, 8, 5}},
        Detail::ElementShapePart<2, 3>{&triangle2, {0, 1, 2, 6, 7, 9}},
        Detail::ElementShapePart<2, 3>{&triangle2, {0, 3, 5, 1, 4, 2}},
        Detail::ElementShapePart<2, 3>{&triangle2, {2, 4, 5, 7, 8, 9}},
    },
    // Edges
    std::vector<Detail::ElementShapePart<1, 3>>{
        Detail::ElementShapePart<1, 3>{&line2, {0, 1, 2}},
        Detail::ElementShapePart<1, 3>{&line2, {0, 3, 5}},
        Detail::ElementShapePart<1, 3>{&line2, {0, 6, 9}},
        Detail::ElementShapePart<1, 3>{&line2, {5, 8, 9}},
        Detail::ElementShapePart<1, 3>{&line2, {2, 7, 9}},
        Detail::ElementShapePart<1, 3>{&line2, {2, 4, 5}},
    },
    // Vertices
    Detail::generateVertices<3>(10)};

namespace Detail {

template <std::size_t dim>
std::vector<ElementShapePart<0, dim>> generateVertices(std::size_t numPoints) {
    std::vector<Detail::ElementShapePart<0, dim>> result(numPoints);
    for (std::size_t i = 0; i < numPoints; ++i) {
        result[i] = ElementShapePart<0, dim>{&point, {i}};
    }
    return result;
}

template <std::size_t dim>
std::vector<ElementShapePart<1, dim>> generateEdges(
    std::initializer_list<std::vector<std::size_t>> linePoints) {
    std::vector<ElementShapePart<1, dim>> result;
    for (auto& linePoint : linePoints) {
        result.emplace_back(&line, linePoint);
    }
    return result;
}
}  // namespace Detail

using Detail::ShapePointerVec;
const TemplateArray<4, ShapePointerVec> hpgemShapes{
    ShapePointerVec<3>{&tetrahedron, &cube, &triangularPrism, &pyramid,
                       &tetrahedron2},
    ShapePointerVec<2>{&triangle, &triangle2, &square},
    ShapePointerVec<1>{&line, &line2},
    ShapePointerVec<0>{&point},
};

}  // namespace Preprocessor