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

#include <Logger.h>
#include <numeric>
#include "elementShape.h"
#include <algorithm>

namespace Preprocessor {


template <std::size_t dimension>
template <std::size_t entityDimension>
std::enable_if_t<(entityDimension == dimension), std::size_t>
    ElementShape<dimension>::getNumberOfEntities() const {
    return 1;
}

template <std::size_t dimension>
template <std::size_t entityDimension>
std::enable_if_t<(entityDimension >= 0 && entityDimension < dimension),
                 std::size_t>
    ElementShape<dimension>::getNumberOfEntities() const {
    return shapeData.template getEntityShapes<entityDimension>().size();
}

template <std::size_t dimension>
template <std::size_t entityDimension>
auto ElementShape<dimension>::getBoundaryShape(std::size_t entityIndex) const
    -> std::enable_if_t<(entityDimension == dimension),
                        const ShapeType<entityDimension>*> {
    logger.assert_debug(entityIndex == 0,
                        "This shape is bounded by only 1 shapes of dimension "
                        "%, but you asked for shape %",
                        dimension, entityIndex);
    return *this;
}

template <std::size_t dimension>
template <std::size_t entityDimension>
auto ElementShape<dimension>::getBoundaryShape(std::size_t entityIndex) const
    -> std::enable_if_t<(entityDimension < dimension),
                        const ShapeType<entityDimension>*> {
    logger.assert_debug(entityIndex < getNumberOfEntities<entityDimension>(),
                        "This shape is bounded by only % shapes of dimension "
                        "%, but you asked for shape %",
                        getNumberOfEntities<entityDimension>(), entityDimension,
                        entityIndex);
    return shapeData.template getEntityShapes<entityDimension>()[entityIndex];
}

template <std::size_t dimension>
template <std::size_t entityDimension, std::size_t targetDimension>
std::enable_if_t<(entityDimension == dimension || targetDimension == dimension),
                 std::vector<std::size_t>>
    ElementShape<dimension>::getAdjacentEntities(
        std::size_t entityIndex) const {
    logger.assert_debug(entityDimension <= dimension,
                        "This shape is not bounded by shapes of dimension %",
                        entityDimension + dimension);
    logger.assert_debug(entityIndex < getNumberOfEntities<entityDimension>(),
                        "This shape is bounded by only % shapes of dimension "
                        "%, but you asked for shape %",
                        getNumberOfEntities<entityDimension>(), dimension,
                        entityIndex);
    if (targetDimension > dimension) return {};
    if (targetDimension == dimension) return {0};
    std::size_t amount = getNumberOfEntities<targetDimension>();
    std::vector<std::size_t> result(amount);
    std::iota(result.begin(), result.end(), 0);
    return result;
}

template <std::size_t dimension>
template <std::size_t entityDimension, std::size_t targetDimension>
std::enable_if_t<(entityDimension < dimension && targetDimension < dimension),
                 std::vector<std::size_t>>
    ElementShape<dimension>::getAdjacentEntities(
        std::size_t entityIndex) const {
    logger.assert_debug(entityIndex < getNumberOfEntities<entityDimension>(),
                        "This shape is bounded by only % shapes of dimension "
                        "%, but you asked for shape %",
                        getNumberOfEntities<entityDimension>(), dimension,
                        entityIndex);
    //    return boundaryShapes.template
    //    get<entityDimension>()[entityIndex].getAdjacentShapes(targetDimension);
    return shapeData.template getAdjacentShapes<
        entityDimension>()[targetDimension][entityIndex];
}

template <std::size_t dimension>
template <int d, std::size_t shapeDimension>
bool ElementShape<dimension>::checkSingleShape(
    const ElementShape<shapeDimension>* boundingShape,
    std::vector<std::size_t> boundaryNodes, std::size_t currentIndex,
    itag<d>) const {
    // We deal here with several shapes (in decreasing dimension)
    // 1. The elementShape E
    // 2. A shape S on the boundary of E, which has index 'currentIndex'
    // 3. A dim-d shape S' on the boundary of S.
    // 4. A dim-d shape E' on the boundary of E that corresponds to S'.
    for (std::size_t i = 0;
         i < boundingShape->template getNumberOfEntities<d>(); ++i) {
        // Compute the number of the nodes of S' with respect to E.
        auto localNodes = boundingShape->template getAdjacentEntities<d, 0>(i);
        std::vector<std::size_t> cellNodes(localNodes.size());
        for (std::size_t j = 0; j < localNodes.size(); ++j) {
            cellNodes[j] = boundaryNodes[localNodes[j]];
        }
        std::sort(cellNodes.begin(), cellNodes.end());
        bool match = false;
        // Find E' by looking for the shape that has the same nodes
        // up to permutation
        for (std::size_t j = 0;
             j < shapeData.template getAdjacentShapes<d>()[0].size(); ++j) {
            auto comparisonNodes =
                shapeData.template getAdjacentShapes<d>()[0][j];
            std::sort(comparisonNodes.begin(), comparisonNodes.end());
            if (cellNodes == comparisonNodes) {
                // Found E'
                match = true;
                auto candindateIndices =
                    getAdjacentEntities<d, shapeDimension>(j);
                // Check that the topological information for E' lists S as a
                // neighbour.
                if (std::find(candindateIndices.begin(),
                              candindateIndices.end(),
                              currentIndex) == candindateIndices.end()) {
                    logger(ERROR,
                           "The %-dimensional shape % bounding %-dimensional "
                           "shape % is not adjacent to it",
                           d, j, shapeDimension, currentIndex);
                    return false;
                }
            }
        }
        if (!match) {
            // No E' found -> error.
            logger(ERROR,
                   "The %-dimensional bounding shape % contains a bounding "
                   "shape that is not a bounding shape of this entity",
                   shapeDimension, currentIndex);
        }
    }
    return checkSingleShape(boundingShape, boundaryNodes, currentIndex,
                            itag<d - 1>{});
}

template <std::size_t dimension>
template <std::size_t shapeDimension>
bool ElementShape<dimension>::checkSingleShape(
    const ElementShape<shapeDimension>* boundingShape,
    std::vector<std::size_t> boundaryNodes, std::size_t currentIndex,
    itag<-1>) const {
    return true;
}

template <std::size_t dimension>
template <int d>
bool ElementShape<dimension>::checkBoundaryShape(itag<d>) const {
    // Check 1: The number of boundary shapes must match the number of boundary
    // shapes for which we have topology information
    if (shapeData.template getEntityShapes<d>().size() !=
        shapeData.template getAdjacentShapes<d>()[0].size()) {
        logger(ERROR,
               "There are % %-dimensional bounding shapes described by their "
               "bounding nodes, but % %-dimensional shapes described by their "
               "topological shape",
               shapeData.template getEntityShapes<d>().size(), d,
               shapeData.template getAdjacentShapes<d>()[0].size(), d);
        return false;
    }
    // Check 2: The individual boundary shapes
    for (std::size_t i = 0; i < getNumberOfEntities<d>(); ++i) {

        auto boundingShape = getBoundaryShape<d>(i);
        // Check 2a: Does the shape itself make sense?
        if (!boundingShape->checkShape()) {
            logger(ERROR,
                   "%-dimensional bounding shape % has an unspecified "
                   "consistency error",
                   d, i);
            return false;
        }
        // Check 2b: The shape has a number of nodes, check that the topological
        // information has enough nodes.
        auto boundaryNodes = shapeData.template getAdjacentShapes<d>()[0][i];
        if (boundaryNodes.size() != boundingShape->getNumberOfNodes()) {
            logger(ERROR,
                   "%-dimensional shape % should be described be % nodes, but "
                   "it has % nodes listed",
                   d, i, boundingShape->getNumberOfNodes(),
                   boundaryNodes.size());
            return false;
        };
        // Check 2c: Check the shape as part of the boundary
        if (!checkSingleShape(boundingShape, boundaryNodes, i, itag<d - 1>{})) {
            logger(ERROR,
                   "There was an unspecified consistency error in the "
                   "boundaries of %-dimensional bounding shape %",
                   d, i);
            return false;
        }
    }
    return checkBoundaryShape(itag<d - 1>{});
}

template <std::size_t dimension>
bool ElementShape<dimension>::checkBoundaryShape(itag<-1>) const {
    return true;
}

template <std::size_t dimension>
template <std::size_t entityDimension, std::size_t targetDimension>
void ElementShape<dimension>::completeSubShapes(tag<entityDimension>,
                                                tag<targetDimension>) {
    // In this function we deal with 3 different shapes:
    // 1. This elementShape: E
    // 2. A boundary shape S for which we need to complete the topological
    //    information (entityDimension)
    // 3. A 'targetDimension' boundary shape S' to which it (could) be connected

    shapeData.template getAdjacentShapes<entityDimension>()[targetDimension]
        .resize(getNumberOfEntities<entityDimension>());
    // Loop over all S to update the dimension
    for (std::size_t entityIndex = 0;
         entityIndex < getNumberOfEntities<entityDimension>(); ++entityIndex) {
        if (targetDimension < entityDimension) {
            // Case 1: The target entities S' is on the boundary of S.
            //
            // An entity S' is on the boundary of S if and only if the vertices
            // of S' are a subset of the vertices of S. The ordering of the
            // vertices is then irrelevant.

            std::vector<std::size_t> result;  // i.e. indices of S'
            std::vector<std::size_t> adjacentNodes =
                getAdjacentEntities<entityDimension, 0>(entityIndex);
            std::sort(adjacentNodes.begin(), adjacentNodes.end());
            // Vector with for each S' the nodes that are used
            const auto& candidates =
                shapeData.template getAdjacentShapes<targetDimension>()[0];
            std::vector<std::size_t> candidate;
            for (std::size_t i = 0; i < candidates.size(); ++i) {
                candidate = candidates[i];
                std::sort(candidate.begin(), candidate.end());
                if (std::includes(adjacentNodes.begin(), adjacentNodes.end(),
                                  candidate.begin(), candidate.end())) {
                    // The vertices of S' are a subset of those of S
                    result.push_back(i);
                }
            }
            shapeData.template getAdjacentShapes<
                entityDimension>()[targetDimension][entityIndex] = result;
        }
        if (targetDimension == entityDimension) {
            // Case 2: S and S' have the same dimension, hence S' == S
            shapeData.template getAdjacentShapes<
                entityDimension>()[targetDimension][entityIndex] = {
                entityIndex};
        }
        if (targetDimension > entityDimension) {
            // Case 3: S is a boundary part of S', the reverse of case 1.
            // Now using that S is part of the boundary of S' if the vertices of
            // S are a subset of the vertices of S'
            std::vector<std::size_t> result;
            std::vector<std::size_t> adjacentNodes =
                getAdjacentEntities<entityDimension, 0>(entityIndex);
            std::sort(adjacentNodes.begin(), adjacentNodes.end());
            const auto& candidates =
                shapeData.template getAdjacentShapes<targetDimension>()[0];
            std::vector<std::size_t> candidate;
            for (std::size_t i = 0; i < candidates.size(); ++i) {
                candidate = candidates[i];
                std::sort(candidate.begin(), candidate.end());
                // this time the candidate is the bigger entity
                if (std::includes(candidate.begin(), candidate.end(),
                                  adjacentNodes.begin(), adjacentNodes.end())) {
                    result.push_back(i);
                }
            }
            shapeData.template getAdjacentShapes<
                entityDimension>()[targetDimension][entityIndex] = result;
        }
    }
    completeSubShapes(tag<entityDimension>{}, tag<targetDimension - 1>{});
}

template <std::size_t dimension>
template <std::size_t entityDimension>
void ElementShape<dimension>::completeSubShapes(tag<entityDimension>, tag<0>) {
    // The connection to 0-dimensional shapes (points/vertices) is used as
    // input. So there should be nothing to do here for this entityDimension.
    completeSubShapes(tag<entityDimension - 1>{}, tag<dimension - 1>{});
}

}  // namespace Preprocessor
