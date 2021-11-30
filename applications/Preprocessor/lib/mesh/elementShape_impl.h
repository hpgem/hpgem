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
std::enable_if_t<(entityDimension <= dimension), std::size_t>
    ElementShape<dimension>::getNumberOfEntities() const {
    return shapeParts.template get<entityDimension>().size();
}

template <std::size_t dimension>
template <std::size_t entityDimension>
auto ElementShape<dimension>::getBoundaryShape(std::size_t entityIndex) const
    -> std::enable_if_t<(entityDimension <= dimension),
                        const ShapeType<entityDimension>*> {
    logger.assert_debug(entityIndex < getNumberOfEntities<entityDimension>(),
                        "This shape is bounded by only % shapes of dimension "
                        "%, but you asked for shape %",
                        getNumberOfEntities<entityDimension>(), entityDimension,
                        entityIndex);
    return shapeParts.template get<entityDimension>()[entityIndex]
        .getReferenceGeometry();
}

template <std::size_t dimension>
template <std::size_t entityDimension, std::size_t targetDimension>
std::enable_if_t<(entityDimension <= dimension && targetDimension <= dimension),
                 std::vector<std::size_t>>
    ElementShape<dimension>::getAdjacentEntities(
        std::size_t entityIndex) const {
    logger.assert_debug(entityIndex < getNumberOfEntities<entityDimension>(),
                        "This shape is bounded by only % shapes of dimension "
                        "%, but you asked for shape %",
                        getNumberOfEntities<entityDimension>(), dimension,
                        entityIndex);
    return shapeParts.template get<entityDimension>()[entityIndex]
        .getAdjacentShapes(targetDimension);
}

template <std::size_t dimension>
template <int d, std::size_t shapeDimension>
bool ElementShape<dimension>::checkSingleShape(
    const ElementShape<shapeDimension>* boundingShape,
    std::vector<std::size_t> boundaryNodes, std::size_t currentIndex,
    itag<d> dimTag) const {
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
        for (std::size_t j = 0; j < shapeParts[dimTag].size(); ++j) {
            auto comparisonNodes = shapeParts[dimTag][j].getNodes();
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
bool ElementShape<dimension>::checkBoundaryShape(itag<d> dimTag) const {
    // Check the individual boundary shapes
    for (std::size_t i = 0; i < getNumberOfEntities<d>(); ++i) {

        auto boundingShape = getBoundaryShape<d>(i);
        // Check 1: Does the shape itself make sense?
        // The case d == dimension is excluded as that implies
        // boundingShape == this, and we would get an infinite loop.
        if (!(d == dimension || boundingShape->checkShape())) {
            logger(ERROR,
                   "%-dimensional bounding shape % has an unspecified "
                   "consistency error",
                   d, i);
            return false;
        }
        // Check 2: The shape has a number of nodes, check that the topological
        // information has enough nodes.
        auto boundaryNodes = shapeParts[dimTag][i].getNodes();
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
void ElementShape<dimension>::completeSubShapes(tag<entityDimension> edimTag,
                                                tag<targetDimension> tdimTag) {
    // In this function we deal with 3 different shapes:
    // 1. This elementShape: E
    // 2. A boundary shape S for which we need to complete the topological
    //    information (entityDimension)
    // 3. A 'targetDimension' boundary shape S' to which it (could) be connected
    //
    // For each entity S it should compute to which S' (of appropriate
    // dimension) it is connected and use that to fill the adjacency
    // information.

    // Loop over all S to update the dimension
    for (std::size_t entityIndex = 0;
         entityIndex < getNumberOfEntities<entityDimension>(); ++entityIndex) {
        if (targetDimension == entityDimension) {
            // Case 1: S and S' have the same dimension, hence S' == S
            shapeParts[edimTag][entityIndex].getAdjacentShapes(
                targetDimension) = {entityDimension};
        } else {
            // Mixed case where either S' is on the boundary of S or vice versa.
            //
            // An entity S' is on the boundary of S if and only if the vertices
            // of S' are a subset of the vertices of S. The ordering of the
            // vertices is then irrelevant.
            //
            // The case S on the boundary of S' follows analogously

            std::vector<std::size_t> result;  // i.e. indices of S'
            std::vector<std::size_t> adjacentNodes =
                getAdjacentEntities<entityDimension, 0>(entityIndex);
            std::sort(adjacentNodes.begin(), adjacentNodes.end());
            // Vector with for each S' the nodes that are used
            const auto& candidates = shapeParts[tdimTag];
            for (std::size_t i = 0; i < candidates.size(); ++i) {
                const Detail::ElementShapePart<targetDimension, dimension>&
                    candidate = candidates[i];
                std::vector<std::size_t> candidateNodes = candidate.getNodes();
                std::sort(candidateNodes.begin(), candidateNodes.end());

                if (targetDimension < entityDimension) {
                    if (std::includes(
                            adjacentNodes.begin(), adjacentNodes.end(),
                            candidateNodes.begin(), candidateNodes.end())) {
                        // The vertices of S' are a subset of those of S
                        result.push_back(i);
                    }
                } else {
                    if (std::includes(
                            candidateNodes.begin(), candidateNodes.end(),
                            adjacentNodes.begin(), adjacentNodes.end())) {
                        // The vertices of S are a subset of those of S'
                        result.push_back(i);
                    }
                }
            }
            // Assign the found connections
            shapeParts[edimTag][entityIndex].getAdjacentShapes(
                targetDimension) = result;
        }
    }
    completeSubShapes(tag<entityDimension>{}, tag<targetDimension - 1>{});
}

template <std::size_t dimension>
template <std::size_t entityDimension>
void ElementShape<dimension>::completeSubShapes(tag<entityDimension>, tag<0>) {
    // The connection to 0-dimensional shapes (points/vertices) is used as
    // input. So there should be nothing to do here for this entityDimension.
    completeSubShapes(tag<entityDimension - 1>{}, tag<dimension>{});
}

}  // namespace Preprocessor
