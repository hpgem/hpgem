/*
 This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
 http://hpgem.org/about-the-code/team
 
 This code is distributed using BSD 3-Clause License. A copy of which can found below.
 
 
 Copyright (c) 2017, University of Twente
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include<Logger.h>
#include<numeric>
#include"elementShape.h"

template<std::size_t dimension>
template<int entityDimension>
std::enable_if_t<(entityDimension < 0), std::size_t> Preprocessor::ElementShape<dimension>::getNumberOfEntities() const {
    if(entityDimension + dimension < 0) return 0;
    return getNumberOfEntities<entityDimension + dimension>();
}

template<std::size_t dimension>
template<int entityDimension>
std::enable_if_t<(entityDimension >= dimension), std::size_t> Preprocessor::ElementShape<dimension>::getNumberOfEntities() const {
    if(entityDimension > dimension) return 0;
    return 1;
}

template<std::size_t dimension>
template<int entityDimension>
std::enable_if_t<(entityDimension >= 0 && entityDimension < dimension), std::size_t> Preprocessor::ElementShape<dimension>::getNumberOfEntities() const {
    return shapeData.template getEntityShapes<entityDimension>().size();
}

template<std::size_t dimension>
template<int entityDimension>
auto Preprocessor::ElementShape<dimension>::getBoundaryShape(std::size_t entityIndex) const -> std::enable_if_t<(entityDimension < 0), const ShapeType<entityDimension>*> {
    logger.assert_debug(entityDimension + dimension >= 0, "This shape is not bounded by shapes of dimension %", entityDimension + dimension);
    return getBoundaryShape<entityDimension + dimension>(entityIndex);
}

template<std::size_t dimension>
template<int entityDimension>
auto Preprocessor::ElementShape<dimension>::getBoundaryShape(std::size_t entityIndex) const -> std::enable_if_t<(entityDimension >= dimension), const ShapeType<entityDimension>*> {
    logger.assert_debug(entityDimension == dimension, "This shape is not bounded by shapes of dimension %");
    logger.assert_debug(entityIndex == 0, "This shape is bounded by only 1 shapes of dimension %, but you asked for shape %", dimension, entityIndex);
    return *this;
}

template<std::size_t dimension>
template<int entityDimension>
auto Preprocessor::ElementShape<dimension>::getBoundaryShape(std::size_t entityIndex) const -> std::enable_if_t<(entityDimension >= 0 && entityDimension < dimension), const ShapeType<entityDimension>*> {
    logger.assert_debug(entityIndex < getNumberOfEntities<entityDimension>(),
                        "This shape is bounded by only % shapes of dimension %, but you asked for shape %",
                        getNumberOfEntities<entityDimension>(), entityDimension, entityIndex);
    return shapeData.template getEntityShapes<entityDimension>()[entityIndex];
}

template<std::size_t dimension>
template<int entityDimension, int targetDimension>
std::enable_if_t<(entityDimension < 0), Preprocessor::stackVector<std::size_t>> Preprocessor::ElementShape<dimension>::getAdjacentEntities(std::size_t entityIndex) const {
    logger.assert_debug(entityDimension + dimension >= 0, "This shape is not bounded by shapes of dimension %", entityDimension + dimension);
    return getAdjacentEntities<entityDimension + dimension, targetDimension>(entityIndex);
}

template<std::size_t dimension>
template<int entityDimension, int targetDimension>
std::enable_if_t<(targetDimension < 0 && entityDimension >= 0), Preprocessor::stackVector<std::size_t>> Preprocessor::ElementShape<dimension>::getAdjacentEntities(std::size_t entityIndex) const {
    if(targetDimension + dimension < 0) return {};
    return getAdjacentEntities<entityDimension, targetDimension + dimension>(entityIndex);
}

template<std::size_t dimension>
template<int entityDimension, int targetDimension>
std::enable_if_t<(entityDimension >= 0 && targetDimension >= 0 && (entityDimension >= dimension || targetDimension >= dimension)), Preprocessor::stackVector<std::size_t>> Preprocessor::ElementShape<dimension>::getAdjacentEntities(std::size_t entityIndex) const {
    logger.assert_debug(entityDimension <= dimension, "This shape is not bounded by shapes of dimension %", entityDimension + dimension);
    logger
    .assert_debug(entityIndex < getNumberOfEntities<entityDimension>(), "This shape is bounded by only % shapes of dimension %, but you asked for shape %",
                  getNumberOfEntities<entityDimension>(), dimension, entityIndex);
    if(targetDimension > dimension) return {};
    if(targetDimension == dimension) return {0};
    std::size_t amount = getNumberOfEntities<targetDimension>();
    stackVector<std::size_t> result(amount);
    std::iota(result.begin(), result.end(), 0);
    return result;
}

template<std::size_t dimension>
template<int entityDimension, int targetDimension>
std::enable_if_t<(entityDimension >= 0 && entityDimension < dimension && targetDimension >= 0 && targetDimension < dimension), Preprocessor::stackVector<std::size_t>>
Preprocessor::ElementShape<dimension>::getAdjacentEntities(std::size_t entityIndex) const {
    logger
    .assert_debug(entityIndex < getNumberOfEntities<entityDimension>(), "This shape is bounded by only % shapes of dimension %, but you asked for shape %",
                  getNumberOfEntities<entityDimension>(), dimension, entityIndex);
    return shapeData.template getAdjacentShapes<entityDimension>()[targetDimension][entityIndex];
}


template<std::size_t dimension>
template<std::size_t d, std::size_t shapeDimension>
bool Preprocessor::ElementShape<dimension>::checkSingleShape(const ElementShape<shapeDimension>* boundingShape, stackVector<std::size_t> boundaryNodes, std::size_t currentIndex, tag<d>) const {
    for(std::size_t i = 0; i < boundingShape->template getNumberOfEntities<d>(); ++i) {
        auto localNodes = boundingShape->template getAdjacentEntities<d, 0>(i);
        stackVector<std::size_t> cellNodes(localNodes.size());
        for(std::size_t j = 0; j < localNodes.size(); ++j) {
            cellNodes[j] = boundaryNodes[localNodes[j]];
        }
        std::sort(cellNodes.begin(), cellNodes.end());
        bool match = false;
        for(std::size_t j = 0; j < shapeData.template getAdjacentShapes<d>()[0].size(); ++j) {
            auto comparisonNodes = shapeData.template getAdjacentShapes<d>()[0][j];
            std::sort(comparisonNodes.begin(), comparisonNodes.end());
            if(cellNodes == comparisonNodes) {
                match = true;
                auto candindateIndices = getAdjacentEntities<d, shapeDimension>(j);
                if(std::find(candindateIndices.begin(), candindateIndices.end(), currentIndex) == candindateIndices.end()) {
                    logger(ERROR, "The %-dimensional shape % bounding %-dimensional shape % is not adjacent to it",
                           d, j, shapeDimension, currentIndex);
                    return false;
                }
            }
        }
        if(!match) {
            logger(ERROR, "The %-dimensional bounding shape % contains a bounding shape that is not a bounding shape of this entity", shapeDimension, currentIndex);
        }
    }
    return checkSingleShape(boundingShape, boundaryNodes, currentIndex, tag<d - 1>{});
}

template<std::size_t dimension>
template<std::size_t shapeDimension>
bool Preprocessor::ElementShape<dimension>::checkSingleShape(const ElementShape <shapeDimension>* boundingShape, stackVector <std::size_t> boundaryNodes, std::size_t currentIndex, tag <0>) const {
    for(std::size_t i = 0; i < boundingShape->template getNumberOfEntities<0>(); ++i) {
        auto localNodes = boundingShape->template getAdjacentEntities<0, 0>(i);
        stackVector<std::size_t> cellNodes(localNodes.size());
        for(std::size_t j = 0; j < localNodes.size(); ++j) {
            cellNodes[j] = boundaryNodes[localNodes[j]];
        }
        std::sort(cellNodes.begin(), cellNodes.end());
        bool match = false;
        for(std::size_t j = 0; j < shapeData.template getAdjacentShapes<0>()[0].size(); ++j) {
            auto comparisonNodes = shapeData.template getAdjacentShapes<0>()[0][j];
            std::sort(comparisonNodes.begin(), comparisonNodes.end());
            if(cellNodes == comparisonNodes) {
                match = true;
                auto candindateIndices = getAdjacentEntities<0, shapeDimension>(j);
                if(std::find(candindateIndices.begin(), candindateIndices.end(), currentIndex) == candindateIndices.end()) {
                    logger(ERROR, "The %-dimensional shape % bounding %-dimensional shape % is not adjacent to it",
                           0, j, shapeDimension, currentIndex);
                    return false;
                }
            }
        }
        if(!match) {
            logger(ERROR, "The %-dimensional bounding shape % contains a bounding shape that is not a bounding shape of this entity", shapeDimension, currentIndex);
        }
    }
    return true;
}

template<std::size_t dimension>
template<std::size_t d>
bool Preprocessor::ElementShape<dimension>::checkBoundaryShape(tag<d>) const {
    if(shapeData.template getEntityShapes<d>().size() != shapeData.template getAdjacentShapes<d>()[0].size()) {
        logger(ERROR, "There are % %-dimensional bounding shapes described by their bounding nodes, but % %-dimensional shapes described by their topological shape",
               shapeData.template getEntityShapes<d>().size(), d, shapeData.template getAdjacentShapes<d>()[0].size(), d);
        return false;
    }
    for(std::size_t i = 0; i < getNumberOfEntities<d>(); ++i) {
        auto boundingShape = getBoundaryShape<d>(i);
        if(!boundingShape->checkShape()) {
            logger(ERROR, "%-dimensional bounding shape % has an unspecified consistency error", d, i);
            return false;
        }
        auto boundaryNodes = shapeData.template getAdjacentShapes<d>()[0][i];
        if(boundaryNodes.size() != boundingShape->getNumberOfNodes()) {
            logger(ERROR, "%-dimensional shape % should be described be % nodes, but it has % nodes listed",
                   d, i, boundingShape->getNumberOfNodes(), boundaryNodes.size());
            return false;
        };
        if(!checkSingleShape(boundingShape, boundaryNodes, i, tag<d - 1>{})) {
            logger(ERROR, "There was an unspecified consistency error in the boundaries of %-dimensional bounding shape %", d, i);
            return false;
        }
    }
    return checkBoundaryShape(tag<d - 1>{});
}

template<std::size_t dimension>
bool Preprocessor::ElementShape<dimension>::checkBoundaryShape(tag<0>) const {
    if(shapeData.template getEntityShapes<0>().size() != shapeData.template getAdjacentShapes<0>()[0].size()) {
        logger(ERROR, "There are % %-dimensional bounding shapes described by their bounding nodes, but % %-dimensional shapes described by their topological shape",
               shapeData.template getEntityShapes<0>().size(), 0, shapeData.template getAdjacentShapes<0>()[0].size(), 0);
        return false;
    }
    for(std::size_t i = 0; i < getNumberOfEntities<0>(); ++i) {
        auto boundingShape = getBoundaryShape<0>(i);
        if(!boundingShape->checkShape()) {
            logger(ERROR, "%-dimensional bounding shape % has an unspecified consistency error", 0, i);
            return false;
        }
        auto boundaryNodes = shapeData.template getAdjacentShapes<0>()[0][i];
        if(boundaryNodes.size() != boundingShape->getNumberOfNodes()) {
            logger(ERROR, "%-dimensional shape % should be described be % nodes, but it has % nodes listed",
                   0, i, boundingShape->getNumberOfNodes(), boundaryNodes.size());
            return false;
        };
    }
    return true;
}


template<std::size_t dimension>
template<std::size_t entityDimension, std::size_t targetDimension>
void Preprocessor::ElementShape<dimension>::completeSubShapes(tag<entityDimension>, tag<targetDimension>) {
    shapeData.template getAdjacentShapes<entityDimension>()[targetDimension].resize(getNumberOfEntities<entityDimension>());
    for(std::size_t entityIndex = 0; entityIndex < getNumberOfEntities<entityDimension>(); ++entityIndex) {
        if(targetDimension < entityDimension) {
            stackVector <std::size_t> result;
            //if the shape is topologically consistent the ordering of the nodes does not matter
            //example: a square has 4 nodes, so it has 6 pairs of nodes; 4 of which correspond to edges and the other 2
            //are diagonals. Since the edges of a cube always are the edges of their bounding squares and never a diagonal
            //we just need to check if an edge of the cube has both bounding nodes as bounding nodes of the square
            //and we don't have to check if they are diagonally across from eachother
            stackVector <std::size_t> adjacentNodes = getAdjacentEntities<entityDimension, 0>(entityIndex);
            std::sort(adjacentNodes.begin(), adjacentNodes.end());
            const auto& candidates = shapeData.template getAdjacentShapes<targetDimension>()[0];
            stackVector <std::size_t> candidate;
            for(std::size_t i = 0; i < candidates.size(); ++i) {
                candidate = candidates[i];
                std::sort(candidate.begin(), candidate.end());
                if(std::includes(adjacentNodes.begin(), adjacentNodes.end(), candidate.begin(), candidate.end())) {
                    result.push_back(i);
                }
            }
            shapeData.template getAdjacentShapes<entityDimension>()[targetDimension][entityIndex] = result;
        }
        if(targetDimension == entityDimension) {
            shapeData.template getAdjacentShapes<entityDimension>()[targetDimension][entityIndex] = {entityIndex};
        }
        if(targetDimension > entityDimension) {
            stackVector <std::size_t> result;
            stackVector <std::size_t> adjacentNodes = getAdjacentEntities<entityDimension, 0>(entityIndex);
            std::sort(adjacentNodes.begin(), adjacentNodes.end());
            const auto& candidates = shapeData.template getAdjacentShapes<targetDimension>()[0];
            stackVector <std::size_t> candidate;
            for(std::size_t i = 0; i < candidates.size(); ++i) {
                candidate = candidates[i];
                std::sort(candidate.begin(), candidate.end());
                //this time the candidate is the bigger entity
                if(std::includes(candidate.begin(), candidate.end(), adjacentNodes.begin(), adjacentNodes.end())) {
                    result.push_back(i);
                }
            }
            shapeData.template getAdjacentShapes<entityDimension>()[targetDimension][entityIndex] = result;
        }
    }
    completeSubShapes(tag<entityDimension>{}, tag<targetDimension - 1>{});
}

template<std::size_t dimension>
template<std::size_t entityDimension>
void Preprocessor::ElementShape<dimension>::completeSubShapes(tag <entityDimension>, tag<0>) {
    completeSubShapes(tag<entityDimension - 1>{}, tag<dimension - 1>{});
}
