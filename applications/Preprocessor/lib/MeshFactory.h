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
#ifndef HPGEM_MESHFACTORY_H
#define HPGEM_MESHFACTORY_H

#include "mesh/Mesh.h"
#include "MergePlan.h"
#include "MeshSource.h"

namespace Preprocessor {

// for use with file readers that can 'guess' the correct numbering of the node
// coordinates file readers that don't do this should overload this function
template <std::size_t dimension>
Mesh<dimension> readFile(MeshSource& file) {
    Mesh<dimension> result;
    logger.assert_always(dimension == file.getDimension(),
                         "Mismatching dimensions");
    for (auto nodeCoordinates : file.getNodeCoordinates()) {
        result.addNode();
        for (auto coordinate : nodeCoordinates.coordinates) {
            logger.assert_debug(
                coordinate.size() == dimension,
                "The coordinates read by this reader have the wrong dimension");
            result.addNodeCoordinate(EntityGId(result.getNumberOfNodes() - 1),
                                     coordinate.data());
        }
    }
    std::vector<CoordId> coords;
    for (auto element : file.getElements()) {
        logger.assert_always(!element.zoneName.empty(),
                             "Element without a zone name");
        // TODO: Move up into MeshSource at a convenient moment
        coords.resize(element.coordinateIds.size());
        for (std::size_t i = 0; i < coords.size(); ++i) {
            coords[i] = CoordId(element.coordinateIds[i]);
        }
        result.addElement(coords, element.zoneName);
    }
    logger.assert_debug(result.isValid(), "Unspecified problem with the mesh");
    return result;
}

template <std::size_t dimension>
Mesh<dimension> fromMeshSource(MeshSource2& file) {
    Mesh<dimension> result;
    logger.assert_always(dimension == file.getDimension(),
                         "Mismatching dimensions");
    // Mapping of the nodeId, source nodeId -> mesh nodeId
    std::map<std::size_t, EntityGId> nodeMapping;
    // Add all the nodes
    for (auto coord : file.getCoordinates()) {
        // Find the node in the mesh, or create a new one if needed
        std::size_t sourceNodeId = coord.nodeId;
        EntityGId meshNodeId;
        auto current = nodeMapping.find(sourceNodeId);
        if (current != nodeMapping.end()) {
            meshNodeId = current->second;
        } else {
            meshNodeId = result.addNode();
            nodeMapping[sourceNodeId] = meshNodeId;
        }
        // Add the actual coordinate
        logger.assert_debug(
            coord.coordinate.size() == dimension,
            "The coordinates read by this reader have the wrong dimension");
        result.addNodeCoordinate(
            meshNodeId,
            LinearAlgebra::SmallVector<dimension>(coord.coordinate.data()));
    }
    // Add all the elements
    std::vector<CoordId> coords;

    for (auto element : file.getElements()) {
        logger.assert_always(!element.zoneName.empty(),
                             "Element without a zone name");
        // TODO: Move up into MeshSource at a convenient moment
        coords.resize(element.coordinateIds.size());
        for (std::size_t i = 0; i < coords.size(); ++i) {
            coords[i] = CoordId(element.coordinateIds[i]);
        }
        result.addElement(coords, element.zoneName);
    }
    logger.assert_debug(result.isValid(), "Unspecified problem with the mesh");

    for (const std::map<std::size_t, std::size_t> rawCoordPairing :
         file.getMerges()) {
        logger(INFO, "Applying coordinate merger with % pairings",
               rawCoordPairing.size());
        // 1-1 translation of the indices
        std::map<CoordId, CoordId> coordPairing;
        for (const auto& rawPair : rawCoordPairing) {
            CoordId first = CoordId(rawPair.first);
            CoordId second = CoordId(rawPair.second);
            coordPairing[first] = second;
            // Print a list of coordinates to be merged for debugging
            logger(DEBUG, "Merging %-% (% -- %)", first, second,
                   result.getCoordinate(first), result.getCoordinate(second));
        }

        MergePlan<dimension> plan =
            MergePlan<dimension>::computeMergePlan(&result, coordPairing);
        plan.executeMerge();
    }
    if (!file.getMerges().empty()) {
        result.removeUnusedEntities();
    }

    return result;
}

}  // namespace Preprocessor

#endif  // HPGEM_MESHFACTORY_H
