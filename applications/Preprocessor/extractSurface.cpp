/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2019, University of Twente
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

// this executable only links against a limited part of the hpGEM kernel, don't
// expect all functionality to be available
#ifdef HPGEM_USE_METIS
#include <metis.h>
#else
#include <cstddef>
using idx_t = std::size_t;
#endif
#include <chrono>
#include <string>
#include "Base/CommandLineOptions.h"
#include "Base/MpiContainer.h"
#include "mesh.h"
#include "hpgem.h"
#include "centaur.h"
#include "meshData.h"
#include "output.h"

using namespace hpgem;

// Describes an axis aligned boundary plane boundary
class Boundary {
   public:
    Boundary() = default;
    Boundary(const char* charDescription) {
        std::string description(charDescription);
        std::transform(description.begin(), description.end(),
                       description.begin(), tolower);
        auto separator = description.find('=');
        logger.assert_always(separator != std::string::npos,
                             "boundary description must contain =");
        std::string coordinate = description.substr(0, separator);
        std::string value = description.substr(separator + 1);
        logger.assert_always(value.size() > 0, "Equations don't end with '='");
        std::size_t processed;
        location = std::stod(value, &processed);
        logger.assert_always(processed == value.size(),
                             "Did not recognize % as a floating point value",
                             value);
        if (coordinate == "x" || coordinate == "x0") {
            coordinateIndex = 0;
        } else if (coordinate == "y" || coordinate == "x1") {
            coordinateIndex = 1;
        } else if (coordinate == "z" || coordinate == "x2") {
            coordinateIndex = 2;
        } else {
            logger(ERROR, "Did not recognize % as a coordinate", description);
        }
    }
    Boundary(const Boundary&) = default;
    Boundary(Boundary&&) = default;

    Boundary& operator=(const Boundary&) = default;
    Boundary& operator=(Boundary&&) = default;

    ///  Check whether a coordinate is on the boundary
    template <std::size_t dimension>
    bool onBoundary(LinearAlgebra::SmallVector<dimension> coordinate) {
        logger.assert_always(coordinateIndex < dimension,
                             "user input assumes a mesh of dimension at least "
                             "%, but the coordinate dimension is %",
                             coordinateIndex + 1, dimension);
        return std::abs(coordinate[coordinateIndex] - location) < 1e-12;
    }

    /// Project a coordinate onto the plane
    template <std::size_t dimension>
    LinearAlgebra::SmallVector<dimension - 1> boundaryCoordinate(
        LinearAlgebra::SmallVector<dimension> coordinate) {
        LinearAlgebra::SmallVector<dimension - 1> result;
        for (std::size_t i = 0, j = 0; i < dimension; ++i) {
            if (i == coordinateIndex) continue;
            result[j++] = coordinate[i];
        }
        return result;
    }

   private:
    std::size_t coordinateIndex;
    double location;
};

auto& inputFileName = Base::register_argument<std::string>(
    '\0', "inFile", "Name of your input file", true);
auto& boundaryData = Base::register_argument<Boundary>(
    '\0', "position", "Location of the boundary to be extracted", true);

template <std::size_t dimension>
void processMesh(Preprocessor::HpgemReader& meshFile) {
    using Preprocessor::CoordId;
    using Preprocessor::EntityGId;

    auto inputMesh = Preprocessor::readFile<dimension>(meshFile);
    Preprocessor::MeshData<std::size_t, dimension, dimension> inputPartitionID(
        &inputMesh);
    auto processorBindings = meshFile.getProcessorBindings();
    {
        auto copyStart = inputPartitionID.data();
        for (auto start = processorBindings.begin();
             start != processorBindings.end(); ++start) {
            *copyStart++ = *start;
        }
    }
    auto position = boundaryData.getValue();
    Preprocessor::Mesh<dimension - 1> resultMesh;
    Preprocessor::MeshData<std::size_t, dimension - 1, dimension - 1>
        resultPartitionID(&resultMesh);
    // Ids of the nodes in the original mesh that are on the boundary
    std::vector<EntityGId> boundaryNodes;
    // Mapping from coordinateId in the old mesh to that in the new mesh. Only
    // has coordinateIds of coordinates on the boundary.
    std::map<CoordId, CoordId> toBoundaryIndex;
    // Coordinates that have been processed, independent of whether they are on
    // the boundary or not.
    std::set<CoordId> processedCoordinates;
    for (auto& node : inputMesh.getNodes()) {
        bool added = false;
        EntityGId newNodeId;  // Id of the node in the processed mesh
        for (std::size_t i = 0; i < node.getNumberOfElements(); ++i) {
            Preprocessor::Element<dimension>& element = node.getElement(i);
            CoordId nodeCoordId =
                element.getCoordinateIndex(node.getLocalIndex(i));
            auto nodeCoordinate = element.getCoordinate(node.getLocalIndex(i));
            if (position.onBoundary(nodeCoordinate) &&
                processedCoordinates.find(nodeCoordId) ==
                    processedCoordinates.end()) {
                if (!added) {
                    // Only add the node to the new mesh if it has at least one
                    // coordinate on the boundary.
                    newNodeId = resultMesh.addNode();
                    boundaryNodes.push_back(node.getGlobalIndex());
                    added = true;
                }
                // Project the coordinate and add it to the new mesh and update
                // the mapping old -> new.
                CoordId newCoordId = resultMesh.addNodeCoordinate(
                    newNodeId, position.boundaryCoordinate(nodeCoordinate));
                toBoundaryIndex[nodeCoordId] = newCoordId;
                logger(DEBUG, "% -> %", nodeCoordId, newCoordId);
            }
            processedCoordinates.insert(nodeCoordId);
        }
    }
    // The faces from the input mesh that are on the boundary will form the
    // elements of the target mesh.
    for (auto& face : inputMesh.getFaces()) {
        // Due to the input mesh restriction we know that a face is on the
        // boundary if all its nodes are on the boundary.
        auto nodeIndices = face.template getIncidenceListAsIndices<0>();
        bool onBoundary = true;
        for (auto index : nodeIndices) {
            if (std::find(boundaryNodes.begin(), boundaryNodes.end(), index) ==
                boundaryNodes.end()) {
                onBoundary = false;
            }
        }
        if (onBoundary) {
            // To get coordinates for the face we need to go through a adjacent
            // element.
            auto inputElement = face.getElement(0);
            auto localNodeIndices =
                inputElement.template getLocalIncidenceListAsIndices<0>(face);
            // Coordinate ids in the original mesh
            std::vector<CoordId> globalCoordinateIndices;
            for (auto localIndex : localNodeIndices) {
                globalCoordinateIndices.push_back(
                    inputElement.getCoordinateIndex(localIndex));
            }
            // Overwrite the coordinates with coordinates from the new mesh
            for (auto& inputIndex : globalCoordinateIndices) {
                inputIndex = toBoundaryIndex[inputIndex];
            }
            resultMesh.addElement(globalCoordinateIndices);
            auto resultElement = resultMesh.getElements().back();
            // Assign the partition of the face based on the element it belonged
            // to.
            resultPartitionID[resultElement] = inputPartitionID[inputElement];
        }
    }
    Preprocessor::outputMesh(resultMesh, resultPartitionID,
                             meshFile.getTargetProcessorCount());
}

/**
 * @brief Tool to generate a surface mesh for one of the sides of generated by
 * the structured mesh generator.
 */
int main(int argc, char** argv) {
    auto start = std::chrono::steady_clock::now();
    Base::parse_options(argc, argv);
    if (Base::MPIContainer::Instance().getNumberOfProcessors() > 1) {
        logger(WARN,
               "This is a sequential code, the target number of processors in "
               "deduced from the input file instead");
    }
    Base::MPIContainer::Instance().onlyOnOneProcessor({[]() {
        std::string fileName = inputFileName.getValue();
        std::transform(fileName.begin(), fileName.end(), fileName.begin(),
                       tolower);
        auto hpgemFile = Preprocessor::HpgemReader(inputFileName.getValue());
        if (hpgemFile.getDimension() == 2) {
            processMesh<2>(hpgemFile);
        } else if (hpgemFile.getDimension() == 3) {
            processMesh<3>(hpgemFile);
        } else {
            logger(ERROR, "Dimension % is not supported",
                   hpgemFile.getDimension());
        }
    }});
    auto end = std::chrono::steady_clock::now();
    logger(INFO, "processing took %ms",
           std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
               .count());
    return 0;
}
