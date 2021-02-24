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
#include <random>
#include "Base/CommandLineOptions.h"
#include "Base/MpiContainer.h"
#include "mesh.h"
#include "hpgem.h"
#include "centaur.h"
#include "meshData.h"
#include "output.h"

using namespace hpgem;

auto& inputFileName = Base::register_argument<std::string>(
    '\0', "inFile", "Name of your input file", true);
auto& fileType = Base::register_argument<std::string>(
    '\0', "type",
    "Type of the file; only needed if this cant be deduced from the extention",
    false);
auto& dimension = Base::register_argument<std::size_t>(
    'd', "dimension", "deprecated - The dimension mentioned in the input file",
    false);
auto& targetMpiCount = Base::register_argument<std::size_t>(
    'n', "MPICount", "Target number of processors", false, 1);

template <std::size_t dimension>
void printMeshStatistics(const Preprocessor::Mesh<dimension>& mesh) {
    logger(INFO, "Mesh counts");
    logger(INFO, "\tElements: %", mesh.getNumberOfElements());
    logger(INFO, "\tFaces: %", mesh.getNumberOfFaces());
    if (dimension > 2) {
        logger(INFO, "\tEdges: %", mesh.getNumberOfEdges());
    }
    if (dimension > 1) {
        logger(INFO, "\tNodes: %", mesh.getNumberOfNodes());
    }

    if (!mesh.getNodeCoordinates().empty()) {
        LinearAlgebra::SmallVector<dimension> minCoord, maxCoord;
        minCoord = mesh.getNodeCoordinates()[0].coordinate;
        maxCoord = minCoord;

        for (const auto& node : mesh.getNodeCoordinates()) {
            for (std::size_t i = 0; i < dimension; ++i) {
                minCoord[i] = std::min(minCoord[i], node.coordinate[i]);
                maxCoord[i] = std::max(maxCoord[i], node.coordinate[i]);
            }
        }
        logger(INFO, "Mesh bounding box min % - max %", minCoord, maxCoord);
    }
}

template <std::size_t dimension>
Preprocessor::MeshData<idx_t, dimension, dimension> partitionMesh(
    Preprocessor::Mesh<dimension> mesh) {
    Preprocessor::MeshData<idx_t, dimension, dimension> partitionID(&mesh);
    idx_t numberOfProcessors = targetMpiCount.getValue();
    if (numberOfProcessors > 1) {
#ifdef HPGEM_USE_METIS
        idx_t numberOfConstraints = 1;
        idx_t numberOfElements = mesh.getNumberOfElements();
        float imbalance = 1.001;
        idx_t totalCutSize;

        std::vector<idx_t> xadj(numberOfElements + 1);
        std::vector<idx_t> adjncy(
            2 * mesh.getNumberOfFaces());  // actually interior faces only
        idx_t connectionsUsed{0};
        for (auto element : mesh.getElements()) {
            xadj[element.getGlobalIndex()] = connectionsUsed;
            for (auto face : element.getFacesList()) {
                if (face.getNumberOfElements() == 2) {
                    if (face.getElement(0) == element) {
                        adjncy[connectionsUsed++] =
                            face.getElement(1).getGlobalIndex();
                    } else {
                        adjncy[connectionsUsed++] =
                            face.getElement(0).getGlobalIndex();
                    }
                }
            }
        }
        xadj.back() = connectionsUsed;
        std::random_device seedGenerator;
        idx_t seed = seedGenerator();

        logger(DEBUG, "Metis ran with the seed %", seed);

        idx_t metisOptions[METIS_NOPTIONS];
        METIS_SetDefaultOptions(metisOptions);
        metisOptions[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM;
        metisOptions[METIS_OPTION_RTYPE] = METIS_RTYPE_FM;
        metisOptions[METIS_OPTION_SEED] = seed;
        METIS_PartGraphKway(&numberOfElements, &numberOfConstraints,
                            xadj.data(), adjncy.data(), NULL, NULL, NULL,
                            &numberOfProcessors, NULL, &imbalance, metisOptions,
                            &totalCutSize, partitionID.data());
#else
        logger(
            ERROR,
            "Please install metis if you want to generate a distributed mesh");
#endif
    }
    return partitionID;
}

template <std::size_t dimension>
void processMesh(Preprocessor::Mesh<dimension> mesh) {
    printMeshStatistics(mesh);

    Preprocessor::MeshData<idx_t, dimension, dimension> partitionID =
        partitionMesh(mesh);
    idx_t numberOfProcessors = targetMpiCount.getValue();
    Preprocessor::outputMesh(mesh, partitionID, numberOfProcessors);
}

/**
 * @brief preprocessing tool to generate distributed mesh files and/or to
 * convert supported file formats to the hpGEM native format
 *
 */
int main(int argc, char** argv) {
    auto start = std::chrono::steady_clock::now();
    Base::parse_options(argc, argv);
    if (Base::MPIContainer::Instance().getNumberOfProcessors() > 1) {
        logger(WARN,
               "This is a sequential code, use the command line options to set "
               "the target number of processors instead");
    }
    Base::MPIContainer::Instance().onlyOnOneProcessor({[]() {
        std::string fileName = inputFileName.getValue();
        std::transform(fileName.begin(), fileName.end(), fileName.begin(),
                       tolower);
        std::size_t nameSize = fileName.size();
        if ((fileType.isUsed() && fileType.getValue() == "hpgem") ||
            (!fileType.isUsed() &&
             fileName.compare(nameSize - 6, 6, ".hpgem") == 0)) {
            auto hpgemFile =
                Preprocessor::HpgemReader(inputFileName.getValue());
            if (dimension.isUsed()) {
                logger.assert_always(
                    dimension.getValue() == hpgemFile.getDimension(),
                    "The input file reports being for dimension %, but the "
                    "code was started for dimension %",
                    hpgemFile.getDimension(), dimension.getValue());
            }
            if (hpgemFile.getDimension() == 1) {
                processMesh(Preprocessor::readFile<1>(hpgemFile));
            } else if (hpgemFile.getDimension() == 2) {
                processMesh(Preprocessor::readFile<2>(hpgemFile));
            } else if (hpgemFile.getDimension() == 3) {
                processMesh(Preprocessor::readFile<3>(hpgemFile));
            } else {
                logger(ERROR, "Dimension % is not supported",
                       hpgemFile.getDimension());
            }
        } else if ((fileType.isUsed() && fileType.getValue() == "centaur") ||
                   (!fileType.isUsed() &&
                    fileName.compare(nameSize - 4, 4, ".hyb") == 0)) {
            auto centaurFile =
                Preprocessor::CentaurReader(inputFileName.getValue());
            if (dimension.isUsed()) {
                logger.assert_always(
                    dimension.getValue() == centaurFile.getDimension(),
                    "The input file reports being for dimension %, but the "
                    "code was started for dimension %",
                    centaurFile.getDimension(), dimension.getValue());
            }
            if (centaurFile.getDimension() == 2) {
                processMesh(Preprocessor::readFile<2>(centaurFile));
            } else if (centaurFile.getDimension() == 3) {
                processMesh(Preprocessor::readFile<3>(centaurFile));
            } else {
                logger(ERROR,
                       "Centaur file should not be able to have dimension %",
                       centaurFile.getDimension());
            }
        } else {
            logger(ERROR, "Don't know what to do with this file");
        }
    }});
    auto end = std::chrono::steady_clock::now();
    logger(INFO, "processing took %ms",
           std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
               .count());
    return 0;
}
