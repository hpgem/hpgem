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
#ifndef HPGEM_TESTMESHES_H
#define HPGEM_TESTMESHES_H

#include "string"
#include "vector"
#include "hpgem-cmake.h"
#include "Logger.h"
#include "ParallelRunTest.h"
#include <limits>

namespace hpgem {

namespace {

const std::size_t ALL_ENTRIES = std::numeric_limits<std::size_t>::max();

/**
 * Helper function to extra a subrange from a vector
 * @param source The source vector
 * @param minLevel The index of the first entry
 * @param maxLevel The index on past the last entry.
 * @return The subrange [minLevel, maxLevel)
 */
std::vector<std::string> limit(const std::vector<std::string>& source,
                               std::size_t minLevel, std::size_t maxLevel) {
    logger.assert_always(
        minLevel < source.size(),
        "Minimum % needs to be smaller than the number of meshes %", minLevel,
        source.size());
    logger.assert_always(
        minLevel <= maxLevel,
        "Maximum % needs to be at least larger than the minimum %", maxLevel,
        minLevel);
    if (maxLevel == ALL_ENTRIES) {
        maxLevel = source.size();
    } else if (maxLevel > source.size()) {
        logger(WARN, "No mesh at level %", maxLevel);
        maxLevel = source.size();
    }
    return std::vector<std::string>(source.begin() + minLevel,
                                    source.begin() + maxLevel);
}
}  // namespace

/**
 * Mesh for the unit line segment. Each mesh is a complete refinement of the
 * previous mesh starting at a single element.
 * @return The file names to the meshes
 */
std::vector<std::string> getUnitSegmentMeshes(
    std::size_t minLevel = 0, std::size_t maxLevel = ALL_ENTRIES) {
    std::string prefix = getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/";
    return limit({prefix + "unitLineN1.hpgem", prefix + "unitLineN2.hpgem",
                  prefix + "unitLineN4.hpgem", prefix + "unitLineN8.hpgem",
                  prefix + "unitLineN16.hpgem", prefix + "unitLineN32.hpgem"},
                 minLevel, maxLevel);
}

/**
 * Same as getUnitSegmentPeriodicMeshes(), with periodic boundary.
 */
std::vector<std::string> getUnitSegmentPeriodicMeshes(
    std::size_t minLevel = 0, std::size_t maxLevel = ALL_ENTRIES) {
    std::string prefix = getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/";
    return limit(
        {prefix + "unitLineN1Per.hpgem", prefix + "unitLineN2Per.hpgem",
         prefix + "unitLineN4Per.hpgem", prefix + "unitLineN8Per.hpgem",
         prefix + "unitLineN16Per.hpgem", prefix + "unitLineN32Per.hpgem"},
        minLevel, maxLevel);
}


/**
 * Meshes for the unit square using triangles. The meshes are subsequent
 * refinements of a mesh with a pair of triangles.
 * @return The file names to the meshes
 */
std::vector<std::string> getSingleProcessorUnitSquareTriangleMeshes(
    std::size_t minLevel = 0, std::size_t maxLevel = ALL_ENTRIES) {
    std::string prefix = getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/";
    return limit({prefix + "unitSquareTrianglesN1.hpgem",
                  prefix + "unitSquareTrianglesN2.hpgem",
                  prefix + "unitSquareTrianglesN4.hpgem",
                  prefix + "unitSquareTrianglesN8.hpgem",
                  prefix + "unitSquareTrianglesN16.hpgem",
                  prefix + "unitSquareTrianglesN32.hpgem",
                  prefix + "unitSquareTrianglesN64.hpgem"},
                 minLevel, maxLevel);
}

/**
 * Meshes for the unit square using triangles. The meshes are subsequent
 * refinements of a mesh with a pair of triangles.
 * @return The file names to the meshes
 */
std::vector<std::string> getDualProcessorUnitSquareTriangleMeshes(
    std::size_t minLevel = 0, std::size_t maxLevel = ALL_ENTRIES) {
    std::string prefix = getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/parallel_meshes/";
    return limit({prefix + "unitSquareTrianglesN1.hpgem",
                  prefix + "unitSquareTrianglesN2.hpgem",
                  prefix + "unitSquareTrianglesN4.hpgem",
                  prefix + "unitSquareTrianglesN8.hpgem",
                  prefix + "unitSquareTrianglesN16.hpgem",
                  prefix + "unitSquareTrianglesN32.hpgem",
                  prefix + "unitSquareTrianglesN64.hpgem"},
                 minLevel, maxLevel);
}

/**
 * Meshes for the unit square using triangles. The meshes are subsequent
 * refinements of a mesh with a pair of triangles.
 * @return The file names to the meshes
 */
std::vector<std::string> getUnitSquareTriangleMeshes(
    std::size_t minLevel = 0, std::size_t maxLevel = ALL_ENTRIES) {
    std::vector<std::string> meshFilenames;

    if (hpgem::isParallelRun()) {
        meshFilenames = getDualProcessorUnitSquareTriangleMeshes();
    } else {
        meshFilenames = getSingleProcessorUnitSquareTriangleMeshes();
    }

    return meshFilenames;
}

    
/**
 * Same as getUnitSquareTriangleMeshes, but with periodic boundaries. As result
 * the smallest mesh starts with 4x4 squares divided into triangles.
 */
std::vector<std::string> getUnitSquarePeriodicTriangleMeshes(
    std::size_t minLevel = 0, std::size_t maxLevel = ALL_ENTRIES) {
    std::string prefix = getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/";
    return limit({prefix + "unitSquareTrianglesN4Per.hpgem",
                  prefix + "unitSquareTrianglesN8Per.hpgem",
                  prefix + "unitSquareTrianglesN16Per.hpgem",
                  prefix + "unitSquareTrianglesN32Per.hpgem",
                  prefix + "unitSquareTrianglesN64Per.hpgem"},
                 minLevel, maxLevel);
}

/**
 * Meshes for the unit cube using cubes. The meshes are subsequent refinements
 * of a mesh with a single cube.
 * @return The file names to the meshes
 */
std::vector<std::string> getUnitCubeCubeMeshes(
    std::size_t minLevel = 0, std::size_t maxLevel = ALL_ENTRIES) {
    std::string prefix = getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/";
    return limit(
        {
            prefix + "unitCubeN1.hpgem",
            prefix + "unitCubeN2.hpgem",
            prefix + "unitCubeN4.hpgem",
            prefix + "unitCubeN8.hpgem",
            prefix + "unitCubeN16.hpgem",
        },
        minLevel, maxLevel);
}

/**
 * Same as getUnitCubeCubeMeshes but with periodic boundaries. Note only level 2
 * and larger are supported.
 */
std::vector<std::string> getUnitCubePeriodicCubeMeshes(
    std::size_t minLevel = 2, std::size_t maxLevel = ALL_ENTRIES) {
    logger.assert_always(
        minLevel >= 2 && maxLevel >= 2,
        "Periodic meshes with 1 or 2 elements are not supported");
    std::string prefix = getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/";
    // Compensate for missing first two meshes
    minLevel -= 2;
    if (maxLevel != ALL_ENTRIES) {
        maxLevel -= 2;
    }
    return limit(
        {
            // These two invalid meshes are commented out because hpgem does not
            // support these meshes due to their topology (at the time of
            // writing) . As example of the strange topology, for the N1 mesh
            // all corners are periodically connected, hence it would have:
            //  - 1 node with 8 coordinates
            //  - 3 edges (one for each x,y,z directions) with the same node on
            //    both sides.
            //  - 3 faces (normal to x,y and z direction), with 4 times the same
            //    node the same edge at top/bottom and left/right and the same
            //    element on both sides.
            //  - 1 Element
            //      prefix + "unitCubeN1Per.hpgem",
            //      prefix + "unitCubeN2Per.hpgem",
            prefix + "unitCubeN4Per.hpgem",
            prefix + "unitCubeN8Per.hpgem",
            prefix + "unitCubeN16Per.hpgem",
        },
        minLevel, maxLevel);
}

/**
 * Meshes for the unit cube using tetrahedra. The meshes are subsequent
 * refinements of a mesh with a single cube.
 * @return The file names to the meshes
 */
std::vector<std::string> getUnitCubeTetMeshes(
    std::size_t minLevel = 0, std::size_t maxLevel = ALL_ENTRIES) {
    std::string prefix = getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/";
    return limit(
        {prefix + "unitCubeN1Tet.hpgem", prefix + "unitCubeN2Tet.hpgem",
         prefix + "unitCubeN4Tet.hpgem", prefix + "unitCubeN8Tet.hpgem"},
        minLevel, maxLevel);
}

/**
 * Meshes for the unit circle using quadratic elements. The meshes are
 * subsequent uniform refinements of the first mesh.
 * @param minLevel
 * @param maxLevel
 * @return
 */
std::vector<std::string> getUnitCircleQuadraticTriangleMeshes(
    std::size_t minLevel = 0, std::size_t maxLevel = ALL_ENTRIES) {
    std::string prefix = getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/";
    return limit(
        {
            prefix + "circle-domain-quadratic-N1.hpgem",
            prefix + "circle-domain-quadratic-N2.hpgem",
            prefix + "circle-domain-quadratic-N3.hpgem",
            prefix + "circle-domain-quadratic-N4.hpgem",
        },
        minLevel, maxLevel);
}

/**
 * Meshes for the unit sphere using quadratic elements. The meshes are
 * subsequent uniform refinements of the first mesh. These meshes are relatively
 * coarse to keep the file size of the meshes down.
 * @param minLevel
 * @param maxLevel
 * @return
 */
std::vector<std::string> getUnitSphereQuadraticTriangleMeshes(
    std::size_t minLevel = 0, std::size_t maxLevel = ALL_ENTRIES) {
    std::string prefix = getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/";
    return limit(
        {
            prefix + "sphere-quadratic-N1.hpgem",
            prefix + "sphere-quadratic-N2.hpgem",
            prefix + "sphere-quadratic-N3.hpgem",
            prefix + "sphere-quadratic-N4.hpgem",
        },
        minLevel, maxLevel);
}

}  // namespace hpgem

#endif  // HPGEM_TESTMESHES_H
