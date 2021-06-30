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

namespace hpgem {

/**
 * Mesh for the unit line segment. Each mesh is a complete refinement of the
 * previous mesh starting at a single element.
 * @return The file names to the meshes
 */
std::vector<std::string> getUnitSegmentMeshes() {
    std::string prefix = getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/";
    return {prefix + "unitLineN1.hpgem",  prefix + "unitLineN2.hpgem",
            prefix + "unitLineN4.hpgem",  prefix + "unitLineN8.hpgem",
            prefix + "unitLineN16.hpgem", prefix + "unitLineN32.hpgem"};
}

/**
 * Meshes for the unit square using triangles. The meshes are subsequent
 * refinements of a mesh with a pair of triangles.
 * @return The file names to the meshes
 */
std::vector<std::string> getUnitSquareTriangleMeshes() {
    std::string prefix = getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/";
    return {prefix + "unitSquareTrianglesN1.hpgem",
            prefix + "unitSquareTrianglesN2.hpgem",
            prefix + "unitSquareTrianglesN4.hpgem",
            prefix + "unitSquareTrianglesN8.hpgem",
            prefix + "unitSquareTrianglesN16.hpgem",
            prefix + "unitSquareTrianglesN32.hpgem",
            prefix + "unitSquareTrianglesN64.hpgem"};
}

/**
 * Meshes for the unit cube using cubes. The meshes are subsequent refinements
 * of a mesh with a single cube.
 * @return The file names to the meshes
 */
std::vector<std::string> getUnitCubeCubeMeshes() {
    std::string prefix = getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/";
    return {
        prefix + "unitCubeN1.hpgem", prefix + "unitCubeN2.hpgem",
        prefix + "unitCubeN4.hpgem", prefix + "unitCubeN8.hpgem",
        prefix + "unitCubeN16.hpgem",
    };
}

/**
 * Meshes for the unit cube using tetrahedra. The meshes are subsequent refinements
 * of a mesh with a single cube.
 * @return The file names to the meshes
 */
std::vector<std::string> getUnitCubeTetMeshes() {
    std::string prefix = getCMAKE_hpGEM_SOURCE_DIR() + "/tests/files/";
    return {
        prefix + "unitCubeN1Tet.hpgem", prefix + "unitCubeN2Tet.hpgem",
        prefix + "unitCubeN4Tet.hpgem", prefix + "unitCubeN8Tet.hpgem",
    };
}

}  // namespace hpgem

#endif  // HPGEM_TESTMESHES_H
