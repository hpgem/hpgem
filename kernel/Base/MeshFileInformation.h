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
#ifndef HPGEM_MESHFILEINFORMATION_H
#define HPGEM_MESHFILEINFORMATION_H

#include <istream>
#include <memory>
#include <vector>

namespace hpgem {
namespace Base {

/**
 * Reader for the header information in the HPGEM mesh format.
 *
 * This will read the meta data of the mesh.
 */
class MeshFileInformation {
   public:
    /**
     * Zone name used as placeholder for version 1 of the mesh format.
     */
    static const std::string MESH_V1_ZONENAME;

    MeshFileInformation() = default;

    /**
     * Read mesh file information from a stream. The stream is expected to be at
     * the start of the mesh file, and all information up to and including the
     * zone names is consumed.
     * @param stream The stream to read from
     */
    void readInformation(std::istream& stream);

    /**
     * Factory method to read the information from file.
     * @param fileName The file name to read from
     * @return The information about the mesh file
     */
    static MeshFileInformation readInformation(const std::string& fileName);

    /**
     * Version of the mesh format used
     */
    std::size_t version = 0;
    /**
     * The dimension of the mesh
     */
    std::size_t dimension = 0;
    /**
     * Count of the number of entities, starting from 0-dimensional (nodes)
     * to elements (d-dimensional)
     */
    std::vector<std::size_t> entityCount;
    /**
     * Number of nodes for each processor
     */
    std::vector<std::size_t> partitionNodeCounts;
    /**
     * Names of the zones. If this is a version 1 mesh it will contain
     * MESH_V1_ZONENAME.
     */
    std::vector<std::string> zoneNames;

    // Helper methods
    std::size_t getProcessorCount() const { return partitionNodeCounts.size(); }
};

}  // namespace Base
}  // namespace hpgem

#endif  // HPGEM_MESHFILEINFORMATION_H
