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
#include "MeshFileInformation.h"

#include "Logger.h"
#include <fstream>
#include <limits>

namespace hpgem {
namespace Base {

void MeshFileInformation::readInformation(std::istream &stream) {
    {
        // File starts with a line 'mesh [version]'
        std::string rawInput;
        std::getline(stream, rawInput);
        logger.assert_always(rawInput.substr(0, 4) == "mesh",
                             "Not a proper hpgem mesh file header");
        std::istringstream temp(rawInput.substr(4));
        temp >> version;
        logger.assert_always(version >= 1 && version <= 2,
                             "Version % is not supported by this reader",
                             version);
    }

    // Header line with the dimension and number of parts
    {
        std::size_t nodes, elements;
        stream >> nodes >> elements >> dimension;
        entityCount.resize(dimension + 1);

        entityCount[0] = nodes;
        entityCount[dimension] = elements;
        // Read codimension
        for (std::size_t i = dimension - 1; i > 0; --i) {
            stream >> entityCount[i];
        }
    }
    // Partition information
    {
        std::size_t partitionCount;
        stream >> partitionCount;
        partitionNodeCounts.resize(partitionCount);
        for (std::size_t i = 0; i < partitionCount; ++i) {
            stream >> partitionNodeCounts[i];
        }
    }
    // Gobble up the remaining part of the line
    stream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // Zones
    if (version == 1) {
        // No zone information was available in version 1.
        zoneNames.resize(1);
        zoneNames[0] = "Main";
    } else {
        std::string line;
        std::getline(stream, line);
        logger.assert_always(line == "zones", "Incorrect zone header");

        std::size_t zoneCount;
        stream >> zoneCount;
        stream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        zoneNames.resize(zoneCount);
        for (std::size_t i = 0; i < zoneCount; ++i) {
            std::getline(stream, zoneNames[i]);
        }
    }
}

MeshFileInformation MeshFileInformation::readInformation(
    const std::string &fileName) {
    std::ifstream input;
    input.open(fileName.c_str());
    logger.assert_always(input.is_open(), "File % could not be opened",
                         fileName);
    logger.assert_always(input.good(), "Opening % resulted in an IO error",
                         fileName);

    MeshFileInformation result;
    result.readInformation(input);

    input.close();

    return result;
}

}  // namespace Base
}  // namespace hpgem