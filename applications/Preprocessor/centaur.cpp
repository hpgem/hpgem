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

#include <algorithm>
#include <iostream>

#include "centaur.h"
#include "Logger.h"

using namespace hpgem;

namespace Preprocessor {

/*
 * General notes.
 * 1. The data is stored in 'Fortran unformatted style'. This means that the
 *   data is stored in data lines, each of the shape:
 *   [start marker] data [end marker]
 *   where both start and end marker are the number of bytes contained in the
 *   data it encloses.
 *
 * 2. As a Fortran dataformat all indices start at 1 (not at 0).
 *
 * 3. Confusingly the hyb files have two 'version' type values. One is a
 *    floating point version that seems to refer to the capabilities (e.g. zone
 *    families are only supported in version >= 5). The second is a fileType
 *    which seem to have 4 categories:
 *     < 0 is used for 2D
 *     4 is 3D without multiline
 *     5 is 3D with multiline
 *     6 is 3D with multiline and int64 indices (for node counts etc.)
 *
 * 4. The data is stored linearly in blocks/groups. For example, after a header
 *    all the node coordinates are stored, followed by all the hexahedra, then
 *    all prisms, etc. Each block/group has a fixed structure:
 *      - 1 data line with the number of entities Ne (e.g. node coordinates) and
 *          for hyb_type 5 & 6 the number of entries per data line
 *      - Nl data lines for the first property
 *      - Nl data lines for the second property
 *      - etc.
 *    The number Nl differs between hyb_type 4 and types 5 & 6. For version 4 Nl
 *    is fixed at 1. For version 5 and 6 it is variable, hence the 'multiline'
 *    capability. The number of lines is such that the first Nl-1 lines contain
 *    exactly Ne entries, and the last line contains the remaining entries.
 *
 *    Note: My guess is that this multiline storage is needed as the start/end
 *    markers are stored as (u?)int32. Hence the number of bytes in each data
 *    line is limited.
 *
 * 5. This reader supports both 2D and 3D hyb meshes, and unlike the examples
 *    from centaur we don't separate them into completely different code blocks.
 */

UnstructuredInputStream<std::istringstream> CentaurReader::readLine() {
    // todo figure out ifstream.get(streambuf)
    logger.assert_always(!!centaurFile,
                         "There is expected to be an extra line, but the file "
                         "is not ready for reading");
    std::uint32_t lineSize, referenceSize;
    centaurFile >> lineSize;
    std::string buffer(lineSize, '\0');
    // workaround: string::data() returns const char* before c++17
    centaurFile.read(&buffer.front(), lineSize);
    centaurFile >> referenceSize;
    if (HPGEM_LOGLEVEL >= Log::DEBUG) {
        // swap the comment if you would like to see doubles instead
        // std::cout << buffer;
        char toHex[] = "0123456789abcdef";
        for (unsigned char input : buffer) {
            std::cout << toHex[input / 16] << toHex[input % 16] << " ";
        }
        // doubleBuffer buffers doubles, it has nothing to do with double
        // buffering note that this might break if your hardware is picky about
        // the alignment of its data
        /*const double* doubleBuffer = reinterpret_cast<const
        double*>(buffer.data()); for(std::size_t i = 0; i < buffer.size() /
        sizeof(double); ++i) { std::cout << doubleBuffer[i] << " ";
        }*/
        std::cout << std::endl;
    }
    logger.assert_always(lineSize == referenceSize,
                         "read error in centaur file");
    return std::istringstream(buffer, std::ios_base::binary);
}

void CentaurReader::readHeader() {
    auto headerLine = readLine();
    headerLine >> version >> centaurFileType;

    logger(INFO, "This mesh is in Centaur version % format", version);
}

CentaurReader::CentaurReader(std::string filename) {
    centaurFile.open(filename, std::ios::binary);
    logger.assert_always(centaurFile.is_open(),
                         "Cannot open Centaur meshfile.");
    logger.assert_always(centaurFile.good(),
                         "Something is not so good about this mesh");

    readHeader();

    /// Scan through the file ///
    /////////////////////////////
    // Scan through the file to find the positions where the interesting data
    // starts.

    // Nodes
    nodeStart = centaurFile.tellg();
    // 1 group with all coordinates
    numberOfNodes = skipGroup();

    // Elements, these are stored by shape
    elementStart = centaurFile.tellg();

    // For each element type there is a group with the 3-8 node indices of the
    // corners.
    numberOfElements = 0;
    if (is3D()) {
        numberOfElements += skipGroup();  // hexahedra
    }
    numberOfElements += skipGroup();  // triangles/prisms
    if (is3D()) {
        numberOfElements += skipGroup();  // pyramids
    }
    numberOfElements += skipGroup();  // quadrilaterals/tetrahedra

    // Additional information

    logger(DEBUG, "done with skipping elements");
    if (is2D()) {
        // In 2D there is the option to specify the boundary condition per node
        // or per edge (=face equivalent).

        // 1 group with pairs of indices (node, group)
        // meaning that the specified node belongs to the boundary group. Nodes
        // can belong to multiple groups
        skipGroup();  // redundant boundary nodes
    }

    // (interface) Boundary faces
    // 3D: Face nodes & Panel id
    // 2D: Boundary edges (2 edge ids) & segment id
    if (is2D() || centaurFileType > 3) {
        skipGroup(2);
    } else {
        // For hybtype < 4 panel id was stored in the face-node table
        skipGroup(1);  // boundary faces
    }

    // Note the following seems to be incorrect:
    // the segment to boundary group mapping should be equivalent to the 3D
    // panel ->PanelGroup mapping. Thus I don't understand what this extra line
    // is doing here.
    if (is2D()) {
        // Segment -> boundary group mapping
        skipGroup();  // segments
    }
    // Panel -> PanelGroup mapping
    skipGroup(1, false);
    // Panel Group -> B.C. type & Name
    skipGroup(2, false);  // group to b.c. type connections
    logger(VERBOSE, "done with skipping information");

    // Periodicity information
    // (<3 only a single periodicity, hence 1 array of periodic node pairs)
    // For hybtype >=4:
    // - 1 data line with the number N of periodic connections
    // - N sets of entries
    //   - 1 data line with 4x4 matrix (xform1) 2D: 3x3
    //   - 1 data line with 4x4 matrix (xform2) 2D: 3x3
    //   - 1 data line with number Nn number of periodic node pairs
    //   - 1 data line with Nn pairs (nod1, nod2), node nod1 is connected to
    //       node nod2. Matrices correspond to the transformation from nod1 to
    //       nod2 and vice versa.
    readPeriodicNodeConnections();

    // Interface panels (3D only)
    // 1 data line with number N of interface panels (could be zero)
    // 1 data line with N entries (panel numbers?)

    // Zones
    // 1 data line with the number N of zones
    // 1? data line with N entries of the form:
    //  - 5 offsets, 4 for the elements, 1 for boundary faces indicating the
    //    start of the elements/boundary faces of this zone (2D: only 2 offsets
    //    for the elements)
    //  - char80 zone name
    // Note 2D no >4.99 seems to allow not having zone information

    // (only version >4.99)
    // 1? data line with N zone family names

    // More internal information is possible, for example polyhedrons.
}

std::function<void(int&)> difference_generator(
    const std::vector<std::size_t>& input) {
    std::size_t i = 1;
    return [=](int& result) mutable {
        result = static_cast<long long>(input[i + 1]) -
                 static_cast<long long>(input[i]);
        i++;
    };
}

Range<std::vector<std::vector<double>>> CentaurReader::getNodeCoordinates() {
    centaurFile.seekg(nodeStart);
    std::size_t dimension = 3;
    if (centaurFileType < 0) dimension = 2;
    std::vector<double> coordinate(dimension);
    std::uint32_t nodesOnLine;
    auto currentLine = readLine();
    currentLine >> nodesOnLine;  // actually available nodes, but for old
                                 // centaur files that is the same
    if (centaurFileType > 4) currentLine >> nodesOnLine;
    logger(VERBOSE, "Processing % nodes per line", nodesOnLine);
    std::uint32_t remainderThisLine = nodesOnLine;
    std::size_t index = 1;
    currentLine = readLine();
    auto increment = [=, currentLine = std::move(currentLine)](
                         std::vector<std::vector<double>>& next) mutable {
        if (remainderThisLine == 0) {
            remainderThisLine = nodesOnLine;
            currentLine = readLine();
        }
        auto position = boundaryConnections.find(index);
        // skip over boundary nodes that are already processed
        while (position != boundaryConnections.end() &&
               position->second[0] != index) {
            position = boundaryConnections.find(++index);
            for (auto& value : coordinate) {
                currentLine >> value;
            }
            if (--remainderThisLine == 0) {
                remainderThisLine = nodesOnLine;
                currentLine = readLine();
            }
        }
        next.resize(1);
        for (auto& value : coordinate) {
            currentLine >> value;
        }
        next[0] = coordinate;
        remainderThisLine--;
        index++;
        if (position != boundaryConnections.end() &&
            position->second.size() > 1) {
            std::istringstream::pos_type currentFilePosition =
                centaurFile.tellg();
            auto oldRemainder = remainderThisLine;
            std::istringstream::pos_type currentLinePosition =
                currentLine.tellg();
            auto oldLine = std::move(currentLine);
            currentLine = std::istringstream(oldLine.str());
            currentLine.seekg(currentLinePosition);
            for (auto offset :
                 Range<int>{static_cast<int>(position->second[1]) -
                                static_cast<int>(position->second[0]),
                            difference_generator(position->second),
                            position->second.size() - 1}) {
                logger.assert_always(offset > 0, "Need to skip % positions",
                                     offset);
                while (static_cast<uint32_t>(offset) > remainderThisLine) {
                    currentLine = readLine();
                    offset -= remainderThisLine;
                    remainderThisLine = nodesOnLine;
                }
                while (offset-- > 0) {
                    --remainderThisLine;
                    for (auto& value : coordinate) {
                        currentLine >> value;
                    }
                }
                next.push_back(coordinate);
            }
            currentLine = std::move(oldLine);
            centaurFile.seekg(currentFilePosition);
            remainderThisLine = oldRemainder;
        }
        if (HPGEM_LOGLEVEL >= Log::VERBOSE) {
            for (auto outputCoordinate : next) {
                std::cout << "(" << outputCoordinate[0];
                for (std::size_t i = 1; i < outputCoordinate.size(); ++i) {
                    std::cout << ", " << outputCoordinate[i];
                }
                std::cout << ")" << std::endl;
            }
        }
        logger(VERBOSE, "");
    };
    std::vector<std::vector<double>> coordinates;
    increment(coordinates);
    return Range<std::vector<std::vector<double>>>{
        coordinates, std::move(increment), numberOfNodes};
}

Range<std::vector<std::size_t>> CentaurReader::getElements() {
    centaurFile.seekg(elementStart);
    std::size_t currentGroupRemainder = 0;
    std::size_t entitiesOnLine = 0;
    std::size_t remainderThisLine = 0;
    std::size_t groupsProcessed = 0;
    UnstructuredInputStream<std::istringstream> currentLine;
    auto increment = [=, currentLine = std::move(currentLine)](
                         std::vector<std::size_t>& next) mutable {
        while (currentGroupRemainder == 0) {
            if (centaurFileType < 2)
                groupsProcessed += 2;
            else
                groupsProcessed++;
            currentLine = readLine();
            currentLine >> currentGroupRemainder;
            entitiesOnLine = currentGroupRemainder;
            if (centaurFileType > 4) currentLine >> entitiesOnLine;
            remainderThisLine = entitiesOnLine;
            if (groupsProcessed == 1)
                next.resize(8);
            else if (groupsProcessed == 3)
                next.resize(5);
            else if (groupsProcessed == 4)
                next.resize(4);
            else if (centaurFileType < 0)
                next.resize(3);
            else
                next.resize(6);
            currentLine = readLine();
        }
        if (remainderThisLine == 0) {
            currentLine = readLine();
            remainderThisLine = entitiesOnLine;
        }
        for (auto& index : next) {
            uint32_t input;
            currentLine >> input;
            index = toHpgemNumbering[input];
        }
        remainderThisLine--;
        currentGroupRemainder--;
    };
    std::vector<std::size_t> first;
    increment(first);
    return {first, std::move(increment), numberOfElements};
}

std::uint32_t CentaurReader::skipGroup(std::size_t linesPerEntity,
                                       bool multiline) {
    std::uint32_t numberOfEntities;
    auto currentLine = readLine();
    currentLine >> numberOfEntities;
    logger(VERBOSE, "skipping % entities", numberOfEntities);
    std::uint32_t numberOfEntitiesPerLine = numberOfEntities;
    if (centaurFileType > 4 && multiline)
        currentLine >> numberOfEntitiesPerLine;
    logger(VERBOSE, "there are % entities per line", numberOfEntitiesPerLine);
    if (numberOfEntities == 0) readLine();
    for (std::size_t i = 0; i < numberOfEntities;
         i += numberOfEntitiesPerLine) {
        for (std::size_t j = 0; j < linesPerEntity; ++j) readLine();
    }
    return numberOfEntities;
}

void CentaurReader::readPeriodicNodeConnections() {
    std::uint32_t numberOfPeriodicTransformations = 1;
    if (centaurFileType > 3) {
        auto currentLine = readLine();
        currentLine >> numberOfPeriodicTransformations;
    }
    for (std::size_t i = 0; i < numberOfPeriodicTransformations; ++i) {
        if (centaurFileType > 3) {
            readLine();
            readLine();  // rotation/scaling information
        }
        std::uint32_t numberOfPeriodicNodes;
        auto currentLine = readLine();
        currentLine >> numberOfPeriodicNodes;
        currentLine = readLine();
        for (std::size_t j = 0; j < numberOfPeriodicNodes; ++j) {
            std::uint32_t matchingNodes[2];
            currentLine >> matchingNodes[0] >> matchingNodes[1];
            std::vector<std::size_t>& first =
                boundaryConnections[matchingNodes[0]];
            std::vector<std::size_t>& second =
                boundaryConnections[matchingNodes[1]];
            std::vector<std::size_t> temp;
            if (matchingNodes[0] > matchingNodes[1])
                std::swap(matchingNodes[0], matchingNodes[1]);
            std::set_union(first.begin(), first.end(), matchingNodes,
                           matchingNodes + 2, std::back_inserter(temp));
            first.clear();
            std::set_union(temp.begin(), temp.end(), second.begin(),
                           second.end(), std::back_inserter(first));
            second = first;
            if (HPGEM_LOGLEVEL >= Log::DEBUG) {
                for (auto k : first) {
                    std::cout << k << " ";
                }
                std::cout << "(" << matchingNodes[0] << ", " << matchingNodes[1]
                          << ")" << std::endl;
            }
        }
    }
    std::size_t hpgemIndex = 0;
    constexpr std::size_t iMax = std::numeric_limits<std::size_t>::max();
    toHpgemNumbering.resize(numberOfNodes + 1, iMax);
    numberOfNodes = 0;
    // the indexing in hpgem has the nodes regrouped such that periodically
    // connected coordinates are together indices in centaur start at 1
    for (std::size_t i = 1; i < toHpgemNumbering.size(); ++i) {
        if (toHpgemNumbering[i] == iMax) {
            numberOfNodes++;
            auto position = boundaryConnections.find(i);
            if (position == boundaryConnections.end()) {
                logger(DEBUG, "% -> %", i, hpgemIndex);
                toHpgemNumbering[i] = hpgemIndex++;
            } else {
                for (auto nodeIndex : position->second) {
                    logger(DEBUG, "% -> %", nodeIndex, hpgemIndex);
                    toHpgemNumbering[nodeIndex] = hpgemIndex++;
                }
            }
        }
    }
}
}  // namespace Preprocessor
