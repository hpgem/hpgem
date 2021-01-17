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
#include "FortranUnformattedFile.h"

namespace Preprocessor {

using namespace hpgem;

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
 *    This is the most common format. For some things (zones, periodic boundary
 *    conditions) the data does not follow this format. For example, Zones in 3D
 *    interleave the zone-element-indices and the zone names (instead of having
 *    two separate arrays). While periodic boundary conditions don't use arrays
 *    to group the storage of the properties, instead storing them seperately
 *    one in each dataline.
 *
 *    Note: My guess is that this multiline storage is needed as the start/end
 *    markers are stored as (u?)int32. Hence the number of bytes in each data
 *    line is limited.
 *
 * 5. This reader supports both 2D and 3D hyb meshes, and unlike the examples
 *    from centaur we don't separate them into completely different code blocks.
 *    This is possible because the basic structure of the format is the same,
 *    though there several tiny differences.
 */

template <std::size_t N>
std::string fromFortranString(const std::array<char, N>& fortranString) {
    // Find the end of the string
    std::size_t length = N;
    while (length > 0 && fortranString[length - 1] == ' ') {
        length--;
    }
    return std::string(fortranString.data(), length);
}

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

CentaurReader::CentaurReader(std::string filename) : centaurFile2(filename) {
    temp(filename);

    centaurFile.open(filename, std::ios::binary);
    logger.assert_always(centaurFile.is_open(),
                         "Cannot open Centaur meshfile.");
    logger.assert_always(centaurFile.good(),
                         "Something is not so good about this mesh");

    readHeader();

    // Don't accept hybtype 1 as this lacks hexahedra & pyramids.
    logger.assert_always(centaurFileType != 1,
                         "File with hyb-type 1 is too old.");

    logger.assert_always(centaurFileType != 6,
                         "File hyb-type 6 uses 64 bit indices for nodes and "
                         "elements, which is not implemented.");
    logger.assert_always(centaurFileType < 6, "File type is too new.");
    if (centaurFileType < 4 && centaurFileType > 0) {
        logger(WARN,
               "Old centaur file format. These files can be upgraded using "
               "hybconvert.");
    }

    /// Scan through the file ///
    /////////////////////////////
    // Scan through the file to find the positions where the interesting data
    // starts.

    readCoordinates();

    // Nodes
    nodeStart = centaurFile.tellg();
    // 1 group with all coordinates
    numberOfNodes = skipGroup();
    // Elements, these are stored by shape
    elementStart = centaurFile.tellg();

    // For each element type there is a group with the 3-8 node indices of the
    // corners.
    std::array<std::uint32_t, 4> elementCounts({0, 0, 0, 0});

    numberOfElements = 0;
    std::size_t numberOfElementTypes = is3D() ? 4 : 2;
    // For 3D the order is hexahedra, prisms, pyramids tetrahedra
    // For 2D the order is triangles, quadrilaterals

    for (std::size_t i = 0; i < numberOfElementTypes; ++i) {
        readElements(i, numNodes(i));
        elementCounts[i] += skipGroup();
        numberOfElements += elementCounts[i];
    }

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

    // Mapping of the boundary
    // 3D: Panel -> PanelGroup mapping
    // 2D: Segment -> boundary group mapping
    skipGroup(1, false);

    // Panel Group -> B.C. type & Name
    // skipGroup(2, false);  // group to b.c. type connections
    readBoundaryGroups();
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
    if (is3D()) {
        skipGroup();
    }

    // Zones
    // 1 data line with the number N of zones
    // 1? data line with N entries of the form:
    //  - 5 endOffsets, 4 for the elements, 1 for boundary faces indicating the
    //    start of the elements/boundary faces of this zone (2D: only 2
    //    endOffsets for the elements)
    //  - char80 zone name
    // (only version >4.99 && 3D)
    //  - 1 data line with N zone family names
    //
    // Note 2D: Version prior <= 4.99 allow not having zones it seems
    readZoneInfo(elementCounts);
    logger(INFO, "Finished scanning through hyb-file");

    // More internal information is possible, for example polyhedrons.
}

void CentaurReader::readHeader() {
    // We do not apriori know the size of the header
    std::uint32_t headerSize = centaurFile2.peekRecordSize();
    std::vector<char> rawInput(headerSize);

    centaurFile2.readRawRecord(headerSize, rawInput.data());
    std::uint32_t offset = 0;
    version = *reinterpret_cast<typeof(version)*>(rawInput.data());
    offset += sizeof(version);
    centaurFileType =
        *reinterpret_cast<typeof(centaurFileType)*>(rawInput.data() + offset);
    // Rest of the data seems to be char[80] with something like
    // "Unstructured hybrid grid format"
    // But can be ignored

    // TODO: Remove
    auto headerLine = readLine();
    headerLine >> version >> centaurFileType;

    logger(INFO, "This mesh is in Centaur version % format", version);
}

void CentaurReader::readCoordinates() {
    readGroup(coordinates, getDimension(), true);
}

void CentaurReader::readElements(std::size_t elemType, std::uint32_t numNodes) {
    logger.assert_always(elemType < 4, "Element type out of range");
    elementCount[elemType] = readGroup(elements[elemType], numNodes, true);
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

Range<MeshSource::Node> CentaurReader::getNodeCoordinates() {
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
                         MeshSource::Node& next) mutable {
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
        next.coordinates.resize(1);
        for (auto& value : coordinate) {
            currentLine >> value;
        }
        next.coordinates[0] = coordinate;
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
                next.coordinates.push_back(coordinate);
            }
            currentLine = std::move(oldLine);
            centaurFile.seekg(currentFilePosition);
            remainderThisLine = oldRemainder;
        }
        if (HPGEM_LOGLEVEL >= Log::VERBOSE) {
            for (auto outputCoordinate : next.coordinates) {
                std::cout << "(" << outputCoordinate[0];
                for (std::size_t i = 1; i < outputCoordinate.size(); ++i) {
                    std::cout << ", " << outputCoordinate[i];
                }
                std::cout << ")" << std::endl;
            }
        }
        logger(VERBOSE, "");
    };
    MeshSource::Node node;
    increment(node);
    return Range<MeshSource::Node>{node, std::move(increment), numberOfNodes};
}

Range<MeshSource::Element> CentaurReader::getElements() {
    centaurFile.seekg(elementStart);
    std::size_t currentGroupRemainder = 0;
    std::size_t currentGroupElements = 0;
    std::size_t currentZoneIndex = 0;
    std::size_t entitiesOnLine = 0;
    std::size_t remainderThisLine = 0;
    std::size_t groupsProcessed = 0;
    UnstructuredInputStream<std::istringstream> currentLine;
    auto increment = [=, currentLine = std::move(currentLine)](
                         MeshSource::Element& next) mutable {
        while (currentGroupRemainder == 0) {
            groupsProcessed++;

            currentGroupElements = 0;
            currentZoneIndex = 0;

            // Read the line with the element count
            currentLine = readLine();
            currentLine >> currentGroupRemainder;
            if (centaurFileType <= 4) {
                // No multiline in versions < 4 (includes 2D)
                entitiesOnLine = currentGroupRemainder;
            } else {
                // Multiline support for hybtype 5/6
                currentLine >> entitiesOnLine;
            }
            remainderThisLine = entitiesOnLine;
            // Resize the node storage
            std::size_t numberOfNodesPerElement = 0;
            if (is2D()) {
                if (groupsProcessed == 1) {
                    // Triangles
                    numberOfNodesPerElement = 3;
                } else if (groupsProcessed == 2) {
                    // Quadrilaterals
                    numberOfNodesPerElement = 4;
                } else {
                    logger.assert_always(false, "Too many groups for 2D");
                }
            } else {
                switch (groupsProcessed) {
                    case 1:
                        // Hexahedra
                        numberOfNodesPerElement = 8;
                        break;
                    case 2:
                        // Prisms
                        numberOfNodesPerElement = 6;
                        break;
                    case 3:
                        // Pyramids
                        numberOfNodesPerElement = 5;
                        break;
                    case 4:
                        // Tetrahedra
                        numberOfNodesPerElement = 4;
                        break;
                    default:
                        logger.assert_always(false, "Too many groups for 3D.");
                        break;
                }
            }
            next.coordinateIds.resize(numberOfNodesPerElement);
            // Read the line with the actual elements
            currentLine = readLine();
        }
        // Increment zone if needed
        while (currentGroupElements >=
               zones[currentZoneIndex].endOffsets[groupsProcessed - 1]) {
            currentZoneIndex++;
            logger.assert_debug(currentZoneIndex < zones.size(),
                                "Zone index out of bounds");
            next.zoneName = zones[currentZoneIndex].getZoneName();
        }

        if (remainderThisLine == 0) {
            // End of a data line in multiline
            currentLine = readLine();
            remainderThisLine = entitiesOnLine;
        }
        // Read the element
        for (auto& index : next.coordinateIds) {
            uint32_t input;
            currentLine >> input;
            index = toHpgemNumbering[input];
        }
        remainderThisLine--;
        currentGroupRemainder--;
    };
    MeshSource::Element first;
    increment(first);
    return {first, std::move(increment), numberOfElements};
}

std::uint32_t CentaurReader::skipGroup(std::size_t linesPerEntity,
                                       bool multiline) {
    GroupSize groupSize = readGroupSize(multiline);
    if (groupSize.totalCount == 0) {
        centaurFile2.skipRecord(0);
    } else {
        for (std::size_t i = 0; i < linesPerEntity; ++i) {
            for (std::size_t read = 0; read < groupSize.totalCount;
                 read += groupSize.perLineCount) {
                // TODO: Remove the need for peeking
                std::uint32_t size = centaurFile2.peekRecordSize();
                centaurFile2.skipRecord(size);
            }
        }
    }
    return groupSize.totalCount;
}

void CentaurReader::readPeriodicNodeConnections() {
    std::uint32_t numberOfPeriodicTransformations = 1;
    if (centaurFileType > 3 || is2D()) {
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

void CentaurReader::readBoundaryGroups() {
    GroupSize groupSize = readGroupSize(false);
    std::uint32_t numberOfGroups = groupSize.totalCount;
    logger(INFO, "Discarding % boundary groups", numberOfGroups);

    std::vector<std::uint32_t> boundaryTypeIds(numberOfGroups);
    std::vector<std::array<char, 80>> rawBoundaryNames(numberOfGroups);

    centaurFile2.readRawRecord(numberOfGroups * sizeof(std::uint32_t),
                               reinterpret_cast<char*>(boundaryTypeIds.data()));
    centaurFile2.readRawRecord(
        numberOfGroups * 80 * sizeof(char),
        reinterpret_cast<char*>(rawBoundaryNames.data()));
    // Output for debugging
    for (std::uint32_t i = 0; i < numberOfGroups; ++i) {
        logger(VERBOSE, "Boundary group \"%\" of type %'", rawBoundaryNames[i],
               boundaryTypeIds[i]);
    }
}

void CentaurReader::readZoneInfo(std::array<std::uint32_t, 4> elementCount) {
    if (is2D() && version <= 4.99 && centaurFile.eof()) {
        // Zone info is only mandatory for 2D, version > 4.99
        zones.resize(1);
        std::fill(zones[0].rawZoneName.begin(), zones[0].rawZoneName.end(),
                  ' ');
        char defaultName[] = "default";
        for (std::size_t i = 0; i < 7; ++i) {
            zones[0].rawZoneName[i] = defaultName[i];
        }
        zones[0].endOffsets[0] = elementCount[0];
        zones[0].endOffsets[1] = elementCount[1];
        return;
    }
    GroupSize groupSize = readGroupSize(false);
    std::uint32_t numberOfZones = groupSize.totalCount;
    zones.resize(numberOfZones);
    logger(INFO, "Reading % zones", numberOfZones);

    // For 2D only the offsets of 2 element types are stored
    // For 3D there are 4 element types and additionally the starting boundary
    // face is stored.
    std::size_t numIndexEntries = is3D() ? 5 : 2;
    // For 2D the zone names are stored separately
    std::size_t byteSize = numIndexEntries * sizeof(std::uint32_t) +
                           (is3D() ? 80 * sizeof(char) : 0);
    byteSize *= numberOfZones;
    std::vector<char> rawZoneData(byteSize);
    centaurFile2.readRawRecord(byteSize, rawZoneData.data());

    // Read all the zones
    std::size_t offset = 0;
    for (std::uint32_t i = 0; i < numberOfZones; ++i) {
        ZoneInformation& zone = zones[i];
        // The centaur file stores the element indices of the first element in
        // the zone. However, we are interested in the last index of each
        // element type in the zone.
        if (i == 0) {
            std::uint32_t entry;
            for (std::size_t j = 0; j < numIndexEntries; ++j) {
                entry = *reinterpret_cast<typeof(entry)*>(rawZoneData.data() +
                                                          offset);
                offset += sizeof(entry);
                logger.assert_always(
                    entry == 1, "Offset for the first zone entry is non zero.");
            }
        } else {
            ZoneInformation& prevZone = zones[i - 1];

            for (std::size_t entryId = 0; entryId < numIndexEntries;
                 ++entryId) {
                prevZone.endOffsets[entryId] =
                    *reinterpret_cast<typeof(std::uint32_t)*>(
                        rawZoneData.data() + offset);
                offset += sizeof(std::uint32_t);
                // Compensate for fortran 1-indexing
                prevZone.endOffsets[entryId]--;
            }
            if (is2D()) {
                // Filling to zero
                prevZone.endOffsets[2] = 0;
                prevZone.endOffsets[3] = 0;
                prevZone.endOffsets[4] = 0;
            }
        }
        if (is3D()) {
            // For 3D the zone name is directly appended to the zone
            auto start = rawZoneData.begin() + offset;
            std::copy(start, start + 80, zone.rawZoneName.begin());
            offset += 80 * sizeof(char);
        }
    }
    // Check whether we exactly read all data.
    logger.assert_always(offset == byteSize, "Incorrect zone reading");

    {
        // Fill the endOffsets of the last zone as there is no next zone to fill
        // them.
        ZoneInformation& lastZone = zones[numberOfZones - 1];
        for (std::size_t i = 0; i < 4; ++i) {
            lastZone.endOffsets[i] = elementCount[i];
        }
        // Not used
        lastZone.endOffsets[4] = 0;
    }

    if (is2D()) {
        // For 2D the zone names are separate
        std::vector<char> names(numberOfZones * 80);
        centaurFile2.readRawRecord(numberOfZones * 80 * sizeof(char),
                                   names.data());
        for (std::size_t i = 0; i < numberOfZones; ++i) {
            auto start = names.begin() + 80 * i;
            std::copy(start, start + 80, zones[i].rawZoneName.begin());
        }
    }

    if (is3D() && version > 4.99) {
        // Read zone families
        std::vector<char> familyNames(numberOfZones * 80);
        centaurFile2.readRawRecord(numberOfZones * 80 * sizeof(char),
                                   familyNames.data());
        for (std::size_t i = 0; i < numberOfZones; ++i) {
            auto start = familyNames.begin() + 80 * i;
            std::copy(start, start + 80, zones[i].rawZoneFamiliyName.begin());
        }
    } else {
        // No zone families: fill the zones with padding
        for (std::uint32_t i = 0; i < numberOfZones; ++i) {
            std::fill(zones[i].rawZoneFamiliyName.begin(),
                      zones[i].rawZoneFamiliyName.end(), ' ');
        }
    }
    for (const ZoneInformation& zone : zones) {
        logger(VERBOSE, "Zone \"%\"", zone.getZoneName());
        logger(VERBOSE, "\tFamily: \"%\"", zone.getZoneFamilyName());
        for (std::size_t i = 0; i < 5; ++i) {
            logger(DEBUG, "\tOffset % : %", i, zone.endOffsets[i]);
        }
    }
}

std::string CentaurReader::ZoneInformation::getZoneName() const {
    return fromFortranString(rawZoneName);
}

std::string CentaurReader::ZoneInformation::getZoneFamilyName() const {
    return fromFortranString(rawZoneFamiliyName);
}

CentaurReader::GroupSize CentaurReader::readGroupSize(bool expectMultiline) {
    // Multiline is only supported in v5 & v6
    bool multiline = expectMultiline && centaurFileType > 4;
    std::array<std::uint32_t, 2> data{};
    centaurFile2.readRawRecord(sizeof(std::uint32_t) * (1 + multiline),
                               reinterpret_cast<char*>(data.data()));
    GroupSize result{};
    result.totalCount = data[0];
    if (multiline) {
        result.perLineCount = data[1];
    } else {
        result.perLineCount = result.totalCount;
    }
    return result;
}

template <typename T>
std::uint32_t CentaurReader::readGroup(std::vector<T>& data,
                                       std::uint32_t numComponents,
                                       bool expectMultiline) {
    GroupSize groupSize = readGroupSize(expectMultiline);
    // Prevent overflow
    std::size_t totalRawEntries =
        static_cast<std::size_t>(groupSize.totalCount) * numComponents;
    data.resize(totalRawEntries);
    if (totalRawEntries == 0) {
        // Read the empty group
        centaurFile2.readRawRecord(0, reinterpret_cast<char*>(data.data()));
    } else {
        std::uint32_t entriesRead = 0;
        while (entriesRead < groupSize.totalCount) {
            // Compute the number of entries on this line
            std::uint32_t lineSize = std::min(
                groupSize.perLineCount, groupSize.totalCount - entriesRead);
            // Read the actual entries
            std::size_t offset = entriesRead * numComponents;
            std::size_t lineByteSize = lineSize * numComponents * sizeof(T);
            centaurFile2.readRawRecord(
                lineByteSize, reinterpret_cast<char*>(data.data() + offset));
            entriesRead += lineSize;
        }
        logger.assert_always(entriesRead == groupSize.totalCount,
                             "Incorrect number of entries read");
    }
    return groupSize.totalCount;
}

void CentaurReader::temp(const std::string& path) {
    FortranUnformattedFile file(path);
    while (!file.eof()) {
        std::uint32_t size = file.peekRecordSize();
        file.skipRecord(size);
        logger(INFO, "Skipped record of size %", size);
    }
    logger(INFO, "Reached end of file");
}

std::uint32_t CentaurReader::numNodes(std::size_t elementType) const {
    if (is2D()) {
        logger.assert_always(elementType < 2, "Invalid element type % for 2D",
                             elementType);
        return elementType == 0 ? 3 : 4;
    } else {
        switch (elementType) {
            case 0:
                return 8;  // Hexahedra;
            case 1:
                return 6;  // Prisms
            case 2:
                return 5;  // Pyramids
            case 3:
                return 4;  // Tetrahedra
            default:
                logger.assert_always(false, "Invalid element  % type for 3D",
                                     elementType);
                return -1;
        }
    }
}

}  // namespace Preprocessor
