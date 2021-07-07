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
#include <cstring>

#include "centaur.h"
#include "Logger.h"
#include "FortranUnformattedFile.h"
#include "SparseUnionFind.h"

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
 *    to group the storage of the properties, instead storing them separately
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

CentaurReader::CentaurReader(std::string filename)
    : centaurFile(filename), elementCount({0, 0, 0, 0}) {
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

    std::size_t numberOfElementTypes = is3D() ? 4 : 2;
    // For 3D the order is hexahedra, prisms, pyramids tetrahedra
    // For 2D the order is triangles, quadrilaterals

    for (std::size_t i = 0; i < numberOfElementTypes; ++i) {
        readElements(i, numNodes(i));
    }

    // Additional information

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
    readZoneInfo();
    logger(INFO, "Finished scanning through hyb-file");

    // More internal information is possible, for example polyhedrons.
}

void CentaurReader::readHeader() {
    // We do not apriori know the size of the header
    std::uint32_t headerSize = centaurFile.peekRecordSize();
    std::vector<char> rawInput(headerSize);

    centaurFile.readRawRecord(headerSize, rawInput.data());

    std::uint32_t offset = 0;
    memcpy(&version, rawInput.data(), sizeof(offset));
    offset += sizeof(version);
    memcpy(&centaurFileType, rawInput.data() + offset, sizeof(centaurFileType));
    // Rest of the data seems to be char[80] with something like
    // "Unstructured hybrid grid format"
    // But can be ignored

    logger(INFO, "This mesh is in Centaur version % format", version);
}

void CentaurReader::readCoordinates() {
    std::vector<double> rawCoordinates;
    readGroup(rawCoordinates, getDimension(), true);
    const std::size_t dim = getDimension();
    std::size_t size = rawCoordinates.size() / dim;
    coordinates.resize(size);
    for (std::size_t i = 0; i < size; ++i) {
        MeshSource2::Coord& coord = coordinates[i];
        coord.nodeId = i;
        coord.coordinate.resize(dim);
        for (std::size_t j = 0; j < dim; ++j) {
            coord.coordinate[j] = rawCoordinates[i * dim + j];
        }
    }
    logger(INFO, "Read % coordinates", coordinates.size());
}

void CentaurReader::readElements(std::size_t elemType, std::uint32_t numNodes) {
    logger.assert_always(elemType < 4, "Element type out of range");
    // Read the continuous stream of coordinate indices
    std::vector<std::uint32_t> coordIds;
    elementCount[elemType] = readGroup(coordIds, numNodes, true);

    // Convert all indices into proper elements
    std::size_t elemOffset = elements.size();
    elements.resize(elemOffset + elementCount[elemType]);
    for (std::size_t elem = 0; elem < elementCount[elemType]; ++elem) {
        elements[elemOffset + elem].zoneName = "UNSET";
        std::vector<std::size_t>& localCoordIds =
            elements[elemOffset + elem].coordinateIds;
        localCoordIds.resize(numNodes);
        for (std::size_t c = 0; c < numNodes; ++c) {
            // -1 due to Fortran indexing
            localCoordIds[c] = coordIds[elem * numNodes + c] - 1;
            // Reorder to hpgem ordering
        }
        reorderElementCoords(localCoordIds, elemType);
    }
    logger(VERBOSE, "Read % elements with % nodes", elementCount[elemType],
           numNodes);
}

std::uint32_t CentaurReader::skipGroup(std::size_t linesPerEntity,
                                       bool multiline) {
    GroupSize groupSize = readGroupSize(multiline);
    if (groupSize.totalCount == 0) {
        centaurFile.skipRecord(0);
    } else {
        for (std::size_t i = 0; i < linesPerEntity; ++i) {
            for (std::size_t read = 0; read < groupSize.totalCount;
                 read += groupSize.perLineCount) {
                // TODO: Remove the need for peeking
                std::uint32_t size = centaurFile.peekRecordSize();
                centaurFile.skipRecord(size);
            }
        }
    }
    return groupSize.totalCount;
}

void CentaurReader::readPeriodicNodeConnections() {
    std::uint32_t numberOfPeriodicTransformations = 0;
    if (centaurFileType > 3 || is2D()) {
        numberOfPeriodicTransformations = readSingleLineGroupSize();
    } else {
        logger.assert_always(false, "Too old format");
    }

    SparseUnionFind connectedCoords;

    for (std::size_t periodic = 0; periodic < numberOfPeriodicTransformations;
         ++periodic) {
        // For each periodic association the following information is stored:
        // - Transformation 'forward' coord1 -> coord2
        // - Transformation 'backward' coord2 -> coord1
        // - Number of paired nodes
        // - Node pairings (coord1, coord2)

        // Each transform is a uint32 with unknown meaning
        // Followed by a (dim+1)^2 matrix describing the rotation and
        // translation. We skip this information.
        std::uint32_t dim = getDimension() + 1;
        std::uint32_t transformSize =
            sizeof(std::uint32_t) + dim * dim * sizeof(double);
        centaurFile.skipRecord(transformSize);
        centaurFile.skipRecord(transformSize);

        // Number of Periodic pairs
        // Periodic pairs
        std::size_t periodicPairs = readSingleLineGroupSize();
        std::vector<std::uint32_t> pairs(2 * periodicPairs);
        centaurFile.readRawRecord(2 * periodicPairs * sizeof(std::uint32_t),
                                  reinterpret_cast<char*>(pairs.data()));
        for (std::size_t i = 0; i < periodicPairs; ++i) {
            // -1 due for Fortran indices
            std::size_t node0 = pairs[2 * i + 0] - 1;
            std::size_t node1 = pairs[2 * i + 1] - 1;
            // Merge the nodes
            connectedCoords.unionSets(node0, node1);
        }
    }
    // Merge the actual node ids.
    for (const std::size_t coord : connectedCoords) {
        std::size_t representative = connectedCoords.findSet(coord);
        coordinates[coord].nodeId = representative;
    }
    logger(INFO, "Read % periodic connections",
           numberOfPeriodicTransformations);
}

void CentaurReader::readBoundaryGroups() {
    GroupSize groupSize = readGroupSize(false);
    std::uint32_t numberOfGroups = groupSize.totalCount;
    logger(INFO, "Discarding % boundary groups", numberOfGroups);

    std::vector<std::uint32_t> boundaryTypeIds(numberOfGroups);
    std::vector<std::array<char, 80>> rawBoundaryNames(numberOfGroups);

    centaurFile.readRawRecord(numberOfGroups * sizeof(std::uint32_t),
                              reinterpret_cast<char*>(boundaryTypeIds.data()));
    centaurFile.readRawRecord(numberOfGroups * 80 * sizeof(char),
                              reinterpret_cast<char*>(rawBoundaryNames.data()));
    // Output for debugging
    for (std::uint32_t i = 0; i < numberOfGroups; ++i) {
        logger(VERBOSE, "Boundary group \"%\" of type %'", rawBoundaryNames[i],
               boundaryTypeIds[i]);
    }
}

void CentaurReader::readZoneInfo() {
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

    // For 2D only the positions of 2 element types are stored
    // For 3D there are 4 element types and additionally the starting boundary
    // face is stored.
    std::size_t numIndexEntries = is3D() ? 5 : 2;
    // For 2D the zone names are stored separately
    std::size_t byteSize = numIndexEntries * sizeof(std::uint32_t) +
                           (is3D() ? 80 * sizeof(char) : 0);
    byteSize *= numberOfZones;
    std::vector<char> rawZoneData(byteSize);
    centaurFile.readRawRecord(byteSize, rawZoneData.data());

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
                memcpy(&entry, rawZoneData.data() + offset, sizeof(entry));
                offset += sizeof(entry);
                logger.assert_always(
                    entry == 1, "Offset for the first zone entry is non zero.");
            }
        } else {
            ZoneInformation& prevZone = zones[i - 1];

            for (std::size_t entryId = 0; entryId < numIndexEntries;
                 ++entryId) {
                memcpy(&(prevZone.endOffsets[entryId]),
                       rawZoneData.data() + offset, sizeof(std::uint32_t));
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
        centaurFile.readRawRecord(numberOfZones * 80 * sizeof(char),
                                  names.data());
        for (std::size_t i = 0; i < numberOfZones; ++i) {
            auto start = names.begin() + 80 * i;
            std::copy(start, start + 80, zones[i].rawZoneName.begin());
        }
    }

    if (is3D() && version > 4.99) {
        // Read zone families
        std::vector<char> familyNames(numberOfZones * 80);
        centaurFile.readRawRecord(numberOfZones * 80 * sizeof(char),
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
    logger(INFO, "Applying % zones", zones.size());
    // For each element type, centaur ensures that the elements are stored in
    // zone order. So first all the cubes for the first zone, then for the
    // second zone etc. and in a separate array all the tetrahedra from the
    // first zone, followed by the tetrahedra from the second zone etc.
    //
    // We do not make any distinction between the different types of elements
    // and just concatenate them into one large vector with elements. This
    // complicates the assignment of zone names, as we need to do it (up to)
    // four different positions in the element vector.

    // Index in the element vector for the four different element types.
    std::array<std::size_t, 4> positions({0, 0, 0, 0});
    // Set at the start of each element type
    for (std::size_t elemType = 1; elemType < 4; ++elemType) {
        positions[elemType] =
            positions[elemType - 1] + elementCount[elemType - 1];
    }
    for (const ZoneInformation& zone : zones) {
        std::string zoneName = zone.getZoneName();
        // Assign the zone for each elementtype
        for (std::size_t elemType = 0; elemType < 4; ++elemType) {
            for (; positions[elemType] < zone.endOffsets[elemType];
                 ++positions[elemType]) {
                elements[positions[elemType]].zoneName = zoneName;
            }
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
    centaurFile.readRawRecord(sizeof(std::uint32_t) * (1 + multiline),
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
        centaurFile.readRawRecord(0, reinterpret_cast<char*>(data.data()));
    } else {
        std::uint32_t entriesRead = 0;
        while (entriesRead < groupSize.totalCount) {
            // Compute the number of entries on this line
            std::uint32_t lineSize = std::min(
                groupSize.perLineCount, groupSize.totalCount - entriesRead);
            // Read the actual entries
            std::size_t offset = entriesRead * numComponents;
            std::size_t lineByteSize = lineSize * numComponents * sizeof(T);
            centaurFile.readRawRecord(
                lineByteSize, reinterpret_cast<char*>(data.data() + offset));
            entriesRead += lineSize;
        }
        logger.assert_always(entriesRead == groupSize.totalCount,
                             "Incorrect number of entries read");
    }
    return groupSize.totalCount;
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

void CentaurReader::reorderElementCoords(std::vector<std::size_t>& coords,
                                         std::size_t elementType) const {
    if (is2D()) {
        logger.assert_always(elementType < 2, "Invalid element type % for 2D",
                             elementType);
        if (elementType == 1) {
            // Usual variation for square faces. Centaur uses the clockwise
            // ordering of coordinates, hpgem the lexicographical.
            std::swap(coords[2], coords[3]);
        }
        // Triangles are the same by definition.
    } else {
        logger.assert_always(elementType < 4, "Invalid element type % for 3D",
                             elementType);
        if (elementType == 0) {
            // Same as for squares, but with the top and bottom face of the cube
            std::swap(coords[2], coords[3]);
            std::swap(coords[6], coords[7]);
        } else if (elementType == 2) {
            // Same as for squares
            std::swap(coords[2], coords[3]);
            // Centaur has the top of the pyramid as last node, hpgem as first
            // one
            std::size_t top = coords[4];
            for (std::size_t i = 4; i > 0; --i) {
                coords[i] = coords[i - 1];
            }
            coords[0] = top;
        }
        // Prism are the same
        // Tetrahedra are the same (by definition)
    }
}

}  // namespace Preprocessor
