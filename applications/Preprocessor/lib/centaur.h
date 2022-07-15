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

#ifndef HPGEM_APP_CENTAUR_H
#define HPGEM_APP_CENTAUR_H

#include "customIterator.h"
#include "FortranUnformattedFile.h"
#include "MeshSource.h"
#include <array>
#include <vector>
#include <fstream>
#include <string>
#include <map>

namespace Preprocessor {

using namespace hpgem;

class CentaurReader : public MeshSource2 {
   public:
    CentaurReader(std::string filename);

    const std::vector<MeshSource2::Coord>& getCoordinates() final {
        return coordinates;
    }
    const std::vector<MeshSource2::Element>& getElements() final {
        return elements;
    }

    const std::vector<std::map<std::size_t, std::size_t>>& getMerges() final {
        return coordinateMerges;
    };

    std::size_t getDimension() const final {
        if (centaurFileType > 0) return 3;

        return 2;
    }

    bool is2D() const { return getDimension() == 2; }
    bool is3D() const { return getDimension() == 3; }

   private:
    /// Helper methods ///
    //////////////////////
    //
    // Helper methods for how centaur uses the fortran file format.

    struct GroupSize {
        /**
         * Total number of entries
         */
        std::uint32_t totalCount;
        /**
         * Entries per line
         */
        std::uint32_t perLineCount;
    };

    // returns the number of skipped entities
    std::uint32_t skipGroup(std::size_t linesPerEntity = 1,
                            bool multiline = true);

    /**
     * Read the size of the next group
     * @param expectMultiline Whether it can be a multiline group (if the file
     * type supports it)
     * @return The size of the next group.
     */
    GroupSize readGroupSize(bool expectMultiline);

    /**
     * Read the size of the next group, assuming it does not support multiline
     * @return The size of the group.
     */
    std::uint32_t readSingleLineGroupSize() {
        GroupSize groupSize = readGroupSize(false);
        return groupSize.totalCount;
    }

    /**
     * Read a group with a uniform data type
     * @tparam T The datatype
     * @param data vector to store the result in
     * @param numComponents The number of components per entry (e.g. 4 corner
     * indices for squares)
     * @param expectMultiline Whether it can be stored as multiline when the
     * version supports it.
     * @return The number of entries that were read
     */
    template <typename T>
    std::uint32_t readGroup(std::vector<T>& data, std::uint32_t numComponents,
                            bool expectMultiline);

    /// Reader Methods ///
    //////////////////////
    //
    // Methods to read parts of the centaur file. Following the same ordering as
    // in the file.

    void readHeader();
    void readCoordinates();
    /**
     * Read a group defining the elements
     * @param elemType the element index
     * @param numNodes The number of nodes for this type of element
     */
    void readElements(std::size_t elemType, std::uint32_t numNodes);

    /**
     * The number of nodes each centaur element type has.
     * @param elementType The element type
     * @return The number of nodes for one element.
     */
    std::uint32_t numNodes(std::size_t elementType) const;

    /**
     * Reorder coordinates to follow hpgem ordering.
     *
     * For some element types the ordering of the coordinates is differs between
     * hpgem and centaur. This function reorders a set of coordinates in centaur
     * ordering to follow hpgem ordering.
     * @param coords The coordinates
     * @param elementType The element type index.
     */
    void reorderElementCoords(std::vector<std::size_t>& coords,
                              std::size_t elementType) const;

    /// Read & discard boundary group information
    /// Note: While the information is discarded, it is very useful in debugging
    /// The boundary group names are human readable and can be easily checked
    /// for correctness.
    void readBoundaryGroups();

    /**
     * Read and apply the periodic connections from the centaur file.
     */
    void readPeriodicNodeConnections();

    /// Read the zone information from the file and update the elements
    void readZoneInfo();

    /// DATA ///
    ////////////
    // All data from the file, following the same ordering as in the file.

    FortranUnformattedFile centaurFile;

    // Header information
    /**
     * Version of the format
     */
    float version;
    /**
     * Type of the hybrid file. Negative values are 2D grids, positive values
     * are 3D.
     */
    std::int32_t centaurFileType;

    /**
     * Coordinates
     */
    std::vector<MeshSource2::Coord> coordinates;
    /**
     * The elements in the mesh. The different element types are placed
     * consecutively.
     */
    std::vector<MeshSource2::Element> elements;
    /**
     * The number of elements of each type.
     */
    std::array<std::uint32_t, 4> elementCount;

    /**
     * A series of coordinate identifications from periodic boundary conditions.
     */
    std::vector<std::map<std::size_t, std::size_t>> coordinateMerges;

    struct ZoneInformation {
        /// The end-offset for element types and (3D only) boundary faces.
        /// Each end offset gives the number of elements/boundary faces in this
        /// zone and all previous zones. Thus when this index is reached the
        /// next zone starts.
        /// For 3D the ordering is:
        /// Hexahedra, Prisms, Pyramids, Tetrahedra, boundary faces
        /// For 2D the ordering is: Triangles, Quadrilaterals
        std::array<std::uint32_t, 5> endOffsets;
        /// Name of the zone, not null terminated, space filled
        std::array<char, 80> rawZoneName;
        /// Name of the zone family (optional)
        /// not null terminated, space filled
        std::array<char, 80> rawZoneFamiliyName;

        std::string getZoneName() const;
        std::string getZoneFamilyName() const;
    };

    /// Storage for the zones.
    std::vector<ZoneInformation> zones;
};
}  // namespace Preprocessor

#endif  // HPGEM_APP_CENTAUR_H
