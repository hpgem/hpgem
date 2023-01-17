/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is azdistributed using BSD 3-Clause License. A copy of which can
 found below.


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

#ifndef HPGEM_APP_GMSH_H
#define HPGEM_APP_GMSH_H

#include <map>
#include <fstream>
#include <string>
#include <vector>

#include "MeshSource.h"
#include "Logger.h"

using namespace hpgem;

namespace Preprocessor {
/**
 * @brief Parser for gmsh files
 * This parser only works with gmsh files of version 2.2.
 * For more info see
 * https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format-version-2-_0028Legacy_0029
 */

class GmshReader final : public MeshSource2 {

   public:
    GmshReader(std::string filename);

    const std::vector<Coord>& getCoordinates() final { return nodes_; }

    const std::vector<Element>& getElements() final { return elements_; }

    const std::vector<std::map<std::size_t, std::size_t>>& getMerges() final {
        return coordinateMerges_;
    };

    std::size_t getDimension() const final { return dimension_; }

   private:
    // Element as used by gmsh
    struct Element {

        static constexpr const std::size_t NO_TAG =
            std::numeric_limits<std::size_t>::max();

        /// Tag of gmsh
        std::size_t gmshTag = NO_TAG;
        /// dimension of this entity/element
        std::size_t dimension_;

        /// Tag corresponding to the physical name
        /// May be overwritten by one from the element data
        std::size_t physicalNameTag = NO_TAG;

        // Reordered to hpgem ordering and using hpgem indices
        std::vector<std::size_t> coordinates_;
    };

    void fillElementTypeMap();
    /// Read the MeshFormat header section
    void readHeader();
    /// Read a PhysicalNames section
    void readPhysicalNames();
    /// Read a Nodes section
    void readNodes(double tol = 1e-12);
    /// Read an Elements section
    void readElements();
    /// Read an ElementData section
    void readElementData();
    /// Read a Periodic section
    void readPBCs();
    /// Post processing read 3D coordinates to the used dimension
    void pruneCoordinatesToDimension();
    /// Skip a section that is not of relevance
    ///
    /// \param sectionName The name of the section
    void skipSection(std::string sectionName);
    /// Read the end of-section marker for a section. Will consume any white
    /// space before it.
    ///
    /// \param sectionName The name of the section
    void readSectionEnd(std::string sectionName);

    /// Convert the raw gmsh elements from rawElements_ into hpgem elements
    /// in elements_. After this call rawElements_ is empty.
    void convertRawElements();

    double version_;
    size_t dimension_;

    // Physical names indexed by dimension and tag
    std::vector<std::map<std::size_t, std::string>> physicalNames_;

    /// Raw elements as read from gmsh. For gmsh every mesh entity is an
    /// element, so depending on the dimension this contains points, lines etc.
    /// Note that 'Entity' in gmsh files refers to model description parts.
    std::vector<Element> rawElements_;
    /// Actual elements for hpgem
    std::vector<MeshSource2::Element> elements_;

    std::vector<MeshSource2::Coord> nodes_;

    std::ifstream Filehandle_;
    std::map<size_t, size_t> nodesPerElementtype_;
    std::map<size_t, size_t> dimensionOfElementtype_;
    /**
     * A series of coordinate identifications from periodic boundary conditions.
     */
    std::vector<std::map<std::size_t, std::size_t>> coordinateMerges_;
};

}  // namespace Preprocessor
#endif  // HPGEM_APP_GMSH_H