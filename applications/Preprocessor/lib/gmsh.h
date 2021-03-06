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

    const std::vector<Coord>& getCoordinates() { return nodes_; }

    const std::vector<Element>& getElements() { return elements_; }

    std::size_t getDimension() const { return dimension_; }

   private:
    void fillElementTypeMap();
    void readHeader();
    void readNodes();
    void readElements();
    void readElementData();
    void readPBCs();
    void purgeLowerDimElements();
    void pruneCoordinatesToDimension();

    size_t determineDimension(double tol = 1e-12) const;

    size_t dimension_;

    std::vector<MeshSource2::Element> elements_;

    std::vector<MeshSource2::Coord> nodes_;

    std::ifstream Filehandle_;
    std::map<size_t, size_t> nodesPerElementtype_;
    std::map<size_t, size_t> dimensionOfElementtype_;
};

}  // namespace Preprocessor
#endif  // HPGEM_APP_GMSH_H