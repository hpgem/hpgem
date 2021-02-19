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

#ifndef HPGEM_APP_GMSH_H
#define HPGEM_APP_GMSH_H

#include "MeshSource.h"
#include <array>
#include <vector>
#include <fstream>
#include <string>
#include <map>

// Most of this was taken from
// https://github.com/Applied-Scientific-Research/gmsh-reader/tree/main/src

using namespace hpgem;

namespace Preprocessor {

class GmshReader : public MeshSource {
   public:
    GmshReader(std::string filename);

    Range<MeshSource::Node> getNodeCoordinates() final;
    Range<MeshSource::Element> getElements() final;

    std::size_t getDimension() const final { return dimension_; }

   private:
    std::ifstream Filehandle_;

    std::vector<MeshSource::Node> nodes_;

    std::vector<MeshSource::Element> elements_;
    std::size_t dimension_;

    void ReadHeader();
    void ReadNodes();
    void ReadElements();
    size_t DetermineDimension(double tol = 1e-12) const;

    void FillElementTypeMap();

    // gmsh elementtype to how many nodes that element has
    std::map<int, int> nodes_per_elementtype_;
};
}  // namespace Preprocessor

#endif  // HPGEM_APP_CENTAUR_H
