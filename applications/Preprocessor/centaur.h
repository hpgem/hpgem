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
#include "unstructuredFile.h"
#include <vector>
#include <fstream>
#include <string>
#include <map>

using namespace hpgem;

namespace Preprocessor {

class CentaurReader {
   public:
    CentaurReader(std::string filename);

    Range<std::vector<std::vector<double>>> getNodeCoordinates();
    Range<std::vector<std::size_t>> getElements();

    std::size_t getDimension() {
        if (centaurFileType > 0) return 3;

        return 2;
    }

    bool is2D() { return getDimension() == 2; }
    bool is3D() { return getDimension() == 3; }

   private:
    UnstructuredInputStream<std::istringstream> readLine();
    // returns the number of skipped entities
    std::uint32_t skipGroup(std::size_t linesPerEntity = 1,
                            bool multiline = true);

    void readHeader();

    void readPeriodicNodeConnections();

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

    UnstructuredInputStream<std::ifstream> centaurFile;
    /// Location in the file where the node/coordinate information starts
    std::ifstream::pos_type nodeStart;
    /// Location in the file where the element definition starts
    std::ifstream::pos_type elementStart;

    /// Number of Nodes
    std::uint32_t numberOfNodes;
    /// Number of Elements (independent of type)
    std::uint32_t numberOfElements;

    std::map<std::uint32_t, std::vector<std::size_t>> boundaryConnections;

    std::vector<std::size_t> toHpgemNumbering;
};
}  // namespace Preprocessor

#endif  // HPGEM_APP_CENTAUR_H
