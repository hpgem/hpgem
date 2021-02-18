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

#include "gmsh.h"
#include "Logger.h"

using namespace hpgem;

namespace Preprocessor {

//
// Non-class helper function for line-reading
// Because getline works only on text files from the OS that it was compiled
// on...
// From:
// https://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf
//
std::istream& safeGetline(std::istream& is, std::string& t) {
    t.clear();

    // The characters in the stream are read one-by-one using a std::streambuf.
    // That is faster than reading them one-by-one using the std::istream.
    // Code that uses streambuf this way must be guarded by a sentry object.
    // The sentry object performs various tasks,
    // such as thread synchronization and updating the stream state.

    std::istream::sentry se(is, true);
    std::streambuf* sb = is.rdbuf();

    for (;;) {
        int c = sb->sbumpc();
        switch (c) {
            case '\n':
                return is;
            case '\r':
                if (sb->sgetc() == '\n') sb->sbumpc();
                return is;
            case std::streambuf::traits_type::eof():
                // Also handle the case when the last line has no line ending
                if (t.empty()) is.setstate(std::ios::eofbit);
                return is;
            default:
                t += (char)c;
        }
    }
}

// find a searchname in the MSH file and return the file stream after that
// line
bool locate_in_file(std::ifstream& filestream, const std::string& searchname) {
    // to find a specific keyword in the MSH file and return the file stream
    std::string temp;
    bool found = false;
    while (!filestream.eof()) {
        safeGetline(filestream, temp);
        found = false;
        if (temp == searchname) {
            found = true;
            break;
        }
    }
    logger.assert_always(found, "The " + searchname + "field was not found");
    return found;
}

GmshReader::GmshReader(std::string filename) {
    // Fileformat is defined in
    // https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
    Filehandle_.open(filename);
    logger.assert_always(Filehandle_.is_open(), "Cannot open msh meshfile.");
    logger.assert_always(Filehandle_.good(),
                         "Something is not so good about this mesh");

    bool found_header = locate_in_file(Filehandle_, "$MeshFormat");
    if (!found_header) {
        logger.assert_always(found_header, "Header of gmsh file % not found",
                             filename);
    } else {
        double version;
        int filetype, datasize;
        Filehandle_ >> version >> filetype >> datasize;

        logger(VERBOSE, "msh version is % and storage format is %", version,
               (filetype == 1 ? "binary" : "ASCII"));
        logger.assert_always(datasize == sizeof(size_t),
                             "Datasize is % vs sizeof(size_t) is %", datasize,
                             sizeof(size_t));
        if (filetype == 1) {
            int endianness;
            Filehandle_ >> endianness;
            if (endianness == 1) {
                logger(VERBOSE, "endianness of machine and file agree");
            }
        }
        logger.assert_always(version >= 4.1,
                             "Error: Cannot read file versions less than 4.1");
        logger.assert_always(filetype == 0,
                             "Error: Cannot read binary format files");
    }
    bool found_nodes = locate_in_file(Filehandle_, "$Nodes");
    if (!found_nodes) {
        logger.assert_always(
            found_header, "Nodes section of gmsh file %  not found", filename);
    } else {
        size_t nodes_total_entities;
        size_t raw_N_nodes;
        size_t nodes_min_index;
        size_t nodes_max_index;
        Filehandle_ >> nodes_total_entities >> raw_N_nodes >> nodes_min_index >>
            nodes_max_index;

        bool is_continuous_data =
            (raw_N_nodes == (nodes_max_index - nodes_min_index + 1));
        logger.assert_always(is_continuous_data,
                             "Error: Cannot read non continous nodes");
        nodes_.reserve(raw_N_nodes);
        for (size_t block_index = 0; block_index < nodes_total_entities;
             block_index++) {
            int entityDim;
            int entityTag;
            int parametric;
            size_t numNodesInBlock;
            Filehandle_ >> entityDim >> entityTag >> parametric >>
                numNodesInBlock;
            std::cout << entityDim << " " << entityTag << " " << parametric
                      << " " << numNodesInBlock << std::endl;

            logger.assert_always(parametric == 0,
                                 "Error: Cannot read parametric coordinates");

            for (size_t nodeindex = 0; nodeindex < numNodesInBlock;
                 nodeindex++) {
                size_t nodeTag;
                Filehandle_ >> nodeTag;
            }
            for (size_t nodeindex = 0; nodeindex < numNodesInBlock;
                 nodeindex++) {
                double x;
                double y;
                double z;
                Filehandle_ >> x >> y >> z;
                MeshSource::Node temp;
                temp.coordinates.push_back(std::vector<double>{x, y, z});
                nodes_.push_back(temp);
            }
        }
    }

    bool found_elements = locate_in_file(Filehandle_, "$Elements");
}

Range<MeshSource::Node> GmshReader::getNodeCoordinates() { ; }
Range<MeshSource::Element> GmshReader::getElements() { ; }

}  // namespace Preprocessor