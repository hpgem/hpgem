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

#include "MeshSource.h"
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

void GmshReader::FillElementTypeMap() {

    nodes_per_elementtype_[1] = 2;  // 2-node line.
    nodes_per_elementtype_[2] = 3;  // 3-node triangle.
    nodes_per_elementtype_[3] = 4;  // 4-node quadrangle.
    nodes_per_elementtype_[4] = 4;  // 4-node tetrahedron.
    nodes_per_elementtype_[5] = 8;  // 8-node hexahedron.
    nodes_per_elementtype_[6] = 6;  // 6-node prism.
    nodes_per_elementtype_[7] = 5;  // 5-node pyramid.
    nodes_per_elementtype_[8] =
        3;  // 3-node second order line (2 nodes associated with the vertices
            // and 1 with the edge).
    nodes_per_elementtype_[9] =
        6;  // 6-node second order triangle (3 nodes associated with the
            // vertices and 3 with the edges).
    nodes_per_elementtype_[10] =
        9;  // 9-node second order quadrangle (4 nodes associated with the
            // vertices, 4 with the edges and 1 with the face).
    nodes_per_elementtype_[11] =
        10;  // 10-node second order tetrahedron (4 nodes associated with the
             // vertices and 6 with the edges).
    nodes_per_elementtype_[12] =
        27;  // 27-node second order hexahedron (8 nodes associated with the
             // vertices, 12 with the edges, 6 with the faces and 1 with the
             // volume).
    nodes_per_elementtype_[13] =
        18;  // 18-node second order prism (6 nodes associated with the
             // vertices, 9 with the edges and 3 with the quadrangular faces).
    nodes_per_elementtype_[14] =
        14;  // 14-node second order pyramid (5 nodes associated with the
             // vertices, 8 with the edges and 1 with the quadrangular face).
    nodes_per_elementtype_[15] = 1;  // 1-node point.
    nodes_per_elementtype_[16] =
        8;  // 8-node second order quadrangle (4 nodes associated with the
            // vertices and 4 with the edges).
    nodes_per_elementtype_[17] =
        20;  // 20-node second order hexahedron (8 nodes associated with the
             // vertices and 12 with the edges).
    nodes_per_elementtype_[18] =
        15;  // 15-node second order prism (6 nodes associated with the vertices
             // and 9 with the edges).
    nodes_per_elementtype_[19] =
        13;  // 13-node second order pyramid (5 nodes associated with the
             // vertices and 8 with the edges).
    nodes_per_elementtype_[20] =
        9;  // 9-node third order incomplete triangle (3 nodes associated with
            // the vertices, 6 with the edges)
    nodes_per_elementtype_[21] =
        10;  // 10-node third order triangle (3 nodes associated with the
             // vertices, 6 with the edges, 1 with the face)
    nodes_per_elementtype_[22] =
        12;  // 12-node fourth order incomplete triangle (3 nodes associated
             // with the vertices, 9 with the edges)
    nodes_per_elementtype_[23] =
        15;  // 15-node fourth order triangle (3 nodes associated with the
             // vertices, 9 with the edges, 3 with the face)
    nodes_per_elementtype_[24] =
        15;  // 15-node fifth order incomplete triangle (3 nodes associated with
             // the vertices, 12 with the edges)
    nodes_per_elementtype_[25] =
        21;  // 21-node fifth order complete triangle (3 nodes associated with
             // the vertices, 12 with the edges, 6 with the face)
    nodes_per_elementtype_[26] =
        4;  // 4-node third order edge (2 nodes associated with the vertices, 2
            // internal to the edge)
    nodes_per_elementtype_[27] =
        5;  // 5-node fourth order edge (2 nodes associated with the vertices, 3
            // internal to the edge)
    nodes_per_elementtype_[28] =
        6;  // 6-node fifth order edge (2 nodes associated with the vertices, 4
            // internal to the edge)
    nodes_per_elementtype_[29] =
        20;  // 20-node third order tetrahedron (4 nodes associated with the
             // vertices, 12 with the edges, 4 with the faces)
    nodes_per_elementtype_[30] =
        35;  // 35-node fourth order tetrahedron (4 nodes associated with the
             // vertices, 18 with the edges, 12 with the faces, 1 in the volume)
    nodes_per_elementtype_[31] =
        56;  // 56-node fifth order tetrahedron (4 nodes associated with the
             // vertices, 24 with the edges, 24 with the faces, 4 in the volume)
    nodes_per_elementtype_[92] =
        64;  // 64-node third order hexahedron (8 nodes associated with the
             // vertices, 24 with the edges, 24 with the faces, 8 in the volume)
    nodes_per_elementtype_[93] =
        125;  // 125-node fourth order hexahedron (8 nodes associated with the
              // vertices, 36 with the edges, 54 with the faces, 27 in the
              // volume)
}

void GmshReader::ReadHeader() {

    bool found_header = locate_in_file(Filehandle_, "$MeshFormat");
    logger.assert_always(found_header, "Header of gmsh file not found");
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

void GmshReader::ReadNodes() {
    bool found_nodes = locate_in_file(Filehandle_, "$Nodes");
    logger.assert_always(found_nodes, "Nodes section of gmsh file not found");
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
        Filehandle_ >> entityDim >> entityTag >> parametric >> numNodesInBlock;

        logger.assert_always(parametric == 0,
                             "Error: Cannot read parametric coordinates");

        for (size_t nodeindex = 0; nodeindex < numNodesInBlock; nodeindex++) {
            size_t nodeTag;
            Filehandle_ >> nodeTag;
        }
        for (size_t nodeindex = 0; nodeindex < numNodesInBlock; nodeindex++) {
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

void GmshReader::ReadElements() {

    bool found_elements = locate_in_file(Filehandle_, "$Elements");
    logger.assert_always(found_elements,
                         "Elements section of gmsh file %  not found");
    size_t numEntityBlocks;
    size_t numElements;
    size_t minElementTag;
    size_t maxElementTag;
    Filehandle_ >> numEntityBlocks >> numElements >> minElementTag >>
        maxElementTag;

    bool is_continuous_data =
        (numElements == (maxElementTag - minElementTag + 1));
    logger.assert_always(is_continuous_data,
                         "Error: Cannot read non continous elements");
    elements_.reserve(numElements);
    for (size_t size_entity_index = 0; size_entity_index < numElements;
         size_entity_index++) {
        int entityDim;
        int entityTag;
        int elementType;
        size_t numElementsInBlock;
        Filehandle_ >> entityDim >> entityTag >> elementType >>
            numElementsInBlock;

        for (size_t size_entity_index = 0;
             size_entity_index < numElementsInBlock; size_entity_index++) {
            int entityDim;
            int entityTag;
            int elementType;
            int numElementsInBlock;
        }
    }
}

size_t GmshReader::DetermineDimension(double tol) const {

    size_t dimension = 3;

    bool two_dim = std::all_of(
        nodes_.begin(), nodes_.end(), [&](const MeshSource::Node& node) {
            return std::all_of(node.coordinates.begin(), node.coordinates.end(),
                               [&](const std::vector<double>& coordinate) {
                                   return (std::abs(coordinate.back()) < tol);
                               });
        });
    if (two_dim) {
        dimension--;
        bool one_dim = std::all_of(
            nodes_.begin(), nodes_.end(), [&](const MeshSource::Node& node) {
                return std::all_of(node.coordinates.begin(),
                                   node.coordinates.end(),
                                   [&](const std::vector<double>& coordinate) {
                                       return (std::abs(coordinate[1]) < tol);
                                   });
            });
        if (one_dim) {
            dimension--;
        }
    }
    logger(VERBOSE, "Dimension of the net is %", dimension);
    return dimension;
}  // namespace Preprocessor

GmshReader::GmshReader(std::string filename) {

    FillElementTypeMap();
    // Fileformat is defined in
    // https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
    Filehandle_.open(filename);
    logger.assert_always(Filehandle_.is_open(), "Cannot open msh meshfile.");
    logger.assert_always(Filehandle_.good(),
                         "Something is not so good about this mesh");

    ReadHeader();
    ReadNodes();
    dimension_ = DetermineDimension();
}

Range<MeshSource::Node> GmshReader::getNodeCoordinates() { ; }
Range<MeshSource::Element> GmshReader::getElements() { ; }

}  // namespace Preprocessor