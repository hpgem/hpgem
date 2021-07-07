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

#include <algorithm>
#include <iostream>
#include <string>
#include <unordered_map>

#include "MeshSource.h"
#include "gmsh.h"
#include "ElementReorder.h"
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
// modified from
// https://stackoverflow.com/questions/33022697/find-word-in-a-text-in-c-and-print-some-next-specific-lines
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
    return found;
}

void GmshReader::fillElementTypeMap() {

    nodesPerElementtype_[1] = 2;  // 2-node line.
    dimensionOfElementtype_[1] = 1;
    nodesPerElementtype_[2] = 3;  // 3-node triangle.
    dimensionOfElementtype_[2] = 2;
    nodesPerElementtype_[3] = 4;  // 4-node quadrangle.
    dimensionOfElementtype_[3] = 2;
    nodesPerElementtype_[4] = 4;  // 4-node tetrahedron.
    dimensionOfElementtype_[4] = 3;
    nodesPerElementtype_[5] = 8;  // 8-node hexahedron.
    dimensionOfElementtype_[5] = 3;
    nodesPerElementtype_[6] = 6;  // 6-node prism.
    dimensionOfElementtype_[6] = 3;
    nodesPerElementtype_[7] = 5;  // 5-node pyramid.
    dimensionOfElementtype_[7] = 3;
    nodesPerElementtype_[8] =
        3;  // 3-node second order line (2 nodes associated with the vertices
            // and 1 with the edge).
    dimensionOfElementtype_[8] = 1;
    nodesPerElementtype_[9] =
        6;  // 6-node second order triangle (3 nodes associated with the
            // vertices and 3 with the edges).
    dimensionOfElementtype_[9] = 2;
    nodesPerElementtype_[10] =
        9;  // 9-node second order quadrangle (4 nodes associated with the
    // vertices, 4 with the edges and 1 with the face).
    dimensionOfElementtype_[10] = 2;
    nodesPerElementtype_[11] =
        10;  // 10-node second order tetrahedron (4 nodes associated with the
             // vertices and 6 with the edges).
    dimensionOfElementtype_[11] = 3;
    nodesPerElementtype_[12] =
        27;  // 27-node second order hexahedron (8 nodes associated with the
             // vertices, 12 with the edges, 6 with the faces and 1 with the
             // volume).
    dimensionOfElementtype_[12] = 3;
    nodesPerElementtype_[13] =
        18;  // 18-node second order prism (6 nodes associated with the
             // vertices, 9 with the edges and 3 with the quadrangular faces).
    dimensionOfElementtype_[13] = 3;
    nodesPerElementtype_[14] =
        14;  // 14-node second order pyramid (5 nodes associated with the
             // vertices, 8 with the edges and 1 with the quadrangular face).
    dimensionOfElementtype_[14] = 3;
    nodesPerElementtype_[15] = 1;  // 1-node point.
    dimensionOfElementtype_[15] = 0;
    nodesPerElementtype_[16] =
        8;  // 8-node second order quadrangle (4 nodes associated with the
            // vertices and 4 with the edges).
    dimensionOfElementtype_[16] = 2;
    nodesPerElementtype_[17] =
        20;  // 20-node second order hexahedron (8 nodes associated with the
             // vertices and 12 with the edges).
    dimensionOfElementtype_[17] = 3;
    nodesPerElementtype_[18] =
        15;  // 15-node second order prism (6 nodes associated with the vertices
             // and 9 with the edges).
    dimensionOfElementtype_[18] = 3;
    nodesPerElementtype_[19] =
        13;  // 13-node second order pyramid (5 nodes associated with the
             // vertices and 8 with the edges).
    dimensionOfElementtype_[19] = 3;
    nodesPerElementtype_[20] =
        9;  // 9-node third order incomplete triangle (3 nodes associated with
            // the vertices, 6 with the edges)
    dimensionOfElementtype_[20] = 2;
    nodesPerElementtype_[21] =
        10;  // 10-node third order triangle (3 nodes associated with the
             // vertices, 6 with the edges, 1 with the face)
    dimensionOfElementtype_[21] = 2;
    nodesPerElementtype_[22] =
        12;  // 12-node fourth order incomplete triangle (3 nodes associated
             // with the vertices, 9 with the edges)
    dimensionOfElementtype_[22] = 2;
    nodesPerElementtype_[23] =
        15;  // 15-node fourth order triangle (3 nodes associated with the
             // vertices, 9 with the edges, 3 with the face)
    dimensionOfElementtype_[23] = 2;
    nodesPerElementtype_[24] =
        15;  // 15-node fifth order incomplete triangle (3 nodes associated with
             // the vertices, 12 with the edges)
    dimensionOfElementtype_[24] = 2;
    nodesPerElementtype_[25] =
        21;  // 21-node fifth order complete triangle (3 nodes associated with
             // the vertices, 12 with the edges, 6 with the face)
    dimensionOfElementtype_[25] = 2;
    nodesPerElementtype_[26] =
        4;  // 4-node third order edge (2 nodes associated with the vertices, 2
            // internal to the edge)
    dimensionOfElementtype_[26] = 1;
    nodesPerElementtype_[27] =
        5;  // 5-node fourth order edge (2 nodes associated with the vertices, 3
    // internal to the edge)
    dimensionOfElementtype_[27] = 1;
    nodesPerElementtype_[28] =
        6;  // 6-node fifth order edge (2 nodes associated with the vertices, 4
            // internal to the edge)
    dimensionOfElementtype_[28] = 1;
    nodesPerElementtype_[29] =
        20;  // 20-node third order tetrahedron (4 nodes associated with the
             // vertices, 12 with the edges, 4 with the faces)
    dimensionOfElementtype_[29] = 3;
    nodesPerElementtype_[30] =
        35;  // 35-node fourth order tetrahedron (4 nodes associated with the
             // vertices, 18 with the edges, 12 with the faces, 1 in the volume)
    dimensionOfElementtype_[30] = 3;
    nodesPerElementtype_[31] =
        56;  // 56-node fifth order tetrahedron (4 nodes associated with the
             // vertices, 24 with the edges, 24 with the faces, 4 in the volume)
    dimensionOfElementtype_[31] = 3;
    nodesPerElementtype_[92] =
        64;  // 64-node third order hexahedron (8 nodes associated with the
             // vertices, 24 with the edges, 24 with the faces, 8 in the volume)
    dimensionOfElementtype_[92] = 3;
    nodesPerElementtype_[93] =
        125;  // 125-node fourth order hexahedron (8 nodes associated with the
              // vertices, 36 with the edges, 54 with the faces, 27 in the
              // volume)
    dimensionOfElementtype_[93] = 3;
}

void GmshReader::readHeader() {

    bool found_header = locate_in_file(Filehandle_, "$MeshFormat");
    logger.assert_always(found_header, "Header of gmsh file not found");
    double version;
    int filetype, datasize;
    Filehandle_ >> version >> filetype >> datasize;

    logger(VERBOSE, "msh version is % and storage format is %", version,
           (filetype == 1 ? "binary" : "ASCII"));

    logger.assert_always(version == 2.2,
                         "Error: Cannot read file versions other than 2.2");
    logger.assert_always(filetype == 0,
                         "Error: Cannot read binary format files");
}

void GmshReader::readNodes() {
    bool found_nodes = locate_in_file(Filehandle_, "$Nodes");
    logger.assert_always(found_nodes, "Nodes section of gmsh file not found");
    size_t number_nodes;
    Filehandle_ >> number_nodes;
    size_t nodeTag;
    double x;
    double y;
    double z;
    nodes_.reserve(number_nodes);
    for (size_t index = 0; index < number_nodes; index++) {
        Filehandle_ >> nodeTag >> x >> y >> z;
        MeshSource2::Coord temp;
        temp.nodeId =
            nodeTag - 1;  // gmsh ids are 1 indexed hpgem is zero indexed
        temp.coordinate = std::vector<double>{x, y, z};
        nodes_.push_back(temp);
    }
}

void GmshReader::readElements() {

    ElementReorder reorder;
    reorder.addElementType(1, "line:", {0, 1});
    reorder.addElementType(2, "triangle:", {0, 1, 2});
    reorder.addElementType(2, "square:", {0, 1, 3, 2});
    reorder.addElementType(3, "tetrahedron:", {0, 3, 1, 2});
    reorder.addElementType(3, "pyramid:", {4, 1, 2, 0, 3});
    reorder.addElementType(3, "cube:", {4, 5, 0, 1, 7, 6, 3, 2});
    reorder.addElementType(3, "prism", {0, 2, 1, 3, 5, 4});

    bool found_elements = locate_in_file(Filehandle_, "$Elements");
    logger.assert_always(found_elements,
                         "Elements section of gmsh file  not found");
    size_t numElements;
    Filehandle_ >> numElements;

    elements_.reserve(numElements);
    for (size_t index = 0; index < numElements; index++) {
        size_t elm_number;
        size_t elm_type;
        size_t no_tags;

        Filehandle_ >> elm_number >> elm_type >> no_tags;
        MeshSource2::Element element;
        element.dimension = dimensionOfElementtype_.at(elm_type);
        element.id =
            elm_number - 1;  // gmsh ids are 1 indexed hpgem is zero indexed
        for (int i = 0; i < no_tags; i++) {
            int dummy;
            Filehandle_ >> dummy;
        }
        int num_points = nodesPerElementtype_.at(elm_type);
        for (int i = 0; i < num_points; i++) {
            size_t nodetag;
            Filehandle_ >> nodetag;
            element.coordinateIds.push_back(
                nodetag - 1);  // gmsh ids are 1 indexed hpgem is zero indexed
        }
        reorder.reorderToHpGem(element.dimension, element.coordinateIds);
        elements_.push_back(element);
    }
}

size_t GmshReader::determineDimension(double tol) const {

    size_t dimension = 3;

    bool two_dim = std::all_of(
        nodes_.begin(), nodes_.end(), [&](const MeshSource2::Coord& node) {
            return (std::abs(node.coordinate.back()) < tol);
        });

    if (two_dim) {
        dimension--;
        bool one_dim = std::all_of(
            nodes_.begin(), nodes_.end(), [&](const MeshSource2::Coord& node) {
                return (std::abs(node.coordinate[1]) < tol);
            });
        if (one_dim) {
            dimension--;
        }
    }
    logger(VERBOSE, "Dimension of the net is %", dimension);
    return dimension;
}

void GmshReader::readElementData() {
    // The zones are determined based on the ElementData section
    bool found_elements = locate_in_file(Filehandle_, "$ElementData");
    if (!found_elements) {
        // Very simple mesh files, like those generated by gmsh itself, may not
        // contain element data. So assign a single zone to all elements.
        for (auto& element : elements_) {
            element.zoneName = "Main";
        }
        return;
    }
    size_t numStringelements;
    Filehandle_ >> numStringelements;
    for (size_t i = 0; i < numStringelements; i++) {
        std::string view_names;
        Filehandle_ >> view_names;
    }
    size_t numrealelements;
    Filehandle_ >> numrealelements;
    for (size_t i = 0; i < numrealelements; i++) {
        double real;
        Filehandle_ >> real;
    }
    size_t numintegerelements;
    Filehandle_ >> numintegerelements;
    logger.assert_always(numintegerelements > 2,
                         "There is no Elementspecific data in gmsh file");
    size_t timestep;
    size_t field_components;
    size_t num_entities;
    Filehandle_ >> timestep >> field_components >> num_entities;

    logger.assert_always(field_components == 1,
                         "We only read in Scalar data as material identifiers");
    for (int i = 3; i < numintegerelements; i++) {
        size_t dummy;
        Filehandle_ >> dummy;
    }
    for (int i = 0; i < num_entities; i++) {
        size_t elementid;
        size_t zoneinfo;
        Filehandle_ >> elementid >> zoneinfo;
        elementid--;  // gmsh is 1 indexed hpgem is zero indexed
        logger.assert_always(elements_[elementid].id == elementid,
                             "Your elements are not indexed continously. This "
                             "parser cannot handle it. % vs %",
                             elements_[elementid].id, elementid);
        elements_[elementid].zoneName = std::to_string(zoneinfo);
    }
}

void GmshReader::readPBCs() {

    bool found_PBCs = locate_in_file(Filehandle_, "$Periodic");
    if (!found_PBCs) {
        Filehandle_.clear();
        Filehandle_.seekg(0);
        return;
    }

    size_t num_periodic_elements;
    Filehandle_ >> num_periodic_elements;

    for (size_t i = 0; i < num_periodic_elements; i++) {
        size_t dimension, entity_tag, master_entity_tag;
        Filehandle_ >> dimension >> entity_tag >> master_entity_tag;
    }

    std::unordered_map<size_t, size_t> pbc_nodes;

    size_t periodic_nodes;
    Filehandle_ >> periodic_nodes;
    logger(VERBOSE, "Found % periodic nodes", periodic_nodes);
    for (size_t i = 0; i < periodic_nodes; i++) {
        size_t node1;
        size_t node2;
        Filehandle_ >> node2 >> node1;
        // gmsh is 1 indexed hpgem is zero indexed
        pbc_nodes[node2 - 1] = node1 - 1;
    }

    // because element ids are indices into the coordinate vector for the pbcs
    // all we have to do is renumber the node ids.

    size_t replaced_nodes = 0;
    for (auto& node : nodes_) {
        auto found = pbc_nodes.find(node.nodeId);
        if (found != pbc_nodes.end()) {
            node.nodeId = found->second;
            replaced_nodes++;
        }
    }
    logger.assert_always(replaced_nodes == pbc_nodes.size(),
                         "Not all pbc pairs were used to replace node ids, so "
                         "probably the indexing is wrong. % pairs unused",
                         pbc_nodes.size() - replaced_nodes);
}

void GmshReader::purgeLowerDimElements() {

    // gmsh elements can be of any dimension smaller than the mesh dimension
    // we sort them to the back of the vector and then resize it
    elements_.erase(std::remove_if(elements_.begin(), elements_.end(),
                                   [&](const Element& e) {
                                       return e.dimension < dimension_;
                                   }),
                    elements_.end());
};

void GmshReader::pruneCoordinatesToDimension() {
    if (dimension_ == 3) {
        return;
    }
    for (auto& node : nodes_) {
        node.coordinate.resize(dimension_);
    }
}

GmshReader::GmshReader(std::string filename) {

    fillElementTypeMap();
    // Fileformat is defined in
    // https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
    Filehandle_.open(filename);
    logger.assert_always(Filehandle_.is_open(), "Cannot open msh meshfile. %",
                         filename);
    logger.assert_always(Filehandle_.good(),
                         "Something is not so good about this mesh in file %",
                         filename);

    readHeader();
    readNodes();
    dimension_ = determineDimension();
    readElements();
    readPBCs();
    readElementData();

    purgeLowerDimElements();
    pruneCoordinatesToDimension();
}

}  // namespace Preprocessor
