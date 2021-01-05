/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2014, University of Twente
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

#include "Logger.h"
#include "VTKSpecificTimeWriter.h"
#include "Base/MpiContainer.h"

#include "Base/MeshManipulator.h"
#include "Geometry/ElementGeometry.h"
#include "Base/Element.h"
#include "Geometry/PointPhysical.h"
#include "Geometry/PhysicalGeometry.h"
#include "Geometry/ReferencePoint.h"
#include "Geometry/ReferenceLine.h"
#include "Geometry/ReferenceTriangle.h"
#include "Geometry/ReferenceSquare.h"
#include "Geometry/ReferenceTetrahedron.h"
#include "Geometry/ReferenceCube.h"
#include "Geometry/ReferenceTriangularPrism.h"
#include "Geometry/ReferencePyramid.h"
#include "Geometry/PointReference.h"
#include "base64.h"
#include "VTKElementOrdering.h"
#include <vector>
#include <unordered_map>

#include <typeindex>

namespace hpgem {

namespace Output {

/////////////////////////////////////
// some VTK specific helper routines//
/////////////////////////////////////

// vtk element types supported by hpGEM
// magic numbers taken from http://www.vtk.org/VTK/img/file-formats.pdf
// smaller underlying type allows for easy conversion to base64
enum class VTKElementName : std::uint8_t {
    VERTEX = 1,
    LINE = 3,
    TRIANGLE = 5,
    QUAD = 9,
    TETRA = 10,
    HEXAHEDRON = 12,
    WEDGE = 13,
    PYRAMID = 14
};

static std::unordered_map<std::type_index, VTKElementName> hpGEMToVTK = {
    {std::type_index(typeid(Geometry::ReferencePoint)), VTKElementName::VERTEX},
    {std::type_index(typeid(Geometry::ReferenceLine)), VTKElementName::LINE},
    {std::type_index(typeid(Geometry::ReferenceTriangle)),
     VTKElementName::TRIANGLE},
    {std::type_index(typeid(Geometry::ReferenceSquare)), VTKElementName::QUAD},
    {std::type_index(typeid(Geometry::ReferenceTetrahedron)),
     VTKElementName::TETRA},
    {std::type_index(typeid(Geometry::ReferenceCube)),
     VTKElementName::HEXAHEDRON},
    {std::type_index(typeid(Geometry::ReferenceTriangularPrism)),
     VTKElementName::WEDGE},
    {std::type_index(typeid(Geometry::ReferencePyramid)),
     VTKElementName::PYRAMID}};

template <std::size_t DIM>
VTKSpecificTimeWriter<DIM>::VTKSpecificTimeWriter(
    const std::string& baseName, const Base::MeshManipulator<DIM>* mesh,
    std::size_t timelevel)
    : totalPoints_(0), totalElements_(0), mesh_(mesh), timelevel_(timelevel) {
    logger.assert_debug(mesh != nullptr, "Invalid mesh passed");
    std::size_t id = Base::MPIContainer::Instance().getProcessorID();
    if (id == 0) {
        writeMasterFileHeader(baseName);
    }
    writeLocalFileHeader(baseName);
}

template <std::size_t DIM>
VTKSpecificTimeWriter<DIM>::~VTKSpecificTimeWriter() {
    std::size_t id = Base::MPIContainer::Instance().getProcessorID();
    if (id == 0) {
        /// Close point data
        masterFile_ << "    </PPointData>" << std::endl;

        /// Describe data format for the points ///
        ///////////////////////////////////////////
        masterFile_ << "    <PPoints>" << std::endl;
        ///\bug assumes all compilers map double to the 64 bit IEEE-754 floating
        /// point data type
        masterFile_
            << "      <PDataArray type=\"Float64\" NumberOfComponents=\"3\"/>"
            << std::endl;
        masterFile_ << "    </PPoints>" << std::endl;
        /// Closure information ///
        ///////////////////////////
        masterFile_ << "  </PUnstructuredGrid>" << std::endl;
        masterFile_ << "</VTKFile>" << std::endl;
        masterFile_.flush();
        masterFile_.close();
    }

    // Close all opened tags in the local file
    localFile_ << "      </PointData>" << std::endl;
    localFile_ << "    </Piece>" << std::endl;
    localFile_ << "  </UnstructuredGrid>" << std::endl;
    localFile_ << "</VTKFile>" << std::endl;
    localFile_.flush();
    localFile_.close();
}

template <std::size_t DIM>
void VTKSpecificTimeWriter<DIM>::write(
    std::function<double(Base::Element*, const Geometry::PointReference<DIM>&,
                         std::size_t)>
        dataCompute,
    const std::string& name) {
    std::size_t id = Base::MPIContainer::Instance().getProcessorID();
    if (id == 0) {
        masterFile_ << "      <PDataArray type=\"Float64\" Name=\"" << name
                    << "\"/>" << std::endl;
    }
    localFile_ << "      <DataArray type=\"Float64\" Name=\"" << name
               << "\" format=\"binary\">" << std::endl;
    std::vector<double> data;
    data.reserve(totalPoints_);
    for (Base::Element* element : mesh_->getElementsList()) {
        for (std::size_t i = 0; i < element->getNumberOfNodes(); ++i) {
            const Geometry::PointReference<DIM>& node =
                element->getReferenceGeometry()->getReferenceNodeCoordinate(
                    tohpGEMOrdering(i, element->getReferenceGeometry()));
            data.push_back(dataCompute(element, node, timelevel_));
        }
    }
    std::uint32_t totalData = sizeof(double) * data.size();
    localFile_ << "        "
               << Detail::toBase64((void*)&totalData, sizeof(totalPoints_));
    if (totalData > 0)
        localFile_ << Detail::toBase64((void*)data.data(), totalData)
                   << std::endl;
    localFile_ << "      </DataArray>" << std::endl;
}

template <std::size_t DIM>
void VTKSpecificTimeWriter<DIM>::write(
    std::function<LinearAlgebra::SmallVector<DIM>(
        Base::Element*, const Geometry::PointReference<DIM>&, std::size_t)>
        dataCompute,
    const std::string& name) {
    std::size_t id = Base::MPIContainer::Instance().getProcessorID();
    if (id == 0) {
        masterFile_ << "      <PDataArray type=\"Float64\" Name=\"" << name
                    << "\" NumberOfComponents=\"3\"/>" << std::endl;
    }
    localFile_ << "      <DataArray type=\"Float64\" Name=\"" << name
               << "\" NumberOfComponents=\"3\" format=\"binary\">" << std::endl;
    std::vector<double> data;
    LinearAlgebra::SmallVector<DIM> newData;
    data.reserve(3 * totalPoints_);
    for (Base::Element* element : mesh_->getElementsList()) {
        for (std::size_t i = 0; i < element->getNumberOfNodes(); ++i) {
            const Geometry::PointReference<DIM>& node =
                element->getReferenceGeometry()->getReferenceNodeCoordinate(i);
            newData = dataCompute(element, node, timelevel_);
            for (std::size_t j = 0; j < newData.size(); ++j) {
                data.push_back(newData[j]);
            }
            for (std::size_t j = newData.size(); j < 3; ++j) {
                data.push_back(0.);
            }
        }
    }
    std::uint32_t totalData = sizeof(double) * data.size();
    localFile_ << "        "
               << Detail::toBase64((void*)&totalData, sizeof(totalPoints_))
               << Detail::toBase64((void*)data.data(), totalData) << std::endl;
    localFile_ << "      </DataArray>" << std::endl;
}

template <std::size_t DIM>
void VTKSpecificTimeWriter<DIM>::write(
    std::function<LinearAlgebra::SmallMatrix<DIM, DIM>(
        Base::Element*, const Geometry::PointReference<DIM>&, std::size_t)>
        dataCompute,
    const std::string& name) {
    std::size_t id = Base::MPIContainer::Instance().getProcessorID();
    if (id == 0) {
        masterFile_ << "      <PDataArray type=\"Float64\" Name=\"" << name
                    << "\" NumberOfComponents=\"3\"/>" << std::endl;
    }
    localFile_ << "      <DataArray type=\"Float64\" Name=\"" << name
               << "\" NumberOfComponents=\"3\" format=\"binary\">" << std::endl;
    std::vector<double> data;
    LinearAlgebra::SmallMatrix<DIM, DIM> newData;
    data.reserve(9 * totalPoints_);
    for (Base::Element* element : mesh_->getElementsList()) {
        for (std::size_t i = 0; i < element->getNumberOfNodes(); ++i) {
            const Geometry::PointReference<DIM>& node =
                element->getReferenceGeometry()->getReferenceNodeCoordinate(i);
            newData = dataCompute(element, node, timelevel_);
            std::size_t j = 0;
            for (; j < newData.getNumberOfRows(); ++j) {
                for (std::size_t k = 0; k < newData.getNumberOfColumns(); ++k) {
                    data.push_back(newData(j, k));
                }
                for (std::size_t k = newData.getNumberOfColumns(); k < 3; ++k) {
                    data.push_back(0.);
                }
            }
            j *= newData.getNumberOfColumns();
            for (; j < 9; ++j) {  // identity matrix
                data.push_back(((j == 0 || j == 4 || j == 8) ? 1. : 0.));
            }
        }
    }
    std::uint32_t totalData = sizeof(double) * data.size();
    localFile_ << "        "
               << Detail::toBase64((void*)&totalData, sizeof(totalPoints_))
               << Detail::toBase64((void*)data.data(), totalData) << std::endl;
    localFile_ << "      </DataArray>" << std::endl;
}

template <std::size_t DIM>
void VTKSpecificTimeWriter<DIM>::writeMasterFileHeader(
    const std::string& baseName) {
    masterFile_.open(baseName + ".pvtu");
    if (!masterFile_.good()) {
        if (baseName.find('/') != std::string::npos) {
            logger(FATAL,
                   "failed to open main paraview output file %.pvtu, does "
                   "the directory % exist?",
                   baseName,
                   baseName.substr(0, baseName.find_last_of('/') + 1));
        } else {
            logger(FATAL, "failed to open main paraview output file %.pvtu",
                   baseName);
        }
    }
    /// Basic header ///
    ////////////////////

    masterFile_ << "<?xml version=\"1.0\"?>" << std::endl;
    masterFile_ << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" "
                   "byte_order=\""
                << (Detail::isBigEndian() ? "BigEndian" : "LittleEndian")
                << "\">" << std::endl;
    masterFile_ << "  <PUnstructuredGrid GhostLevel=\"0\">" << std::endl;
    std::size_t numberOfProcs =
        Base::MPIContainer::Instance().getNumberOfProcessors();

    /// Serial file locations ///
    /////////////////////////////

    for (std::size_t i = 0; i < numberOfProcs; ++i) {
        std::string fileName = baseName;
        if (fileName.find('/') != std::string::npos) {
            fileName = fileName.substr(fileName.find_last_of('/') + 1);
        }
        masterFile_ << "    <Piece Source=\"" << fileName << "." << i
                    << ".vtu\"/>" << std::endl;
    }
    /// Start of the point data ///
    ///////////////////////////////

    // This starts the point data that will be added later
    // tag will be closed in the destructor
    masterFile_ << "    <PPointData>" << std::endl;
}

template <std::size_t DIM>
void VTKSpecificTimeWriter<DIM>::writeLocalFileHeader(
    const std::string& baseName) {
    std::size_t id = Base::MPIContainer::Instance().getProcessorID();
    using namespace std::string_literals;
    localFile_.open(baseName + "."s + std::to_string(id) + ".vtu");
    if (!localFile_.good()) {
        logger(ERROR,
               "failed to open local paraview output file %.vtu, part of the "
               "output will not be written",
               baseName);
    }

    /// Basic Header ///
    ////////////////////
    localFile_ << "<?xml version=\"1.0\"?>" << std::endl;
    std::string endianness =
        Detail::isBigEndian() ? "BigEndian" : "LittleEndian";
    localFile_
        << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\""
        << endianness << "\">" << std::endl;
    localFile_ << "  <UnstructuredGrid>" << std::endl;

    /// Prepare data ///
    ////////////////////

    // first pass compute sizes
    for (Base::Element* element : mesh_->getElementsList()) {
        totalPoints_ += element->getNumberOfNodes();
        ++totalElements_;
    }
    // Second pass
    // We need 4 DataArray's to describe the mesh
    // - array with all the coordinates of the points
    // - array with the VTK-cell-types (line, triangle, square, etc.)
    // - array with connectivity, for each cell giving the indices in the point
    //   data array for the corner points of the cell (in VTK defined order)
    // - array with offsets in the connectivity array, the i-th is the
    //   cummulative number of points of the first i elements. (Alternatively it
    //   is the index of the first entry after the i-th element in the
    //   connectivity array).

    // Point coordinates, VTK needs exactly three coordinates per point
    std::vector<double> pointCoordinates(totalPoints_ * 3, 0.0);
    // Types
    std::vector<VTKElementName> types(totalElements_);
    // Connectivity
    std::vector<std::uint32_t> connectivity(totalPoints_);
    // Offsets
    std::vector<std::uint32_t> offsets(totalElements_);
    // Number of points in previous elements
    std::uint32_t pointCount(0);
    std::size_t elementId = 0;

    for (Base::Element* element : mesh_->getElementsList()) {
        std::uint32_t localNumberOfPoints = element->getNumberOfNodes();
        const Geometry::ReferenceGeometry& referenceGeometry =
            *element->getReferenceGeometry();

        types[elementId] =
            hpGEMToVTK.at(std::type_index(typeid(referenceGeometry)));
        // Individual nodes
        for (std::size_t i = 0; i < localNumberOfPoints; ++i) {
            // Use the VTK ordering for the points. This makes the connectivity
            // trivial, but we need to reorder the coordinates
            connectivity[pointCount + i] = pointCount + i;
            Geometry::PointPhysical<DIM> coordinate =
                element->getPhysicalGeometry()->getLocalNodeCoordinates(
                    tohpGEMOrdering(i, element->getReferenceGeometry()));
            // Copy coordinate
            for (std::size_t j = 0; j < std::min(DIM, 3ul); ++j) {
                pointCoordinates[3 * (pointCount + i) + j] = coordinate[j];
            }
        }

        // Offset array is including the current element
        pointCount += localNumberOfPoints;
        offsets[elementId] = pointCount;
        elementId++;
    }

    /// Write the local part of the mesh ///
    ////////////////////////////////////////

    // clang-format off
    // Manual formatting to show the XML structure
    localFile_ << "    <Piece NumberOfPoints=\"" << totalPoints_ << "\" NumberOfCells=\"" << totalElements_ << "\">" << std::endl;
    localFile_ << "      <Points>" << std::endl;
    localFile_ << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"binary\">" << std::endl;
    localFile_ << "          ";
    writeBinaryDataArrayData(pointCoordinates);
    localFile_ << "        </DataArray>" << std::endl;
    localFile_ << "      </Points>" << std::endl;
    localFile_ << "      <Cells>" << std::endl;
    localFile_ << "        <DataArray type=\"UInt32\" Name=\"connectivity\" format=\"binary\">" << std::endl;
    localFile_ << "          ";
    writeBinaryDataArrayData(connectivity);
    localFile_ << "        </DataArray>" << std::endl;

    localFile_ << "        <DataArray type=\"UInt32\" Name=\"offsets\" format=\"binary\">" << std::endl;
    localFile_ << "          ";
    writeBinaryDataArrayData(offsets);
    localFile_ << "        </DataArray>" << std::endl;
    localFile_ << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"binary\">" << std::endl;
    localFile_ << "          ";
    writeBinaryDataArrayData(types);
    localFile_ << "        </DataArray>" << std::endl;
    localFile_ << "      </Cells>" << std::endl;
    // clang-format on

    /// Prepare for following point data ///
    ////////////////////////////////////////
    localFile_ << "      <PointData>" << std::endl;
}

template <std::size_t DIM>
template <class T>
void VTKSpecificTimeWriter<DIM>::writeBinaryDataArrayData(std::vector<T> data) {
    // Note: according to https://vtk.org/Wiki/VTK_XML_Formats we need a header
    // with the number of bytes. We use the default of using a uint32 as header
    // type.
    std::uint32_t size = data.size() * sizeof(T);
    localFile_ << Detail::toBase64(static_cast<void*>(&size), sizeof(size));
    localFile_ << Detail::toBase64(static_cast<void*>(data.data()), size);
    // Optional newline, but for a nice layout
    localFile_ << std::endl;
}

}  // namespace Output
}  // namespace hpgem