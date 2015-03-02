/*
 This file forms part of hpGEM. This package has been developed over a number of years by various people at the University of Twente and a full list of contributors can be found at
 http://hpgem.org/about-the-code/team
 
 This code is distributed using BSD 3-Clause License. A copy of which can found below.
 
 
 Copyright (c) 2014, University of Twente
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "Base/Logger.h"
#include "VTKSpecificTimeWriter.hpp"
#include "Base/MpiContainer.hpp"

#include "Base/MeshManipulator.hpp"
#include "Geometry/ElementGeometry.hpp"
#include "Base/Element.hpp"
#include "Geometry/PointPhysical.hpp"
#include "Geometry/PhysicalGeometry.hpp"
#include "Geometry/ReferencePoint.hpp"
#include "Geometry/ReferenceLine.hpp"
#include "Geometry/ReferenceTriangle.hpp"
#include "Geometry/ReferenceSquare.hpp"
#include "Geometry/ReferenceTetrahedron.hpp"
#include "Geometry/ReferenceCube.hpp"
#include "Geometry/ReferenceTriangularPrism.hpp"
#include "Geometry/ReferencePyramid.hpp"
#include "Geometry/PointReference.hpp"
#include "Base/ElementCacheData.hpp"
#include "base64.hpp"
#include "VTKElementOrdering.hpp"
#include <vector>
#include <unordered_map>
#include <typeindex>

static Logger<LOG_LEVEL> logger("VTK discontinuous solution writer (single time)");

/////////////////////////////////////
//some VTK specific helper routines//
/////////////////////////////////////

//vtk element types supported by hpGEM
//magic numbers taken from http://www.vtk.org/VTK/img/file-formats.pdf
//smaller underlying type allows for easy conversion to base64
enum class VTKElementName : uint8_t 
{
    VERTEX=1,
    LINE=3,
    TRIANGLE=5,
    QUAD=9,
    TETRA=10,
    HEXAHEDRON=12,
    WEDGE=13,
    PYRAMID=14
};

static std::unordered_map<std::type_index,VTKElementName> hpGEMToVTK = 
{
    {std::type_index(typeid(Geometry::ReferencePoint))          , VTKElementName::VERTEX},
    {std::type_index(typeid(Geometry::ReferenceLine))           , VTKElementName::LINE},
    {std::type_index(typeid(Geometry::ReferenceTriangle))       , VTKElementName::TRIANGLE},
    {std::type_index(typeid(Geometry::ReferenceSquare))         , VTKElementName::QUAD},
    {std::type_index(typeid(Geometry::ReferenceTetrahedron))    , VTKElementName::TETRA},
    {std::type_index(typeid(Geometry::ReferenceCube))           , VTKElementName::HEXAHEDRON},
    {std::type_index(typeid(Geometry::ReferenceTriangularPrism)), VTKElementName::WEDGE},
    {std::type_index(typeid(Geometry::ReferencePyramid))        , VTKElementName::PYRAMID}
};

Output::VTKSpecificTimeWriter::VTKSpecificTimeWriter(const std::string& baseName, const Base::MeshManipulator* mesh, std::size_t timelevel) : mesh_(mesh), totalPoints_(0), timelevel_(timelevel)
{
    std::size_t id = Base::MPIContainer::Instance().getProcessorID();
    std::uint32_t totalData;
    if(id == 0)
    {
        masterFile_.open(baseName + ".pvtu");
        if(!masterFile_.good()){
            logger(FATAL, "failed to open main paraview output file %.pvtu",baseName);
            exit(1);
        }
        masterFile_ << "<?xml version=\"1.0\"?>" << std::endl;
        masterFile_ << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"" << (Detail::isBigEndian()?"BigEndian":"LittleEndian") << "\">" << std::endl;
        masterFile_ << "  <PUnstructuredGrid GhostLevel=\"0\">" << std::endl;
        std::size_t numProcs = Base::MPIContainer::Instance().getNumProcessors();
        for (std::size_t i = 0; i < numProcs; ++i)
        {
            std::string fileName=baseName;
            if(fileName.find('/')!=std::string::npos)
            {
                fileName=fileName.substr(fileName.find_last_of('/')+1);
            }
            masterFile_ << "    <Piece Source=\"" << fileName << i << ".vtu\"/>" <<std::endl;
        }
        masterFile_ << "    <PPointData>" << std::endl;
    }
    localFile_.open(baseName + std::to_string(id) + ".vtu");
    if(!localFile_.good()){
        logger(ERROR, "failed to open local paraview output file %.vtu, part of the output will not be written",baseName);
    }
    localFile_ << "<?xml version=\"1.0\"?>" << std::endl;
    localFile_ << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"" << (Detail::isBigEndian()?"BigEndian":"LittleEndian") << "\">" << std::endl;
    localFile_ << "  <UnstructuredGrid>" << std::endl;
    //the number of points is not an inherent quantity of the mesh, because we have to repeat nodes to allow discontinuous data
    std::uint32_t totalElements{(std::uint32_t)mesh_->getElementsList().size()};
    for(Base::Element* element:mesh_->getElementsList())
    {
        totalPoints_+=element->getNrOfNodes();
    }
    localFile_ << "    <Piece NumberOfPoints=\"" << totalPoints_ << "\" NumberOfCells=\"" << totalElements << "\">" << std::endl;
    localFile_ << "      <Points>" << std::endl;
    localFile_ << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"binary\">" << std::endl << "          ";
    totalData=3*totalPoints_*sizeof(double);
    localFile_ << Detail::toBase64((void*)&totalData,sizeof(totalData));
    Geometry::PointPhysical usefullNode(3);
    const int DIM = mesh->dimension();
    Geometry::PointPhysical actualNode(DIM);
    std::vector<std::uint32_t> cumulativeNodesPerElement;
    cumulativeNodesPerElement.reserve(totalElements+1);
    cumulativeNodesPerElement.push_back(0);
    std::vector<VTKElementName> elementTypes;
    elementTypes.reserve(totalElements);
    for(Base::Element* element:mesh_->getElementsList())
    {
        cumulativeNodesPerElement.push_back(element->getNrOfNodes()+cumulativeNodesPerElement.back());
        elementTypes.push_back(hpGEMToVTK.at(std::type_index(typeid(*element->getReferenceGeometry()))));
        for (std::size_t i = 0; i < element->getNrOfNodes(); ++i)
        {
            element->getPhysicalGeometry()->getLocalNodeCoordinates(tohpGEMOrdering(i,element->getReferenceGeometry()),actualNode);
            for (std::size_t j = 0; j < DIM; ++j)
            {
                usefullNode[j]=actualNode[j];
            }
            //this bit will only work correctly if the data size is a multiple of 3, but VTK requires 3D coordinates anyway
            localFile_ << Detail::toBase64((void*)usefullNode.data(),3*sizeof(double));
        }
    }
    localFile_ << std::endl << "        </DataArray>" << std::endl;
    localFile_ << "      </Points>" << std::endl;
    localFile_ << "      <Cells>" << std::endl;
    localFile_ << "        <DataArray type=\"UInt32\" Name=\"connectivity\" format=\"binary\">" << std::endl << "          ";
    totalData=totalPoints_*sizeof(totalPoints_);
    localFile_ << Detail::toBase64((void*)&totalData,sizeof(totalPoints_));
    std::vector<std::uint32_t> index(totalPoints_);
    for(std::size_t i=0;i<totalPoints_;++i)
    {
        index[i]=i;
        
    }
    localFile_ << Detail::toBase64((void*)index.data(),totalData) << std::endl;
    localFile_ << "        </DataArray>" << std::endl;
    localFile_ << "        <DataArray type=\"UInt32\" Name=\"offsets\" format=\"binary\">" << std::endl;
    totalData=totalElements*sizeof(totalElements);
    localFile_ << "          " << Detail::toBase64((void*)&totalData, sizeof(totalElements)) << Detail::toBase64((void*)(cumulativeNodesPerElement.data()+1), totalData) << std::endl;
    localFile_ << "        </DataArray>" << std::endl;
    localFile_ << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"binary\">" << std::endl;
    totalData = totalElements * sizeof(VTKElementName);
    localFile_ << "          " << Detail::toBase64((void*)&totalData, sizeof(totalElements)) << Detail::toBase64((void*)elementTypes.data(), totalData) << std::endl;
    localFile_ << "        </DataArray>" << std::endl;
    localFile_ << "      </Cells>" << std::endl;
    localFile_ << "      <PointData>" << std::endl;
}

Output::VTKSpecificTimeWriter::~VTKSpecificTimeWriter()
{
    std::size_t id = Base::MPIContainer::Instance().getProcessorID();
    if(id == 0)
    {
        masterFile_ << "    </PPointData>" << std::endl;
        masterFile_ << "    <PPoints>" << std::endl;
        ///\bug assumes all compilers map double to the 64 bit IEEE-754 floating point data type
        masterFile_ << "      <PDataArray type=\"Float64\" NumberOfComponents=\"3\"/>" << std::endl;
        masterFile_ << "    </PPoints>" << std::endl;
        masterFile_ << "  </PUnstructuredGrid>" << std::endl;
        masterFile_ << "</VTKFile>" << std::endl;
        masterFile_.flush();
        masterFile_.close();
    }
    localFile_ << "      </PointData>" << std::endl;
    localFile_ << "    </Piece>" << std::endl;
    localFile_ << "  </UnstructuredGrid>" << std::endl;
    localFile_ << "</VTKFile>" << std::endl;
    localFile_.flush();
    localFile_.close();
}

void Output::VTKSpecificTimeWriter::write(std::function<double(Base::Element*, const Geometry::PointReference&, std::size_t) > dataCompute, const std::string& name)
{
    std::size_t id = Base::MPIContainer::Instance().getProcessorID();
    if(id == 0)
    {
        masterFile_ << "      <PDataArray type=\"Float64\" Name=\"" << name << "\"/>" << std::endl;
    }
    localFile_ << "      <DataArray type=\"Float64\" Name=\"" << name << "\" format=\"binary\">" << std::endl;
    std::vector<double> data;
    data.reserve(totalPoints_);
    Geometry::PointReference node(mesh_->dimension());
    for(Base::Element* element:mesh_->getElementsList())
    {
        for (std::size_t i = 0; i < element->getNrOfNodes(); ++i)
        {
            element->getReferenceGeometry()->getNode(tohpGEMOrdering(i,element->getReferenceGeometry()),node);
            data.push_back(dataCompute(element, node, timelevel_));
        }
    }
    std::uint32_t totalData = sizeof(double)*data.size();
    localFile_ << "        " << Detail::toBase64((void*)&totalData,sizeof(totalPoints_)) << Detail::toBase64((void*)data.data(), totalData) << std::endl;
    localFile_ << "      </DataArray>" << std::endl;
}

void Output::VTKSpecificTimeWriter::write(std::function<LinearAlgebra::NumericalVector(Base::Element*, const Geometry::PointReference&, std::size_t) > dataCompute, const std::string& name)
{
    std::size_t id = Base::MPIContainer::Instance().getProcessorID();
    if(id == 0)
    {
        masterFile_ << "      <PDataArray type=\"Float64\" Name=\"" << name << "\" NumberOfComponents=\"3\"/>" << std::endl;
    }
    localFile_ << "      <DataArray type=\"Float64\" Name=\"" << name << "\" NumberOfComponents=\"3\" format=\"binary\">" << std::endl;
    std::vector<double> data;
    LinearAlgebra::NumericalVector newData;
    data.reserve(3*totalPoints_);
    Geometry::PointReference node(mesh_->dimension());
    for(Base::Element* element:mesh_->getElementsList())
    {
        for (std::size_t i = 0; i < element->getNrOfNodes(); ++i)
        {
            element->getReferenceGeometry()->getNode(i,node);
            newData = dataCompute(element,node, timelevel_);
            for (std::size_t j = 0; j < newData.size(); ++j)
            {
                data.push_back(newData[j]);
            }
            for (std::size_t j = newData.size(); j < 3; ++j)
            {
                data.push_back(0.);
            }
        }
    }
    std::uint32_t totalData = sizeof(double)*data.size();
    localFile_ << "        " << Detail::toBase64((void*)&totalData,sizeof(totalPoints_)) << Detail::toBase64((void*)data.data(), totalData) << std::endl;
    localFile_ << "      </DataArray>" << std::endl;
}

void Output::VTKSpecificTimeWriter::write(std::function<LinearAlgebra::Matrix(Base::Element*, const Geometry::PointReference&, std::size_t) > dataCompute, const std::string& name)
{
    std::size_t id = Base::MPIContainer::Instance().getProcessorID();
    if(id == 0)
    {
        masterFile_ << "      <PDataArray type=\"Float64\" Name=\"" << name << "\" NumberOfComponents=\"3\"/>" << std::endl;
    }
    localFile_ << "      <DataArray type=\"Float64\" Name=\"" << name << "\" NumberOfComponents=\"3\" format=\"binary\">" << std::endl;
    std::vector<double> data;
    LinearAlgebra::Matrix newData;
    data.reserve(3*totalPoints_);
    Geometry::PointReference node(mesh_->dimension());
    for(Base::Element* element:mesh_->getElementsList())
    {
        for (std::size_t i = 0; i < element->getNrOfNodes(); ++i)
        {
            element->getReferenceGeometry()->getNode(i,node);
            newData = dataCompute(element,node, timelevel_);
            std::size_t j = 0;
            for(;j<newData.getNRows();++j)
            {
                for (std::size_t k = 0; k < newData.getNCols(); ++k)
                {
                    data.push_back(newData(j,k));
                }
                for (std::size_t k = newData.getNCols(); k < 3; ++k)
                {
                    data.push_back(0.);
                }
            }
            j*=newData.getNCols();
            for(;j<9;++j)
            {//identity matrix
                data.push_back(((j==0||j==4||j==8)?1.:0.));
            }
        }
    }
    std::uint32_t totalData = sizeof(double)*data.size();
    localFile_ << "        " << Detail::toBase64((void*)&totalData,sizeof(totalPoints_)) << Detail::toBase64((void*)data.data(), totalData) << std::endl;
    localFile_ << "      </DataArray>" << std::endl;
}

