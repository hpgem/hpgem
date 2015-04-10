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

#ifndef BASE_HPP
#define BASE_HPP

#include <vector>

#include "MeshManipulator.h"
#include "GlobalNamespaceBase.h"
#include "CommandLineOptions.h"
#include "GlobalData.h"
#include "RectangularMeshDescriptor.h"

namespace Base
{
    class MeshMoverBase;
    struct RectangularMeshDescriptor;
    struct GlobalData;
    struct ConfigurationData;
    
    /// Basic interface that can create meshes. This class can be used to construct more advanced interfaces.
    class HpgemAPIBase
    {
    public:
        using ConstElementIterator = MeshManipulator::ConstElementIterator;
        using ElementIterator = MeshManipulator::ElementIterator;
        using ConstFaceIterator = MeshManipulator::ConstFaceIterator;
        using FaceIterator = MeshManipulator::FaceIterator;
        using VectorOfMeshManipulatorT = std::vector<MeshManipulator*>;
        using PointPhysicalT = Geometry::PointPhysical;
        using PointReferenceT = Geometry::PointReference;

        using MeshId = std::size_t;
        using VectorOfUIntegers = std::vector<std::size_t>;
        using String = std::string;
        


    public:
        
        HpgemAPIBase(GlobalData * const global, const ConfigurationData* config);

        /// \brief Destructor, destructs the meshes, configData_ and globalData_
        virtual ~HpgemAPIBase();
        
        //If you want a copy-constructor, please make sure to make deep copies of
        //meshes_, configData_ and globalData_.
        HpgemAPIBase(const HpgemAPIBase &other) = delete;
        HpgemAPIBase& operator=(const HpgemAPIBase &other) = delete;

        /// \brief Gives the pointer of meshMoverBase class to mesh.
        virtual bool initialiseMeshMover(const MeshMoverBase* meshMoverBase, std::size_t meshID);

        /// Creating a mesh with in-house remesher.
        MeshId addMesh(const RectangularMeshDescriptor& meshDescriptor, const MeshType& meshType = MeshType::RECTANGULAR, std::size_t nrOfElementMatrixes = 0, std::size_t nrOfElementVectors = 0, std::size_t nrOfFaceMatrixes = 0, std::size_t nrOfFaceVectors = 0);
        
        /// Reading a mesh from a file, currently only Centaur is supported.
        MeshId addMesh(const String& fileName, std::size_t nrOfElementMatrixes = 0, std::size_t nrOfElementVectors = 0, std::size_t nrOfFaceMatrixes = 0, std::size_t nrOfFaceVectors = 0);

        std::size_t getNumberOfElements(MeshId id) const
        {
            return meshes_[id]->getNumberOfElements();
        }
        
        ConstElementIterator elementColBegin(MeshId mId = 0) const;
        ConstElementIterator elementColEnd(MeshId mId = 0) const;

        ElementIterator elementColBegin(MeshId mId = 0);
        ElementIterator elementColEnd(MeshId mId = 0);

        ConstFaceIterator faceColBegin(MeshId mId = 0) const;
        ConstFaceIterator faceColEnd(MeshId mId = 0) const;

        FaceIterator faceColBegin(MeshId mId = 0);
        FaceIterator faceColEnd(MeshId mId = 0);
        
    protected:
        VectorOfMeshManipulatorT meshes_;

        GlobalData * const globalData_;
        const ConfigurationData * const configData_;
    };
}
#endif

