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

#include "HpgemAPIBase.h"

#include "Base/Element.h"

#include "RectangularMeshDescriptor.h"
#include "Geometry/PointPhysical.h"
#include "LinearAlgebra/NumericalVector.h"
#include "ElementCacheData.h"
#include "FaceCacheData.h"
#include "ConfigurationData.h"

namespace Base
{
    HpgemAPIBase::HpgemAPIBase(GlobalData * const global, const ConfigurationData* config)
            : meshes_(), globalData_(global), configData_(config)
    {
        if (!parse_isDone())
        {
            std::cerr << "Warning: Command line arguments have not been parsed.\n" << "  Please call Base::parse_options(argc, argv); first.\n" << "  This application may not behave as intended." << std::endl;
        }
    }
    
    //Destructor, destructs the meshes, configData_ and globalData_
    HpgemAPIBase::~HpgemAPIBase()
    {
        for (std::size_t i = 0; i < meshes_.size(); ++i)
            delete meshes_[i];
        delete configData_;
        delete globalData_;
    }
    
    bool HpgemAPIBase::initialiseMeshMover(const MeshMoverBaseT* meshMoverBase, std::size_t meshID)
    {
        meshes_[meshID]->setMeshMover(meshMoverBase);
        return true;
    }
    
    HpgemAPIBase::MeshId HpgemAPIBase::addMesh(const RectangularMeshDescriptorT& meshDscr, const MeshType& meshType, std::size_t nrOfElementMatrixes, std::size_t nrOfElementVectors, std::size_t nrOfFaceMatrixes, std::size_t nrOfFaceVectors)
    {
        std::size_t numOfMeshes = meshes_.size();
        MeshManipulator* mesh = new MeshManipulator(configData_, meshDscr.boundaryConditions_[0], (configData_->dimension_ > 1) ? meshDscr.boundaryConditions_[1] : false, (configData_->dimension_ > 2) ? meshDscr.boundaryConditions_[2] : false, configData_->polynomialOrder_, 0, nrOfElementMatrixes, nrOfElementVectors, nrOfFaceMatrixes, nrOfFaceVectors);
        
        if (meshType == RECTANGULAR)
        {
            mesh->createRectangularMesh(meshDscr.bottomLeft_, meshDscr.topRight_, meshDscr.numElementsInDIM_);
            mesh->getElementsList();
            meshes_.push_back(mesh);
        }
        else if (meshType == TRIANGULAR)
        {
            mesh->createTriangularMesh(meshDscr.bottomLeft_, meshDscr.topRight_, meshDscr.numElementsInDIM_);
            mesh->getElementsList();
            meshes_.push_back(mesh);
        }
        else
        {
            std::cerr << "Other types are yet to be implemented! " << std::endl;
        }
        std::cout << "I just created a mesh!!!" << std::endl;
        return numOfMeshes;
    }
    
    HpgemAPIBase::MeshId HpgemAPIBase::addMesh(const HpgemAPIBase::String& fileName, std::size_t nrOfElementMatrixes, std::size_t nrOfElementVectors, std::size_t nrOfFaceMatrixes, std::size_t nrOfFaceVectors)
    {
        std::size_t numOfMeshes = meshes_.size();
        MeshManipulator* mesh = new MeshManipulator(configData_, false, false, false, configData_->polynomialOrder_, 0, nrOfElementMatrixes, nrOfElementVectors, nrOfFaceMatrixes, nrOfFaceVectors);
        mesh->readCentaurMesh(fileName); //boundary information (^) is ignored
        mesh->getElementsList();
        meshes_.push_back(mesh);
        std::cout << "I just read a mesh!!!" << std::endl;
        return numOfMeshes;
    }
    
    HpgemAPIBase::ConstElementIterator HpgemAPIBase::elementColBegin(MeshId mId) const
    {
        return meshes_[mId]->elementColBegin();
    }
    
    HpgemAPIBase::ConstElementIterator HpgemAPIBase::elementColEnd(MeshId mId) const
    {
        return meshes_[mId]->elementColEnd();
    }
    
    HpgemAPIBase::ElementIterator HpgemAPIBase::elementColBegin(MeshId mId)
    {
        return meshes_[mId]->elementColBegin();
    }
    
    HpgemAPIBase::ElementIterator HpgemAPIBase::elementColEnd(MeshId mId)
    {
        return meshes_[mId]->elementColEnd();
    }
    
    HpgemAPIBase::ConstFaceIterator HpgemAPIBase::faceColBegin(MeshId mId) const
    {
        return meshes_[mId]->faceColBegin();
    }
    
    HpgemAPIBase::ConstFaceIterator HpgemAPIBase::faceColEnd(MeshId mId) const
    {
        return meshes_[mId]->faceColEnd();
    }
    
    HpgemAPIBase::FaceIterator HpgemAPIBase::faceColBegin(MeshId mId)
    {
        return meshes_[mId]->faceColBegin();
    }
    
    HpgemAPIBase::FaceIterator HpgemAPIBase::faceColEnd(MeshId mId)
    {
        return meshes_[mId]->faceColEnd();
    }

}

