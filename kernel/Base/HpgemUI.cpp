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

#include "HpgemUI.hpp"

#include "Base/Element.hpp"
#include "Integration/ReturnTrait1.hpp"

#include "RectangularMeshDescriptor.hpp"
#include "Geometry/PointPhysical.hpp"
#include "LinearAlgebra/NumericalVector.hpp"
#include "ElementCacheData.hpp"
#include "FaceCacheData.hpp"
#include "ConfigurationData.hpp"

namespace Base
{

    HpgemUI::HpgemUI(GlobalData* const global, const ConfigurationData* config):
        meshes_(),
        globalData_(global),
        configData_(config)
    {
    }
    
    HpgemUI::~HpgemUI()
    {
    	std::cout << "is this called yet?"<<std::endl;
            for(int i = 0; i < meshes_.size() ; ++i)
                delete meshes_[i];
    }
    
    bool
    HpgemUI::initialiseMeshMover(const MeshMoverBaseT* meshMoverBase, unsigned int meshID)
    {
        meshes_[meshID]->setMeshMover(meshMoverBase);
        return true;
    }
    
    
    typename HpgemUI::MeshId
    HpgemUI::addMesh(const RectangularMeshDescriptorT& meshDscr, const MeshType& meshType,int nrOfElementMatrixes, int nrOfElementVectors, int nrOfFaceMatrixes, int nrOfFaceVectors)
    {
        unsigned int numOfMeshes=meshes_.size();
        MeshManipulator* mesh = new MeshManipulator(configData_,meshDscr.boundaryConditions_[0],meshDscr.boundaryConditions_[1],meshDscr.boundaryConditions_[2],configData_->polynomialOrder_,0,nrOfElementMatrixes,nrOfElementVectors,nrOfFaceMatrixes,nrOfFaceVectors);//
        
        if (meshType== RECTANGULAR)
        {   
            mesh->createRectangularMesh(meshDscr.bottomLeft_, meshDscr.topRight_, meshDscr.numElementsInDIM_);
            meshes_.push_back(mesh);
        }
        else if (meshType==TRIANGULAR)
        {
        	mesh->createTriangularMesh(meshDscr.bottomLeft_,meshDscr.topRight_,meshDscr.numElementsInDIM_);
        	meshes_.push_back(mesh);
        }
        else
        {
            std::cerr << "Other types are yet to be implemented! " << std::endl;
        }
        std::cout<<"I just created a mesh!!!"<<std::endl;
            //mesh->outputMesh(std::cout);
        return numOfMeshes;
    }
    
    
    HpgemUI::MeshId HpgemUI::addMesh(const HpgemUI::String& fileName, int nrOfElementMatrixes, int nrOfElementVectors,int nrOfFaceMatrixes, int nrOfFaceVectors)
    {
        unsigned int numOfMeshes=meshes_.size();
        MeshManipulator* mesh = new MeshManipulator(configData_,false,false,false,configData_->polynomialOrder_,0,nrOfElementMatrixes,nrOfElementVectors,nrOfFaceMatrixes,nrOfFaceVectors);
        mesh->readCentaurMesh(fileName);  //boundary information (^) is ignored
        meshes_.push_back(mesh);
        std::cout<<"I just read a mesh!!!"<<std::endl;
            //mesh->outputMesh(std::cout);
        return numOfMeshes;
    }
    
    typename HpgemUI::ConstElementIterator
    HpgemUI::elementColBegin(MeshId mId)const
    {
        return meshes_[mId]->elementColBegin();
    }

    typename HpgemUI::ConstElementIterator
    HpgemUI::elementColEnd(MeshId mId)const
    {
        return meshes_[mId]->elementColEnd();
    }

    typename HpgemUI::ElementIterator
    HpgemUI::elementColBegin(MeshId mId)
    {
        return meshes_[mId]->elementColBegin();
    }

    typename HpgemUI::ElementIterator
    HpgemUI::elementColEnd(MeshId mId)
    {
        return meshes_[mId]->elementColEnd();
    }
    
    typename HpgemUI::ConstFaceIterator
    HpgemUI::faceColBegin(MeshId mId)const
    {
        return meshes_[mId]->faceColBegin();
    }

    typename HpgemUI::ConstFaceIterator
    HpgemUI::faceColEnd(MeshId mId)const
    {
        return meshes_[mId]->faceColEnd();
    }

    typename HpgemUI::FaceIterator
    HpgemUI::faceColBegin(MeshId mId)
    {
        return meshes_[mId]->faceColBegin();
    }

    typename HpgemUI::FaceIterator
    HpgemUI::faceColEnd(MeshId mId)
    {
        return meshes_[mId]->faceColEnd();
    }

}
