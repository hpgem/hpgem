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
#include "LinearAlgebra/MiddleSizeVector.h"
#include "ElementCacheData.h"
#include "FaceCacheData.h"
#include "ConfigurationData.h"

namespace Base
{
    template<std::size_t DIM>
    HpgemAPIBase<DIM>::HpgemAPIBase(GlobalData * const global, const ConfigurationData* config)
            : meshes_(), globalData_(global), configData_(config)
    {
        if (!parse_isDone())
        {
            logger(WARN, "Warning: Command line arguments have not been parsed.\n" 
                    "  Please call Base::parse_options(argc, argv); first.\n"
                    "  This application may not behave as intended.");
        }
    }
    
    //Destructor, destructs the meshes, configData_ and globalData_
    template<std::size_t DIM>
    HpgemAPIBase<DIM>::~HpgemAPIBase()
    {
        for (std::size_t i = 0; i < meshes_.size(); ++i)
            delete meshes_[i];
        delete configData_;
        delete globalData_;
    }

    template<std::size_t DIM>
    bool HpgemAPIBase<DIM>::initialiseMeshMover(const MeshMoverBase<DIM>* meshMoverBase, std::size_t meshID)
    {
        meshes_[meshID]->setMeshMover(meshMoverBase);
        return true;
    }

    template<std::size_t DIM>
    typename HpgemAPIBase<DIM>::MeshId HpgemAPIBase<DIM>::addMesh(const RectangularMeshDescriptor<DIM>& meshDscr, const MeshType& meshType, std::size_t nrOfElementMatrixes, std::size_t nrOfElementVectors, std::size_t nrOfFaceMatrixes, std::size_t nrOfFaceVectors)
    {
        std::size_t numOfMeshes = meshes_.size();
        MeshManipulator<DIM>* mesh = new MeshManipulator<DIM>(configData_, meshDscr.boundaryConditions_[0],
                                                    (configData_->dimension_ > 1) ? meshDscr.boundaryConditions_[1] : BoundaryType::SOLID_WALL, (configData_->dimension_ > 2) ? meshDscr.boundaryConditions_[2] : BoundaryType::SOLID_WALL, configData_->polynomialOrder_, 0, nrOfElementMatrixes, nrOfElementVectors, nrOfFaceMatrixes, nrOfFaceVectors);
        
        if (meshType == MeshType::RECTANGULAR)
        {
            mesh->createRectangularMesh(meshDscr.bottomLeft_, meshDscr.topRight_, meshDscr.numElementsInDIM_);
            mesh->getElementsList();
            meshes_.push_back(mesh);
        }
        else if (meshType == MeshType::TRIANGULAR)
        {
            mesh->createTriangularMesh(meshDscr.bottomLeft_, meshDscr.topRight_, meshDscr.numElementsInDIM_);
            mesh->getElementsList();
            meshes_.push_back(mesh);
        }
        /*else
        {
            logger(ERROR, "The only mesh types that are implemented are RECTANGULAR and TRIANGULAR. % is not implemented.", meshType);
        }*/
        logger(INFO, "HpgemAPIBase::addMesh created a mesh.");
        return numOfMeshes;
    }

    template<std::size_t DIM>
    typename HpgemAPIBase<DIM>::MeshId HpgemAPIBase<DIM>::addMesh(const std::string& fileName, std::size_t nrOfElementMatrixes, std::size_t nrOfElementVectors, std::size_t nrOfFaceMatrixes, std::size_t nrOfFaceVectors)
    {
        std::size_t numOfMeshes = meshes_.size();
        MeshManipulator<DIM>* mesh = new MeshManipulator<DIM>(configData_, BoundaryType::SOLID_WALL, BoundaryType::SOLID_WALL, BoundaryType::SOLID_WALL, configData_->polynomialOrder_, 0, nrOfElementMatrixes, nrOfElementVectors, nrOfFaceMatrixes, nrOfFaceVectors);
        mesh->readCentaurMesh(fileName);                             //boundary information (^) is ignored
        mesh->getElementsList();
        meshes_.push_back(mesh);
        logger(INFO, "HpgemAPIBase::addMesh read a mesh.");
        return numOfMeshes;
    }

    template<std::size_t DIM>
    void HpgemAPIBase<DIM>::synchronize(const std::size_t timeLevel)
    {
#ifdef HPGEM_USE_MPI
        //Now, set it up.
        Base::MeshManipulator<DIM> * meshManipulator = this->meshes_[0];
        Base::Submesh& mesh = meshManipulator->getMesh().getSubmesh();

        const auto& pushes = mesh.getPushElements();
        const auto& pulls = mesh.getPullElements();

        //receive first for lower overhead
        for (const auto& it : pulls)
        {
            for (Base::Element *ptrElement : it.second)
            {
                logger.assert(ptrElement->getTimeLevelDataVector(timeLevel).size() == ptrElement->getNrOfBasisFunctions() * this->configData_->numberOfUnknowns_ , "Size of time level % data vector is wrong: % instead of %.", timeLevel, ptrElement->getTimeLevelDataVector(timeLevel).size(), this->configData_->numberOfBasisFunctions_ * this->configData_->numberOfUnknowns_);

                Base::MPIContainer::Instance().receive(ptrElement->getTimeLevelDataVector(timeLevel), it.first, ptrElement->getID());
            }
        }
        for (const auto& it : pushes)
        {
            for (Base::Element *ptrElement : it.second)
            {
                logger.assert(ptrElement->getTimeLevelDataVector(timeLevel).size() == ptrElement->getNrOfBasisFunctions() * this->configData_->numberOfUnknowns_, "Size of time level % data vector is wrong: % instead of %.", timeLevel, ptrElement->getTimeLevelDataVector(timeLevel).size(), this->configData_->numberOfBasisFunctions_ * this->configData_->numberOfUnknowns_);

                Base::MPIContainer::Instance().send(ptrElement->getTimeLevelDataVector(timeLevel), it.first, ptrElement->getID());
            }
        }
        Base::MPIContainer::Instance().sync();
#endif
    }

    template<std::size_t DIM>
    typename HpgemAPIBase<DIM>::ConstElementIterator HpgemAPIBase<DIM>::elementColBegin(MeshId mId) const
    {
        return meshes_[mId]->elementColBegin();
    }

    template<std::size_t DIM>
    typename HpgemAPIBase<DIM>::ConstElementIterator HpgemAPIBase<DIM>::elementColEnd(MeshId mId) const
    {
        return meshes_[mId]->elementColEnd();
    }

    template<std::size_t DIM>
    typename HpgemAPIBase<DIM>::ElementIterator HpgemAPIBase<DIM>::elementColBegin(MeshId mId)
    {
        return meshes_[mId]->elementColBegin();
    }

    template<std::size_t DIM>
    typename HpgemAPIBase<DIM>::ElementIterator HpgemAPIBase<DIM>::elementColEnd(MeshId mId)
    {
        return meshes_[mId]->elementColEnd();
    }

    template<std::size_t DIM>
    typename HpgemAPIBase<DIM>::ConstFaceIterator HpgemAPIBase<DIM>::faceColBegin(MeshId mId) const
    {
        return meshes_[mId]->faceColBegin();
    }

    template<std::size_t DIM>
    typename HpgemAPIBase<DIM>::ConstFaceIterator HpgemAPIBase<DIM>::faceColEnd(MeshId mId) const
    {
        return meshes_[mId]->faceColEnd();
    }

    template<std::size_t DIM>
    typename HpgemAPIBase<DIM>::FaceIterator HpgemAPIBase<DIM>::faceColBegin(MeshId mId)
    {
        return meshes_[mId]->faceColBegin();
    }

    template<std::size_t DIM>
    typename HpgemAPIBase<DIM>::FaceIterator HpgemAPIBase<DIM>::faceColEnd(MeshId mId)
    {
        return meshes_[mId]->faceColEnd();
    }

}

