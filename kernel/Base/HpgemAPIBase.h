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
    template<std::size_t DIM>
    class MeshMoverBase;
    template<std::size_t DIM>
    struct RectangularMeshDescriptor;
    struct GlobalData;
    struct ConfigurationData;


    /// Basic interface that can create meshes. This class can be used to construct more advanced interfaces.
    template<std::size_t DIM>
    class HpgemAPIBase
    {
    public:
        using ConstElementIterator = typename MeshManipulator<DIM>::ConstElementIterator;
        using ElementIterator = typename MeshManipulator<DIM>::ElementIterator;
        using ConstFaceIterator = typename MeshManipulator<DIM>::ConstFaceIterator;
        using FaceIterator = typename MeshManipulator<DIM>::FaceIterator;
        using PointPhysicalT = Geometry::PointPhysical<DIM>;
        using PointReferenceT = Geometry::PointReference<DIM>;
        using PointReferenceOnFaceT = Geometry::PointReference<DIM - 1>;
        
        HpgemAPIBase(GlobalData * const global, const ConfigurationData* config);

        /// \brief Destructor, destructs the meshes, configData_ and globalData_
        virtual ~HpgemAPIBase();
        
        //If you want a copy-constructor, please make sure to make deep copies of
        //meshes_, configData_ and globalData_.
        HpgemAPIBase(const HpgemAPIBase &other) = delete;
        HpgemAPIBase& operator=(const HpgemAPIBase &other) = delete;

        /// \brief Gives the pointer of meshMoverBase class to mesh.
        virtual bool initialiseMeshMover(const MeshMoverBase<DIM>* meshMoverBase, std::size_t meshID);
        
        /// Reading a mesh from a file, currently only Centaur is supported.
        std::size_t addMesh(const std::string& fileName, std::size_t numberOfElementMatrixes = 0, std::size_t numberOfElementVectors = 0, std::size_t numberOfFaceMatrixes = 0, std::size_t numberOfFaceVectors = 0);

        /// \brief Synchronize between the different submeshes (when using MPI)
        /// You should call this function after you update a timeIntegrationVector, but before
        /// you need information from neighboring elements, including when you integrate over a Face to compute a flux
        virtual void synchronize(const std::size_t timeIntegrationVectorId);

        /// \brief Set the number of time integration vectors for every element.
        void setNumberOfTimeIntegrationVectorsGlobally(std::size_t numberOfTimeIntegrationVectors);
        
    
        /// \brief Copy the data of a time integration vector to the data of a certain time level.
        void copyTimeIntegrationToTimeLevelData(std::size_t timeIntegrationVectorId, std::size_t timeLevel, std::size_t meshId = 0);
        
        /// \brief Copy the data of a time level to a time integration vector.
        void copyTimeLevelToTimeIntegrationData(std::size_t timeLevel, std::size_t timeIntegrationVectorId, std::size_t meshId = 0);
        
        
        std::size_t getNumberOfElements(std::size_t id) const
        {
            return meshes_[id]->getNumberOfElements();
        }
        
        ConstElementIterator elementColBegin(std::size_t mId = 0) const;
        ConstElementIterator elementColEnd(std::size_t mId = 0) const;

        ElementIterator elementColBegin(std::size_t mId = 0);
        ElementIterator elementColEnd(std::size_t mId = 0);

        ConstFaceIterator faceColBegin(std::size_t mId = 0) const;
        ConstFaceIterator faceColEnd(std::size_t mId = 0) const;

        FaceIterator faceColBegin(std::size_t mId = 0);
        FaceIterator faceColEnd(std::size_t mId = 0);
        
    protected:
        std::vector<MeshManipulator<DIM>*> meshes_;

        GlobalData * const globalData_;
        const ConfigurationData * const configData_;
        
        /// \brief Number of time integration vectors for every element
        /// \details Use this constant if the number of time integration vectors is the same for all elements. This is the case for many standard time integation methods. However, for some methods, like local time-stepping methods or hybrid time integration methods, the required number of time integration vectors can differ per element.
        std::size_t globalNumberOfTimeIntegrationVectors_;
    };
}

#include "HpgemAPIBase_Impl.h"

#endif

