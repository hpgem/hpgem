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

#ifndef MESHMANIPULATORBASE_H_
#define MESHMANIPULATORBASE_H_

#include <vector>
#include <fstream>

#include "Geometry/FaceGeometry.h"
#include "FaceCacheData.h"
#include "Mesh.h"
#include "GlobalNamespaceBase.h"
#include "BasisFunctionSet.h"

namespace Base
{
    template<std::size_t DIM>
    class MeshManipulator;
}

template<std::size_t DIM>
std::ostream& operator<<(std::ostream&, const Base::MeshManipulator<DIM>&);
namespace Geometry
{
    template<std::size_t DIM>
    class PointPhysical;
    template<std::size_t DIM>
    class PointReference;
}

namespace Base
{
    class BasisFunctionSet;
    class OrientedBasisFunctionSet;
    class Face;
    template<std::size_t DIM>
    class MeshMoverBase;
    template<class V>
    class LevelTree;
    class Element;
    struct ConfigurationData;
    class Edge;
    
    //class is made final so we don't have to create a v-table specifically for the destructor
    class MeshManipulatorBase
    {
    public:
        
        using PointIndexT = std::size_t;
        using ElementT = Element;
        using FaceT = Face;
        using BasisFunctionSetT = Base::BasisFunctionSet;

        using ElementLevelTreeT = LevelTree<ElementT>;
        using FaceLevelTreeT = LevelTree<FaceT>;

        using ListOfFacesT = std::vector<FaceT*>;
        using ListOfElementsT = std::vector<ElementT*>;
        using VectorOfElementPtrT = std::vector<ElementT* >;
        using VectorOfPointIndicesT = std::vector<PointIndexT>;
        using CollectionOfBasisFunctionSets = std::vector<std::shared_ptr<const BasisFunctionSetT>>;
        using VecOfElementLevelTreePtrT = std::vector<ElementLevelTreeT*>;
        using VecOfFaceLevelTreePtrT = std::vector<FaceLevelTreeT*>;
        
        using ConstElementIterator = ListOfElementsT::const_iterator;
        using ElementIterator = ListOfElementsT::iterator;

        using ConstFaceIterator = ListOfFacesT::const_iterator;
        using FaceIterator = ListOfFacesT::iterator;

    public:
        /// idRangeBegin is the beginning of the range, from where the Element's ids should be assigned.
        /// In case of multiple meshes, one has to take care of empty intersection of those ranges!!!
        MeshManipulatorBase(const ConfigurationData* configData, BoundaryType xPer = BoundaryType::SOLID_WALL, BoundaryType yPer = BoundaryType::SOLID_WALL, BoundaryType zPer = BoundaryType::SOLID_WALL, std::size_t orderOfFEM = 1, std::size_t idRangeBegin = 0, std::size_t nrOfElementMatrixes = 0, std::size_t nrOfElementVectors = 0, std::size_t nrOfFaceMatrixes = 0, std::size_t nrOfFaceVectors = 0);

        MeshManipulatorBase(const MeshManipulatorBase& other);

        virtual ~MeshManipulatorBase() = default;

        virtual std::size_t getNumberOfElements(IteratorType part = IteratorType::LOCAL) const = 0;
        
        virtual std::size_t getNumberOfFaces(IteratorType part = IteratorType::LOCAL) const = 0;
        
        virtual std::size_t getNumberOfEdges(IteratorType part = IteratorType::LOCAL) const = 0;
        
        virtual std::size_t getNumberOfNodes(IteratorType part = IteratorType::LOCAL) const = 0;
        
        /// *****************Iteration through the Elements*******************
        
        virtual ConstElementIterator elementColBegin(IteratorType part = IteratorType::LOCAL) const = 0;
        
        virtual ConstElementIterator elementColEnd(IteratorType part = IteratorType::LOCAL) const = 0;
        
        virtual ElementIterator elementColBegin(IteratorType part = IteratorType::LOCAL) = 0;
        
        virtual ElementIterator elementColEnd(IteratorType part = IteratorType::LOCAL) = 0;
        
        virtual ConstFaceIterator faceColBegin(IteratorType part = IteratorType::LOCAL) const = 0;
        
        virtual ConstFaceIterator faceColEnd(IteratorType part = IteratorType::LOCAL) const = 0;
        
        virtual FaceIterator faceColBegin(IteratorType part = IteratorType::LOCAL) = 0;
        
        virtual FaceIterator faceColEnd(IteratorType part = IteratorType::LOCAL) = 0;
        
        virtual std::vector<Edge*>::const_iterator edgeColBegin(IteratorType part = IteratorType::LOCAL) const = 0;
        
        virtual std::vector<Edge*>::const_iterator edgeColEnd(IteratorType part = IteratorType::LOCAL) const = 0;
        
        virtual std::vector<Edge*>::iterator edgeColBegin(IteratorType part = IteratorType::LOCAL) = 0;
        
        virtual std::vector<Edge*>::iterator edgeColEnd(IteratorType part = IteratorType::LOCAL) = 0;
        
        virtual std::vector<Node*>::const_iterator nodeColBegin(IteratorType part = IteratorType::LOCAL) const = 0;
        
        virtual std::vector<Node*>::const_iterator nodeColEnd(IteratorType part = IteratorType::LOCAL) const = 0;
        
        virtual std::vector<Node*>::iterator nodeColBegin(IteratorType part = IteratorType::LOCAL) = 0;
        
        virtual std::vector<Node*>::iterator nodeColEnd(IteratorType part = IteratorType::LOCAL) = 0;
        //  *****************Iteration through the Elements*******************
        
        //! Get const list of elements
        virtual const ListOfElementsT& getElementsList(IteratorType part = IteratorType::LOCAL) const = 0;
        
        //! Get non-const list of elements
        virtual ListOfElementsT& getElementsList(IteratorType part = IteratorType::LOCAL) = 0;
        
        //! Get const list of faces
        virtual const ListOfFacesT& getFacesList(IteratorType part = IteratorType::LOCAL) const = 0;
        
        //! Get non-const list of faces
        virtual ListOfFacesT& getFacesList(IteratorType part = IteratorType::LOCAL) = 0;
        
        virtual const std::vector<Edge*>& getEdgesList(IteratorType part = IteratorType::LOCAL) const = 0;
        
        virtual std::vector<Edge*>& getEdgesList(IteratorType part = IteratorType::LOCAL) = 0;
        
        virtual const std::vector<Node*>& getNodesList(IteratorType part = IteratorType::LOCAL) const = 0;
        
        virtual std::vector<Node*>& getNodesList(IteratorType part = IteratorType::LOCAL) = 0;
                
        std::size_t dimension() const;
        
    protected:

        const ConfigurationData* configData_;
        //! Periodicity in x-direction.
        bool periodicX_;

        //! Periodicity in y-direction.
        bool periodicY_;

        //! Periodicity in z-direction.
        bool periodicZ_;
        
        std::size_t numberOfElementMatrixes_;
        std::size_t numberOfFaceMatrixes_;
        std::size_t numberOfElementVectors_;
        std::size_t numberOfFaceVectors_;
    };
    

}

#endif /* MESHMANIPULATORBASE_H_ */
