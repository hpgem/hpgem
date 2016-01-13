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

#include "Base/MpiContainer.h"
#include "GlobalVector.h"
#include <vector>
#include "Base/MeshManipulatorBase.h"
#include "Base/Edge.h"
#include "Base/Face.h"
#include "Base/Element.h"
#include "Base/ElementCacheData.h"
#include "LinearAlgebra/MiddleSizeVector.h"
#include "Base/FaceCacheData.h"
#include "Geometry/PointPhysical.h"
#include "Geometry/PhysicalGeometry.h"
#include "Geometry/ReferenceGeometry.h"
#include "Geometry/PointReference.h"
#include <Logger.h>
#include <numeric>
#if defined(HPGEM_USE_PETSC) || defined(HPGEM_USE_COMPLEX_PETSC)
    #include "petscis.h"
#endif
#if defined(HPGEM_USE_SUNDIALS)
	#include "nvector/nvector_serial.h"
#endif

namespace Utilities
{
    
    GlobalVector::GlobalVector(Base::MeshManipulatorBase* theMesh, int elementVectorID, int faceVectorID)
            : meshLevel_(-2), elementVectorID_(elementVectorID), faceVectorID_(faceVectorID), startPositionsOfElementsInTheVector_(), theMesh_(theMesh)
    {
        logger.assert(theMesh!=nullptr, "Invalid mesh passed");
    }
    
#if defined(HPGEM_USE_PETSC) || defined(HPGEM_USE_COMPLEX_PETSC)
    
    GlobalPetscVector::GlobalPetscVector(Base::MeshManipulatorBase* theMesh, int elementVectorID, int faceVectorID)
            : GlobalVector(theMesh, elementVectorID, faceVectorID)
    {
        logger.assert(theMesh!=nullptr, "Invalid mesh passed");
        PetscBool petscRuns;
        PetscInitialized(&petscRuns);
        logger.assert(petscRuns == PETSC_TRUE, "Early call, firstly the command line arguments should be parsed");
        VecCreateSeq(PETSC_COMM_SELF, 1, &b_);
        
        reset();
    }
    
    GlobalPetscVector::~GlobalPetscVector()
    {
        int ierr = VecDestroy(&b_);
        CHKERRV(ierr);
    }
    
    GlobalPetscVector::operator Vec()
    {
        if(HPGEM_LOGLEVEL==Log::DEBUG)
        {
            VecView(b_, PETSC_VIEWER_STDOUT_WORLD);
        }
        return b_;
    }
    
    std::vector<PetscInt> GlobalPetscVector::makePositionsInVector(const Base::Element* element)
    {
        logger.assert(element!=nullptr, "invalid element passed");
        //we need storage for the amount of basis functions to return
        std::vector<PetscInt> positions(element->getNumberOfBasisFunctions() * element->getNumberOfUnknowns());
        
        auto pos = positions.begin();
        for(std::size_t index = 0; index < element->getNumberOfUnknowns(); ++index)
        {
            std::size_t numberOfElementBasisFunctions = element->getLocalNumberOfBasisFunctions();
            //First step: construct ids for the functions of the current element itself
            for (std::size_t i = 0; i < numberOfElementBasisFunctions; ++i)
            {
                *pos = i + startPositionsOfElementsInTheVector_[element->getID()] + index * numberOfElementBasisFunctions;
                pos++;
            }

            //Push forward our iterator
            std::size_t numberOfFaces = element->getPhysicalGeometry()->getNumberOfFaces();
            for (std::size_t i = 0; i < numberOfFaces; ++i)
            {
                std::size_t numberOfFaceBasisFunctions = element->getFace(i)->getLocalNumberOfBasisFunctions();
                for (std::size_t j = 0; j < numberOfFaceBasisFunctions; ++j)
                {
                    *pos = j + startPositionsOfFacesInTheVector_[element->getFace(i)->getID()] + index * numberOfFaceBasisFunctions;
                    pos++;
                }
            }

            std::size_t numberOfEdges = element->getNumberOfEdges();
            for (std::size_t i = 0; i < numberOfEdges; ++i)
            {
                std::size_t numberOfEdgeBasisFunctions = element->getEdge(i)->getLocalNumberOfBasisFunctions();
                for (std::size_t j = 0; j < numberOfEdgeBasisFunctions; ++j)
                {
                    *pos = j + startPositionsOfEdgesInTheVector_[element->getEdge(i)->getID()] + index * numberOfEdgeBasisFunctions;
                    pos++;
                }
            }

            std::size_t numberOfNodes = element->getNumberOfNodes();
            for (std::size_t i = 0; i < numberOfNodes; ++i)
            {
                std::size_t numberOfNodeBasisFunctions = element->getNode(i)->getLocalNumberOfBasisFunctions();
                for (std::size_t j = 0; j < numberOfNodeBasisFunctions; ++j)
                {
                    *pos = j + startPositionsOfNodesInTheVector_[element->getNode(i)->getID()] + index * numberOfNodeBasisFunctions;
                    pos++;
                }
            }
        }
        logger.assert(pos == positions.end(), "GlobalVector: did not process all elements correctly");
        return positions;
    }
    
    void GlobalPetscVector::reset()
    {
        int ierr = VecDestroy(&b_);
#ifdef HPGEM_USE_MPI
        std::size_t n = Base::MPIContainer::Instance().getNumberOfProcessors();
        //offset by one to put a 0 in front
        std::vector<int> MPISendElementCounts(n+1,0), MPISendFaceCounts(n+1,0), MPISendEdgeCounts(n+1,0), MPISendNodeCounts(n+1,0);

        int rank = Base::MPIContainer::Instance().getProcessorID();

        MPISendElementCounts[rank+1] = theMesh_->getNumberOfElements();
        MPISendFaceCounts[rank+1] = theMesh_->getNumberOfFaces();
        MPISendEdgeCounts[rank+1] = theMesh_->getNumberOfEdges();
        MPISendNodeCounts[rank+1] = theMesh_->getNumberOfNodes();

        //tell the rest of the processes how much info the others have
        MPI::Intracomm& comm = Base::MPIContainer::Instance().getComm();
        comm.Allgather(MPI_IN_PLACE,0,Base::Detail::toMPIType(rank),MPISendElementCounts.data()+1,1,Base::Detail::toMPIType(rank));
        comm.Allgather(MPI_IN_PLACE,0,Base::Detail::toMPIType(rank),MPISendFaceCounts.data()+1,1,Base::Detail::toMPIType(rank));
        comm.Allgather(MPI_IN_PLACE,0,Base::Detail::toMPIType(rank),MPISendEdgeCounts.data()+1,1,Base::Detail::toMPIType(rank));
        comm.Allgather(MPI_IN_PLACE,0,Base::Detail::toMPIType(rank),MPISendNodeCounts.data()+1,1,Base::Detail::toMPIType(rank));

        std::vector<int> MPISendElementStarts(n+1), MPISendFaceStarts(n+1), MPISendEdgeStarts(n+1), MPISendNodeStarts(n+1);

        std::partial_sum(MPISendElementCounts.begin(), MPISendElementCounts.end(), MPISendElementStarts.begin());
        std::partial_sum(MPISendFaceCounts.begin(), MPISendFaceCounts.end(), MPISendFaceStarts.begin());
        std::partial_sum(MPISendEdgeCounts.begin(), MPISendEdgeCounts.end(), MPISendEdgeStarts.begin());
        std::partial_sum(MPISendNodeCounts.begin(), MPISendNodeCounts.end(), MPISendNodeStarts.begin());

        //pack the computed data to send it using MPI
        std::vector<std::size_t> MPISendElementNumbers(MPISendElementStarts.back(), std::numeric_limits<std::size_t>::max());
        std::vector<std::size_t> MPISendFaceNumbers(MPISendFaceStarts.back(), std::numeric_limits<std::size_t>::max());
        std::vector<std::size_t> MPISendEdgeNumbers(MPISendEdgeStarts.back(), std::numeric_limits<std::size_t>::max());
        std::vector<std::size_t> MPISendNodeNumbers(MPISendNodeStarts.back(), std::numeric_limits<std::size_t>::max());

        std::vector<std::size_t> MPISendElementPositions(MPISendElementStarts.back(), std::numeric_limits<std::size_t>::max());
        std::vector<std::size_t> MPISendFacePositions(MPISendFaceStarts.back(), std::numeric_limits<std::size_t>::max());
        std::vector<std::size_t> MPISendEdgePositions(MPISendEdgeStarts.back(), std::numeric_limits<std::size_t>::max());
        std::vector<std::size_t> MPISendNodePositions(MPISendNodeStarts.back(), std::numeric_limits<std::size_t>::max());

        auto currentElementNumber = MPISendElementNumbers.begin() + MPISendElementStarts[rank];
        auto currentFaceNumber = MPISendFaceNumbers.begin() + MPISendFaceStarts[rank];
        auto currentEdgeNumber = MPISendEdgeNumbers.begin() + MPISendEdgeStarts[rank];
        auto currentNodeNumber = MPISendNodeNumbers.begin() + MPISendNodeStarts[rank];

        auto currentElementPosition = MPISendElementPositions.begin() + MPISendElementStarts[rank];
        auto currentFacePosition = MPISendFacePositions.begin() + MPISendFaceStarts[rank];
        auto currentEdgePosition = MPISendEdgePositions.begin() + MPISendEdgeStarts[rank];
        auto currentNodePosition = MPISendNodePositions.begin() + MPISendNodeStarts[rank];
#endif
        std::size_t totalNumberOfDOF(0), DIM(theMesh_->dimension());
        startPositionsOfElementsInTheVector_.resize(theMesh_->getNumberOfElements(Base::IteratorType::GLOBAL));
        startPositionsOfFacesInTheVector_.resize(theMesh_->getNumberOfFaces(Base::IteratorType::GLOBAL));
        startPositionsOfEdgesInTheVector_.resize(theMesh_->getNumberOfEdges(Base::IteratorType::GLOBAL));
        startPositionsOfNodesInTheVector_.resize(theMesh_->getNumberOfNodes(Base::IteratorType::GLOBAL));
        for (Base::Element* element : theMesh_->getElementsList())
        {
#ifdef HPGEM_USE_MPI
            *currentElementNumber=element->getID();
            *currentElementPosition=totalNumberOfDOF;
            ++currentElementNumber;
            ++currentElementPosition;
#else
            startPositionsOfElementsInTheVector_[element->getID()] = totalNumberOfDOF;
#endif
            totalNumberOfDOF += element->getLocalNumberOfBasisFunctions() * element->getNumberOfUnknowns();
            for (std::size_t i = 0; i < element->getNumberOfFaces(); ++i)
            {
                //faces at the boundary of the subdomain should also be added only once, so add them here is the left element of the face is in the subdomain
                if ((element->getFace(i)->getFaceType() == Geometry::FaceType::SUBDOMAIN_BOUNDARY || element->getFace(i)->getFaceType() == Geometry::FaceType::PERIODIC_SUBDOMAIN_BC)
                        && element->getFace(i)->getPtrElementLeft() == element)
                {
#ifdef HPGEM_USE_MPI
                    *currentFaceNumber=element->getFace(i)->getID();
                    *currentFacePosition=totalNumberOfDOF;
                    ++currentFaceNumber;
                    ++currentFacePosition;
#else
                    startPositionsOfFacesInTheVector_[element->getFace(i)->getID()] = totalNumberOfDOF;
#endif
                    totalNumberOfDOF += element->getFace(i)->getLocalNumberOfBasisFunctions() * element->getNumberOfUnknowns();
                }
                
            }
        }
        for (Base::Face* face : theMesh_->getFacesList())
        {
            //skip faces at the subdomain boundary because we already treated them
            if (face->getFaceType() != Geometry::FaceType::SUBDOMAIN_BOUNDARY && face->getFaceType() != Geometry::FaceType::PERIODIC_SUBDOMAIN_BC)
            {
#ifdef HPGEM_USE_MPI
                *currentFaceNumber=face->getID();
                *currentFacePosition=totalNumberOfDOF;
                ++currentFaceNumber;
                ++currentFacePosition;
#else
                startPositionsOfFacesInTheVector_[face->getID()] = totalNumberOfDOF;
#endif
                totalNumberOfDOF += face->getLocalNumberOfBasisFunctions() * face->getPtrElementLeft()->getNumberOfUnknowns();
            }
        }
        for (Base::Edge* edge : theMesh_->getEdgesList())
        {
#ifdef HPGEM_USE_MPI
            *currentEdgeNumber=edge->getID();
            *currentEdgePosition=totalNumberOfDOF;
            ++currentEdgeNumber;
            ++currentEdgePosition;
#else
            startPositionsOfEdgesInTheVector_[edge->getID()] = totalNumberOfDOF;
#endif
            totalNumberOfDOF += edge->getLocalNumberOfBasisFunctions() * edge->getElement(0)->getNumberOfUnknowns();
        }
        //DIM == 1 faces and nodes are the same entities, skip one of them
        if (DIM > 1)
        {
            for (Base::Node* node : theMesh_->getNodesList())
            {
#ifdef HPGEM_USE_MPI
                *currentNodeNumber=node->getID();
                *currentNodePosition=totalNumberOfDOF;
                ++currentNodeNumber;
                ++currentNodePosition;
#else
                startPositionsOfNodesInTheVector_[node->getID()] = totalNumberOfDOF;
#endif
                totalNumberOfDOF += node->getLocalNumberOfBasisFunctions() * node->getElement(0)->getNumberOfUnknowns();
            }
        }
        
#ifdef HPGEM_USE_MPI
        //offset by one to insert a zero at the front
        std::vector<std::size_t> cumulativeDOF(n + 1);
        cumulativeDOF[rank+1]=totalNumberOfDOF;

        //communicate...
        comm.Allgather(MPI_IN_PLACE,0,Base::Detail::toMPIType(n),cumulativeDOF.data()+1,1,Base::Detail::toMPIType(n));

        comm.Allgatherv(MPI_IN_PLACE,0,Base::Detail::toMPIType(n),MPISendElementNumbers.data(),MPISendElementCounts.data()+1,MPISendElementStarts.data(),Base::Detail::toMPIType(n));
        comm.Allgatherv(MPI_IN_PLACE,0,Base::Detail::toMPIType(n),MPISendFaceNumbers.data(),MPISendFaceCounts.data()+1,MPISendFaceStarts.data(),Base::Detail::toMPIType(n));
        comm.Allgatherv(MPI_IN_PLACE,0,Base::Detail::toMPIType(n),MPISendEdgeNumbers.data(),MPISendEdgeCounts.data()+1,MPISendEdgeStarts.data(),Base::Detail::toMPIType(n));
        comm.Allgatherv(MPI_IN_PLACE,0,Base::Detail::toMPIType(n),MPISendNodeNumbers.data(),MPISendNodeCounts.data()+1,MPISendNodeStarts.data(),Base::Detail::toMPIType(n));

        comm.Allgatherv(MPI_IN_PLACE,0,Base::Detail::toMPIType(n),MPISendElementPositions.data(),MPISendElementCounts.data()+1,MPISendElementStarts.data(),Base::Detail::toMPIType(n));
        comm.Allgatherv(MPI_IN_PLACE,0,Base::Detail::toMPIType(n),MPISendFacePositions.data(),MPISendFaceCounts.data()+1,MPISendFaceStarts.data(),Base::Detail::toMPIType(n));
        comm.Allgatherv(MPI_IN_PLACE,0,Base::Detail::toMPIType(n),MPISendEdgePositions.data(),MPISendEdgeCounts.data()+1,MPISendEdgeStarts.data(),Base::Detail::toMPIType(n));
        comm.Allgatherv(MPI_IN_PLACE,0,Base::Detail::toMPIType(n),MPISendNodePositions.data(),MPISendNodeCounts.data()+1,MPISendNodeStarts.data(),Base::Detail::toMPIType(n));

        //and unpack the information
        std::partial_sum(cumulativeDOF.begin(),cumulativeDOF.end(),cumulativeDOF.begin());

        currentElementNumber = MPISendElementNumbers.begin();
        currentFaceNumber = MPISendFaceNumbers.begin();
        currentEdgeNumber = MPISendEdgeNumbers.begin();
        currentNodeNumber = MPISendNodeNumbers.begin();

        currentElementPosition = MPISendElementPositions.begin();
        currentFacePosition = MPISendFacePositions.begin();
        currentEdgePosition = MPISendEdgePositions.begin();
        currentNodePosition = MPISendNodePositions.begin();

        std::size_t currentDomain = 0;
        auto startOFNextDomain = MPISendElementNumbers.begin() + MPISendElementCounts[currentDomain+1];
        std::size_t offset = cumulativeDOF[currentDomain];
        for(;currentElementNumber!=MPISendElementNumbers.end();++currentElementNumber,++currentElementPosition)
        {   
            if(currentElementNumber==startOFNextDomain)
            {   
                currentDomain++;
                startOFNextDomain += MPISendElementCounts[currentDomain+1];
                offset = cumulativeDOF[currentDomain];
            }
            logger.assert(*currentElementNumber != std::numeric_limits<std::size_t>::max(), "currentElementNumber = -1");
            startPositionsOfElementsInTheVector_[*currentElementNumber]=*currentElementPosition+offset;
        }

        currentDomain = 0;
        startOFNextDomain = MPISendFaceNumbers.begin() + MPISendFaceCounts[currentDomain+1];
        offset = cumulativeDOF[currentDomain];
        for(;currentFaceNumber!=MPISendFaceNumbers.end();++currentFaceNumber,++currentFacePosition)
        {   
            if(currentFaceNumber==startOFNextDomain)
            {   
                currentDomain++;
                startOFNextDomain += MPISendFaceCounts[currentDomain+1];
                offset = cumulativeDOF[currentDomain];
            }
            if (*currentFaceNumber != std::numeric_limits<std::size_t>::max())
            startPositionsOfFacesInTheVector_[*currentFaceNumber]=*currentFacePosition+offset;
        }

        currentDomain = 0;
        startOFNextDomain = MPISendEdgeNumbers.begin() + MPISendEdgeCounts[currentDomain+1];
        offset = cumulativeDOF[currentDomain];
        for(;currentEdgeNumber!=MPISendEdgeNumbers.end();++currentEdgeNumber,++currentEdgePosition)
        {   
            if(currentEdgeNumber==startOFNextDomain)
            {   
                currentDomain++;
                startOFNextDomain += MPISendEdgeCounts[currentDomain+1];
                offset = cumulativeDOF[currentDomain];
            }
            logger.assert(*currentEdgeNumber != std::numeric_limits<std::size_t>::max(), "currentEdgeNumber=-1");
            startPositionsOfEdgesInTheVector_[*currentEdgeNumber]=*currentEdgePosition+offset;
        }

        currentDomain = 0;
        startOFNextDomain = MPISendNodeNumbers.begin() + MPISendNodeCounts[currentDomain+1];
        offset = cumulativeDOF[currentDomain];
        for(;currentNodeNumber!=MPISendNodeNumbers.end();++currentNodeNumber,++currentNodePosition)
        {   
            if(currentNodeNumber==startOFNextDomain)
            {   
                currentDomain++;
                startOFNextDomain += MPISendNodeCounts[currentDomain+1];
                offset = cumulativeDOF[currentDomain];
            }
            if (*currentNodeNumber != std::numeric_limits<std::size_t>::max())
            startPositionsOfNodesInTheVector_[*currentNodeNumber]=*currentNodePosition+offset;
        }
#endif
        
        ierr = VecCreateMPI(PETSC_COMM_WORLD, totalNumberOfDOF, PETSC_DETERMINE, &b_);
        CHKERRV(ierr);
    }
    
    void GlobalPetscVector::assemble()
    {
        reset();
        
        LinearAlgebra::MiddleSizeVector elementVector;
        
        if (elementVectorID_ >= 0)
        {
            for (Base::Element* element : theMesh_->getElementsList())
            {
                std::vector<PetscInt> positions = makePositionsInVector(element);
                elementVector = element->getElementVector(elementVectorID_);
                int ierr = VecSetValues(b_, positions.size(), positions.data(), elementVector.data(), ADD_VALUES);
                CHKERRV(ierr);
                
            }
        }
        
        LinearAlgebra::MiddleSizeVector faceVector;
        
        if (faceVectorID_ >= 0)
        {
            for (Base::Face* face : theMesh_->getFacesList())
            {
                std::vector<PetscInt> positions = makePositionsInVector(face->getPtrElementLeft());
                if (face->isInternal())
                {
                    std::vector<PetscInt> rightPositions = makePositionsInVector(face->getPtrElementRight());
                    positions.reserve(positions.size() + rightPositions.size());
                    for (auto& a : rightPositions)
                    {
                        positions.push_back(a);
                    }
                }
                faceVector = face->getFaceVector(faceVectorID_);
                int ierr = VecSetValues(b_, positions.size(), positions.data(), faceVector.data(), ADD_VALUES);
                CHKERRV(ierr);
            }
        }
        
        int ierr = VecAssemblyBegin(b_);
        ierr = VecAssemblyEnd(b_);
        CHKERRV(ierr);
    }
    
    void GlobalPetscVector::constructFromTimeIntegrationVector(std::size_t timeIntegrationVectorId, std::size_t solutionVar)
    {
        reset();
        
        LinearAlgebra::MiddleSizeVector elementData;
        for (Base::Element* element : theMesh_->getElementsList())
        {
            std::size_t numberOfBasisFunctions = element->getNumberOfBasisFunctions();
            std::vector<PetscInt> positions = makePositionsInVector(element);
            elementData.resize(numberOfBasisFunctions * element->getNumberOfUnknowns());
            for (std::size_t i = 0; i < numberOfBasisFunctions; ++i)
            {
                elementData[element->convertToSingleIndex(i, solutionVar)] = element->getTimeIntegrationData(timeIntegrationVectorId, solutionVar, i);
                for(std::size_t j = 0; j < element->getNumberOfUnknowns(); ++j)
                {
                    if(j != solutionVar)
                    {
                        positions[element->convertToSingleIndex(i, j)] = -1;
                    }
                }
            }
            int ierr = VecSetValues(b_, numberOfBasisFunctions * element->getNumberOfUnknowns(), positions.data(), elementData.data(), INSERT_VALUES);
            CHKERRV(ierr);
        }
        
        int ierr = VecAssemblyBegin(b_);
        ierr = VecAssemblyEnd(b_);
        CHKERRV(ierr);
    }
    
    void GlobalPetscVector::constructFromTimeIntegrationVector(std::size_t timeIntegrationVectorId)
    {
        reset();

        LinearAlgebra::MiddleSizeVector elementData;
        for (Base::Element* element : theMesh_->getElementsList())
        {
            std::size_t numberOfBasisFunctions = element->getNumberOfBasisFunctions();
            std::vector<PetscInt> positions = makePositionsInVector(element);
            int ierr = VecSetValues(b_, numberOfBasisFunctions * element->getNumberOfUnknowns(), positions.data(), element->getTimeIntegrationVector(timeIntegrationVectorId).data(), INSERT_VALUES);
            CHKERRV(ierr);
        }

        int ierr = VecAssemblyBegin(b_);
        ierr = VecAssemblyEnd(b_);
        CHKERRV(ierr);
    }

    void GlobalPetscVector::writeTimeIntegrationVector(std::size_t timeIntegrationVectorId)
    {
        PetscScalar *data;
        
        VecScatter scatter;
        Vec localB;
        //create a local vector...
        VecScatterCreateToAll(b_, &scatter, &localB);
        //we dont need the scatter context, just the vector
        VecScatterDestroy(&scatter);
        
        std::vector<PetscInt> positions;
        std::size_t totalPositions = 0;
        for(Base::Element* element : theMesh_->getElementsList())
        {
            totalPositions += element->getNumberOfBasisFunctions() * element->getNumberOfUnknowns();
        }
        positions.reserve(totalPositions);
        for (Base::Element* element : theMesh_->getElementsList())
        {
            std::vector<PetscInt> newPositions = makePositionsInVector(element);
            for (auto& a : newPositions)
            {
                positions.push_back(a);
            }
        }
        
        IS scatterIS;
        ISCreateGeneral(PETSC_COMM_SELF, positions.size(), positions.data(), PETSC_COPY_VALUES, &scatterIS);
        ISSortRemoveDups(scatterIS);
        VecScatterCreate(b_, scatterIS, localB, scatterIS, &scatter);
        VecScatterBegin(scatter, b_, localB, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(scatter, b_, localB, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterDestroy(&scatter);
        ISDestroy(&scatterIS);
        
        int ierr = VecGetArray(localB, &data);
        CHKERRV(ierr);
        //Question: wat doet deze lijn
        for (Base::MeshManipulatorBase::ElementIterator it = theMesh_->elementColBegin(); it != theMesh_->elementColEnd(); ++it)
        {
        	//Create vector for the local data that has to be written
            LinearAlgebra::MiddleSizeVector localData((*it)->getNumberOfBasisFunctions() * (*it)->getNumberOfUnknowns());
            std::size_t runningTotal = 0;
            for(std::size_t index = 0; index < (*it)->getNumberOfUnknowns(); ++index) // for every variable iV
            {
                for(std::size_t i = 0; i < (*it)->getLocalNumberOfBasisFunctions(); ++i) // Get the local basis functions of the element
                {
                	//Copy the values from data to the localData and update the running number
                    localData[runningTotal] = std::real(data[startPositionsOfElementsInTheVector_[(*it)->getID()] + i + index * (*it)->getLocalNumberOfBasisFunctions()]);
                    ++runningTotal;
                }
                for (std::size_t i = 0; i < (*it)->getPhysicalGeometry()->getNumberOfFaces(); ++i) // for all faces of the element
                {
                    for (std::size_t j = 0; j < (*it)->getFace(i)->getLocalNumberOfBasisFunctions(); ++j) // get local basis functions of a face
                    {
                    	//Copy the values from data to the localData and update the running number
                        localData[runningTotal] = std::real(data[startPositionsOfFacesInTheVector_[(*it)->getFace(i)->getID()] + j + index * (*it)->getFace(i)->getLocalNumberOfBasisFunctions()]);
                        ++runningTotal;
                    }
                }
                for (std::size_t i = 0; i < (*it)->getNumberOfEdges(); ++i) // For all edges of the element
                {
                    for (std::size_t j = 0; j < (*it)->getEdge(i)->getLocalNumberOfBasisFunctions(); ++j) //Get the local basis function of an edge
                    {
                    	//Copy the values from data to the localData and update the running number
                        localData[runningTotal] = std::real(data[startPositionsOfEdgesInTheVector_[(*it)->getEdge(i)->getID()] + j + index * (*it)->getEdge(i)->getLocalNumberOfBasisFunctions()]);
                        ++runningTotal;
                    }
                }
                if (theMesh_->dimension() > 1) // There are no nodes in a 1D problem
                {
                    for (std::size_t i = 0; i < (*it)->getNumberOfNodes(); ++i) //For all nodes
                    {
                        for (std::size_t j = 0; j < (*it)->getNode(i)->getLocalNumberOfBasisFunctions(); ++j) //Get the local number of basis function of a node
                        {
                        	//Copy the values from data to the localData and update the running number
                            localData[runningTotal] = std::real(data[startPositionsOfNodesInTheVector_[(*it)->getNode(i)->getID()] + j + index * (*it)->getNode(i)->getLocalNumberOfBasisFunctions()]);
                            ++runningTotal;
                        }
                    }
                }
            }
            //Put the localData in the element
            logger.assert(localData.size() == runningTotal, "not enough info to fill the vector");
            (*it)->setTimeIntegrationVector(timeIntegrationVectorId, localData);
        }
        ierr = VecRestoreArray(localB, &data);
        VecDestroy(&localB);
        CHKERRV(ierr);
    }

    void GlobalPetscVector::writeTimeIntegrationVector(std::size_t timeIntegrationVectorId, std::size_t unknown)
    {
        PetscScalar *data;

        VecScatter scatter;
        Vec localB;
        //create a local vector...
        VecScatterCreateToAll(b_, &scatter, &localB);
        //we dont need the scatter context, just the vector
        VecScatterDestroy(&scatter);

        std::vector<PetscInt> positions;
        std::size_t totalPositions = 0;
        for(Base::Element* element : theMesh_->getElementsList())
        {
            totalPositions += element->getNumberOfBasisFunctions() * element->getNumberOfUnknowns();
        }
        positions.reserve(totalPositions);
        for (Base::Element* element : theMesh_->getElementsList())
        {
            std::vector<PetscInt> newPositions = makePositionsInVector(element);
            for (auto& a : newPositions)
            {
                positions.push_back(a);
            }
        }

        IS scatterIS;
        ISCreateGeneral(PETSC_COMM_SELF, positions.size(), positions.data(), PETSC_COPY_VALUES, &scatterIS);
        ISSortRemoveDups(scatterIS);
        VecScatterCreate(b_, scatterIS, localB, scatterIS, &scatter);
        VecScatterBegin(scatter, b_, localB, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(scatter, b_, localB, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterDestroy(&scatter);
        ISDestroy(&scatterIS);

        int ierr = VecGetArray(localB, &data);
        CHKERRV(ierr);
        for (Base::MeshManipulatorBase::ElementIterator it = theMesh_->elementColBegin(); it != theMesh_->elementColEnd(); ++it)
        {
            LinearAlgebra::MiddleSizeVector localData((*it)->getNumberOfBasisFunctions() * (*it)->getNumberOfUnknowns());
            std::size_t runningTotal = 0;
            for(std::size_t index = 0; index < (*it)->getNumberOfUnknowns(); ++index)
            {
                for(std::size_t i = 0; i < (*it)->getLocalNumberOfBasisFunctions(); ++i)
                {
                    localData[runningTotal] = std::real(data[startPositionsOfElementsInTheVector_[(*it)->getID()] + i + index * (*it)->getLocalNumberOfBasisFunctions()]);
                    ++runningTotal;
                }
                for (std::size_t i = 0; i < (*it)->getPhysicalGeometry()->getNumberOfFaces(); ++i)
                {
                    for (std::size_t j = 0; j < (*it)->getFace(i)->getLocalNumberOfBasisFunctions(); ++j)
                    {
                        localData[runningTotal] = std::real(data[startPositionsOfFacesInTheVector_[(*it)->getFace(i)->getID()] + j + index * (*it)->getFace(i)->getLocalNumberOfBasisFunctions()]);
                        ++runningTotal;
                    }
                }
                for (std::size_t i = 0; i < (*it)->getNumberOfEdges(); ++i)
                {
                    for (std::size_t j = 0; j < (*it)->getEdge(i)->getLocalNumberOfBasisFunctions(); ++j)
                    {
                        localData[runningTotal] = std::real(data[startPositionsOfEdgesInTheVector_[(*it)->getEdge(i)->getID()] + j + index * (*it)->getEdge(i)->getLocalNumberOfBasisFunctions()]);
                        ++runningTotal;
                    }
                }
                if (theMesh_->dimension() > 1)
                {
                    for (std::size_t i = 0; i < (*it)->getNumberOfNodes(); ++i)
                    {
                        for (std::size_t j = 0; j < (*it)->getNode(i)->getLocalNumberOfBasisFunctions(); ++j)
                        {
                            localData[runningTotal] = std::real(data[startPositionsOfNodesInTheVector_[(*it)->getNode(i)->getID()] + j + index * (*it)->getNode(i)->getLocalNumberOfBasisFunctions()]);
                            ++runningTotal;
                        }
                    }
                }
            }
            logger.assert(localData.size() == runningTotal, "not enough info to fill the vector");
            LinearAlgebra::MiddleSizeVector singleUnknownData((*it)->getNumberOfBasisFunctions());
            for(std::size_t i = 0; i < (*it)->getNumberOfBasisFunctions(); ++i)
            {
                singleUnknownData[i] = localData[(*it)->convertToSingleIndex(i, unknown)];
            }
            (*it)->setTimeIntegrationSubvector(timeIntegrationVectorId, unknown, singleUnknownData);
        }
        ierr = VecRestoreArray(localB, &data);
        VecDestroy(&localB);
        CHKERRV(ierr);
    }
#endif


#if defined(HPGEM_USE_SUNDIALS)

    //Constructor that will initialise the vector with zero entries but the correct length of the problem, nDOF
    GlobalSundialsVector::GlobalSundialsVector(Base::MeshManipulatorBase* theMesh, int elementVectorID, int faceVectorID)
            : GlobalVector(theMesh, elementVectorID, faceVectorID)
    {
        logger.assert(theMesh!=nullptr, "Invalid mesh passed");
        //b_ = N_VNew_Serial(1); //Create a temporary placeholder vector
        reset();
    }

    GlobalSundialsVector::~GlobalSundialsVector()
    {
    	//N_VDestroy_Serial(b_);
    }

    GlobalSundialsVector::operator N_Vector()
    {
        return b_;
    }

    /// \brief Resets the N_Vector
    void GlobalSundialsVector::reset()
    {
        //Destroy old vector
    	//N_VDestroy_Serial(b_);

    	//Assign private variables to the global vector
    	//Initialise totalNrOfDOF to zero, this will be calculated as we go
        std::size_t totalNumberOfDOF(0), DIM(theMesh_->dimension());
        //Resize the startPosition vectors
        startPositionsOfElementsInTheVector_.resize(theMesh_->getNumberOfElements(Base::IteratorType::GLOBAL));
        startPositionsOfFacesInTheVector_.resize(theMesh_->getNumberOfFaces(Base::IteratorType::GLOBAL));
        startPositionsOfEdgesInTheVector_.resize(theMesh_->getNumberOfEdges(Base::IteratorType::GLOBAL));
        startPositionsOfNodesInTheVector_.resize(theMesh_->getNumberOfNodes(Base::IteratorType::GLOBAL));

        for (Base::Element* element : theMesh_->getElementsList())
        {
        	//Append the start of a new element to the vector, Then update the new total number of DOF
            startPositionsOfElementsInTheVector_[element->getID()] = totalNumberOfDOF; //New element position is appended to the vector
            totalNumberOfDOF += element->getLocalNumberOfBasisFunctions() * element->getNumberOfUnknowns(); //Calculate the new total DOF

            for (std::size_t i = 0; i < element->getNumberOfFaces(); ++i)
            {
                //faces at the boundary of the subdomain should also be added only once, so add them here is the left element of the face is in the subdomain
                if ((element->getFace(i)->getFaceType() == Geometry::FaceType::SUBDOMAIN_BOUNDARY || element->getFace(i)->getFaceType() == Geometry::FaceType::PERIODIC_SUBDOMAIN_BC)
                        && element->getFace(i)->getPtrElementLeft() == element)
                {
                    startPositionsOfFacesInTheVector_[element->getFace(i)->getID()] = totalNumberOfDOF;
                    totalNumberOfDOF += element->getFace(i)->getLocalNumberOfBasisFunctions() * element->getNumberOfUnknowns();
                    std::cout << "face local basis functions: " << element->getFace(i)->getLocalNumberOfBasisFunctions() << std::endl;
                }

            }
        }
        for (Base::Face* face : theMesh_->getFacesList())
        {
            //skip faces at the subdomain boundary because we already treated them
            if (face->getFaceType() != Geometry::FaceType::SUBDOMAIN_BOUNDARY && face->getFaceType() != Geometry::FaceType::PERIODIC_SUBDOMAIN_BC)
            {
                startPositionsOfFacesInTheVector_[face->getID()] = totalNumberOfDOF;
                totalNumberOfDOF += face->getLocalNumberOfBasisFunctions() * face->getPtrElementLeft()->getNumberOfUnknowns();
            }
        }
        for (Base::Edge* edge : theMesh_->getEdgesList())
        {
            startPositionsOfEdgesInTheVector_[edge->getID()] = totalNumberOfDOF;
            totalNumberOfDOF += edge->getLocalNumberOfBasisFunctions() * edge->getElement(0)->getNumberOfUnknowns();
        }
        //DIM == 1 faces and nodes are the same entities, skip one of them
        if (DIM > 1)
        {
            for (Base::Node* node : theMesh_->getNodesList())
            {
                startPositionsOfNodesInTheVector_[node->getID()] = totalNumberOfDOF;
                totalNumberOfDOF += node->getLocalNumberOfBasisFunctions() * node->getElement(0)->getNumberOfUnknowns();
            }
        }

        //b_ = N_VNew_Serial(totalNumberOfDOF);
        totalNumberOfDOF_ = totalNumberOfDOF;
    }

    void GlobalSundialsVector::constructFromTimeIntegrationVector(std::size_t timeIntegrationVectorId)
   {
       	//Extract data pointer from the input vector
    	double *data = NV_DATA_S(b_);

    		//For all elements
            for (Base::Element* element : theMesh_->getElementsList())
            {
            	//Create vector for the local data that has to be read
                LinearAlgebra::MiddleSizeVector localData = element->getTimeIntegrationVector(timeIntegrationVectorId);
                std::size_t runningTotal = 0;
                for(std::size_t index = 0; index < element->getNumberOfUnknowns(); ++index) // for every variable iV
                {
                    for(std::size_t i = 0; i < element->getLocalNumberOfBasisFunctions(); ++i) // Get the local basis functions of the element
                    {
                    	//Copy the values from localData to data and update the running number
                        data[startPositionsOfElementsInTheVector_[element->getID()] + i + index * element->getLocalNumberOfBasisFunctions()] = localData[runningTotal];
                        ++runningTotal;

                    }
                    // The code below does not return anything in case of DG, however in case of conforming FEM they will return data
                    for (std::size_t i = 0; i < element->getPhysicalGeometry()->getNumberOfFaces(); ++i) // for all faces of the element
                    {
                        for (std::size_t j = 0; j < element->getFace(i)->getLocalNumberOfBasisFunctions(); ++j) // get local basis functions of a face
                        {
                            data[startPositionsOfFacesInTheVector_[element->getFace(i)->getID()] + j + index * element->getFace(i)->getLocalNumberOfBasisFunctions()] = localData[runningTotal];
                            ++runningTotal;
                        }
                    }
                    for (std::size_t i = 0; i < element->getNumberOfEdges(); ++i) // For all edges of the element
                    {
                        for (std::size_t j = 0; j < element->getEdge(i)->getLocalNumberOfBasisFunctions(); ++j) //Get the local basis function of an edge
                        {
                        	//Copy the values from localData to data and update the running number
                        	data[startPositionsOfEdgesInTheVector_[element->getEdge(i)->getID()] + j + index * element->getEdge(i)->getLocalNumberOfBasisFunctions()] = localData[runningTotal];
                            ++runningTotal;
                        }
                    }
                    if (theMesh_->dimension() > 1) // There are no nodes in a 1D problem
                    {
                        for (std::size_t i = 0; i < element->getNumberOfNodes(); ++i) //For all nodes
                        {
                            for (std::size_t j = 0; j < element->getNode(i)->getLocalNumberOfBasisFunctions(); ++j) //Get the local number of basis function of a node
                            {
                            	//Copy the values from localData to data and update the running number
                            	data[startPositionsOfNodesInTheVector_[element->getNode(i)->getID()] + j + index * element->getNode(i)->getLocalNumberOfBasisFunctions()] = localData[runningTotal];
                                ++runningTotal;
                            }
                        }
                    }
                }
                logger.assert(localData.size() == runningTotal, "not enough info to fill the vector");
            }
    }

    void GlobalSundialsVector::setScale(LinearAlgebra::MiddleSizeVector scaleFactor)
    {
       	//Extract data pointer from the input vector
    	double *data = NV_DATA_S(b_);

    	//For all elements
    	for (Base::Element* element : theMesh_->getElementsList())
    	{
    		std::size_t runningTotal = 0;
    		for(std::size_t index = 0; index < element->getNumberOfUnknowns(); ++index) // for every variable iV
    		{
    			for(std::size_t i = 0; i < element->getLocalNumberOfBasisFunctions(); ++i) // Get the local basis functions of the element
    			{
    				//Copy the scale value
    				data[startPositionsOfElementsInTheVector_[element->getID()] + i + index * element->getLocalNumberOfBasisFunctions()] = scaleFactor(index);
    				++runningTotal;

    			}
    			// The code below does not return anything in case of DG, however in case of conforming FEM they will return data
    			for (std::size_t i = 0; i < element->getPhysicalGeometry()->getNumberOfFaces(); ++i) // for all faces of the element
    			{
    				for (std::size_t j = 0; j < element->getFace(i)->getLocalNumberOfBasisFunctions(); ++j) // get local basis functions of a face
    				{
    					data[startPositionsOfFacesInTheVector_[element->getFace(i)->getID()] + j + index * element->getFace(i)->getLocalNumberOfBasisFunctions()] = scaleFactor(index);
    					++runningTotal;
    				}
    			}
    			for (std::size_t i = 0; i < element->getNumberOfEdges(); ++i) // For all edges of the element
    			{
    				for (std::size_t j = 0; j < element->getEdge(i)->getLocalNumberOfBasisFunctions(); ++j) //Get the local basis function of an edge
    				{
    					//Copy the scale value
    					data[startPositionsOfEdgesInTheVector_[element->getEdge(i)->getID()] + j + index * element->getEdge(i)->getLocalNumberOfBasisFunctions()] = scaleFactor(index);
    					++runningTotal;
    				}
    			}
    			if (theMesh_->dimension() > 1) // There are no nodes in a 1D problem
    			{
    				for (std::size_t i = 0; i < element->getNumberOfNodes(); ++i) //For all nodes
    				{
    					for (std::size_t j = 0; j < element->getNode(i)->getLocalNumberOfBasisFunctions(); ++j) //Get the local number of basis function of a node
    					{
    						//Copy the scale value
    						data[startPositionsOfNodesInTheVector_[element->getNode(i)->getID()] + j + index * element->getNode(i)->getLocalNumberOfBasisFunctions()] = scaleFactor(index);
    						++runningTotal;
    					}
    				}
    			}
    		}
    	}
    }

    void GlobalSundialsVector::writeTimeIntegrationVector(std::size_t timeIntegrationVectorId)
    {
    	//do something here
    	double *data = NV_DATA_S(b_);

        for (Base::Element* element : theMesh_->getElementsList())
        {
        	//Create vector for the local data that has to be written
            LinearAlgebra::MiddleSizeVector localData(element->getNumberOfBasisFunctions() * element->getNumberOfUnknowns());
            std::size_t runningTotal = 0;
            for(std::size_t index = 0; index < element->getNumberOfUnknowns(); ++index) // for every variable iV
            {
                for(std::size_t i = 0; i < element->getLocalNumberOfBasisFunctions(); ++i) // Get the local basis functions of the element
                {
                	//Copy the values from data to the localData and update the running number
                    localData[runningTotal] = std::real(data[startPositionsOfElementsInTheVector_[element->getID()] + i + index * element->getLocalNumberOfBasisFunctions()]);
                    ++runningTotal;
                }
                // The code below does not return anything in case of DG, however in case of conforming FEM they will return data
                for (std::size_t i = 0; i < element->getPhysicalGeometry()->getNumberOfFaces(); ++i) // for all faces of the element
                {
                    for (std::size_t j = 0; j < element->getFace(i)->getLocalNumberOfBasisFunctions(); ++j) // get local basis functions of a face
                    {
                    	//Copy the values from data to the localData and update the running number
                        localData[runningTotal] = std::real(data[startPositionsOfFacesInTheVector_[element->getFace(i)->getID()] + j + index * element->getFace(i)->getLocalNumberOfBasisFunctions()]);
                        ++runningTotal;
                    }
                }
                for (std::size_t i = 0; i < element->getNumberOfEdges(); ++i) // For all edges of the element
                {
                    for (std::size_t j = 0; j < element->getEdge(i)->getLocalNumberOfBasisFunctions(); ++j) //Get the local basis function of an edge
                    {
                    	//Copy the values from data to the localData and update the running number
                        localData[runningTotal] = std::real(data[startPositionsOfEdgesInTheVector_[element->getEdge(i)->getID()] + j + index * element->getEdge(i)->getLocalNumberOfBasisFunctions()]);
                        ++runningTotal;
                    }
                }
                if (theMesh_->dimension() > 1) // There are no nodes in a 1D problem
                {
                    for (std::size_t i = 0; i < element->getNumberOfNodes(); ++i) //For all nodes
                    {
                        for (std::size_t j = 0; j < element->getNode(i)->getLocalNumberOfBasisFunctions(); ++j) //Get the local number of basis function of a node
                        {
                        	//Copy the values from data to the localData and update the running number
                            localData[runningTotal] = std::real(data[startPositionsOfNodesInTheVector_[element->getNode(i)->getID()] + j + index * element->getNode(i)->getLocalNumberOfBasisFunctions()]);
                            ++runningTotal;
                        }
                    }
                }
            }
            //Put the localData in the element
            logger.assert(localData.size() == runningTotal, "not enough info to fill the vector");
            element->setTimeIntegrationVector(timeIntegrationVectorId, localData);
        }
    }


    void GlobalSundialsVector::print()
    {
    	std::size_t length = NV_LENGTH_S(b_);
    	double *data = NV_DATA_S(b_);

		std::cout << "[" << data[0];
    	for (std::size_t i = 1; i < length; i++)
    	{
    		std::cout << ", " << data[i];
    	}
    	std::cout << ']' << std::endl;
    }

    //Possibly deprecated
/*    LinearAlgebra::MiddleSizeVector GlobalSundialsVector::getLocalVector(const Base::Element* ptrElement)
    {
    	//Create the local vector
    	std::size_t numberOfDOF = ptrElement->getNumberOfBasisFunctions() * ptrElement->getNumberOfUnknowns();
    	LinearAlgebra::MiddleSizeVector localVector(numberOfDOF);

    	//Fill the local vector
    	for (std::size_t i = 0; i < numberOfDOF; i++)
    	{
    		localVector(i) = startPositionsOfElementsInTheVector_[ptrElement->getID()] + i;
    	}

    	return localVector;
    }*/

    void GlobalSundialsVector::setVector(N_Vector b)
    {
    	b_ = b;
    }

#endif

}

