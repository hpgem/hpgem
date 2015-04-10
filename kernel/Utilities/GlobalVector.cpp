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
#include "Base/MeshManipulator.h"
#include "Base/Edge.h"
#include "Base/Face.h"
#include "Base/Element.h"
#include "Base/ElementCacheData.h"
#include "LinearAlgebra/NumericalVector.h"
#include "Base/FaceCacheData.h"
#include "Geometry/PointPhysical.h"
#include "Geometry/PhysicalGeometry.h"
#include "Geometry/ReferenceGeometry.h"
#include "Geometry/PointReference.h"
#include <Logger.h>
#include <numeric>

namespace Utilities
{
    
    GlobalVector::GlobalVector(Base::MeshManipulator* theMesh, int elementVectorID, int faceVectorID)
            : meshLevel_(-2), elementVectorID_(elementVectorID), faceVectorID_(faceVectorID), startPositionsOfElementsInTheVector_(), theMesh_(theMesh)
    {
        logger.assert(theMesh!=nullptr, "Invalid mesh passed");
    }
    
#if defined(HPGEM_USE_PETSC) || defined(HPGEM_USE_COMPLEX_PETSC)
    
    GlobalPetscVector::GlobalPetscVector(Base::MeshManipulator* theMesh, int elementVectorID, int faceVectorID)
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
        return b_;
    }
    
    std::vector<PetscInt> GlobalPetscVector::makePositionsInVector(const Base::Element* element)
    {
        logger.assert(element!=nullptr, "invalid element passed");
        //we need storage for the amount of basis functions to return
        std::vector<PetscInt> positions(element->getNrOfBasisFunctions());
        
        auto pos = positions.begin();
        
        std::size_t numElementBasisFuncs = element->getLocalNrOfBasisFunctions();
        //First step: construct ids for the functions of the current element itself
        for (std::size_t i = 0; i < numElementBasisFuncs; ++i)
        {
            *pos = i + startPositionsOfElementsInTheVector_[element->getID()];
            pos++;
        }
        
        //Push forward our iterator
        std::size_t numFaces = element->getPhysicalGeometry()->getNrOfFaces();
        for (std::size_t i = 0; i < numFaces; ++i)
        {
            std::size_t numFaceBasisFuncs = element->getFace(i)->getLocalNrOfBasisFunctions();
            for (std::size_t j = 0; j < numFaceBasisFuncs; ++j)
            {
                *pos = j + startPositionsOfFacesInTheVector_[element->getFace(i)->getID()];
                pos++;
            }
        }
        
        std::size_t numEdges = element->getNrOfEdges();
        for (std::size_t i = 0; i < numEdges; ++i)
        {
            std::size_t numEdgeBasisFuncs = element->getEdge(i)->getLocalNrOfBasisFunctions();
            for (std::size_t j = 0; j < numEdgeBasisFuncs; ++j)
            {
                *pos = j + startPositionsOfEdgesInTheVector_[element->getEdge(i)->getID()];
                pos++;
            }
        }
        
        std::size_t numNodes = element->getNrOfNodes();
        for (std::size_t i = 0; i < numNodes; ++i)
        {
            std::size_t numNodeBasisFuncs = element->getNode(i)->getLocalNrOfBasisFunctions();
            for (std::size_t j = 0; j < numNodeBasisFuncs; ++j)
            {
                *pos = j + startPositionsOfVerticesInTheVector_[element->getNode(i)->getID()];
                pos++;
            }
        }
        logger.assert(pos == positions.end(), "GlobalVector: did not process all elements correctly");
        return positions;
    }
    
    void GlobalPetscVector::reset()
    {
        int ierr = VecDestroy(&b_);
#ifdef HPGEM_USE_MPI
        std::size_t n = Base::MPIContainer::Instance().getNumProcessors();
        //offset by one to put a 0 in front
        std::vector<int> MPISendElementCounts(n+1,0), MPISendFaceCounts(n+1,0), MPISendEdgeCounts(n+1,0), MPISendNodeCounts(n+1,0);

        int rank = Base::MPIContainer::Instance().getProcessorID();

        MPISendElementCounts[rank+1] = theMesh_->getNumberOfElements();
        MPISendFaceCounts[rank+1] = theMesh_->getNumberOfFaces();
        MPISendEdgeCounts[rank+1] = theMesh_->getNumberOfEdges();
        MPISendNodeCounts[rank+1] = theMesh_->getNumberOfVertices();

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
        std::size_t totalNrOfDOF(0), DIM(theMesh_->dimension());
        startPositionsOfElementsInTheVector_.resize(theMesh_->getNumberOfElements(Base::IteratorType::GLOBAL));
        startPositionsOfFacesInTheVector_.resize(theMesh_->getNumberOfFaces(Base::IteratorType::GLOBAL));
        startPositionsOfEdgesInTheVector_.resize(theMesh_->getNumberOfEdges(Base::IteratorType::GLOBAL));
        startPositionsOfVerticesInTheVector_.resize(theMesh_->getNumberOfVertices(Base::IteratorType::GLOBAL));
        for (Base::Element* element : theMesh_->getElementsList())
        {
#ifdef HPGEM_USE_MPI
            *currentElementNumber=element->getID();
            *currentElementPosition=totalNrOfDOF;
            ++currentElementNumber;
            ++currentElementPosition;
#else
            startPositionsOfElementsInTheVector_[element->getID()] = totalNrOfDOF;
#endif
            totalNrOfDOF += element->getLocalNrOfBasisFunctions();
            for (std::size_t i = 0; i < element->getNrOfFaces(); ++i)
            {
                //faces at the boundary of the subdomain should also be added only once, so add them here is the left element of the face is in the subdomain
                if (element->getFace(i)->getFaceType() == Geometry::FaceType::SUBDOMAIN_BOUNDARY && element->getFace(i)->getPtrElementLeft() == element)
                {
#ifdef HPGEM_USE_MPI
                    *currentFaceNumber=element->getFace(i)->getID();
                    *currentFacePosition=totalNrOfDOF;
                    ++currentFaceNumber;
                    ++currentFacePosition;
#else
                    startPositionsOfFacesInTheVector_[element->getFace(i)->getID()] = totalNrOfDOF;
#endif
                    totalNrOfDOF += element->getFace(i)->getLocalNrOfBasisFunctions();
                }
                
            }
        }
        //DIM == 1 faces and nodes are the same entities, skip one of them
        if (DIM > 1)
        {
            for (Base::Face* face : theMesh_->getFacesList())
            {
                //skip faces at the subdomain boundary because we already treated them
                if (face->getFaceType() != Geometry::FaceType::SUBDOMAIN_BOUNDARY)
                {
#ifdef HPGEM_USE_MPI
                    *currentFaceNumber=face->getID();
                    *currentFacePosition=totalNrOfDOF;
                    ++currentFaceNumber;
                    ++currentFacePosition;
#else
                    startPositionsOfFacesInTheVector_[face->getID()] = totalNrOfDOF;
#endif
                    totalNrOfDOF += face->getLocalNrOfBasisFunctions();
                }
            }
        }
        for (Base::Edge* edge : theMesh_->getEdgesList())
        {
#ifdef HPGEM_USE_MPI
            *currentEdgeNumber=edge->getID();
            *currentEdgePosition=totalNrOfDOF;
            ++currentEdgeNumber;
            ++currentEdgePosition;
#else
            startPositionsOfEdgesInTheVector_[edge->getID()] = totalNrOfDOF;
#endif
            totalNrOfDOF += edge->getLocalNrOfBasisFunctions();
        }
        for (Base::Node* node : theMesh_->getVerticesList())
        {
#ifdef HPGEM_USE_MPI
            *currentNodeNumber=node->getID();
            *currentNodePosition=totalNrOfDOF;
            ++currentNodeNumber;
            ++currentNodePosition;
#else
            startPositionsOfVerticesInTheVector_[node->getID()] = totalNrOfDOF;
#endif
            totalNrOfDOF += node->getLocalNrOfBasisFunctions();
        }
        
#ifdef HPGEM_USE_MPI
        //offset by one to insert a zero at the front
        std::vector<std::size_t> cumulativeDOF(n + 1);
        cumulativeDOF[rank+1]=totalNrOfDOF;

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
            logger.assert(*currentNodeNumber != std::numeric_limits<std::size_t>::max(), "currentNodeNumber = -1");
            startPositionsOfVerticesInTheVector_[*currentNodeNumber]=*currentNodePosition+offset;
        }
#endif
        
        ierr = VecCreateMPI(PETSC_COMM_WORLD, totalNrOfDOF, PETSC_DETERMINE, &b_);
        CHKERRV(ierr);
    }
    
    void GlobalPetscVector::assemble()
    {
        reset();
        
        LinearAlgebra::NumericalVector elementVector;
        
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
        
        LinearAlgebra::NumericalVector faceVector;
        
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
    
    void GlobalPetscVector::constructFromTimeLevelData(std::size_t timelevel, std::size_t solutionVar)
    {
        reset();
        
        LinearAlgebra::NumericalVector elementData;
        for (Base::Element* element : theMesh_->getElementsList())
        {
            std::size_t numBasisFuncs = element->getNrOfBasisFunctions();
            std::vector<PetscInt> positions = makePositionsInVector(element);
            elementData.resize(numBasisFuncs);
            for (std::size_t i = 0; i < numBasisFuncs; ++i)
            {
                elementData[i] = element->getData(timelevel, solutionVar, i);
            }
            int ierr = VecSetValues(b_, numBasisFuncs, positions.data(), elementData.data(), INSERT_VALUES);
            CHKERRV(ierr);
        }
        
        int ierr = VecAssemblyBegin(b_);
        ierr = VecAssemblyEnd(b_);
        CHKERRV(ierr);
    }
    
    void GlobalPetscVector::writeTimeLevelData(std::size_t timeLevel, std::size_t variable)
    {
        PetscScalar *data;
        
        VecScatter scatter;
        Vec localB;
        //create a local vector...
        VecScatterCreateToAll(b_, &scatter, &localB);
        //we dont need the scatter context, just the vector
        VecScatterDestroy(&scatter);
        
        std::vector<PetscInt> positions;
        for (Base::Element* element : theMesh_->getElementsList())
        {
            std::vector<PetscInt> newPositions = makePositionsInVector(element);
            positions.reserve(positions.size() + newPositions.size());
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
        for (Base::MeshManipulator::ElementIterator it = theMesh_->elementColBegin(); it != theMesh_->elementColEnd(); ++it)
        {
            std::size_t numBasisFuns = (*it)->getNrOfBasisFunctions();
            LinearAlgebra::NumericalVector localData(&data[startPositionsOfElementsInTheVector_[(*it)->getID()]], (*it)->getLocalNrOfBasisFunctions());
            localData.resize(numBasisFuns);
            std::size_t runningTotal((*it)->getLocalNrOfBasisFunctions());
            if (theMesh_->dimension() > 1)
                for (std::size_t i = 0; i < (*it)->getPhysicalGeometry()->getNrOfFaces(); ++i)
                {
                    for (std::size_t j = 0; j < (*it)->getFace(i)->getLocalNrOfBasisFunctions(); ++j)
                    {
                        localData[runningTotal] = std::real(data[startPositionsOfFacesInTheVector_[(*it)->getFace(i)->getID()] + j]);
                        ++runningTotal;
                    }
                }
            for (std::size_t i = 0; i < (*it)->getNrOfEdges(); ++i)
            {
                for (std::size_t j = 0; j < (*it)->getEdge(i)->getLocalNrOfBasisFunctions(); ++j)
                {
                    localData[runningTotal] = std::real(data[startPositionsOfEdgesInTheVector_[(*it)->getEdge(i)->getID()] + j]);
                    ++runningTotal;
                }
            }
            for (std::size_t i = 0; i < (*it)->getNrOfNodes(); ++i)
            {
                for (std::size_t j = 0; j < (*it)->getNode(i)->getLocalNrOfBasisFunctions(); ++j)
                {
                    localData[runningTotal] = std::real(data[startPositionsOfVerticesInTheVector_[(*it)->getNode(i)->getID()] + j]);
                    ++runningTotal;
                }
            }
            (*it)->setTimeLevelData(timeLevel, variable, localData);
        }
        ierr = VecRestoreArray(localB, &data);
        VecDestroy(&localB);
        CHKERRV(ierr);
    }
#endif
}

