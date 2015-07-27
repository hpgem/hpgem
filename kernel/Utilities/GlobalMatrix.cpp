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
#include "GlobalMatrix.h"
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
#include "Base/Mesh.h"
#include "Logger.h"
#include <numeric>

namespace Utilities
{
    
    GlobalMatrix::GlobalMatrix(Base::MeshManipulatorBase* theMesh, int elementMatrixID, int faceMatrixID)
            : meshLevel_(-2), elementMatrixID_(elementMatrixID), faceMatrixID_(faceMatrixID), theMesh_(theMesh)
    {
        logger.assert(theMesh!=nullptr,"Invalid mesh passed");
    }
    
    void GlobalMatrix::getMatrixBCEntries(const Base::Face* face, std::size_t& numberOfEntries, std::vector<int>& entries)
    {
        logger.assert(face!=nullptr, "Invalid face passed");
        for(std::size_t index = 0; index < face->getPtrElementLeft()->getNumberOfUnknowns(); ++index)
        {
            std::size_t number = face->getLocalNumberOfBasisFunctions();
            numberOfEntries += number;
            for (std::size_t i = 0; i < number; ++i)
            {
                entries.push_back(startPositionsOfFacesInTheMatrix_[face->getID()] + i + index * number);
            }
            std::vector<std::size_t> nodeEntries = face->getPtrElementLeft()->getPhysicalGeometry()->getGlobalFaceNodeIndices(face->localFaceNumberLeft());
            std::vector<std::size_t> edgeIndex(2);
            for (std::size_t i = 0; i < face->getPtrElementLeft()->getNumberOfEdges(); ++i)
            {
                edgeIndex = face->getPtrElementLeft()->getReferenceGeometry()->getCodim2EntityLocalIndices(i);
                edgeIndex[0] = face->getPtrElementLeft()->getPhysicalGeometry()->getNodeIndex(edgeIndex[0]);
                edgeIndex[1] = face->getPtrElementLeft()->getPhysicalGeometry()->getNodeIndex(edgeIndex[1]);
                bool firstFound(false), secondFound(false);
                for (std::size_t j = 0; j < nodeEntries.size(); ++j)
                {
                    if (nodeEntries[j] == edgeIndex[0])
                        firstFound = true;
                    if (nodeEntries[j] == edgeIndex[1])
                        secondFound = true;
                }
                if (firstFound && secondFound)
                {
                    number = face->getPtrElementLeft()->getEdge(i)->getLocalNrOfBasisFunctions();
                    numberOfEntries += number;
                    for (std::size_t j = 0; j < number; ++j)
                    {
                        entries.push_back(startPositionsOfEdgesInTheMatrix_[face->getPtrElementLeft()->getEdge(i)->getID()] + j + index * number);
                    }
                }
            }
            nodeEntries = face->getPtrElementLeft()->getPhysicalGeometry()->getLocalFaceNodeIndices(face->localFaceNumberLeft());
            for (std::size_t i : nodeEntries)
            {
                number = face->getPtrElementLeft()->getNode(i)->getLocalNrOfBasisFunctions();
                numberOfEntries += number;
                for (std::size_t j = 0; j < number; ++j)
                    entries.push_back(startPositionsOfNodesInTheMatrix_[face->getPtrElementLeft()->getNode(i)->getID()] + j + index * number);
            }
        }
    }
#if defined(HPGEM_USE_PETSC) || defined(HPGEM_USE_COMPLEX_PETSC)
    
    GlobalPetscMatrix::GlobalPetscMatrix(Base::MeshManipulatorBase* theMesh, int elementMatrixID, int faceMatrixID)
            : GlobalMatrix(theMesh, elementMatrixID, faceMatrixID)
    {
        logger.assert(theMesh!=nullptr, "Invalid mesh passed");
        PetscBool petscRuns;
        PetscInitialized(&petscRuns);
        logger.assert(petscRuns == PETSC_TRUE, "Early call, firstly the command line arguments should be parsed");
        //temporary
        MatCreateSeqAIJ(PETSC_COMM_SELF, 1, 1, 1, PETSC_NULL, &A_);
        
        reAssemble();
    }
    
    GlobalPetscMatrix::~GlobalPetscMatrix()
    {
        int ierr = MatDestroy(&A_);
        //giving error about Petsc has generated inconsistent data and likely memory corruption in heap
        CHKERRV(ierr);
    }
    
    GlobalPetscMatrix::operator Mat()
    {
        if(HPGEM_LOGLEVEL==Log::DEBUG)
        {
            MatView(A_, PETSC_VIEWER_STDOUT_WORLD);
        }
        return A_;
    }
    
    std::vector<PetscInt> GlobalPetscMatrix::makePositionsInMatrix(const Base::Element* element)
    {
        logger.assert(element!=nullptr, "Invalid element passed");
        //we need storage for the amount of basis functions to return
        std::vector<PetscInt> positions(element->getNrOfBasisFunctions() * element->getNrOfUnknowns());
        
        auto pos = positions.begin();
        for(std::size_t index = 0; index < element->getNrOfUnknowns(); ++index)
        {
            std::size_t numElementBasisFuncs = element->getLocalNrOfBasisFunctions();
            //First step: construct ids for the functions of the current element itself
            for (std::size_t i = 0; i < numElementBasisFuncs; ++i)
            {
                *pos = i + startPositionsOfElementsInTheMatrix_[element->getID()] + index * numElementBasisFuncs;
                pos++;
            }

            //Push forward our iterator
            std::size_t numFaces = element->getPhysicalGeometry()->getNrOfFaces();
            for (std::size_t i = 0; i < numFaces; ++i)
            {
                std::size_t numFaceBasisFuncs = element->getFace(i)->getLocalNrOfBasisFunctions();
                for (std::size_t j = 0; j < numFaceBasisFuncs; ++j)
                {
                    *pos = j + startPositionsOfFacesInTheMatrix_[element->getFace(i)->getID()] + index * numFaceBasisFuncs;
                    pos++;
                }
            }

            std::size_t numEdges = element->getNrOfEdges();
            for (std::size_t i = 0; i < numEdges; ++i)
            {
                std::size_t numEdgeBasisFuncs = element->getEdge(i)->getLocalNrOfBasisFunctions();
                for (std::size_t j = 0; j < numEdgeBasisFuncs; ++j)
                {
                    *pos = j + startPositionsOfEdgesInTheMatrix_[element->getEdge(i)->getID()] + index * numEdgeBasisFuncs;
                    pos++;
                }
            }

            std::size_t numNodes = element->getNrOfNodes();
            for (std::size_t i = 0; i < numNodes; ++i)
            {
                std::size_t numNodeBasisFuncs = element->getNode(i)->getLocalNrOfBasisFunctions();
                for (std::size_t j = 0; j < numNodeBasisFuncs; ++j)
                {
                    *pos = j + startPositionsOfNodesInTheMatrix_[element->getNode(i)->getID()] + index * numNodeBasisFuncs;
                    pos++;
                }
            }
        }
        logger.assert(pos == positions.end(), "Not all positions are processed.");
        return positions;
    }
    
    void GlobalPetscMatrix::reset()
    {
        int ierr = MatZeroEntries(A_);
        CHKERRV(ierr);
        
        LinearAlgebra::MiddleSizeMatrix elementMatrix;
        
        if (elementMatrixID_ >= 0)
        {
            for (Base::Element* element : theMesh_->getElementsList())
            {
                std::vector<PetscInt> positions = makePositionsInMatrix(element);
                elementMatrix = element->getElementMatrix(elementMatrixID_);
                ierr = MatSetValues(A_, positions.size(), positions.data(), positions.size(), positions.data(), elementMatrix.data(), ADD_VALUES);
                CHKERRV(ierr);
            }
        }
        
        LinearAlgebra::MiddleSizeMatrix faceMatrix;
        
        if (faceMatrixID_ >= 0)
        {
            for (Base::Face* face : theMesh_->getFacesList())
            {
                std::vector<PetscInt> positions = makePositionsInMatrix(face->getPtrElementLeft());
                if (face->isInternal())
                {
                    std::vector<PetscInt> rightPositions = makePositionsInMatrix(face->getPtrElementRight());
                    positions.reserve(positions.size() + rightPositions.size());
                    for (auto& a : rightPositions)
                    {
                        positions.push_back(a);
                    }
                }
                faceMatrix = face->getFaceMatrixMatrix(faceMatrixID_);
                ierr = MatSetValues(A_, positions.size(), positions.data(), positions.size(), positions.data(), faceMatrix.data(), ADD_VALUES);
                CHKERRV(ierr);
            }
        }
        
        ierr = MatAssemblyBegin(A_, MAT_FINAL_ASSEMBLY);
        ierr = MatAssemblyEnd(A_, MAT_FINAL_ASSEMBLY);
        
        CHKERRV(ierr);
    }
    
    ///\todo figure out a nice way to keep local data local
    void GlobalPetscMatrix::reAssemble()
    {
        MatDestroy(&A_);
#ifdef HPGEM_USE_MPI
        std::size_t n = Base::MPIContainer::Instance().getNumProcessors();
        //offset by one to put a 0 in front (type int, because MPI wants it to be int)
        std::vector<int> MPISendElementCounts(n+1,0), MPISendFaceCounts(n+1,0), MPISendEdgeCounts(n+1,0), MPISendNodeCounts(n+1,0);

        int rank = Base::MPIContainer::Instance().getProcessorID();

        MPISendElementCounts[rank+1] = theMesh_->getNumberOfElements();
        MPISendFaceCounts[rank+1] = theMesh_->getNumberOfFaces();
        MPISendEdgeCounts[rank+1] = theMesh_->getNumberOfEdges();
        MPISendNodeCounts[rank+1] = theMesh_->getNumberOfNodes();

        //tell the rest of the processes how much info the others have
        MPI::Intracomm& comm = Base::MPIContainer::Instance().getComm();
        comm.Allgather(MPI_IN_PLACE,1,Base::Detail::toMPIType(*MPISendElementCounts.data()),MPISendElementCounts.data()+1,1,Base::Detail::toMPIType(*MPISendElementCounts.data()));
        comm.Allgather(MPI_IN_PLACE,1,Base::Detail::toMPIType(*MPISendElementCounts.data()),MPISendFaceCounts.data()+1,1,Base::Detail::toMPIType(*MPISendElementCounts.data()));
        comm.Allgather(MPI_IN_PLACE,1,Base::Detail::toMPIType(*MPISendElementCounts.data()),MPISendEdgeCounts.data()+1,1,Base::Detail::toMPIType(*MPISendElementCounts.data()));
        comm.Allgather(MPI_IN_PLACE,1,Base::Detail::toMPIType(*MPISendElementCounts.data()),MPISendNodeCounts.data()+1,1,Base::Detail::toMPIType(*MPISendElementCounts.data()));

        std::vector<int> MPISendElementStarts(n+1), MPISendFaceStarts(n+1), MPISendEdgeStarts(n+1), MPISendNodeStarts(n+1);

        std::partial_sum(MPISendElementCounts.begin(), MPISendElementCounts.end(), MPISendElementStarts.begin());
        std::partial_sum(MPISendFaceCounts.begin(), MPISendFaceCounts.end(), MPISendFaceStarts.begin());
        std::partial_sum(MPISendEdgeCounts.begin(), MPISendEdgeCounts.end(), MPISendEdgeStarts.begin());
        std::partial_sum(MPISendNodeCounts.begin(), MPISendNodeCounts.end(), MPISendNodeStarts.begin());

        logger.assert(MPISendElementStarts.back()==theMesh_->getNumberOfElements(Base::IteratorType::GLOBAL),
                "MPI does not see the right amount of elements");
        logger.assert(MPISendFaceStarts.back()>=theMesh_->getNumberOfFaces(Base::IteratorType::GLOBAL),
                "MPI does not see the right amount of faces");
        logger.assert(MPISendEdgeStarts.back()==theMesh_->getNumberOfEdges(Base::IteratorType::GLOBAL),
                "MPI does not see the right amount of edges");
        logger.assert(MPISendNodeStarts.back()==theMesh_->getNumberOfNodes(Base::IteratorType::GLOBAL),
                "MPI does not see the right amount of nodes");

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
        startPositionsOfElementsInTheMatrix_.assign(theMesh_->getNumberOfElements(Base::IteratorType::GLOBAL), -1);
        startPositionsOfFacesInTheMatrix_.assign(theMesh_->getNumberOfFaces(Base::IteratorType::GLOBAL), -1);
        startPositionsOfEdgesInTheMatrix_.assign(theMesh_->getNumberOfEdges(Base::IteratorType::GLOBAL), -1);
        startPositionsOfNodesInTheMatrix_.assign(theMesh_->getNumberOfNodes(Base::IteratorType::GLOBAL), -1);
        for (Base::Element* element : theMesh_->getElementsList())
        {
#ifdef HPGEM_USE_MPI
            *currentElementNumber=element->getID();
            *currentElementPosition=totalNrOfDOF;
            ++currentElementNumber;
            ++currentElementPosition;
#else
            logger.assert(startPositionsOfElementsInTheMatrix_[element->getID()] == -1, "Duplicate element detected");
            startPositionsOfElementsInTheMatrix_[element->getID()] = totalNrOfDOF;
#endif
            totalNrOfDOF += element->getLocalNrOfBasisFunctions() * element->getNrOfUnknowns();
            for (std::size_t i = 0; i < element->getNrOfFaces(); ++i)
            {
                //faces at the boundary of the subdomain should also be added only once, so add them here if the left element of the face is in the subdomain
                if ((element->getFace(i)->getFaceType() == Geometry::FaceType::SUBDOMAIN_BOUNDARY || element->getFace(i)->getFaceType() == Geometry::FaceType::PERIODIC_SUBDOMAIN_BC)
                        && element->getFace(i)->getPtrElementLeft() == element)
                {
#ifdef HPGEM_USE_MPI
                    *currentFaceNumber=element->getFace(i)->getID();
                    *currentFacePosition=totalNrOfDOF;
                    ++currentFaceNumber;
                    ++currentFacePosition;
#else
                    logger.assert(startPositionsOfFacesInTheMatrix_[element->getFace(i)->getID()] == -1, "Duplicate face detected");
                    startPositionsOfFacesInTheMatrix_[element->getFace(i)->getID()] = totalNrOfDOF;
#endif
                    totalNrOfDOF += element->getFace(i)->getLocalNrOfBasisFunctions() * element->getNrOfUnknowns();
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
                *currentFacePosition=totalNrOfDOF;
                ++currentFaceNumber;
                ++currentFacePosition;
#else
                logger.assert(startPositionsOfFacesInTheMatrix_[face->getID()] == -1, "Duplicate face detected");
                startPositionsOfFacesInTheMatrix_[face->getID()] = totalNrOfDOF;
#endif
                totalNrOfDOF += face->getLocalNrOfBasisFunctions() * face->getPtrElementLeft()->getNrOfUnknowns();
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
            logger.assert(startPositionsOfEdgesInTheMatrix_[edge->getID()] == -1, "Duplicate edge detected");
            startPositionsOfEdgesInTheMatrix_[edge->getID()] = totalNrOfDOF;
#endif
            totalNrOfDOF += edge->getLocalNrOfBasisFunctions() * edge->getElement(0)->getNrOfUnknowns();
        }
        //DIM == 1 faces and nodes are the same entities, skip one of them
        if (DIM > 1)
        {
            for (Base::Node* node : theMesh_->getNodesList())
            {
#ifdef HPGEM_USE_MPI
                *currentNodeNumber=node->getID();
                *currentNodePosition=totalNrOfDOF;
                ++currentNodeNumber;
                ++currentNodePosition;
#else
                logger.assert(startPositionsOfNodesInTheMatrix_[node->getID()] == -1, "Duplicate node detected");
                startPositionsOfNodesInTheMatrix_[node->getID()] = totalNrOfDOF;
#endif
                totalNrOfDOF += node->getLocalNrOfBasisFunctions() * node->getElement(0)->getNrOfUnknowns();
            }
        }
        
#ifdef HPGEM_USE_MPI
        
        logger.assert(currentElementNumber == MPISendElementNumbers.begin() + MPISendElementStarts[rank+1],
                "MPI is not at the correct element.");
        logger.assert(currentFaceNumber <= MPISendFaceNumbers.begin() + MPISendFaceStarts[rank+1],
                "MPI is not at the correct face.");
        logger.assert(currentEdgeNumber == MPISendEdgeNumbers.begin() + MPISendEdgeStarts[rank+1],
                "MPI is not at the correct edge.");
        logger.assert(currentNodeNumber == MPISendNodeNumbers.begin() + MPISendNodeStarts[rank+1] || DIM == 1,
                "MPI is not at the correct node.");

        logger.assert(currentElementPosition == MPISendElementPositions.begin() + MPISendElementStarts[rank+1],
                "MPI is not at the correct element.");
        logger.assert(currentFacePosition <= MPISendFacePositions.begin() + MPISendFaceStarts[rank+1],
                "MPI is not at the correct face.");
        logger.assert(currentEdgePosition == MPISendEdgePositions.begin() + MPISendEdgeStarts[rank+1],
                "MPI is not at the correct edge.");
        logger.assert(currentNodePosition == MPISendNodePositions.begin() + MPISendNodeStarts[rank+1] || DIM == 1,
                "MPI is not at the correct node.");

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

        std::partial_sum(cumulativeDOF.begin(),cumulativeDOF.end(),cumulativeDOF.begin());

        //and unpack the information
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
            logger.assert(startPositionsOfElementsInTheMatrix_[*currentElementNumber] == -1, "Duplicate element detected");
            startPositionsOfElementsInTheMatrix_[*currentElementNumber]=*currentElementPosition+offset;
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
            {
                logger.assert(startPositionsOfFacesInTheMatrix_[*currentFaceNumber] == -1, "Duplicate face detected");
                startPositionsOfFacesInTheMatrix_[*currentFaceNumber]=*currentFacePosition+offset;
            }
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
            logger.assert(*currentEdgeNumber != std::numeric_limits<std::size_t>::max(), "currentEdgeNumber = -1");
            logger.assert(startPositionsOfEdgesInTheMatrix_[*currentEdgeNumber] == -1, "Duplicate edge detected");
            startPositionsOfEdgesInTheMatrix_[*currentEdgeNumber]=*currentEdgePosition+offset;
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
            {
                logger.assert(startPositionsOfNodesInTheMatrix_[*currentNodeNumber] == -1, "Duplicate node detected");
                startPositionsOfNodesInTheMatrix_[*currentNodeNumber]=*currentNodePosition+offset;
            }
        }
        std::size_t MPIOffset = cumulativeDOF[rank];
        std::size_t end = cumulativeDOF[n] + 1;
#else
        std::size_t MPIOffset = 0;
        std::size_t end = totalNrOfDOF + 1;
#endif
        logger.assert(*std::min_element(startPositionsOfElementsInTheMatrix_.begin(), startPositionsOfElementsInTheMatrix_.end()) > -1, "Missing an element index");
        logger.assert(*std::min_element(startPositionsOfFacesInTheMatrix_.begin(), startPositionsOfFacesInTheMatrix_.end()) > -1, "Missing a face index");
        logger.assert(DIM < 3 || *std::min_element(startPositionsOfEdgesInTheMatrix_.begin(), startPositionsOfEdgesInTheMatrix_.end()) > -1, "Missing an edge index");
        logger.assert(DIM < 2 || *std::min_element(startPositionsOfNodesInTheMatrix_.begin(), startPositionsOfNodesInTheMatrix_.end()) > -1, "Missing node index %", std::min_element(startPositionsOfElementsInTheMatrix_.begin(), startPositionsOfElementsInTheMatrix_.end()) - startPositionsOfElementsInTheMatrix_.begin());
        //could be equal to totalNrOfDOF if there are no degrees of freedom accosiated with the element/face/w.e.
        logger.assert(*std::max_element(startPositionsOfElementsInTheMatrix_.begin(), startPositionsOfElementsInTheMatrix_.end()) < end, "Start of entries for element % out of bounds (%, with only % total entries)", std::max_element(startPositionsOfElementsInTheMatrix_.begin(), startPositionsOfElementsInTheMatrix_.end()) - startPositionsOfElementsInTheMatrix_.begin(), *std::max_element(startPositionsOfElementsInTheMatrix_.begin(), startPositionsOfElementsInTheMatrix_.end()), end - 1);
        logger.assert(*std::max_element(startPositionsOfFacesInTheMatrix_.begin(), startPositionsOfFacesInTheMatrix_.end()) < end, "Missing a face index");
        logger.assert(DIM < 3 || *std::max_element(startPositionsOfEdgesInTheMatrix_.begin(), startPositionsOfEdgesInTheMatrix_.end()) < end, "Missing an edge index");
        logger.assert(DIM < 2 || *std::max_element(startPositionsOfNodesInTheMatrix_.begin(), startPositionsOfNodesInTheMatrix_.end()) < end, "Missing a node index");
        
        //now construct the only bit of data where PETSc expects a local numbering...
        std::vector<PetscInt> numberOfPositionsPerRow(totalNrOfDOF, 0);
        std::vector<PetscInt> offDiagonalPositionsPerRow(totalNrOfDOF, 0);
        
        for (Base::Element* element : theMesh_->getElementsList())
        {
            for (int j = 0; j < element->getLocalNrOfBasisFunctions() * element->getNrOfUnknowns(); ++j)
            {
                numberOfPositionsPerRow[startPositionsOfElementsInTheMatrix_[element->getID()] + j - MPIOffset] += element->getNrOfBasisFunctions() * element->getNrOfUnknowns();
            }
            for (int i = 0; i < element->getReferenceGeometry()->getNrOfCodim1Entities(); ++i)
            {
                //conforming contributions
                for (int j = 0; j < element->getFace(i)->getLocalNrOfBasisFunctions() * element->getNrOfUnknowns(); ++j)
                {
                    numberOfPositionsPerRow[startPositionsOfFacesInTheMatrix_[element->getFace(i)->getID()] + j - MPIOffset] += element->getNrOfBasisFunctions() * element->getNrOfUnknowns();
                }
                
                //flux term for DG
                for (int j = 0; j < element->getLocalNrOfBasisFunctions() * element->getNrOfUnknowns(); ++j)
                {
                    if (element->getFace(i)->getFaceType() == Geometry::FaceType::SUBDOMAIN_BOUNDARY || element->getFace(i)->getFaceType() == Geometry::FaceType::PERIODIC_SUBDOMAIN_BC)
                    {
                        offDiagonalPositionsPerRow[startPositionsOfElementsInTheMatrix_[element->getID()] + j - MPIOffset] += element->getNrOfBasisFunctions() * element->getNrOfUnknowns();
                    }
                    else if (element->getFace(i)->isInternal())
                    {
                        numberOfPositionsPerRow[startPositionsOfElementsInTheMatrix_[element->getID()] + j - MPIOffset] += element->getNrOfBasisFunctions() * element->getNrOfUnknowns();
                    }
                }
            }
            for (int i = 0; i < element->getNrOfEdges(); ++i)
            {
                for (int j = 0; j < element->getEdge(i)->getLocalNrOfBasisFunctions() * element->getNrOfUnknowns(); ++j)
                {
                    if(startPositionsOfEdgesInTheMatrix_[element->getEdge(i)->getID()] + j - MPIOffset < totalNrOfDOF && -1 < startPositionsOfEdgesInTheMatrix_[element->getEdge(i)->getID()] + j - MPIOffset)
                    {
                        numberOfPositionsPerRow[startPositionsOfEdgesInTheMatrix_[element->getEdge(i)->getID()] + j - MPIOffset] += element->getNrOfBasisFunctions() * element->getNrOfUnknowns();
                    }
                }
            }
            for (int i = 0; i < element->getNrOfNodes(); ++i)
            {
                for (int j = 0; j < element->getNode(i)->getLocalNrOfBasisFunctions() * element->getNrOfUnknowns(); ++j)
                {
                    if(startPositionsOfNodesInTheMatrix_[element->getNode(i)->getID()] + j - MPIOffset < totalNrOfDOF && -1 < startPositionsOfNodesInTheMatrix_[element->getNode(i)->getID()] + j - MPIOffset)
                    numberOfPositionsPerRow[startPositionsOfNodesInTheMatrix_[element->getNode(i)->getID()] + j - MPIOffset] += element->getNrOfBasisFunctions() * element->getNrOfUnknowns();
                }
            }
        }
        
        for (int i = 0; i < totalNrOfDOF; ++i)
        {
            if (numberOfPositionsPerRow[i] > totalNrOfDOF)
            {
                numberOfPositionsPerRow[i] = totalNrOfDOF; //a row cant have more nonzero entries than the number of columns
            }
        }
        
        int ierr = MatCreateAIJ(PETSC_COMM_WORLD, totalNrOfDOF, totalNrOfDOF, PETSC_DETERMINE, PETSC_DETERMINE, -1, numberOfPositionsPerRow.data(), 0, offDiagonalPositionsPerRow.data(), &A_);
        CHKERRV(ierr);
        MatSetOption(A_, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE); //the estimate is known to be wrong for conforming cases
        ierr = MatSetUp(A_);
        CHKERRV(ierr);
        reset();
    }
#endif
}

