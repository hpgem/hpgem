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
                    number = face->getPtrElementLeft()->getEdge(i)->getLocalNumberOfBasisFunctions();
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
                number = face->getPtrElementLeft()->getNode(i)->getLocalNumberOfBasisFunctions();
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
        std::vector<PetscInt> positions(element->getNumberOfBasisFunctions() * element->getNumberOfUnknowns());
        
        auto pos = positions.begin();
        for(std::size_t index = 0; index < element->getNumberOfUnknowns(); ++index)
        {
            std::size_t numberOfElementBasisFunctions = element->getLocalNumberOfBasisFunctions();
            //First step: construct ids for the functions of the current element itself
            for (std::size_t i = 0; i < numberOfElementBasisFunctions; ++i)
            {
                *pos = i + startPositionsOfElementsInTheMatrix_[element->getID()] + index * numberOfElementBasisFunctions;
                pos++;
            }

            //Push forward our iterator
            std::size_t numberOfFaces = element->getPhysicalGeometry()->getNumberOfFaces();
            for (std::size_t i = 0; i < numberOfFaces; ++i)
            {
                std::size_t numberOfFaceBasisFunctions = element->getFace(i)->getLocalNumberOfBasisFunctions();
                for (std::size_t j = 0; j < numberOfFaceBasisFunctions; ++j)
                {
                    *pos = j + startPositionsOfFacesInTheMatrix_[element->getFace(i)->getID()] + index * numberOfFaceBasisFunctions;
                    pos++;
                }
            }

            std::size_t numberOfEdges = element->getNumberOfEdges();
            for (std::size_t i = 0; i < numberOfEdges; ++i)
            {
                std::size_t numberOfEdgeBasisFunctions = element->getEdge(i)->getLocalNumberOfBasisFunctions();
                for (std::size_t j = 0; j < numberOfEdgeBasisFunctions; ++j)
                {
                    *pos = j + startPositionsOfEdgesInTheMatrix_[element->getEdge(i)->getID()] + index * numberOfEdgeBasisFunctions;
                    pos++;
                }
            }

            std::size_t numberOfNodes = element->getNumberOfNodes();
            for (std::size_t i = 0; i < numberOfNodes; ++i)
            {
                std::size_t numberOfNodeBasisFunctions = element->getNode(i)->getLocalNumberOfBasisFunctions();
                for (std::size_t j = 0; j < numberOfNodeBasisFunctions; ++j)
                {
                    *pos = j + startPositionsOfNodesInTheMatrix_[element->getNode(i)->getID()] + index * numberOfNodeBasisFunctions;
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
                //work-around: both subdomains have the boundary face so by default it is added twice, but it should only be added once
                if(face->getFaceType() == Geometry::FaceType::SUBDOMAIN_BOUNDARY || face->getFaceType() == Geometry::FaceType::PERIODIC_SUBDOMAIN_BC)
                {
                    faceMatrix *= 0.5;
                }
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
        std::size_t n = Base::MPIContainer::Instance().getNumberOfProcessors();
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
        std::size_t totalNumberOfDOF(0), DIM(theMesh_->dimension());
        for (Base::Element* element : theMesh_->getElementsList())
        {
#ifdef HPGEM_USE_MPI
            *currentElementNumber=element->getID();
            *currentElementPosition=totalNumberOfDOF;
            ++currentElementNumber;
            ++currentElementPosition;
#else
            startPositionsOfElementsInTheMatrix_[element->getID()] = totalNumberOfDOF;
#endif
            totalNumberOfDOF += element->getLocalNumberOfBasisFunctions() * element->getNumberOfUnknowns();
            for (std::size_t i = 0; i < element->getNumberOfFaces(); ++i)
            {
                //faces at the boundary of the subdomain should also be added only once, so add them here if the left element of the face is in the subdomain
                if ((element->getFace(i)->getFaceType() == Geometry::FaceType::SUBDOMAIN_BOUNDARY || element->getFace(i)->getFaceType() == Geometry::FaceType::PERIODIC_SUBDOMAIN_BC)
                        && element->getFace(i)->getPtrElementLeft() == element)
                {
#ifdef HPGEM_USE_MPI
                    *currentFaceNumber=element->getFace(i)->getID();
                    *currentFacePosition=totalNumberOfDOF;
                    ++currentFaceNumber;
                    ++currentFacePosition;
#else
                    startPositionsOfFacesInTheMatrix_[element->getFace(i)->getID()] = totalNumberOfDOF;
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
                startPositionsOfFacesInTheMatrix_[face->getID()] = totalNumberOfDOF;
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
            startPositionsOfEdgesInTheMatrix_[edge->getID()] = totalNumberOfDOF;
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
                startPositionsOfNodesInTheMatrix_[node->getID()] = totalNumberOfDOF;
#endif
                totalNumberOfDOF += node->getLocalNumberOfBasisFunctions() * node->getElement(0)->getNumberOfUnknowns();
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
                startPositionsOfNodesInTheMatrix_[*currentNodeNumber]=*currentNodePosition+offset;
            }
        }
        std::size_t MPIOffset = cumulativeDOF[rank];
        std::size_t end = cumulativeDOF[n] + 1;
#else
        std::size_t MPIOffset = 0;
        std::size_t end = totalNumberOfDOF + 1;
#endif
        logger.assert(startPositionsOfElementsInTheMatrix_.size() == theMesh_->getNumberOfElements(Base::IteratorType::GLOBAL), "Missing an element");
        logger.assert(startPositionsOfFacesInTheMatrix_.size() == theMesh_->getNumberOfFaces(Base::IteratorType::GLOBAL), "Missing a face");
        logger.assert(DIM < 3 || startPositionsOfEdgesInTheMatrix_.size() == theMesh_->getNumberOfEdges(Base::IteratorType::GLOBAL), "Missing a edge");
        logger.assert(DIM < 2 || startPositionsOfNodesInTheMatrix_.size() == theMesh_->getNumberOfNodes(Base::IteratorType::GLOBAL), "Missing a node");
        
        //now construct the only bit of data where PETSc expects a local numbering...
        std::vector<PetscInt> numberOfPositionsPerRow(totalNumberOfDOF, 0);
        std::vector<PetscInt> offDiagonalPositionsPerRow(totalNumberOfDOF, 0);
        
        for (Base::Element* element : theMesh_->getElementsList())
        {
            for (int j = 0; j < element->getLocalNumberOfBasisFunctions() * element->getNumberOfUnknowns(); ++j)
            {
                numberOfPositionsPerRow[startPositionsOfElementsInTheMatrix_[element->getID()] + j - MPIOffset] += element->getNumberOfBasisFunctions() * element->getNumberOfUnknowns();
            }
            for (int i = 0; i < element->getReferenceGeometry()->getNumberOfCodim1Entities(); ++i)
            {
                //conforming contributions
                for (int j = 0; j < element->getFace(i)->getLocalNumberOfBasisFunctions() * element->getNumberOfUnknowns(); ++j)
                {
                    if(startPositionsOfFacesInTheMatrix_[element->getFace(i)->getID()] + j < totalNumberOfDOF + MPIOffset && MPIOffset < startPositionsOfFacesInTheMatrix_[element->getFace(i)->getID()] + j + 1)
                    {
                        numberOfPositionsPerRow[startPositionsOfFacesInTheMatrix_[element->getFace(i)->getID()] + j - MPIOffset] += element->getNumberOfBasisFunctions() * element->getNumberOfUnknowns();
                    }
                }
            }
            for (int i = 0; i < element->getNumberOfEdges(); ++i)
            {
                for (int j = 0; j < element->getEdge(i)->getLocalNumberOfBasisFunctions() * element->getNumberOfUnknowns(); ++j)
                {
                    if(startPositionsOfEdgesInTheMatrix_[element->getEdge(i)->getID()] + j < totalNumberOfDOF + MPIOffset && MPIOffset < startPositionsOfEdgesInTheMatrix_[element->getEdge(i)->getID()] + j + 1)
                    {
                        numberOfPositionsPerRow[startPositionsOfEdgesInTheMatrix_[element->getEdge(i)->getID()] + j - MPIOffset] += element->getNumberOfBasisFunctions() * element->getNumberOfUnknowns();
                    }
                }
            }
            for (int i = 0; i < element->getNumberOfNodes(); ++i)
            {
                for (int j = 0; j < element->getNode(i)->getLocalNumberOfBasisFunctions() * element->getNumberOfUnknowns(); ++j)
                {
                    if(startPositionsOfNodesInTheMatrix_[element->getNode(i)->getID()] + j < totalNumberOfDOF + MPIOffset && MPIOffset < startPositionsOfNodesInTheMatrix_[element->getNode(i)->getID()] + j + 1)
                    {
                        numberOfPositionsPerRow[startPositionsOfNodesInTheMatrix_[element->getNode(i)->getID()] + j - MPIOffset] += element->getNumberOfBasisFunctions() * element->getNumberOfUnknowns();
                    }
                }
            }
        }
        
        for(Base::Face* face : theMesh_->getFacesList())
        {
            if(face->isInternal())
            {
                std::vector<int>& changeVec = (face->getFaceType() == Geometry::FaceType::SUBDOMAIN_BOUNDARY || face->getFaceType() == Geometry::FaceType::PERIODIC_SUBDOMAIN_BC) ? offDiagonalPositionsPerRow : numberOfPositionsPerRow;
                std::size_t nDuplicates = 0;
                std::vector<int> duplicates;
                getMatrixBCEntries(face, nDuplicates, duplicates);
                for(Base::Element* element : {face->getPtrElementLeft(), face->getPtrElementRight()})
                {
                    for (int j = 0; j < element->getLocalNumberOfBasisFunctions() * element->getNumberOfUnknowns(); ++j)
                    {
                        if(startPositionsOfElementsInTheMatrix_[element->getID()] + j < totalNumberOfDOF + MPIOffset && MPIOffset < startPositionsOfElementsInTheMatrix_[element->getID()] + j + 1)
                        {
                            changeVec[startPositionsOfElementsInTheMatrix_[element->getID()] + j - MPIOffset] += element->getNumberOfBasisFunctions() * element->getNumberOfUnknowns() - nDuplicates;
                        }
                    }
                    for (int i = 0; i < element->getReferenceGeometry()->getNumberOfCodim1Entities(); ++i)
                    {
                        //conforming contributions
                        for (int j = 0; j < element->getFace(i)->getLocalNumberOfBasisFunctions() * element->getNumberOfUnknowns(); ++j)
                        {
                            if(startPositionsOfFacesInTheMatrix_[element->getFace(i)->getID()] + j < totalNumberOfDOF + MPIOffset && MPIOffset < startPositionsOfFacesInTheMatrix_[element->getFace(i)->getID()] + j + 1)
                            {
                                changeVec[startPositionsOfFacesInTheMatrix_[element->getFace(i)->getID()] + j - MPIOffset] += element->getNumberOfBasisFunctions() * element->getNumberOfUnknowns() - nDuplicates;
                            }
                        }
                    }
                    for (int i = 0; i < element->getNumberOfEdges(); ++i)
                    {
                        for (int j = 0; j < element->getEdge(i)->getLocalNumberOfBasisFunctions() * element->getNumberOfUnknowns(); ++j)
                        {
                            if(startPositionsOfEdgesInTheMatrix_[element->getEdge(i)->getID()] + j < totalNumberOfDOF + MPIOffset && MPIOffset < startPositionsOfEdgesInTheMatrix_[element->getEdge(i)->getID()] + j + 1)
                            {
                                changeVec[startPositionsOfEdgesInTheMatrix_[element->getEdge(i)->getID()] + j - MPIOffset] += element->getNumberOfBasisFunctions() * element->getNumberOfUnknowns() - nDuplicates;
                            }
                        }
                    }
                    for (int i = 0; i < element->getNumberOfNodes(); ++i)
                    {
                        for (int j = 0; j < element->getNode(i)->getLocalNumberOfBasisFunctions() * element->getNumberOfUnknowns(); ++j)
                        {
                            if(startPositionsOfNodesInTheMatrix_[element->getNode(i)->getID()] + j < totalNumberOfDOF + MPIOffset && MPIOffset < startPositionsOfNodesInTheMatrix_[element->getNode(i)->getID()] + j + 1)
                            {
                                changeVec[startPositionsOfNodesInTheMatrix_[element->getNode(i)->getID()] + j - MPIOffset] += element->getNumberOfBasisFunctions() * element->getNumberOfUnknowns() - nDuplicates;
                            }
                        }
                    }
                }
                for(std::size_t i : duplicates)
                {
                    if(i < totalNumberOfDOF + MPIOffset && MPIOffset < i + 1)
                    {
                        changeVec[i-MPIOffset] -= face->getPtrElementLeft()->getNumberOfBasisFunctions() * face->getPtrElementLeft()->getNumberOfUnknowns() - nDuplicates;
                    }
                }
            }
        }

        for (int i = 0; i < totalNumberOfDOF; ++i)
        {
            if (numberOfPositionsPerRow[i] > totalNumberOfDOF)
            {
                numberOfPositionsPerRow[i] = totalNumberOfDOF; //a row cant have more nonzero entries than the number of columns
            }
        }
        
        int ierr = MatCreateAIJ(PETSC_COMM_WORLD, totalNumberOfDOF, totalNumberOfDOF, PETSC_DETERMINE, PETSC_DETERMINE, -1, numberOfPositionsPerRow.data(), 0, offDiagonalPositionsPerRow.data(), &A_);
        CHKERRV(ierr);
        MatSetOption(A_, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE); //the estimate is known to be wrong for mixed element cases and conforming parallel cases
        ierr = MatSetUp(A_);
        CHKERRV(ierr);
        reset();
    }
#endif
}

