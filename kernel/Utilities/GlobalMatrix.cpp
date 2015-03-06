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

#include "Base/MpiContainer.hpp"
#include "GlobalMatrix.hpp"
#include <vector>
#include "Base/MeshManipulator.hpp"
#include "Base/Edge.hpp"
#include "Base/Face.hpp"
#include "Base/Element.hpp"
#include "Base/ElementCacheData.hpp"
#include "LinearAlgebra/NumericalVector.hpp"
#include "Base/FaceCacheData.hpp"
#include "Geometry/PointPhysical.hpp"
#include "Geometry/PhysicalGeometry.hpp"
#include "Geometry/ReferenceGeometry.hpp"
#include "Geometry/PointReference.hpp"
#include "Base/Norm2.hpp"
#include "Base/Mesh.hpp"
#include "cassert"
#include <numeric>

namespace Utilities {

    GlobalMatrix::GlobalMatrix(Base::MeshManipulator* theMesh, int elementMatrixID, int faceMatrixID) :
    meshLevel_(-2), theMesh_(theMesh), elementMatrixID_(elementMatrixID), faceMatrixID_(faceMatrixID) {
    }

    void GlobalMatrix::getMatrixBCEntries(const Base::Face* face, int& numberOfEntries, std::vector<int>& entries) {
        int number = face->getLocalNrOfBasisFunctions();
        numberOfEntries += number;
        for (int i = 0; i < number; ++i) {
            entries.push_back(startPositionsOfFacesInTheMatrix_[face->getID()] + i);
        }
        std::vector<std::size_t> nodeEntries = face->getPtrElementLeft()->getPhysicalGeometry()->getGlobalFaceNodeIndices(face->localFaceNumberLeft());
        //for (int i = 0; i < face->getReferenceGeometry()->getNumberOfVertices(); ++i) {
        //    for (int j = 0; j < face->getPtrElementLeft()->getLocalNrOfBasisFunctionsVertex(); ++j) {
        //        entries.push_back(startPositionsOfVerticesInTheMatrix_[nodeEntries[i]] + j);
        //    }
        //}
        std::vector<std::size_t> edgeIndex(2);
        for (int i = 0; i < face->getPtrElementLeft()->getNrOfEdges(); ++i) {
            edgeIndex = face->getPtrElementLeft()->getReferenceGeometry()->getCodim2EntityLocalIndices(i);
            edgeIndex[0] = face->getPtrElementLeft()->getPhysicalGeometry()->getNodeIndex(edgeIndex[0]);
            edgeIndex[1] = face->getPtrElementLeft()->getPhysicalGeometry()->getNodeIndex(edgeIndex[1]);
            bool firstFound(false), secondFound(false);
            for (int j = 0; j < nodeEntries.size(); ++j) {
                if (nodeEntries[j] == edgeIndex[0])
                    firstFound = true;
                if (nodeEntries[j] == edgeIndex[1])
                    secondFound = true;
            }
            if (firstFound && secondFound) {
                numberOfEntries += face->getPtrElementLeft()->getEdge(i)->getLocalNrOfBasisFunctions();
                for (int j = 0; j < face->getPtrElementLeft()->getEdge(i)->getLocalNrOfBasisFunctions(); ++j) {
                    entries.push_back(startPositionsOfEdgesInTheMatrix_[face->getPtrElementLeft()->getEdge(i)->getID()] + j);
                }
            }
        }
        nodeEntries = face->getPtrElementLeft()->getPhysicalGeometry()->getLocalFaceNodeIndices(face->localFaceNumberLeft());
        for(std::size_t i:nodeEntries) {
            numberOfEntries += face->getPtrElementLeft()->getNode(i)->getLocalNrOfBasisFunctions();
            for (std::size_t j = 0; j < face->getPtrElementLeft()->getNode(i)->getLocalNrOfBasisFunctions(); ++j)
            entries.push_back(startPositionsOfVerticesInTheMatrix_[face->getPtrElementLeft()->getNode(i)->getID()]+j);
        }
    }
#if defined(HPGEM_USE_PETSC) || defined(HPGEM_USE_COMPLEX_PETSC)

    GlobalPetscMatrix::GlobalPetscMatrix(Base::MeshManipulator* theMesh, int elementMatrixID, int faceMatrixID) :
    GlobalMatrix(theMesh, elementMatrixID, faceMatrixID) {
        PetscBool petscRuns;
        PetscInitialized(&petscRuns);
        if (petscRuns == PETSC_FALSE) {//PETSc thinks real bools are troublesome...
            throw "this happened way too early; please parse the command line arguments as the first thing in your code";
        }

        //temporary
        MatCreateSeqAIJ(PETSC_COMM_SELF, 1, 1, 1, PETSC_NULL, &A_);

        reAssemble();
    }

    GlobalPetscMatrix::~GlobalPetscMatrix() {
        int ierr = MatDestroy(&A_);
        //giving error about Petsc has generated inconsistent data and likely memory corruption in heap
        CHKERRV(ierr);
    }

    GlobalPetscMatrix::operator Mat() {
        //if(meshLevel_!=theMesh_->getActiveLevel(0)){
        //	std::cout<<"Warning: global matrix does not match currently active refinement level!";
        //}
        //MatView(A_,PETSC_VIEWER_STDOUT_WORLD);
        return A_;
    }

    std::vector<PetscInt> GlobalPetscMatrix::makePositionsInMatrix(const Base::Element* element) {
        //we need storage for the amount of basis functions to return
        std::vector<PetscInt> positions(element->getNrOfBasisFunctions());

        auto pos = positions.begin();
        
        std::size_t numElementBasisFuncs = element->getLocalNrOfBasisFunctions();
        //First step: construct ids for the functions of the current element itself
        for (std::size_t i = 0; i < numElementBasisFuncs; ++i) {
            *pos = i + startPositionsOfElementsInTheMatrix_[element->getID()];
            pos++;
            //positions[i] = i + startPositionsOfElementsInTheMatrix_[element->getID()];
        }
        
        //Push forward our iterator
        std::size_t numFaces = element->getPhysicalGeometry()->getNrOfFaces();
        for (std::size_t i = 0; i < numFaces; ++i) {
            std::size_t numFaceBasisFuncs = element->getFace(i)->getLocalNrOfBasisFunctions();
            for (std::size_t j = 0; j < numFaceBasisFuncs; ++j) {
                *pos = j + startPositionsOfFacesInTheMatrix_[element->getFace(i)->getID()];
                pos++;
                //positions[j + usedEntries] = j + startPositionsOfFacesInTheMatrix_[element->getFace(i)->getID()];
            }
        }

        std::size_t numEdges = element->getNrOfEdges();
        for (std::size_t i = 0; i < numEdges; ++i) {
            std::size_t numEdgeBasisFuncs = element->getEdge(i)->getLocalNrOfBasisFunctions();
            for (std::size_t j = 0; j < numEdgeBasisFuncs; ++j) {
                *pos = j + startPositionsOfEdgesInTheMatrix_[element->getEdge(i)->getID()];
                pos++;
                //positions[j + usedEntries] = j + startPositionsOfEdgesInTheMatrix_[element->getEdge(i)->getID()];
            }
        }

        std::size_t numNodes = element->getNrOfNodes();
        for (std::size_t i = 0; i < numNodes; ++i) {
            std::size_t numNodeBasisFuncs = element->getNode(i)->getLocalNrOfBasisFunctions();
            for (std::size_t j = 0; j < numNodeBasisFuncs; ++j) {
                *pos = j + startPositionsOfVerticesInTheMatrix_[element->getNode(i)->getID()];
                pos++;
//                positions[j + usedEntries] = j + startPositionsOfVerticesInTheMatrix_[element->getNode(i)->getID()];
            }
        }
        assert( pos == positions.end() );
        return positions;
    }

    void GlobalPetscMatrix::reset() {
        int ierr = MatZeroEntries(A_);
        CHKERRV(ierr);

        LinearAlgebra::Matrix elementMatrix;

        if (elementMatrixID_ >= 0) {
            for (Base::Element* element : theMesh_->getElementsList() ) {
                std::vector<PetscInt> positions = makePositionsInMatrix(element);
                elementMatrix = element->getElementMatrix(elementMatrixID_);
                ierr = MatSetValues(A_, positions.size(), positions.data(), positions.size(), positions.data(), elementMatrix.data(), ADD_VALUES);
                CHKERRV(ierr);
            } 
        }

        LinearAlgebra::Matrix faceMatrix;
        
        if (faceMatrixID_ >= 0) 
        {
            for(Base::Face* face : theMesh_->getFacesList() ) 
            {
                std::vector<PetscInt> positions = makePositionsInMatrix(face->getPtrElementLeft());
                if(face->isInternal())
                {
                    std::vector<PetscInt> rightPositions = makePositionsInMatrix(face->getPtrElementRight());
                    positions.reserve(positions.size() + rightPositions.size());
                    for(auto& a : rightPositions)
                    {
                        positions.push_back(a);
                    }
                }
                faceMatrix = face->getFaceMatrixMatrix (faceMatrixID_);
                ierr = MatSetValues(A_, positions.size(), positions.data(), positions.size(), positions.data(), faceMatrix.data(), ADD_VALUES);
                CHKERRV(ierr);
            }
        }

        ierr = MatAssemblyBegin(A_, MAT_FINAL_ASSEMBLY);
        ierr = MatAssemblyEnd(A_, MAT_FINAL_ASSEMBLY);
        
        CHKERRV(ierr);
    }

    ///\todo figure out a nice way to keep local data local
    void GlobalPetscMatrix::reAssemble() {
        //if(meshLevel_!=theMesh_->getActiveLevel(0)){
        //meshLevel_=theMesh_->getActiveLevel(0);
        MatDestroy(&A_);
#ifdef HPGEM_USE_MPI
        std::size_t n = Base::MPIContainer::Instance().getNumProcessors();
        //offset by one to put a 0 in front (type int, because MPI wants it to be int)
        std::vector<int> MPISendElementCounts(n+1,0), MPISendFaceCounts(n+1,0), MPISendEdgeCounts(n+1,0), MPISendNodeCounts(n+1,0);
        
        int rank = Base::MPIContainer::Instance().getProcessorID();
        
        MPISendElementCounts[rank+1] = theMesh_->getNumberOfElements();
        MPISendFaceCounts[rank+1] = theMesh_->getNumberOfFaces();
        MPISendEdgeCounts[rank+1] = theMesh_->getNumberOfEdges();
        MPISendNodeCounts[rank+1] = theMesh_->getNumberOfVertices();
        
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
        
        assert(MPISendElementStarts.back()==theMesh_->getNumberOfElements(Base::IteratorType::GLOBAL));
        assert(MPISendFaceStarts.back()>=theMesh_->getNumberOfFaces(Base::IteratorType::GLOBAL));
        assert(MPISendEdgeStarts.back()==theMesh_->getNumberOfEdges(Base::IteratorType::GLOBAL));
        assert(MPISendNodeStarts.back()==theMesh_->getNumberOfVertices(Base::IteratorType::GLOBAL));
        
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
        startPositionsOfElementsInTheMatrix_.resize(theMesh_->getNumberOfElements(Base::IteratorType::GLOBAL));
        startPositionsOfFacesInTheMatrix_.resize(theMesh_->getNumberOfFaces(Base::IteratorType::GLOBAL));
        startPositionsOfEdgesInTheMatrix_.resize(theMesh_->getNumberOfEdges(Base::IteratorType::GLOBAL));
        startPositionsOfVerticesInTheMatrix_.resize(theMesh_->getNumberOfVertices(Base::IteratorType::GLOBAL));
        for(Base::Element* element : theMesh_->getElementsList())
        {
#ifdef HPGEM_USE_MPI
            *currentElementNumber=element->getID();
            *currentElementPosition=totalNrOfDOF;
            ++currentElementNumber;
            ++currentElementPosition;
#else
            startPositionsOfElementsInTheMatrix_[element->getID()] = totalNrOfDOF;
#endif
            totalNrOfDOF += element->getLocalNrOfBasisFunctions();
            for (std::size_t i = 0; i < element->getNrOfFaces(); ++i)
            {
                //faces at the boundary of the subdomain should also be added only once, so add them here is the left element of the face is in the subdomain
                if(element->getFace(i)->getFaceType() == Geometry::FaceType::SUBDOMAIN_BOUNDARY && element->getFace(i)->getPtrElementLeft() == element)
                {
#ifdef HPGEM_USE_MPI
                    *currentFaceNumber=element->getFace(i)->getID();
                    *currentFacePosition=totalNrOfDOF;
                    ++currentFaceNumber;
                    ++currentFacePosition;
#else
                    startPositionsOfFacesInTheMatrix_[element->getFace(i)->getID()] = totalNrOfDOF;
#endif
                    totalNrOfDOF += element->getFace(i)->getLocalNrOfBasisFunctions();
                }
                    
            }
        }
        //DIM == 1 faces and nodes are the same entities, skip one of them
        if(DIM > 1)
        {
            for(Base::Face* face : theMesh_->getFacesList())
            {
                //skip faces at the subdomain boundary because we already treated them
                if(face->getFaceType() != Geometry::FaceType::SUBDOMAIN_BOUNDARY)
                {
#ifdef HPGEM_USE_MPI
                    *currentFaceNumber=face->getID();
                    *currentFacePosition=totalNrOfDOF;
                    ++currentFaceNumber;
                    ++currentFacePosition;
#else
                    startPositionsOfFacesInTheMatrix_[face->getID()] = totalNrOfDOF;
#endif
                    totalNrOfDOF += face->getLocalNrOfBasisFunctions();
                }
            }
        }
        for(Base::Edge* edge : theMesh_->getEdgesList())
        {
#ifdef HPGEM_USE_MPI
            *currentEdgeNumber=edge->getID();
            *currentEdgePosition=totalNrOfDOF;
            ++currentEdgeNumber;
            ++currentEdgePosition;
#else
            startPositionsOfEdgesInTheMatrix_[edge->getID()] = totalNrOfDOF;
#endif
            totalNrOfDOF += edge->getLocalNrOfBasisFunctions();
        }
        for(Base::Node* node : theMesh_->getVerticesList())
        {
#ifdef HPGEM_USE_MPI
            *currentNodeNumber=node->getID();
            *currentNodePosition=totalNrOfDOF;
            ++currentNodeNumber;
            ++currentNodePosition;
#else
            startPositionsOfVerticesInTheMatrix_[node->getID()] = totalNrOfDOF;
#endif
            totalNrOfDOF += node->getLocalNrOfBasisFunctions();
        }

#ifdef HPGEM_USE_MPI
        
        assert(currentElementNumber == MPISendElementNumbers.begin() + MPISendElementStarts[rank+1]);
        assert(currentFaceNumber <= MPISendFaceNumbers.begin() + MPISendFaceStarts[rank+1]);
        assert(currentEdgeNumber == MPISendEdgeNumbers.begin() + MPISendEdgeStarts[rank+1]);
        assert(currentNodeNumber == MPISendNodeNumbers.begin() + MPISendNodeStarts[rank+1]);
        
        assert(currentElementPosition == MPISendElementPositions.begin() + MPISendElementStarts[rank+1]);
        assert(currentFacePosition <= MPISendFacePositions.begin() + MPISendFaceStarts[rank+1]);
        assert(currentEdgePosition == MPISendEdgePositions.begin() + MPISendEdgeStarts[rank+1]);
        assert(currentNodePosition == MPISendNodePositions.begin() + MPISendNodeStarts[rank+1]);
        
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
            assert(*currentElementNumber != std::numeric_limits<std::size_t>::max());
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
                startPositionsOfFacesInTheMatrix_[*currentFaceNumber]=*currentFacePosition+offset;
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
            assert(*currentEdgeNumber != std::numeric_limits<std::size_t>::max());
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
            assert(*currentNodeNumber != std::numeric_limits<std::size_t>::max());
            startPositionsOfVerticesInTheMatrix_[*currentNodeNumber]=*currentNodePosition+offset;
        }
        std::size_t MPIOffset = cumulativeDOF[rank];
#else
        std::size_t MPIOffset = 0;
#endif
       
        //now construct the only bit of data where PETSc expects a local numbering...
        std::vector<PetscInt> numberOfPositionsPerRow(totalNrOfDOF,0);
        std::vector<PetscInt> offDiagonalPositionsPerRow(totalNrOfDOF,0);

        for (Base::Element* element : theMesh_->getElementsList()) {
            for (int j = 0; j < element->getLocalNrOfBasisFunctions(); ++j) {
                numberOfPositionsPerRow[startPositionsOfElementsInTheMatrix_[element->getID()] + j - MPIOffset] += element->getNrOfBasisFunctions();
            }
            for (int i = 0; i < element->getReferenceGeometry()->getNrOfCodim1Entities(); ++i) {
                //conforming contributions
                for (int j = 0; j < element->getFace(i)->getLocalNrOfBasisFunctions(); ++j) {
                    numberOfPositionsPerRow[startPositionsOfFacesInTheMatrix_[element->getFace(i)->getID()] + j - MPIOffset] += element->getNrOfBasisFunctions();
                }
                
                //flux term for DG
                for (int j = 0; j < element->getLocalNrOfBasisFunctions(); ++j) {
                    if(element->getFace(i)->getFaceType()==Geometry::FaceType::SUBDOMAIN_BOUNDARY) {
                        offDiagonalPositionsPerRow[startPositionsOfElementsInTheMatrix_[element->getID()] + j - MPIOffset] += element->getNrOfBasisFunctions();
                    } else if (element->getFace(i)->isInternal()){
                        numberOfPositionsPerRow[startPositionsOfElementsInTheMatrix_[element->getID()] + j - MPIOffset] += element->getNrOfBasisFunctions();
                    }
                }
            }
            for (int i = 0; i < element->getNrOfEdges(); ++i) {
                for (int j = 0; j < element->getEdge(i)->getLocalNrOfBasisFunctions(); ++j) {
                    numberOfPositionsPerRow[startPositionsOfEdgesInTheMatrix_[element->getEdge(i)->getID()] + j - MPIOffset] += element->getNrOfBasisFunctions();
                }
            }
            for (int i = 0; i < element->getNrOfNodes(); ++i) {
                for (int j = 0; j < element->getNode(i)->getLocalNrOfBasisFunctions(); ++j) {
                    numberOfPositionsPerRow[startPositionsOfVerticesInTheMatrix_[element->getNode(i)->getID()] + j - MPIOffset] += element->getNrOfBasisFunctions();
                }
            }
        }
        
        for (int i = 0; i < totalNrOfDOF; ++i) {
            if (numberOfPositionsPerRow[i] > totalNrOfDOF) {
                numberOfPositionsPerRow[i] = totalNrOfDOF; //a row cant have more nonzero entries than the number of columns
            }
        }
        
        int ierr = MatCreateAIJ(PETSC_COMM_WORLD, totalNrOfDOF, totalNrOfDOF, PETSC_DETERMINE, PETSC_DETERMINE, -1, numberOfPositionsPerRow.data(), 0, offDiagonalPositionsPerRow.data(), &A_);
        MatSetOption(A_, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE); //the estimate is known to be wrong for conforming cases
        ierr = MatSetUp(A_);
        CHKERRV(ierr);
        //}
        reset();
    }
#endif
}



