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
#include "GlobalVector.hpp"
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
#include <cassert>
#include <numeric>

namespace Utilities {

    GlobalVector::GlobalVector(Base::MeshManipulator* theMesh, int elementVectorID, int faceVectorID) :
    theMesh_(theMesh), startPositionsOfElementsInTheVector_(), meshLevel_(-2), elementVectorID_(elementVectorID), faceVectorID_(faceVectorID) {
    }

#ifdef HPGEM_USE_PETSC

    GlobalPetscVector::GlobalPetscVector(Base::MeshManipulator* theMesh, int elementVectorID, int faceVectorID) :
    GlobalVector(theMesh, elementVectorID, faceVectorID) {
        PetscBool petscRuns;
        PetscInitialized(&petscRuns);
        if (petscRuns == PETSC_FALSE) {//PETSc thinks real bools are troublesome...

            throw "this happened way too early; please parse the command line arguments as the first thing in your code";
        }

        VecCreateSeq(PETSC_COMM_SELF, 1, &b_);

        reset();
    }

    GlobalPetscVector::~GlobalPetscVector() {
        int ierr = VecDestroy(&b_);
        CHKERRV(ierr);
    }

    GlobalPetscVector::operator Vec() {
        //if(meshLevel_!=theMesh_->getActiveLevel(0)){
        //std::cout<<"Warning: global vector does not match currently active refinement level!";
        //}
        //VecView(b_,PETSC_VIEWER_DRAW_WORLD);
        return b_;
    }

    std::vector<PetscInt> GlobalPetscVector::makePositionsInVector(const Base::Element* element) {
        //we need storage for the amount of basis functions to return
        std::vector<PetscInt> positions(element->getNrOfBasisFunctions());

        auto pos = positions.begin();
        
        std::size_t numElementBasisFuncs = element->getLocalNrOfBasisFunctions();
        //First step: construct ids for the functions of the current element itself
        for (std::size_t i = 0; i < numElementBasisFuncs; ++i) {
            *pos = i + startPositionsOfElementsInTheVector_[element->getID()];
            pos++;
            //positions[i] = i + startPositionsOfElementsInTheMatrix_[element->getID()];
        }
        
        //Push forward our iterator
        std::size_t numFaces = element->getPhysicalGeometry()->getNrOfFaces();
        for (std::size_t i = 0; i < numFaces; ++i) {
            std::size_t numFaceBasisFuncs = element->getFace(i)->getLocalNrOfBasisFunctions();
            for (std::size_t j = 0; j < numFaceBasisFuncs; ++j) {
                *pos = j + startPositionsOfFacesInTheVector_[element->getFace(i)->getID()];
                pos++;
                //positions[j + usedEntries] = j + startPositionsOfFacesInTheMatrix_[element->getFace(i)->getID()];
            }
        }

        std::size_t numEdges = element->getNrOfEdges();
        for (std::size_t i = 0; i < numEdges; ++i) {
            std::size_t numEdgeBasisFuncs = element->getEdge(i)->getLocalNrOfBasisFunctions();
            for (std::size_t j = 0; j < numEdgeBasisFuncs; ++j) {
                *pos = j + startPositionsOfEdgesInTheVector_[element->getEdge(i)->getID()];
                pos++;
                //positions[j + usedEntries] = j + startPositionsOfEdgesInTheMatrix_[element->getEdge(i)->getID()];
            }
        }

        std::size_t numNodes = element->getNrOfNodes();
        for (std::size_t i = 0; i < numNodes; ++i) {
            std::size_t numNodeBasisFuncs = element->getNode(i)->getLocalNrOfBasisFunctions();
            for (std::size_t j = 0; j < numNodeBasisFuncs; ++j) {
                *pos = j + startPositionsOfVerticesInTheVector_[element->getNode(i)->getID()];
                pos++;
//                positions[j + usedEntries] = j + startPositionsOfVerticesInTheMatrix_[element->getNode(i)->getID()];
            }
        }
        assert( pos == positions.end() );
        return positions;
    }

    void GlobalPetscVector::reset() {
        //if(meshLevel_!=theMesh_->getActiveLevel(0)){
        //meshLevel_=theMesh_->getActiveLevel(0);
        int ierr = VecDestroy(&b_);
#ifdef HPGEM_USE_MPI
        size_t n = Base::MPIContainer::Instance().getNumProcessors();
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
        std::vector<size_t> MPISendElementNumbers(MPISendElementStarts.back(),std::numeric_limits<size_t>::max());
        std::vector<size_t> MPISendFaceNumbers(MPISendFaceStarts.back(),std::numeric_limits<size_t>::max());
        std::vector<size_t> MPISendEdgeNumbers(MPISendEdgeStarts.back(),std::numeric_limits<size_t>::max());
        std::vector<size_t> MPISendNodeNumbers(MPISendNodeStarts.back(),std::numeric_limits<size_t>::max());
        
        std::vector<size_t> MPISendElementPositions(MPISendElementStarts.back(),std::numeric_limits<size_t>::max());
        std::vector<size_t> MPISendFacePositions(MPISendFaceStarts.back(),std::numeric_limits<size_t>::max());
        std::vector<size_t> MPISendEdgePositions(MPISendEdgeStarts.back(),std::numeric_limits<size_t>::max());
        std::vector<size_t> MPISendNodePositions(MPISendNodeStarts.back(),std::numeric_limits<size_t>::max());
        
        auto currentElementNumber = MPISendElementNumbers.begin() + MPISendElementStarts[rank];
        auto currentFaceNumber = MPISendFaceNumbers.begin() + MPISendFaceStarts[rank];
        auto currentEdgeNumber = MPISendEdgeNumbers.begin() + MPISendEdgeStarts[rank];
        auto currentNodeNumber = MPISendNodeNumbers.begin() + MPISendNodeStarts[rank];
        
        auto currentElementPosition = MPISendElementPositions.begin() + MPISendElementStarts[rank];
        auto currentFacePosition = MPISendFacePositions.begin() + MPISendFaceStarts[rank];
        auto currentEdgePosition = MPISendEdgePositions.begin() + MPISendEdgeStarts[rank];
        auto currentNodePosition = MPISendNodePositions.begin() + MPISendNodeStarts[rank];
#endif
        size_t totalNrOfDOF(0), DIM(theMesh_->dimension());
        startPositionsOfElementsInTheVector_.resize(theMesh_->getNumberOfElements(Base::IteratorType::GLOBAL));
        startPositionsOfFacesInTheVector_.resize(theMesh_->getNumberOfFaces(Base::IteratorType::GLOBAL));
        startPositionsOfEdgesInTheVector_.resize(theMesh_->getNumberOfEdges(Base::IteratorType::GLOBAL));
        startPositionsOfVerticesInTheVector_.resize(theMesh_->getNumberOfVertices(Base::IteratorType::GLOBAL));
        for(Base::Element* element : theMesh_->getElementsList())
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
            for(size_t i = 0; i < element->getNrOfFaces(); ++i)
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
                    startPositionsOfFacesInTheVector_[element->getFace(i)->getID()] = totalNrOfDOF;
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
                    startPositionsOfFacesInTheVector_[face->getID()] = totalNrOfDOF;
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
            startPositionsOfEdgesInTheVector_[edge->getID()] = totalNrOfDOF;
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
            startPositionsOfVerticesInTheVector_[node->getID()] = totalNrOfDOF;
#endif
            totalNrOfDOF += node->getLocalNrOfBasisFunctions();
        }

#ifdef HPGEM_USE_MPI
        //offset by one to insert a zero at the front
        std::vector<size_t> cumulativeDOF(n+1);
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
        
        size_t currentDomain = 0;
        auto startOFNextDomain = MPISendElementNumbers.begin() + MPISendElementCounts[currentDomain+1];
        size_t offset = cumulativeDOF[currentDomain];
        for(;currentElementNumber!=MPISendElementNumbers.end();++currentElementNumber,++currentElementPosition)
        {
            if(currentElementNumber==startOFNextDomain)
            {
                currentDomain++;
                startOFNextDomain += MPISendElementCounts[currentDomain+1];
                offset = cumulativeDOF[currentDomain];
            }
            assert(*currentElementNumber!=std::numeric_limits<size_t>::max());
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
            if(*currentFaceNumber!=std::numeric_limits<size_t>::max())
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
            assert(*currentEdgeNumber!=std::numeric_limits<size_t>::max());
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
            assert(*currentNodeNumber!=std::numeric_limits<size_t>::max());
            startPositionsOfVerticesInTheVector_[*currentNodeNumber]=*currentNodePosition+offset;
        }
#endif
       
        
        ierr = VecCreateMPI(PETSC_COMM_WORLD, totalNrOfDOF, PETSC_DETERMINE, &b_);
        CHKERRV(ierr);
        //}else{
        //	int ierr=VecZeroEntries(b_);
        //	CHKERRV(ierr);
        //}
    }

    void GlobalPetscVector::assemble() {
        reset();

        LinearAlgebra::NumericalVector elementVector;

        if (elementVectorID_ >= 0) 
        {
            for(Base::Element* element : theMesh_->getElementsList())
            {
                std::vector<PetscInt> positions = makePositionsInVector(element);
                element->getElementVector(elementVector, elementVectorID_);
                int ierr = VecSetValues(b_, positions.size(), positions.data(), elementVector.data(), ADD_VALUES);
                CHKERRV(ierr);
                
            }
        }

        LinearAlgebra::NumericalVector faceVector;
        
        if (faceVectorID_ >= 0) 
        {
            for(Base::Face* face : theMesh_->getFacesList()) 
            {
                std::vector<PetscInt> positions = makePositionsInVector(face->getPtrElementLeft());
                if(face->isInternal())
                {
                    std::vector<PetscInt> rightPositions = makePositionsInVector(face->getPtrElementRight());
                    positions.reserve(positions.size() + rightPositions.size());
                    for (auto& a : rightPositions)
                    {
                        positions.push_back(a);
                    }
                }
                face->getFaceVector(faceVector,faceVectorID_);
                int ierr = VecSetValues(b_, positions.size(), positions.data(), faceVector.data(), ADD_VALUES);
                CHKERRV(ierr);
            }
        }

        int ierr = VecAssemblyBegin(b_);
        ierr = VecAssemblyEnd(b_);
        CHKERRV(ierr);
    }

    void GlobalPetscVector::constructFromTimeLevelData(int timelevel, int solutionVar) {
        reset();

        LinearAlgebra::NumericalVector elementData;
        for(Base::Element* element : theMesh_->getElementsList())
        {
            size_t numBasisFuncs = element->getNrOfBasisFunctions();
            std::vector<PetscInt> positions = makePositionsInVector(element);
            elementData.resize(numBasisFuncs);
            for(size_t i=0; i < numBasisFuncs; ++i)
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

    void GlobalPetscVector::writeTimeLevelData(int timeLevel, int variable) {
        double *data;
        
        VecScatter scatter;
        Vec localB;
        //create a local vector...
        VecScatterCreateToAll(b_, &scatter, &localB);
        //we dont need the scatter context, just the vector
        VecScatterDestroy(&scatter);
        
        std::vector<PetscInt> positions;
        for(Base::Element* element : theMesh_->getElementsList())
        {
            std::vector<PetscInt> newPositions = makePositionsInVector(element);
            positions.reserve(positions.size() + newPositions.size());
            for(auto& a : newPositions)
            {
                positions.push_back(a);
            }
        }
        
        IS scatterIS;
        ISCreateGeneral(PETSC_COMM_SELF, positions.size(), positions.data(), PETSC_COPY_VALUES, &scatterIS);
        ISSortRemoveDups(scatterIS);
        VecScatterCreate(b_,scatterIS, localB, scatterIS, &scatter);
        VecScatterBegin(scatter, b_, localB, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(scatter, b_, localB, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterDestroy(&scatter);
        ISDestroy(&scatterIS);
        
        int ierr = VecGetArray(localB, &data);
        CHKERRV(ierr);
        for (Base::MeshManipulator::ElementIterator it = theMesh_->elementColBegin(); it != theMesh_->elementColEnd(); ++it) {
            int numBasisFuns = (*it)->getNrOfBasisFunctions();
            LinearAlgebra::NumericalVector localData(&data[startPositionsOfElementsInTheVector_[(*it)->getID()]], numBasisFuns);
            int runningTotal((*it)->getLocalNrOfBasisFunctions());
            if (theMesh_->dimension() > 1)
                for (int i = 0; i < (*it)->getPhysicalGeometry()->getNrOfFaces(); ++i) {
                    for (int j = 0; j < (*it)->getFace(i)->getLocalNrOfBasisFunctions(); ++j) {
                        localData[runningTotal] = data[startPositionsOfFacesInTheVector_[(*it)->getFace(i)->getID()] + j];
                        ++runningTotal;
                    }
                }
            for (int i = 0; i < (*it)->getNrOfEdges(); ++i) {
                for (int j = 0; j < (*it)->getEdge(i)->getLocalNrOfBasisFunctions(); ++j) {
                    localData[runningTotal] = data[startPositionsOfEdgesInTheVector_[(*it)->getEdge(i)->getID()] + j];
                    ++runningTotal;
                }
            }
            for (int i = 0; i < (*it)->getNrOfNodes(); ++i) {
                for (int j = 0; j < (*it)->getNode(i)->getLocalNrOfBasisFunctions(); ++j) {
                    localData[runningTotal] = data[startPositionsOfVerticesInTheVector_[(*it)->getNode(i)->getID()] + j];
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


