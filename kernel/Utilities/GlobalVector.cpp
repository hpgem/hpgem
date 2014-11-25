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

//#define HPGEM_USE_PETSC//temporarily activating this definition makes development easier on some IDEs

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

        int maxNrOfDOF(0), totalNrOfDOF(0), DOFForAVertex(0), DIM(theMesh_->dimension());

        startPositionsOfElementsInTheVector_.resize(theMesh_->getNumberOfElements(Base::IteratorType::GLOBAL));
        startPositionsOfFacesInTheVector_.resize(theMesh_->getNumberOfFaces(Base::IteratorType::GLOBAL));
        startPositionsOfEdgesInTheVector_.resize(theMesh_->getNumberOfEdges(Base::IteratorType::GLOBAL));
        startPositionsOfVerticesInTheVector_.resize(theMesh_->getNumberOfVertices(Base::IteratorType::GLOBAL));
        for (Base::MeshManipulator::ElementIterator it = theMesh_->elementColBegin(Base::IteratorType::GLOBAL); it != theMesh_->elementColEnd(Base::IteratorType::GLOBAL); ++it) {
            GlobalVector::startPositionsOfElementsInTheVector_[(*it)->getID()] = totalNrOfDOF;
            if ((*it)->getNrOfBasisFunctions() > maxNrOfDOF) {
                maxNrOfDOF = (*it)->getNrOfBasisFunctions();
            }
            totalNrOfDOF += (*it)->getLocalNrOfBasisFunctions();
        }
        if (DIM > 1) {
            for (Base::MeshManipulator::FaceIterator it = theMesh_->faceColBegin(Base::IteratorType::GLOBAL); it != theMesh_->faceColEnd(Base::IteratorType::GLOBAL); ++it) {
                startPositionsOfFacesInTheVector_[(*it)->getID()] = totalNrOfDOF;
                totalNrOfDOF += (*it)->getLocalNrOfBasisFunctions();
            }
        }
        for (std::vector< Base::Edge*>::iterator it = theMesh_->edgeColBegin(Base::IteratorType::GLOBAL); it != theMesh_->edgeColEnd(Base::IteratorType::GLOBAL); ++it) {
            startPositionsOfEdgesInTheVector_[(*it)->getID()] = totalNrOfDOF;
            totalNrOfDOF += (*it)->getLocalNrOfBasisFunctions();
        }
        for(Base::Node* node : theMesh_->getVerticesList(Base::IteratorType::GLOBAL))
        {
            startPositionsOfVerticesInTheVector_[node->getID()] = totalNrOfDOF;
            totalNrOfDOF += node->getLocalNrOfBasisFunctions();
        }

        ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, totalNrOfDOF, &b_);
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
        
        ///\TODO highly suboptimal...
        VecScatter scatter;
        Vec localB;
        VecScatterCreateToAll(b_, &scatter, &localB);
        VecScatterBegin(scatter, b_, localB, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(scatter, b_, localB, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterDestroy(&scatter);
        
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


