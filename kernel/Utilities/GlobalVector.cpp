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

    void GlobalPetscVector::makePositionsInVector(int amountOfPositions, const Base::Element* element, int positions[]) {
        if (amountOfPositions > 0) {
            int n = element->getLocalNrOfBasisFunctions();
            int usedEntries(0);
            for (int i = 0; i < n; ++i) {
                positions[i] = i + startPositionsOfElementsInTheVector_[element->getID()];
            }
            usedEntries += n;
            int m = element->getPhysicalGeometry()->getNrOfFaces();
            for (int i = 0; i < m; ++i) {
                n = element->getFace(i)->getLocalNrOfBasisFunctions();
                for (int j = 0; j < n; ++j) {
                    positions[j + usedEntries] = j + startPositionsOfFacesInTheVector_[element->getFace(i)->getID()];
                }
                usedEntries += n;
            }
            m = element->getNrOfEdges();
            for (int i = 0; i < m; ++i) {
                n = element->getEdge(i)->getLocalNrOfBasisFunctions();
                for (int j = 0; j < n; ++j) {
                    positions[j + usedEntries] = j + startPositionsOfEdgesInTheVector_[element->getEdge(i)->getID()];
                }
                usedEntries += n;
            }
            m = element->getNrOfNodes();
            for (int i = 0; i < m; ++i) {
                n = element->getLocalNrOfBasisFunctionsVertex();
                for (int j = 0; j < n; ++j) {
                    positions[j + usedEntries] = j + startPositionsOfVerticesInTheVector_[element->getPhysicalGeometry()->getNodeIndex(i)];
                }
                usedEntries += n;
            }
        }


        /*for(unsigned int i=0;i<amountOfPositions;++i){
                int usedEntries(0);
                if(i<element->getLocalNrOfBasisFunctions()){
                        positions[i]=i+startPositionsOfElementsInTheVector_[element->getID()];
                }
                usedEntries+=element->getLocalNrOfBasisFunctions();
                for(int j=0;j<element->getPhysicalGeometry()->getNrOfFaces();++j){
                        if(i-usedEntries<element->getFace(j)->getLocalNrOfBasisFunctions()){
                                positions[i]=i-usedEntries+startPositionsOfFacesInTheVector_[element->getFace(j)->getID()];
                        }
                        usedEntries+=element->getFace(j)->getLocalNrOfBasisFunctions();
                }
                for(int j=0;j<element->getNrOfEdges();++j){
                        if(i-usedEntries<element->getEdge(j)->getLocalNrOfBasisFunctions()){
                                positions[i]=i-usedEntries+startPositionsOfEdgesInTheVector_[element->getEdge(j)->getID()];
                        }
                        usedEntries+=element->getEdge(j)->getLocalNrOfBasisFunctions();
                }
                for(int j=0;j<element->getNrOfNodes();++j){
                        if(i-usedEntries<element->getLocalNrOfBasisFunctionsVertex()){
                                positions[i]=i-usedEntries+startPositionsOfVerticesInTheVector_[element->getPhysicalGeometry()->getNodeIndex(j)];
                        }
                        usedEntries+=element->getLocalNrOfBasisFunctionsVertex();
                }
        }*/
    }

    void GlobalPetscVector::reset() {
        //if(meshLevel_!=theMesh_->getActiveLevel(0)){
        //meshLevel_=theMesh_->getActiveLevel(0);
        int ierr = VecDestroy(&b_);

        int maxNrOfDOF(0), totalNrOfDOF(0), DOFForAVertex(0), DIM(theMesh_->dimension());

        startPositionsOfElementsInTheVector_.resize(theMesh_->getNumberOfElements(Base::IteratorType::GLOBAL));
        startPositionsOfFacesInTheVector_.resize(theMesh_->getNumberOfFaces(Base::IteratorType::GLOBAL));
        startPositionsOfEdgesInTheVector_.resize(theMesh_->getNumberOfEdges(Base::IteratorType::GLOBAL));
        startPositionsOfVerticesInTheVector_.resize(theMesh_->getNumberOfNodes());
        for (Base::MeshManipulator::ElementIterator it = theMesh_->elementColBegin(Base::IteratorType::GLOBAL); it != theMesh_->elementColEnd(Base::IteratorType::GLOBAL); ++it) {
            GlobalVector::startPositionsOfElementsInTheVector_[(*it)->getID()] = totalNrOfDOF;
            if ((*it)->getNrOfBasisFunctions() > maxNrOfDOF) {
                maxNrOfDOF = (*it)->getNrOfBasisFunctions();
                DOFForAVertex = (*it)->getLocalNrOfBasisFunctionsVertex(); //just needs to be set at some point, no need to collect this over and over again
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
        int size = theMesh_->getNodes().size();
        for (int i = 0; i < size; ++i) {//vertices
            startPositionsOfVerticesInTheVector_[i] = totalNrOfDOF;
            totalNrOfDOF += DOFForAVertex;
        }


        //make edges and vertixes periodic
        std::vector<unsigned int> leftIndexes, rightIndexes;
        Geometry::PointReference centre(DIM - 1), leftRef(DIM), rightRef(DIM);
        Geometry::PointPhysical leftPhys(DIM), rightPhys(DIM), displacement(DIM);
        for (Base::Face* face : theMesh_->getFacesList(Base::IteratorType::GLOBAL)) {
            if (face->isInternal()) {
                face->getReferenceGeometry()->getCenter(centre);
                face->mapRefFaceToRefElemL(centre, leftRef);
                face->mapRefFaceToRefElemR(centre, rightRef);
                face->getPtrElementLeft()->referenceToPhysical(leftRef, leftPhys);
                face->getPtrElementRight()->referenceToPhysical(rightRef, rightPhys);
                displacement = leftPhys;
                displacement -= rightPhys;
                if (Utilities::norm2(displacement) > 1e-9) {//this is a periodic boundary
                    face->getPtrElementLeft()->getPhysicalGeometry()->getGlobalFaceNodeIndices(face->localFaceNumberLeft(), leftIndexes);
                    face->getPtrElementRight()->getPhysicalGeometry()->getGlobalFaceNodeIndices(face->localFaceNumberRight(), rightIndexes);
                    const std::vector<Geometry::PointPhysical> locations = face->getPtrElementLeft()->getPhysicalGeometry()->getNodes();
                    for (int i = 0; i < leftIndexes.size(); ++i) {
                        for (int j = 0; j < rightIndexes.size(); ++j) {
                            if (Utilities::norm2(displacement - locations[leftIndexes[i]] + locations[rightIndexes[j]]) < 1e-9) {
                                if (startPositionsOfVerticesInTheVector_[leftIndexes[i]] != startPositionsOfVerticesInTheVector_[rightIndexes[j]]) {
                                    int oldPosition = startPositionsOfVerticesInTheVector_[rightIndexes[j]];
                                    for (int k = 0; k < leftIndexes[i]; ++k) {//before leftIndices[i] connect and shift
                                        if (startPositionsOfVerticesInTheVector_[k] == oldPosition) {
                                            startPositionsOfVerticesInTheVector_[k] = startPositionsOfVerticesInTheVector_[leftIndexes[i]];
                                        }
                                        if (startPositionsOfVerticesInTheVector_[k] > oldPosition) {
                                            startPositionsOfVerticesInTheVector_[k] -= DOFForAVertex;
                                        }
                                    }
                                    for (int k = leftIndexes[i]; k < startPositionsOfVerticesInTheVector_.size(); ++k) {
                                        if (startPositionsOfVerticesInTheVector_[k] > oldPosition) {
                                            startPositionsOfVerticesInTheVector_[k] -= DOFForAVertex; //shift all startposition later then index(j) back a bit
                                        } else if (startPositionsOfVerticesInTheVector_[k] == oldPosition) {//and connect the nodes (also connect previous connections))
                                            startPositionsOfVerticesInTheVector_[k] = startPositionsOfVerticesInTheVector_[leftIndexes[i]];
                                        }
                                    }
                                    totalNrOfDOF -= DOFForAVertex;
                                }
                            }
                        }
                    }
                }
            }
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

        if (elementVectorID_ >= 0) {
            for (Base::MeshManipulator::ElementIterator it = theMesh_->elementColBegin(); it != theMesh_->elementColEnd(); ++it) {
                int n((*it)->getNrOfBasisFunctions()*(*it)->getNrOfUnknows()), positions[n];
                makePositionsInVector(n, *it, positions);
                elementVector.resize(n);
                (*it)->getElementVector(elementVector, elementVectorID_);

                int ierr = VecSetValues(b_, n, positions, &elementVector[0], ADD_VALUES);
                CHKERRV(ierr);
            }
        }

        if (faceVectorID_ >= 0) {
            for (Base::MeshManipulator::FaceIterator it = theMesh_->faceColBegin(); it != theMesh_->faceColEnd(); ++it) {
                const Base::Element * elLeft((*it)->getPtrElementLeft()), *elRight((*it)->getPtrElementRight());
                int nLeft(elLeft->getNrOfBasisFunctions()), n(nLeft);
                if (elRight != NULL)
                    n += elRight->getNrOfBasisFunctions();
                int positions[n];
                makePositionsInVector(nLeft, elLeft, positions);
                makePositionsInVector(n - nLeft, elRight, positions + nLeft);
                elementVector.resize(n);
                (*it)->getFaceVector(elementVector, faceVectorID_);
                int ierr = VecSetValues(b_, n, positions, &elementVector[0], ADD_VALUES);
                CHKERRV(ierr);
            }
        }

        int ierr = VecAssemblyBegin(b_);
        ierr = VecAssemblyEnd(b_);
        CHKERRV(ierr);
    }

    void GlobalPetscVector::constructFromTimeLevelData(int timelevel, int variable) {
        reset();

        LinearAlgebra::NumericalVector elementData;
        for (Base::MeshManipulator::ElementIterator it = theMesh_->elementColBegin(); it != theMesh_->elementColEnd(); ++it) {
            int n((*it)->getNrOfBasisFunctions()), positions[n];
            makePositionsInVector(n, (*it), positions);
            elementData.resize(n);
            for (int i = 0; i < n; ++i) {
                elementData[i] = (*it)->getData(timelevel, variable, i);
            }
            int ierr = VecSetValues(b_, n, positions, &elementData[0], INSERT_VALUES);
        }

        int ierr = VecAssemblyBegin(b_);
        ierr = VecAssemblyEnd(b_);
        CHKERRV(ierr);
    }

    void GlobalPetscVector::writeTimeLevelData(int timeLevel, int variable) {
        double *data;
        int ierr = VecGetArray(b_, &data);
        CHKERRV(ierr);
        for (Base::MeshManipulator::ElementIterator it = theMesh_->elementColBegin(); it != theMesh_->elementColEnd(); ++it) {
            int n = (*it)->getNrOfBasisFunctions();
            LinearAlgebra::NumericalVector localData(&data[startPositionsOfElementsInTheVector_[(*it)->getID()]], n);
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
                for (int j = 0; j < (*it)->getLocalNrOfBasisFunctionsVertex(); ++j) {
                    localData[runningTotal] = data[startPositionsOfVerticesInTheVector_[(*it)->getPhysicalGeometry()->getNodeIndex(i)] + j];
                    ++runningTotal;
                }
            }
            (*it)->setTimeLevelData(timeLevel, variable, localData);
        }
        ierr = VecRestoreArray(b_, &data);
        CHKERRV(ierr);
    }
#endif
}


