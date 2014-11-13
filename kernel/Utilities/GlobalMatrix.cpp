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

///\TODO optimize placement of matrix rows to minimize communication

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
        std::vector<unsigned int> nodeEntries(face->getReferenceGeometry()->getNumberOfNodes());
        face->getPtrElementLeft()->getPhysicalGeometry()->getGlobalFaceNodeIndices(face->localFaceNumberLeft(), nodeEntries);
        for (int i = 0; i < face->getReferenceGeometry()->getNumberOfNodes(); ++i) {
            numberOfEntries += face->getPtrElementLeft()->getLocalNrOfBasisFunctionsVertex();
            for (int j = 0; j < face->getPtrElementLeft()->getLocalNrOfBasisFunctionsVertex(); ++j) {
                entries.push_back(startPositionsOfVerticesInTheMatrix_[nodeEntries[i]] + j);
            }
        }
        std::vector<unsigned int> edgeIndex(2);
        for (int i = 0; i < face->getPtrElementLeft()->getNrOfEdges(); ++i) {
            face->getPtrElementLeft()->getReferenceGeometry()->getCodim2EntityLocalIndices(i, edgeIndex);
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
    }
#ifdef HPGEM_USE_PETSC

    GlobalPetscMatrix::GlobalPetscMatrix(Base::MeshManipulator* theMesh, int elementMatrixID, int faceMatrixID) :
    GlobalMatrix(theMesh, elementMatrixID, faceMatrixID) {
        PetscBool petscRuns;
        PetscInitialized(&petscRuns);
        if (petscRuns == PETSC_FALSE) {//PETSc thinks real bools are troublesome...
            throw "this happened way too early; please parse the command line arguments as the first thing in your code";
        }

        //temporary
        MatCreateSeqAIJ(PETSC_COMM_SELF, 1, 1, 1, NULL, &A_);

        reAssemble();
    }

    GlobalPetscMatrix::~GlobalPetscMatrix() {
        int ierr = MatDestroy(&A_);
        CHKERRV(ierr);
    }

    GlobalPetscMatrix::operator Mat() {
        //if(meshLevel_!=theMesh_->getActiveLevel(0)){
        //	std::cout<<"Warning: global matrix does not match currently active refinement level!";
        //}
        //MatView(A_,PETSC_VIEWER_DRAW_WORLD);
        return A_;
    }

    void GlobalPetscMatrix::makePositionsInMatrix(int amountOfPositions, const Base::Element* element, int* positions) {
        if (amountOfPositions > 0) {
            int n = element->getLocalNrOfBasisFunctions();
            int usedEntries(0);
            for (int i = 0; i < n; ++i) {
                positions[i] = i + startPositionsOfElementsInTheMatrix_[element->getID()];
            }
            usedEntries += n;
            int m = element->getPhysicalGeometry()->getNrOfFaces();
            for (int i = 0; i < m; ++i) {
                n = element->getFace(i)->getLocalNrOfBasisFunctions();
                for (int j = 0; j < n; ++j) {
                    positions[j + usedEntries] = j + startPositionsOfFacesInTheMatrix_[element->getFace(i)->getID()];
                }
                usedEntries += n;
            }
            m = element->getNrOfEdges();
            for (int i = 0; i < m; ++i) {
                n = element->getEdge(i)->getLocalNrOfBasisFunctions();
                for (int j = 0; j < n; ++j) {
                    positions[j + usedEntries] = j + startPositionsOfEdgesInTheMatrix_[element->getEdge(i)->getID()];
                }
                usedEntries += n;
            }
            m = element->getNrOfNodes();
            for (int i = 0; i < m; ++i) {
                n = element->getLocalNrOfBasisFunctionsVertex();
                for (int j = 0; j < n; ++j) {
                    positions[j + usedEntries] = j + startPositionsOfVerticesInTheMatrix_[element->getPhysicalGeometry()->getNodeIndex(i)];
                }
                usedEntries += n;
            }
            assert(usedEntries == amountOfPositions);
        }
    }

    void GlobalPetscMatrix::reset() {
        int ierr = MatZeroEntries(A_);
        CHKERRV(ierr);

        LinearAlgebra::Matrix elementMatrix;

        if (elementMatrixID_ >= 0) {
            for (Base::MeshManipulator::ElementIterator it = theMesh_->elementColBegin(); it != theMesh_->elementColEnd(); ++it) {
                int n((*it)->getNrOfBasisFunctions()), positions[n];
                makePositionsInMatrix(n, *it, positions);
                elementMatrix.resize(n, n);
                (*it)->getElementMatrix(elementMatrix, elementMatrixID_);

                ierr = MatSetValues(A_, n, positions, n, positions, &elementMatrix[0], ADD_VALUES);
                CHKERRV(ierr);
            }
        }

        if (faceMatrixID_ >= 0) {
            for (Base::MeshManipulator::FaceIterator it = theMesh_->faceColBegin(); it != theMesh_->faceColEnd(); ++it) {
                const Base::Element * elLeft((*it)->getPtrElementLeft()), *elRight((*it)->getPtrElementRight());
                int nLeft(elLeft->getNrOfBasisFunctions()), n(nLeft);
                if ((*it)->isInternal()) {
                    n += elRight->getNrOfBasisFunctions();
                }
                int positions[n];
                makePositionsInMatrix(nLeft, elLeft, positions);
                makePositionsInMatrix(n - nLeft, elRight, positions + nLeft);
                elementMatrix.resize(n, n);
                (*it)->getFaceMatrix(elementMatrix, faceMatrixID_);
                ierr = MatSetValues(A_, n, positions, n, positions, &elementMatrix[0], ADD_VALUES);
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

        int maxNrOfDOF(0), totalNrOfDOF(0), DIM(theMesh_->dimension()), DOFForAVertex(0);
        startPositionsOfElementsInTheMatrix_.resize(theMesh_->getNumberOfElements(Base::IteratorType::GLOBAL));
        startPositionsOfFacesInTheMatrix_.resize(theMesh_->getNumberOfFaces(Base::IteratorType::GLOBAL));
        startPositionsOfEdgesInTheMatrix_.resize(theMesh_->getNumberOfEdges(Base::IteratorType::GLOBAL));
        startPositionsOfVerticesInTheMatrix_.resize(theMesh_->getNumberOfNodes());
        for (Base::MeshManipulator::ElementIterator it = theMesh_->elementColBegin(Base::IteratorType::GLOBAL); it != theMesh_->elementColEnd(Base::IteratorType::GLOBAL); ++it) {
            if ((*it)->getNrOfBasisFunctions() > maxNrOfDOF) {
                maxNrOfDOF = (*it)->getNrOfBasisFunctions();
                DOFForAVertex = (*it)->getLocalNrOfBasisFunctionsVertex(); //just needs to be set at some point, no need to collect this over and over again
            }
            startPositionsOfElementsInTheMatrix_[(*it)->getID()] = totalNrOfDOF;
            totalNrOfDOF += (*it)->getLocalNrOfBasisFunctions();
        }
        if (DIM > 1) {
            for (Base::MeshManipulator::FaceIterator it = theMesh_->faceColBegin(Base::IteratorType::GLOBAL); it != theMesh_->faceColEnd(Base::IteratorType::GLOBAL); ++it) {
                startPositionsOfFacesInTheMatrix_[(*it)->getID()] = totalNrOfDOF;
                totalNrOfDOF += (*it)->getLocalNrOfBasisFunctions();
            }
        }
        for (std::vector< Base::Edge*>::iterator it = theMesh_->edgeColBegin(Base::IteratorType::GLOBAL); it != theMesh_->edgeColEnd(Base::IteratorType::GLOBAL); ++it) {
            startPositionsOfEdgesInTheMatrix_[(*it)->getID()] = totalNrOfDOF;
            totalNrOfDOF += (*it)->getLocalNrOfBasisFunctions();
        }
        int size = theMesh_->getNodes().size();
        for (int i = 0; i < size; ++i) {//vertices
            startPositionsOfVerticesInTheMatrix_[i] = totalNrOfDOF;
            totalNrOfDOF += DOFForAVertex;
        }

        //make vertices periodic ///\bug periodic boundary conditions not set for edges
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
                            if (Utilities::norm2(displacement - locations[leftIndexes[i]] + locations[rightIndexes[j]]) < 1e-9) {//if the nodes have the same displacement as the centre
                                if (startPositionsOfVerticesInTheMatrix_[leftIndexes[i]] != startPositionsOfVerticesInTheMatrix_[rightIndexes[j]]) {//and they do not map into the same location in the matrix yet
                                    int oldPosition = startPositionsOfVerticesInTheMatrix_[rightIndexes[j]];
                                    for (int k = 0; k < leftIndexes[i]; ++k) {//before leftIndices[i] connect and shift
                                        if (startPositionsOfVerticesInTheMatrix_[k] == oldPosition) {
                                            startPositionsOfVerticesInTheMatrix_[k] = startPositionsOfVerticesInTheMatrix_[leftIndexes[i]];
                                        }
                                        if (startPositionsOfVerticesInTheMatrix_[k] > oldPosition) {
                                            startPositionsOfVerticesInTheMatrix_[k] -= DOFForAVertex;
                                        }
                                    }
                                    for (int k = leftIndexes[i]; k < startPositionsOfVerticesInTheMatrix_.size(); ++k) {
                                        if (startPositionsOfVerticesInTheMatrix_[k] > oldPosition) {
                                            startPositionsOfVerticesInTheMatrix_[k] -= DOFForAVertex; //shift all startposition later then index(j) back a bit
                                        } else if (startPositionsOfVerticesInTheMatrix_[k] == oldPosition) {//and connect the nodes (also connect previous connections))
                                            startPositionsOfVerticesInTheMatrix_[k] = startPositionsOfVerticesInTheMatrix_[leftIndexes[i]];
                                        }
                                    }
                                    totalNrOfDOF -= DOFForAVertex;
                                    //std::cout<<leftIndexes[i]<<" "<<rightIndexes[j]<<std::endl;
                                }
                            }
                        }
                    }
                }
            }
        }


        int *numberOfPositionsPerRow=new int[totalNrOfDOF];
        for (int i = 0; i < totalNrOfDOF; ++i) {
            numberOfPositionsPerRow[i] = 0;
        }

        for (Base::Element* element : theMesh_->getElementsList()) {
            for (int j = 0; j < element->getLocalNrOfBasisFunctions(); ++j) {
                numberOfPositionsPerRow[startPositionsOfElementsInTheMatrix_[element->getID()] + j] += element->getNrOfBasisFunctions();
            }
            for (int i = 0; i < element->getReferenceGeometry()->getNrOfCodim1Entities(); ++i) {
                for (int j = 0; j < element->getFace(i)->getLocalNrOfBasisFunctions(); ++j) {
                    numberOfPositionsPerRow[startPositionsOfFacesInTheMatrix_[element->getFace(i)->getID()] + j] += element->getNrOfBasisFunctions();
                }
                for (int j = 0; j < element->getLocalNrOfBasisFunctions(); ++j) {
                    numberOfPositionsPerRow[startPositionsOfElementsInTheMatrix_[element->getID()] + j] += element->getNrOfBasisFunctions();
                }
            }
            for (int i = 0; i < element->getNrOfEdges(); ++i) {
                for (int j = 0; j < element->getEdge(i)->getLocalNrOfBasisFunctions(); ++j) {
                    numberOfPositionsPerRow[startPositionsOfEdgesInTheMatrix_[element->getEdge(i)->getID()] + j] += element->getNrOfBasisFunctions();
                }
            }
            for (int i = 0; i < element->getNrOfNodes(); ++i) {
                for (int j = 0; j < DOFForAVertex; ++j) {
                    numberOfPositionsPerRow[startPositionsOfVerticesInTheMatrix_[element->getPhysicalGeometry()->getNodeIndex(i)] + j] += element->getNrOfBasisFunctions();
                }
            }
        }

        
        
        for (int i = 0; i < totalNrOfDOF; ++i) {
            if (numberOfPositionsPerRow[i] > totalNrOfDOF) {
                numberOfPositionsPerRow[i] = totalNrOfDOF; //a row cant have more nonzero entries than the number of columns
            }
        }
        
        int ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, totalNrOfDOF, totalNrOfDOF, -1, numberOfPositionsPerRow, 0, NULL, &A_);
        MatSetOption(A_, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE); //the estimate is accurate for most matrices, but if you want to do something conforming with face matrices you need more
        ierr = MatSetUp(A_);
        CHKERRV(ierr);
        delete[] numberOfPositionsPerRow;
        //}
        reset();
    }
#endif
}



