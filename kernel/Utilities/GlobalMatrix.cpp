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
        //for (int i = 0; i < face->getReferenceGeometry()->getNumberOfVertices(); ++i) {
        //    for (int j = 0; j < face->getPtrElementLeft()->getLocalNrOfBasisFunctionsVertex(); ++j) {
        //        entries.push_back(startPositionsOfVerticesInTheMatrix_[nodeEntries[i]] + j);
        //    }
        //}
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
        face->getPtrElementLeft()->getPhysicalGeometry()->getLocalFaceNodeIndices(face->localFaceNumberLeft(),nodeEntries);
        for(unsigned int i:nodeEntries) {
            numberOfEntries += face->getPtrElementLeft()->getNode(i)->getLocalNrOfBasisFunctions();
            for(size_t j=0;j<face->getPtrElementLeft()->getNode(i)->getLocalNrOfBasisFunctions();++j)
            entries.push_back(startPositionsOfVerticesInTheMatrix_[face->getPtrElementLeft()->getNode(i)->getID()]+j);
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
                element->getElementMatrix(elementMatrix, elementMatrixID_);
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
                face->getFaceMatrix(faceMatrix, faceMatrixID_);
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

        int maxNrOfDOF(0), totalNrOfDOF(0), DIM(theMesh_->dimension()), DOFForAVertex(0);
        startPositionsOfElementsInTheMatrix_.resize(theMesh_->getNumberOfElements(Base::IteratorType::GLOBAL));
        startPositionsOfFacesInTheMatrix_.resize(theMesh_->getNumberOfFaces(Base::IteratorType::GLOBAL));
        startPositionsOfEdgesInTheMatrix_.resize(theMesh_->getNumberOfEdges(Base::IteratorType::GLOBAL));
        startPositionsOfVerticesInTheMatrix_.resize(theMesh_->getNumberOfVertices(Base::IteratorType::GLOBAL));
        for (Base::MeshManipulator::ElementIterator it = theMesh_->elementColBegin(Base::IteratorType::GLOBAL); it != theMesh_->elementColEnd(Base::IteratorType::GLOBAL); ++it) {
            if ((*it)->getNrOfBasisFunctions() > maxNrOfDOF) {
                maxNrOfDOF = (*it)->getNrOfBasisFunctions();
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
        for(Base::Node* node : theMesh_->getVerticesList(Base::IteratorType::GLOBAL))
        {
            startPositionsOfVerticesInTheMatrix_[node->getID()] = totalNrOfDOF;
            totalNrOfDOF += node->getLocalNrOfBasisFunctions();
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
                for (int j = 0; j < element->getNode(i)->getLocalNrOfBasisFunctions(); ++j) {
                    numberOfPositionsPerRow[startPositionsOfVerticesInTheMatrix_[element->getNode(i)->getID()] + j] += element->getNrOfBasisFunctions();
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



