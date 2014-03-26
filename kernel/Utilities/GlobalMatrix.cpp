/*
 * GlobalMatrix.cpp
 *
 *  Created on: Mar 20, 2014
 *      Author: brinkf
 */

#define hpGEM_INCLUDE_PETSC_SUPPORT

#include "GlobalMatrix.hpp"
#include <list>
#include "Base/MeshManipulator.hpp"
#include "Base/Edge.hpp"

namespace Utilities{

	GlobalMatrix::GlobalMatrix(Base::MeshManipulator* theMesh, int elementMatrixID,  int faceMatrixID):
		meshLevel_(-2),theMesh_(theMesh),elementMatrixID_(elementMatrixID),faceMatrixID_(faceMatrixID){}

	void GlobalMatrix::getMatrixBCEntries(Base::Face* face, int& numberOfEntries, std::vector<int>& entries){
		int number=face->getLocalNrOfBasisFunctions();
		numberOfEntries+=number;
		for(int i=0;i<number;++i){
			entries.push_back(startPositionsOfFacesInTheMatrix_[face->getID()]+i);
		}
		std::vector<unsigned int> nodeEntries(face->getReferenceGeometry()->getNumberOfNodes());
		face->getPtrElementLeft()->getPhysicalGeometry()->getGlobalFaceNodeIndices(face->localFaceNumberLeft(),nodeEntries);
		for(int i=0;i<face->getReferenceGeometry()->getNumberOfNodes();++i){
			numberOfEntries+=face->getPtrElementLeft()->getLocalNrOfBasisFunctionsVertex();
			for(int j=0;j<face->getPtrElementLeft()->getLocalNrOfBasisFunctionsVertex();++j){
				entries.push_back(startPositionsOfVerticesInTheMatrix_[nodeEntries[i]]+j);
			}
		}
		std::vector<unsigned int> edgeIndex(2);
		for(int i=0;i<face->getPtrElementLeft()->getNrOfEdges();++i){
			face->getPtrElementLeft()->getReferenceGeometry()->getCodim2EntityLocalIndices(i,edgeIndex);
			edgeIndex[0]=face->getPtrElementLeft()->getPhysicalGeometry()->getNodeIndex(edgeIndex[0]);
			edgeIndex[1]=face->getPtrElementLeft()->getPhysicalGeometry()->getNodeIndex(edgeIndex[1]);
			bool firstFound(false),secondFound(false);
			for(int j=0;j<nodeEntries.size();++j){
				if(nodeEntries[j]==edgeIndex[0])
					firstFound=true;
				if(nodeEntries[j]==edgeIndex[1])
					secondFound=true;
			}
			if(firstFound&&secondFound){
				numberOfEntries+=face->getPtrElementLeft()->getEdge(i)->getLocalNrOfBasisFunctions();
				for(int j=0;j<face->getPtrElementLeft()->getEdge(i)->getLocalNrOfBasisFunctions();++j){
					entries.push_back(startPositionsOfEdgesInTheMatrix_[face->getPtrElementLeft()->getEdge(i)->getID()]+j);
				}
			}
		}
	}
#ifdef hpGEM_INCLUDE_PETSC_SUPPORT
	GlobalPetscMatrix::GlobalPetscMatrix(Base::MeshManipulator* theMesh, int elementMatrixID,  int faceMatrixID):
		GlobalMatrix(theMesh,elementMatrixID,faceMatrixID){
		PetscBool petscRuns;
		PetscInitialized(&petscRuns);
		if(petscRuns==PETSC_FALSE){//PETSc thinks real bools are troublesome...
			int ierr=PetscInitializeNoArguments();//If a user knows to use command line arguments, they know how to properly start PETSc themselves
			CHKERRV(ierr);
		}

		//temporary
		MatCreateSeqAIJ(MPI_COMM_SELF,1,1,1,NULL,&A_);

		reAssemble();
	}

	GlobalPetscMatrix::~GlobalPetscMatrix(){
		int ierr=MatDestroy(&A_);
		CHKERRV(ierr);
	}

	GlobalPetscMatrix::operator Mat(){
		if(meshLevel_!=theMesh_->getActiveLevel(0)){
			std::cout<<"Warning: global matrix does not match currently active refinement level!";
		}
		//MatView(A_,PETSC_VIEWER_DRAW_WORLD);
		return A_;
	}

	void GlobalPetscMatrix::makePositionsInMatrix(int amountOfPositions, const Base::Element* element, int positions[]){
		for(unsigned int i=0;i<amountOfPositions;++i){
			int usedEntries(0);
			if(i<element->getLocalNrOfBasisFunctions()){
				positions[i]=i+startPositionsOfElementsInTheMatrix_[element->getID()];
			}
			usedEntries+=element->getLocalNrOfBasisFunctions();
			for(int j=0;j<element->getPhysicalGeometry()->getNrOfFaces();++j){
				if(i-usedEntries<element->getFace(j)->getLocalNrOfBasisFunctions()){
					positions[i]=i-usedEntries+startPositionsOfFacesInTheMatrix_[element->getFace(j)->getID()];
				}
				usedEntries+=element->getFace(j)->getLocalNrOfBasisFunctions();
			}
			for(int j=0;j<element->getNrOfEdges();++j){
				if(i-usedEntries<element->getEdge(j)->getLocalNrOfBasisFunctions()){
					positions[i]=i-usedEntries+startPositionsOfEdgesInTheMatrix_[element->getEdge(j)->getID()];
				}
				usedEntries+=element->getEdge(j)->getLocalNrOfBasisFunctions();
			}
			for(int j=0;j<element->getNrOfNodes();++j){
				if(i-usedEntries<element->getLocalNrOfBasisFunctionsVertex()){
					positions[i]=i-usedEntries+startPositionsOfVerticesInTheMatrix_[element->getPhysicalGeometry()->getNodeIndex(j)];
				}
				usedEntries+=element->getLocalNrOfBasisFunctionsVertex();
			}
		}
	}

	void GlobalPetscMatrix::reset(){
		int ierr=MatZeroEntries(A_);
		CHKERRV(ierr);

		LinearAlgebra::Matrix elementMatrix;

		if(elementMatrixID_>=0){
			for(Base::MeshManipulator::ElementIterator it=theMesh_->elementColBegin();it!=theMesh_->elementColEnd();++it){
				int n((*it)->getNrOfBasisFunctions()),positions[n];
				makePositionsInMatrix(n,*it,positions);
				elementMatrix.resize(n,n);
				(*it)->getElementMatrix(elementMatrix,elementMatrixID_);

				ierr=MatSetValues(A_,n,positions,n,positions,&elementMatrix[0],ADD_VALUES);
				CHKERRV(ierr);
			}
		}

		if(faceMatrixID_>=0){
			for(Base::MeshManipulator::FaceIterator it=theMesh_->faceColBegin();it!=theMesh_->faceColEnd();++it){
				const Base::Element* elLeft((*it)->getPtrElementLeft()),*elRight((*it)->getPtrElementRight());
				int nLeft(elLeft->getNrOfBasisFunctions()),n(nLeft);
				if((*it)->isInternal())
					n+=elRight->getNrOfBasisFunctions();
				int positions[n];
				makePositionsInMatrix(n-nLeft,elRight,positions+nLeft);
				makePositionsInMatrix(nLeft,elLeft,positions);
				elementMatrix.resize(n,n);
				(*it)->getFaceMatrix(elementMatrix,faceMatrixID_);
				ierr=MatSetValues(A_,n,positions,n,positions,&elementMatrix[0],ADD_VALUES);
				CHKERRV(ierr);
			}
		}

		ierr=MatAssemblyBegin(A_,MAT_FINAL_ASSEMBLY);
		ierr=MatAssemblyEnd(A_,MAT_FINAL_ASSEMBLY);
		CHKERRV(ierr);
	}

	void GlobalPetscMatrix::reAssemble(){
		//if(meshLevel_!=theMesh_->getActiveLevel(0)){
		meshLevel_=theMesh_->getActiveLevel(0);
		MatDestroy(&A_);

		int maxNrOfDOF(0),totalNrOfDOF(0),DIM(theMesh_->dimension()),DOFForAVertex(0);
		startPositionsOfElementsInTheMatrix_.resize(theMesh_->getNumberOfElements());
		startPositionsOfFacesInTheMatrix_.resize(theMesh_->getNumberOfFaces());
		startPositionsOfEdgesInTheMatrix_.resize(theMesh_->getNumberOfEdges());
		startPositionsOfVerticesInTheMatrix_.resize(theMesh_->getNumberOfNodes());
		for(Base::MeshManipulator::ElementIterator it=theMesh_->elementColBegin();it!=theMesh_->elementColEnd();++it){
			if((*it)->getNrOfBasisFunctions()>maxNrOfDOF){
				maxNrOfDOF=(*it)->getNrOfBasisFunctions();
				DOFForAVertex=(*it)->getLocalNrOfBasisFunctionsVertex();//just needs to be set at some point, no need to collect this over and over again
			}
			startPositionsOfElementsInTheMatrix_[(*it)->getID()]=totalNrOfDOF;
			totalNrOfDOF+=(*it)->getLocalNrOfBasisFunctions();
		}
		if(DIM>1){
			for(Base::MeshManipulator::FaceIterator it=theMesh_->faceColBegin();it!=theMesh_->faceColEnd();++it){
				startPositionsOfFacesInTheMatrix_[(*it)->getID()]=totalNrOfDOF;
				totalNrOfDOF+=(*it)->getLocalNrOfBasisFunctions();
			}
		}
		for(std::list< Base::Edge*>::iterator it=theMesh_->edgeColBegin();it!=theMesh_->edgeColEnd();++it){
			startPositionsOfEdgesInTheMatrix_[(*it)->getID()]=totalNrOfDOF;
			totalNrOfDOF+=(*it)->getLocalNrOfBasisFunctions();
		}
		int size=theMesh_->getNodes().size();
		for(int i=0;i<size;++i){//vertices
			startPositionsOfVerticesInTheMatrix_[i]=totalNrOfDOF;
			totalNrOfDOF+=DOFForAVertex;
		}

		int numberOfPositionsPerRow[totalNrOfDOF];
		for(int i=0;i<totalNrOfDOF;++i){
			numberOfPositionsPerRow[i]=0;
		}

		for(Base::Element* element:theMesh_->getElementsList()){
			for(int j=0;j<element->getLocalNrOfBasisFunctions();++j){
				numberOfPositionsPerRow[startPositionsOfElementsInTheMatrix_[element->getID()]+j]+=element->getNrOfBasisFunctions();
			}
			for(int i=0;i<element->getReferenceGeometry()->getNrOfCodim1Entities();++i){
				for(int j=0;j<element->getFace(i)->getLocalNrOfBasisFunctions();++j){
					numberOfPositionsPerRow[startPositionsOfFacesInTheMatrix_[element->getFace(i)->getID()]+j]+=element->getNrOfBasisFunctions();
				}
				for(int j=0;j<element->getLocalNrOfBasisFunctions();++j){//for DG
					numberOfPositionsPerRow[startPositionsOfElementsInTheMatrix_[element->getID()]+j]+=element->getNrOfBasisFunctions();
				}
			}
			for(int i=0;i<element->getNrOfEdges();++i){
				for(int j=0;j<element->getEdge(i)->getLocalNrOfBasisFunctions();++j){
					numberOfPositionsPerRow[startPositionsOfEdgesInTheMatrix_[element->getEdge(i)->getID()]+j]+=element->getNrOfBasisFunctions();
				}
			}
			for(int i=0;i<element->getNrOfNodes();++i){
				for(int j=0;j<DOFForAVertex;++j){
					numberOfPositionsPerRow[startPositionsOfVerticesInTheMatrix_[element->getPhysicalGeometry()->getNodeIndex(i)]+j]+=element->getNrOfBasisFunctions();
				}
			}
		}

		for(int i=0;i<totalNrOfDOF;++i){
			if(totalNrOfDOF<numberOfPositionsPerRow[i])
				numberOfPositionsPerRow[i]=totalNrOfDOF;
		}


		int ierr=MatCreateAIJ(MPI_COMM_WORLD,totalNrOfDOF,totalNrOfDOF,PETSC_DETERMINE,PETSC_DETERMINE,-1,numberOfPositionsPerRow,0,NULL,&A_);
		ierr=MatSetUp(A_);
		CHKERRV(ierr);
		//}
		reset();
	}
#endif
}



