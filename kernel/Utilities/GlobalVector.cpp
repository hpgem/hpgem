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

#define hpGEM_INCLUDE_PETSC_SUPPORT

#include "GlobalVector.hpp"
#include <list>
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

namespace Utilities{

	GlobalVector::GlobalVector(Base::MeshManipulator* theMesh,int elementVectorID, int faceVectorID):
		theMesh_(theMesh),startPositionsOfElementsInTheVector_(),meshLevel_(-2),elementVectorID_(elementVectorID),faceVectorID_(faceVectorID){}

#ifdef hpGEM_INCLUDE_PETSC_SUPPORT
	GlobalPetscVector::GlobalPetscVector(Base::MeshManipulator* theMesh,int elementVectorID,int faceVectorID):
		GlobalVector(theMesh,elementVectorID,faceVectorID){
		PetscBool petscRuns;
		PetscInitialized(&petscRuns);
		if(petscRuns==PETSC_FALSE){//PETSc thinks real bools are troublesome...
			int ierr=PetscInitializeNoArguments();//If a user knows to use command line arguments, they know how to properly start PETSc themselves
			CHKERRV(ierr);
		}

		VecCreateSeq(MPI_COMM_SELF,1,&b_);

		reset();
	}

	GlobalPetscVector::~GlobalPetscVector(){
		int ierr=VecDestroy(&b_);
		CHKERRV(ierr);
	}

	GlobalPetscVector::operator Vec(){
		//if(meshLevel_!=theMesh_->getActiveLevel(0)){
			//std::cout<<"Warning: global vector does not match currently active refinement level!";
		//}
		//VecView(b_,PETSC_VIEWER_DRAW_WORLD);
		return b_;
	}

	void GlobalPetscVector::makePositionsInVector(int amountOfPositions, const Base::Element* element, int positions[]){
		for(unsigned int i=0;i<amountOfPositions;++i){
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
		}
	}

	void GlobalPetscVector::reset(){
		//if(meshLevel_!=theMesh_->getActiveLevel(0)){
			//meshLevel_=theMesh_->getActiveLevel(0);
			int ierr=VecDestroy(&b_);

			int maxNrOfDOF(0),totalNrOfDOF(0),DOFForAVertex(0),DIM(theMesh_->dimension());

			startPositionsOfElementsInTheVector_.resize(theMesh_->getNumberOfElements());
			startPositionsOfFacesInTheVector_.resize(theMesh_->getNumberOfFaces());
			startPositionsOfEdgesInTheVector_.resize(theMesh_->getNumberOfEdges());
			startPositionsOfVerticesInTheVector_.resize(theMesh_->getNumberOfNodes());
			for(Base::MeshManipulator::ElementIterator it=theMesh_->elementColBegin();it!=theMesh_->elementColEnd();++it){
				GlobalVector::startPositionsOfElementsInTheVector_[(*it)->getID()]=totalNrOfDOF;
				if((*it)->getNrOfBasisFunctions()>maxNrOfDOF){
					maxNrOfDOF=(*it)->getNrOfBasisFunctions();
					DOFForAVertex=(*it)->getLocalNrOfBasisFunctionsVertex();//just needs to be set at some point, no need to collect this over and over again
				}
				totalNrOfDOF+=(*it)->getLocalNrOfBasisFunctions();
			}
			if(DIM>1){
				for(Base::MeshManipulator::FaceIterator it=theMesh_->faceColBegin();it!=theMesh_->faceColEnd();++it){
					startPositionsOfFacesInTheVector_[(*it)->getID()]=totalNrOfDOF;
					totalNrOfDOF+=(*it)->getLocalNrOfBasisFunctions();
				}
			}
			for(std::list< Base::Edge*>::iterator it=theMesh_->edgeColBegin();it!=theMesh_->edgeColEnd();++it){
				startPositionsOfEdgesInTheVector_[(*it)->getID()]=totalNrOfDOF;
				totalNrOfDOF+=(*it)->getLocalNrOfBasisFunctions();
			}
			int size=theMesh_->getNodes().size();
			for(int i=0;i<size;++i){//vertices
				startPositionsOfVerticesInTheVector_[i]=totalNrOfDOF;
				totalNrOfDOF+=DOFForAVertex;
			}

			//at the moment assumes fully local vector in the d_nz/o_nz part of the arguments
			ierr=VecCreateMPI(MPI_COMM_WORLD,totalNrOfDOF,PETSC_DETERMINE,&b_);
			CHKERRV(ierr);

		//}else{
		//	int ierr=VecZeroEntries(b_);
		//	CHKERRV(ierr);
		//}
	}

	void GlobalPetscVector::assemble(){
		reset();

		LinearAlgebra::NumericalVector elementVector;

		if(elementVectorID_>=0){
			for(Base::MeshManipulator::ElementIterator it=theMesh_->elementColBegin();it!=theMesh_->elementColEnd();++it){
				int n((*it)->getNrOfBasisFunctions()*(*it)->getNrOfUnknows()),positions[n];
				makePositionsInVector(n,*it,positions);
				elementVector.resize(n);
				(*it)->getElementVector(elementVector,elementVectorID_);

				int ierr=VecSetValues(b_,n,positions,&elementVector[0],ADD_VALUES);
				CHKERRV(ierr);
			}
		}

		if(faceVectorID_>=0){
			for(Base::MeshManipulator::FaceIterator it=theMesh_->faceColBegin();it!=theMesh_->faceColEnd();++it){
				const Base::Element* elLeft((*it)->getPtrElementLeft()),*elRight((*it)->getPtrElementRight());
				int nLeft(elLeft->getNrOfBasisFunctions()),n(nLeft);
				if(elRight!=NULL)
					n+=elRight->getNrOfBasisFunctions();
				int positions[n];
				makePositionsInVector(nLeft,elLeft,positions);
				makePositionsInVector(n-nLeft,elRight,positions+nLeft);
				elementVector.resize(n);
				(*it)->getFaceVector(elementVector,faceVectorID_);
				int ierr=VecSetValues(b_,n,positions,&elementVector[0],ADD_VALUES);
				CHKERRV(ierr);
			}
		}

		int ierr=VecAssemblyBegin(b_);
		ierr=VecAssemblyEnd(b_);
		CHKERRV(ierr);
	}

	void GlobalPetscVector::constructFromTimeLevelData(int timelevel){
		reset();

		LinearAlgebra::Matrix elementData;
		for(Base::MeshManipulator::ElementIterator it=theMesh_->elementColBegin();it!=theMesh_->elementColEnd();++it){
			int n((*it)->getNrOfBasisFunctions()),positions[n];
			makePositionsInVector(n,(*it),positions);
			elementData.resize(n,1);
			elementData=(*it)->getTimeLevelData(timelevel);
			//coefficients belonging to conforming basisfunctions appear in multiple elements so don't add
			int ierr=VecSetValues(b_,n,positions,&elementData[0],INSERT_VALUES);
		}

		int ierr=VecAssemblyBegin(b_);
		ierr=VecAssemblyEnd(b_);
		CHKERRV(ierr);
	}

	void GlobalPetscVector::writeTimeLevelData(int timeLevel){
		double *data;
		int ierr=VecGetArray(b_,&data);
		CHKERRV(ierr);
		for(Base::MeshManipulator::ElementIterator it=theMesh_->elementColBegin();it!=theMesh_->elementColEnd();++it){
			int n=(*it)->getNrOfBasisFunctions();
			LinearAlgebra::NumericalVector localData(&data[startPositionsOfElementsInTheVector_[(*it)->getID()]],n);
			int runningTotal((*it)->getLocalNrOfBasisFunctions());
			if(theMesh_->dimension()>1)
			for(int i=0;i<(*it)->getPhysicalGeometry()->getNrOfFaces();++i){
				for(int j=0;j<(*it)->getFace(i)->getLocalNrOfBasisFunctions();++j){
					localData[runningTotal]=data[startPositionsOfFacesInTheVector_[(*it)->getFace(i)->getID()]+j];
					++runningTotal;
				}
			}
			for(int i=0;i<(*it)->getNrOfEdges();++i){
				for(int j=0;j<(*it)->getEdge(i)->getLocalNrOfBasisFunctions();++j){
					localData[runningTotal]=data[startPositionsOfEdgesInTheVector_[(*it)->getEdge(i)->getID()]+j];
					++runningTotal;
				}
			}
			for(int i=0;i<(*it)->getNrOfNodes();++i){
				for(int j=0;j<(*it)->getLocalNrOfBasisFunctionsVertex();++j){
					localData[runningTotal]=data[startPositionsOfVerticesInTheVector_[(*it)->getPhysicalGeometry()->getNodeIndex(i)]+j];
					++runningTotal;
				}
			}
			(*it)->setTimeLevelData(timeLevel,0,localData);
		}
		ierr=VecRestoreArray(b_,&data);
		CHKERRV(ierr);
	}
#endif
}


