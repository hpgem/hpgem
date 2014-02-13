/*
 * GlobalAssembly.cpp
 *
 *  Created on: Jan 17, 2014
 *      Author: brinkf
 */

#include "GlobalAssembly.hpp"

namespace Utilities{

	GlobalMatrix::GlobalMatrix(Base::MeshManipulator* theMesh,unsigned int elementMatrixID, unsigned int faceMatrixID):
		meshLevel_(-2),theMesh_(theMesh),elementMatrixID_(elementMatrixID),faceMatrixID_(faceMatrixID){}

	GlobalPetscMatrix::GlobalPetscMatrix(Base::MeshManipulator* theMesh,unsigned int elementMatrixID, unsigned int faceMatrixID):
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
		//MatView(A_,PETSC_VIEWER_STDOUT_WORLD);
		return A_;
	}

	void GlobalPetscMatrix::reset(){
		int ierr=MatZeroEntries(A_);
		CHKERRV(ierr);

		LinearAlgebra::Matrix elementMatrix;

		for(Base::MeshManipulator::ElementIterator it=theMesh_->elementColBegin();it!=theMesh_->elementColEnd();++it){
			int n((*it)->getNrOfBasisFunctions()*(*it)->getNrOfUnknows()),positions[n];
			for(int i=0;i<n;++i){
				positions[i]=i+startPositionsOfElementsInTheMatrix_[(*it)->getID()];
			}
			elementMatrix.resize(n,n);
			(*it)->getElementMatrix(elementMatrix,elementMatrixID_);

			ierr=MatSetValues(A_,n,positions,n,positions,&elementMatrix[0],ADD_VALUES);
			CHKERRV(ierr);
		}

		for(Base::MeshManipulator::FaceIterator it=theMesh_->faceColBegin();it!=theMesh_->faceColEnd();++it){
			const Base::Element* elLeft((*it)->getPtrElementLeft()),*elRight((*it)->getPtrElementRight());
			int nLeft(elLeft->getNrOfBasisFunctions()),n(nLeft);
			if(elRight!=NULL)
				n+=elRight->getNrOfBasisFunctions();
			int positions[n];
			for(int i=0;i<n;++i){
				if(i<nLeft)
					positions[i]=i+startPositionsOfElementsInTheMatrix_[elLeft->getID()];
				else
					positions[i]=i-nLeft+startPositionsOfElementsInTheMatrix_[elRight->getID()];
			}
			elementMatrix.resize(n,n);
			(*it)->getFaceMatrix(elementMatrix,faceMatrixID_);
			ierr=MatSetValues(A_,n,positions,n,positions,&elementMatrix[0],ADD_VALUES);
			CHKERRV(ierr);
		}

		ierr=MatAssemblyBegin(A_,MAT_FINAL_ASSEMBLY);
		ierr=MatAssemblyEnd(A_,MAT_FINAL_ASSEMBLY);
		CHKERRV(ierr);
	}

	void GlobalPetscMatrix::reAssemble(){
		if(meshLevel_!=theMesh_->getActiveLevel(0)){
			meshLevel_=theMesh_->getActiveLevel(0);
			MatDestroy(&A_);

			int maxNrOfDOF(0),totalNrOfDOF(0),DIM(0);
			startPositionsOfElementsInTheMatrix_.resize(theMesh_->getNumberOfElements());
			for(Base::MeshManipulator::ElementIterator it=theMesh_->elementColBegin();it!=theMesh_->elementColEnd();++it){
				if((*it)->getNrOfBasisFunctions()>maxNrOfDOF){
					maxNrOfDOF=(*it)->getNrOfBasisFunctions();
					DIM=(*it)->getReferenceGeometry()->getGaussQuadratureRule(0)->dimension();
				}
				startPositionsOfElementsInTheMatrix_[(*it)->getID()]=totalNrOfDOF;
				totalNrOfDOF+=(*it)->getNrOfBasisFunctions();
			}

			int ierr=MatCreateAIJ(MPI_COMM_WORLD,totalNrOfDOF,totalNrOfDOF,PETSC_DETERMINE,PETSC_DETERMINE,(1+2*DIM)*maxNrOfDOF,NULL,0,NULL,&A_);
			ierr=MatSetUp(A_);
			CHKERRV(ierr);
		}
		reset();
	}


	GlobalVector::GlobalVector(Base::MeshManipulator* theMesh,int elementVectorID, int faceVectorID):
		theMesh_(theMesh),startPositionsOfElementsInTheVector_(),meshLevel_(-2),elementVectorID_(elementVectorID),faceVectorID_(faceVectorID){}

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
		if(meshLevel_!=theMesh_->getActiveLevel(0)){
			std::cout<<"Warning: global vector does not match currently active refinement level!";
		}
		//VecView(b_,PETSC_VIEWER_STDOUT_WORLD);
		return b_;
	}

	void GlobalPetscVector::reset(){
		if(meshLevel_!=theMesh_->getActiveLevel(0)){
			meshLevel_=theMesh_->getActiveLevel(0);
			int ierr=VecDestroy(&b_);

			int maxNrOfDOF(0),totalNrOfDOF(0);

			startPositionsOfElementsInTheVector_.resize(theMesh_->getNumberOfElements());
			for(Base::MeshManipulator::ElementIterator it=theMesh_->elementColBegin();it!=theMesh_->elementColEnd();++it){
				GlobalVector::startPositionsOfElementsInTheVector_[(*it)->getID()]=totalNrOfDOF;
				if((*it)->getNrOfBasisFunctions()>maxNrOfDOF){
					maxNrOfDOF=(*it)->getNrOfBasisFunctions();
				}
				totalNrOfDOF+=(*it)->getNrOfBasisFunctions();
			}

			//at the moment assumes fully local vector in the d_nz/o_nz part of the arguments
			ierr=VecCreateMPI(MPI_COMM_WORLD,totalNrOfDOF,PETSC_DETERMINE,&b_);
			CHKERRV(ierr);

		}else{
			int ierr=VecZeroEntries(b_);
			CHKERRV(ierr);
		}
	}

	void GlobalPetscVector::assemble(){
		reset();

		LinearAlgebra::NumericalVector elementVector;

		for(Base::MeshManipulator::ElementIterator it=theMesh_->elementColBegin();it!=theMesh_->elementColEnd();++it){
			int n((*it)->getNrOfBasisFunctions()*(*it)->getNrOfUnknows()),positions[n];
			for(int i=0;i<n;++i){
				positions[i]=i+startPositionsOfElementsInTheVector_[(*it)->getID()];
			}
			elementVector.resize(n);
			(*it)->getElementVector(elementVector,elementVectorID_);

			int ierr=VecSetValues(b_,n,positions,&elementVector[0],ADD_VALUES);
			CHKERRV(ierr);
		}

		for(Base::MeshManipulator::FaceIterator it=theMesh_->faceColBegin();it!=theMesh_->faceColEnd();++it){
			const Base::Element* elLeft((*it)->getPtrElementLeft()),*elRight((*it)->getPtrElementRight());
			int nLeft(elLeft->getNrOfBasisFunctions()),n(nLeft);
			if(elRight!=NULL)
				n+=elRight->getNrOfBasisFunctions();
			int positions[n];
			for(int i=0;i<n;++i){
				if(i<nLeft)
					positions[i]=i+startPositionsOfElementsInTheVector_[elLeft->getID()];
				else
					positions[i]=i-nLeft+startPositionsOfElementsInTheVector_[elRight->getID()];
			}
			elementVector.resize(n);
			(*it)->getFaceVector(elementVector,faceVectorID_);
			int ierr=VecSetValues(b_,n,positions,&elementVector[0],ADD_VALUES);
			CHKERRV(ierr);
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
			(*it)->setTimeLevelData(timeLevel,0,localData);
		}
		ierr=VecRestoreArray(b_,&data);
		CHKERRV(ierr);
	}
};
