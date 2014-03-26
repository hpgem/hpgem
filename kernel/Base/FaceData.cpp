/*
 * faceData.cpp
 *
 *  Created on: Feb 5, 2014
 *      Author: brinkf
 */

#include "FaceData.hpp"
#include "TestErrorDebug.hpp"

Base::FaceData::FaceData(unsigned int numberOfDOF, unsigned int numberOfFaceMatrices, unsigned int numberOfFaceVectors):
	faceMatrix_(numberOfFaceMatrices),faceVector_(numberOfFaceVectors)
{
	for(std::vector<LinearAlgebra::Matrix>::iterator it=faceMatrix_.begin();it!=faceMatrix_.end();++it){
		it->resize(numberOfDOF,numberOfDOF);
	}
	for(std::vector<LinearAlgebra::NumericalVector>::iterator it=faceVector_.begin();it!=faceVector_.end();++it){
		it->resize(numberOfDOF);
	}
}

void
Base::FaceData::setFaceMatrix(const LinearAlgebra::Matrix& matrix, unsigned int matrixID)
{
	if(matrixID>=faceMatrix_.size()){
		std::cout<<"Warning: Setting a face matrix that was not preallocated. If this is expected, please allocate more face matrixes in the constructor";
		faceMatrix_.resize(matrixID+1);
	}
	faceMatrix_[matrixID].resize(matrix.getNRows(),matrix.getNCols());
	faceMatrix_[matrixID]=matrix;
}

void
Base::FaceData::getFaceMatrix(LinearAlgebra::Matrix& matrix, unsigned int matrixID) const
{
	TestErrorDebug(matrixID<faceMatrix_.size(),"insufficient face matrixes stored");
	matrix=faceMatrix_[matrixID];
}

void
Base::FaceData::setFaceVector(const LinearAlgebra::NumericalVector& vector, unsigned int vectorID)
{
	if(vectorID>=faceVector_.size()){
		std::cout<<"Warning: Setting a face vector that was not preallocated. If this is expected, please allocate more face vectors in the constructor";
		faceVector_.resize(vectorID+1);
	}
	faceVector_[vectorID].resize(vector.size());
	faceVector_[vectorID]=vector;
}

void
Base::FaceData::getFaceVector(LinearAlgebra::NumericalVector& vector, unsigned int vectorID) const
{
	TestErrorDebug(vectorID<faceVector_.size(),"insufficient face vectors stored");
	vector=faceVector_[vectorID];
}

