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

