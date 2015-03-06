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
#include <iostream>

#include "FaceCacheData.hpp"
#include "LinearAlgebra/Matrix.hpp"
#include "LinearAlgebra/NumericalVector.hpp"
#include "Logger.h"


Base::FaceData::FaceData(std::size_t numberOfDOF, std::size_t numberOfFaceMatrices, std::size_t numberOfFaceVectors) :
faceMatrix_(numberOfFaceMatrices), faceVector_(numberOfFaceVectors)
{
    //std::cout<<"numberOfFaceMatrices "<<numberOfFaceMatrices<<std::endl;
    //std::cout<<"FaceMatrix_ size "<<faceMatrix_.size()<<std::endl;
    //std::cout<<"numberOfFaceVectors "<<numberOfFaceVectors<<std::endl;
    //std::cout<<"faceVector_ size = "<<faceVector_.size()<<std::endl;

}

/// \param[in] matrix The standard matrix used to set the FaceMatrix.
/// \param[in] matrixID The index to specify which FaceMatrix should be set.
/// \details To set a FaceMatrix using a standard matrix we must also know the number of basis functions corresponding to the left and right element. Since this class has no access to these numbers, we will assume the number of basis functions are the same on both sides and the input matrix should be a square matrix. For safety and also efficiency it is advised to use the other version of this function instead, which takes a FaceMatrix as input. This is actually a dated function and should be removed.
void Base::FaceData::setFaceMatrix(const LinearAlgebra::Matrix& matrix, std::size_t matrixID)
{
    if (matrixID >= faceMatrix_.size())
    {
        std::cout << "Warning: Setting a face matrix that was not preallocated. If this is expected, please allocate more face matrixes in the mesh generator" << std::endl;
        faceMatrix_.resize(matrixID + 1);
    }
    
    // Check if the input matrix is square.
    logger.assert(matrix.getNRows() == matrix.getNCols(), "FaceMatrix is not square.");
    
    std::size_t nDOFLeft = std::size_t(matrix.getNRows() / 2);
    std::size_t nDOFRight = matrix.getNRows() - nDOFLeft;
    faceMatrix_[matrixID].resize(nDOFLeft, nDOFRight);
    faceMatrix_[matrixID].setEntireMatrix(matrix);
}

/// \param[in] matrix The standard matrix used to set the FaceMatrix.
/// \param[in] matrixID The index to specify which FaceMatrix should be set.
void Base::FaceData::setFaceMatrix(const Base::FaceMatrix &faceMatrix, std::size_t matrixID)
{
    if (matrixID >= faceMatrix_.size())
    {
        std::cout << "Warning: Setting a face matrix that was not preallocated. If this is expected, please allocate more face matrixes in the mesh generator" << std::endl;
        faceMatrix_.resize(matrixID + 1);
    }
    
    faceMatrix_[matrixID].resize(faceMatrix.getNrOfDegreesOfFreedom(Base::Side::LEFT), faceMatrix.getNrOfDegreesOfFreedom(Base::Side::RIGHT));
    faceMatrix_[matrixID] = faceMatrix;
}

/// \param[in] matrix The standard matrix which will be used to get the face matrix as one entire matrix.
/// \param[in] matrixID The index to specify which FaceMatrix to get.
/// \details To convert a face matrix into a standard matrix is slow and inefficient. It is advised to use the other version of this function that returns a FaceMatrix. This is actually a dated function and should be removed.
void Base::FaceData::getFaceMatrix(LinearAlgebra::Matrix& matrix, std::size_t matrixID) const
{
    // Check if there are enough faces matrices stored.
    logger.assert(matrixID < faceMatrix_.size(), "Not enough face matrices stored.");
    
    matrix = faceMatrix_[matrixID].getEntireMatrix();
}

/// \param[in] matrixID The index to specify which FaceMatrix to get.
const Base::FaceMatrix & Base::FaceData::getFaceMatrix(std::size_t matrixID) const
{
    // Check if there are enough faces matrices stored.
    logger.assert(matrixID < faceMatrix_.size(), "Not enough face matrices stored.");
    
    return faceMatrix_[matrixID];
}

void Base::FaceData::setFaceVector(const LinearAlgebra::NumericalVector& vector, std::size_t vectorID)
{
    if (vectorID >= faceVector_.size())
    {
        std::cout << "Warning: Setting a face vector that was not preallocated. If this is expected, please allocate more face vectors in the mesh generator" << std::endl;
        faceVector_.resize(vectorID + 1);
    }
    faceVector_[vectorID].resize(vector.size());
    faceVector_[vectorID] = vector;
}

void Base::FaceData::getFaceVector(LinearAlgebra::NumericalVector& vector, std::size_t vectorID) const
{
    logger.assert(vectorID < faceVector_.size(), "insufficient face vectors stored");
    vector = faceVector_[vectorID];
}

const LinearAlgebra::NumericalVector& Base::FaceData::getResidue() const
{
    return residual_;
}

void Base::FaceData::setResidue(LinearAlgebra::NumericalVector& residue)
{
    residual_ = residue;
}

