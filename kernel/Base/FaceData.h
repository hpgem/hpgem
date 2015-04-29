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

//----------------------------------------------------------------
#ifndef FaceData_h
#define FaceData_h
//----------------------------------------------------------------
#include <vector>
#include "Base/FaceMatrix.h"
#include "LinearAlgebra/Matrix.h"
#include "LinearAlgebra/NumericalVector.h"
#include "FaceCacheData.h"

namespace LinearAlgebra
{
    class NumericalVector;
}

namespace Base
{
    
    struct FaceCacheData;
    class UserFaceData;
    
    class FaceData
    {
    public:
        using CacheT = FaceCacheData;
        using VecCacheT = std::vector<CacheT>;

    public:
        FaceData(std::size_t numberOfDOF, std::size_t numberOfFaceMatrices = 0, std::size_t numberOfFaceVactors = 0);        
        
        FaceData(const FaceData& other);
        
        virtual ~FaceData()
        {
        }

        /// \brief Sets face matrix number 'matrixID' using a standard matrix.
        /// \deprecated For safety and also efficiency it is advised to use the other version
        /// of this function instead, which takes a FaceMatrix as input. This is actually
        /// a dated function and should be removed.
        void setFaceMatrix(const LinearAlgebra::Matrix& matrix, std::size_t matrixID = 0);

        /// \brief Sets face matrix number 'matrixID' using a standard matrix.
        void setFaceMatrix(const FaceMatrix &faceMatrix, std::size_t matrixID = 0);

        /// \brief Gets face matrix number 'matrixID' and return it as a standard matrix. It is advised to use the other version instead, which returns a FaceMatrix.
        virtual LinearAlgebra::Matrix getFaceMatrixMatrix(std::size_t matrixID = 0) const;

        /// \brief Returns face matrix number 'matrixID'.
        const FaceMatrix & getFaceMatrix(std::size_t matrixID = 0) const;

        void setFaceVector(const LinearAlgebra::NumericalVector& vector, std::size_t vectorID = 0);

        virtual LinearAlgebra::NumericalVector getFaceVector(std::size_t vectorID = 0) const;

        virtual VecCacheT& getVecCacheData()
        {
            return vecCacheData_;
        }        
        
        virtual UserFaceData* getUserData() const
        {
            return userData_;
        }
        
        void setUserData(UserFaceData* data)
        {
            //the user may pass any kind of data he/she wants (including nullptr) even if this does not seem to make sense
            userData_ = data;
        }
        
        virtual const LinearAlgebra::NumericalVector& getResidue() const;

        void setResidue(LinearAlgebra::NumericalVector& residue);
        
        std::size_t getNumberFaceMatrices() const;
        std::size_t getNumberFaceVectors() const;
        
    protected:
        VecCacheT vecCacheData_;

    private:
        
        UserFaceData* userData_;
        std::vector<FaceMatrix> faceMatrix_;
        std::vector<LinearAlgebra::NumericalVector> faceVector_;

        //a concatenation of the flux contributions to the residuals in the left and the right elements
        LinearAlgebra::NumericalVector residual_;
    };
}
#endif
