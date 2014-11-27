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

#include "ElementData.hpp"
#include "TestErrorDebug.hpp"

#include "LinearAlgebra/Matrix.hpp"
#include "LinearAlgebra/NumericalVector.hpp"
#include <iostream>

namespace Base
{
    ElementData::ElementData(unsigned int timeLevels,
                             unsigned int nrOfUnkowns,
                             unsigned int nrOfBasisFunctions,
                             unsigned int nrOfElementMatrixes,
                             unsigned int nrOfElementVectors) :
    timeLevels_(timeLevels),
    nrOfUnkowns_(nrOfUnkowns),
    nrOfBasisFunctions_(nrOfBasisFunctions),
    expansionCoefficients_(timeLevels_),
    userData_(nullptr),
    elementMatrix_(nrOfElementMatrixes),
    elementVector_(nrOfElementVectors) { }
    
    
    void ElementData::setElementMatrix(const LinearAlgebra::Matrix& matrix, int matrixID)
    {
        if (matrixID >= elementMatrix_.size())
        {
            std::cout << "Warning: Setting an element matrix that was not preallocated. If this is expected, please allocate more element matrixes in the mesh generator" << std::endl;
            elementMatrix_.resize(matrixID + 1);
        }
        elementMatrix_[matrixID] = matrix;
    }
    
    void ElementData::getElementMatrix(LinearAlgebra::Matrix& matrix, int matrixID) const
    {
        TestErrorDebug(matrixID < elementMatrix_.size(), "insufficient element matrixes stored");
        matrix = elementMatrix_[matrixID];
    }
    
    void ElementData::setElementVector(const LinearAlgebra::NumericalVector& vector, int vectorID)
    {
        if (vectorID >= elementVector_.size())
        {
            std::cout << "Warning: Setting an element vector that was not preallocated. If this is expected, please allocate more element vectors in the mesh generator" << std::endl;
            elementVector_.resize(vectorID + 1);
        }
        elementVector_[vectorID] = vector;
    }
    
    void ElementData::getElementVector(LinearAlgebra::NumericalVector& vector, int vectorID) const
    {
        TestErrorDebug(vectorID < elementVector_.size(), "insufficient element vectors stored");
        vector = elementVector_[vectorID];
    }
    
    void ElementData::setNumberOfBasisFunctions(unsigned int number)
    {
        nrOfBasisFunctions_ = number;
    }
    
    const LinearAlgebra::Matrix& ElementData::getTimeLevelData(size_t timeLevel) const
    {
        if (timeLevel < timeLevels_)
        {
            if (expansionCoefficients_[timeLevel].size() != nrOfUnkowns_ * nrOfBasisFunctions_)
            {
                const_cast<LinearAlgebra::Matrix *> (&expansionCoefficients_[timeLevel])->resize(nrOfUnkowns_, nrOfBasisFunctions_);
            }
            return expansionCoefficients_[timeLevel];
        }
        else
        {
            throw "Error: Asked for a time level greater than the amount of time levels";
        }
    }
    
    double ElementData::getData(unsigned int timeLevel, unsigned int unknown, unsigned int basisFunction) const
    {
        if (timeLevel < timeLevels_ && unknown < nrOfUnkowns_ * nrOfBasisFunctions_)
        {
            if (expansionCoefficients_[timeLevel].size() != nrOfUnkowns_ * nrOfBasisFunctions_)
            {
                const_cast<LinearAlgebra::Matrix *> (&expansionCoefficients_[timeLevel])->resize(nrOfUnkowns_, nrOfBasisFunctions_);
            }
            return expansionCoefficients_[timeLevel](unknown, basisFunction);
        }
        else
        {
            throw "Error: Asked for a time level, or unknown, greater than the amount of time levels";
        }
    }
    
    void ElementData::setData(unsigned int timeLevel, unsigned int unknown, unsigned int basisFunction, double val)
    {
        if (timeLevel < timeLevels_ && unknown < nrOfUnkowns_ * nrOfBasisFunctions_)
        {
            if (expansionCoefficients_[timeLevel].size() != nrOfUnkowns_ * nrOfBasisFunctions_)
            {
                expansionCoefficients_[timeLevel].resize(nrOfUnkowns_, nrOfBasisFunctions_);
            }
            expansionCoefficients_[timeLevel](unknown, basisFunction) = val;
        }
        else
        {
            throw "Error: Asked for a time level, or unknown, greater than the amount of time levels";
        }
    }

    ///Rewrite with swap!!! and for all variables immediately
    void ElementData::setTimeLevelData(unsigned int timeLevel, unsigned int solutionId, const LinearAlgebra::NumericalVector& unknown)
    {
        if (timeLevel < timeLevels_ && solutionId < nrOfUnkowns_)
        {
            if (expansionCoefficients_[timeLevel].size() != nrOfUnkowns_ * nrOfBasisFunctions_)
            {
                expansionCoefficients_[timeLevel].resize(nrOfUnkowns_, nrOfBasisFunctions_);
            }
            LinearAlgebra::Matrix& mat = expansionCoefficients_[timeLevel];

            for (size_t i = 0; i < unknown.size(); ++i)
            {
                mat(solutionId, i) = unknown[i];
            }
        }
        else
        {
            throw "Error: Asked for a time level, or unknown, greater than the amount of time levels";
        }
    }
    
    void ElementData::setTimeLevelData(unsigned int timeLevel, const LinearAlgebra::Matrix& unknown)
    {
        if (timeLevel < timeLevels_)
        {
            expansionCoefficients_[timeLevel] = LinearAlgebra::Matrix(unknown);
        }
        else
        {
            throw "Error: Asked for a time level, or unknown, greater than the amount of time levels";
        }
    }
    
    int ElementData::getNrOfUnknows() const
    {
        return nrOfUnkowns_;
    }
    
    int ElementData::getNrOfBasisFunctions() const
    {
        return nrOfBasisFunctions_;
    }
    
    const typename ElementData::VectorOfDoubles& ElementData::getResidue() const
    {
        return residue_;
    }
    
    void ElementData::setResidue(VectorOfDoubles& residue)
    {
        residue_ = residue;
    }
    
    void ElementData::setUserData(UserElementData* data)
    {
        userData_ = data;
    }
}
