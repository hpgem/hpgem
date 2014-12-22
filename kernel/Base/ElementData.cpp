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
#include "Element.hpp"

#include "LinearAlgebra/Matrix.hpp"
#include "LinearAlgebra/NumericalVector.hpp"
#include <iostream>
#include <functional>

namespace Base
{
    ElementData::ElementData(unsigned int timeLevels,
                             unsigned int nrOfUnknowns,
                             unsigned int nrOfBasisFunctions,
                             unsigned int nrOfElementMatrixes,
                             unsigned int nrOfElementVectors) :
    timeLevels_(timeLevels),
    nrOfUnknowns_(nrOfUnknowns),
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
    
    const LinearAlgebra::NumericalVector ElementData::getTimeLevelData(size_t timeLevel) const
    {
        if (timeLevel < timeLevels_)
        {
            LinearAlgebra::Matrix thisLevelMatrix = expansionCoefficients_[timeLevel];
            LinearAlgebra::NumericalVector coeffVec(nrOfBasisFunctions_);
            for (size_t i = 0; i < nrOfBasisFunctions_; ++i)
            {
                coeffVec[i] = thisLevelMatrix[i];
            }
            return coeffVec;
        }
        else
        {
            throw "Error: Asked for a time level greater than the amount of time levels";
        }
    }
    
    double ElementData::getData(unsigned int timeLevel, unsigned int unknown, unsigned int basisFunction) const
    {
        if (timeLevel < timeLevels_ && unknown < nrOfUnknowns_ && basisFunction < nrOfBasisFunctions_)
        {
            return expansionCoefficients_[timeLevel](unknown, basisFunction);
        }
        else
        {
            throw "Error: Asked for a time level, or unknown, greater than the amount of time levels";
        }
    }
    
    void ElementData::setData(unsigned int timeLevel, unsigned int unknown, unsigned int basisFunction, double val)
    {
        if (timeLevel < timeLevels_ && unknown < nrOfUnknowns_ && basisFunction < nrOfBasisFunctions_)
        {
            if (expansionCoefficients_[timeLevel].size() != nrOfUnknowns_ * nrOfBasisFunctions_)
            {
                expansionCoefficients_[timeLevel].resize(nrOfUnknowns_, nrOfBasisFunctions_);
            }
            expansionCoefficients_[timeLevel](unknown, basisFunction) = val;
        }
        else
        {
            throw "Error: Asked for a time level, or unknown, greater than the amount of time levels";
        }
    }
    
    /*
    ** -- Date Friday 19 Dec 2014
    ** -- Developer @dducks
    **
    ** Guys - can we have a small inline debate about the return type of getTimeLevelData?
    ** 
    ** So, getTimeLevelData / setTimeLevelData convert their elements into a NumericalVector from a Matrix
    ** Why the hell do we do that? This makes additional copies, and causes issues with the MPI synchronisation, since
    ** everything is done asynchronous. This makes me sad and stuff.
    **
    ** Also, there seems to be an asymmetry for solutionId / etc. 
    **
    ** For now I've reimplemented these methods which return the direct matrix. PloxFix? <3
    **
    */
    /**
      \brief Returns (and creates if unavailable) the time level data matrix for this element
      
      This method returns the TimeLevelData matrix present in this Element for the given timeLevel
      If this matrix does not exist yet (or better said, is of dimension 0x0), it will be initialised
      with the proper dimensions.
      
      \arg timeLevel the corresponding timeLevel
      \return a reference to the actual matrix
    */
    LinearAlgebra::Matrix& ElementData::getTimeLevelDataMatrix(std::size_t timeLevel) 
    {
      if (timeLevel < timeLevels_)
      {
        // The matrix can be of dimension 0x0 if it hasn't been used before.
        // So lets resize it first!
        if (expansionCoefficients_[timeLevel].size() != nrOfUnknowns_ * nrOfBasisFunctions_)
        {
          expansionCoefficients_[timeLevel].resize(nrOfUnknowns_, nrOfBasisFunctions_);
        }
        
        return expansionCoefficients_[timeLevel];
      }
      else
      {
        throw "Error: Asked for a time level greater than the amount of time levels!";
      }
      
    }
    

    ///Rewrite with swap!!! and for all variables immediately
    void ElementData::setTimeLevelData(unsigned int timeLevel, unsigned int solutionId, const LinearAlgebra::NumericalVector& unknown)
    {
        if (timeLevel < timeLevels_ && solutionId < nrOfUnknowns_)
        {
            if (expansionCoefficients_[timeLevel].size() != nrOfUnknowns_ * nrOfBasisFunctions_)
            {
                expansionCoefficients_[timeLevel].resize(nrOfUnknowns_, nrOfBasisFunctions_);
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
    
    void ElementData::setTimeLevelData(unsigned int timeLevel, const LinearAlgebra::NumericalVector& unknown)
    {
        if (timeLevel < timeLevels_)
        {
            if (expansionCoefficients_[timeLevel].size() != nrOfUnknowns_ * nrOfBasisFunctions_)
            {
                expansionCoefficients_[timeLevel].resize(nrOfUnknowns_, nrOfBasisFunctions_);
            }            
            
            LinearAlgebra::Matrix& mat = expansionCoefficients_[timeLevel];

            for (std::size_t i = 0; i < unknown.size(); ++i)
            {
                mat[i] = unknown[i];
            }
        }
        else
        {
            throw "Error: Asked for a time level, or unknown, greater than the amount of time levels";
        }
    }
    
    int ElementData::getNrOfUnknows() const
    {
        return nrOfUnknowns_;
    }
    
    int ElementData::getNrOfBasisFunctions() const
    {
        return nrOfBasisFunctions_;
    }
    
    const typename LinearAlgebra::NumericalVector& ElementData::getResidue() const
    {
        return residue_;
    }
    
    void ElementData::setResidue(LinearAlgebra::NumericalVector& residue)
    {
        residue_ = residue;
    }
    
    void ElementData::setUserData(UserElementData* data)
    {
        userData_ = data;
    }
}
