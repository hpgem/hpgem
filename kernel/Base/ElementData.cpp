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
#include "Logger.h"
#include "Element.hpp"

#include "LinearAlgebra/Matrix.hpp"
#include "LinearAlgebra/NumericalVector.hpp"
#include <iostream>
#include <functional>

namespace Base
{
    ElementData::ElementData(std::size_t timeLevels,
                             std::size_t nrOfUnknowns,
                             std::size_t nrOfBasisFunctions,
                             std::size_t nrOfElementMatrixes,
                             std::size_t nrOfElementVectors) :
    timeLevels_(timeLevels),
    nrOfUnknowns_(nrOfUnknowns),
    nrOfBasisFunctions_(nrOfBasisFunctions),
    expansionCoefficients_(timeLevels_),
    userData_(nullptr),
    elementMatrix_(nrOfElementMatrixes),
    elementVector_(nrOfElementVectors) {
        //std::cout<<"nrOfElementMatrixes "<<nrOfElementMatrixes<<std::endl;
        //std::cout<<"elementMatrix_ size "<<elementMatrix_.size()<<std::endl;
        //std::cout<<"nrOfElementVectors "<<nrOfElementVectors<<std::endl;
        //std::cout<<"elementVector_ size = "<<elementVector_.size()<<std::endl;
    }
    
    
    void ElementData::setElementMatrix(const LinearAlgebra::Matrix& matrix, std::size_t matrixID)
    {
        //std::cout<<"matrix ID = "<<matrixID<<std::endl;
        //std::cout<<"elementMatrix_ size = "<<elementMatrix_.size()<<std::endl;
        if (matrixID >= elementMatrix_.size())
        {
            std::cout << "Warning: Setting an element matrix that was not preallocated. If this is expected, please allocate more element matrixes in the mesh generator" << std::endl;
            elementMatrix_.resize(matrixID + 1);
        }
        elementMatrix_[matrixID] = matrix;
    }
    
    const LinearAlgebra::Matrix & ElementData::getElementMatrix(std::size_t matrixID) const{
        logger.assert(matrixID < elementMatrix_.size(), "Requested matrix %, "
                "while there are only % matrices for this element.", matrixID, 
                      elementMatrix_.size());
        return elementMatrix_[matrixID];
    }
    
    void ElementData::setElementVector(const LinearAlgebra::NumericalVector& vector, std::size_t vectorID)
    {
        //std::cout<<"VectorID : "<<vectorID<<std::endl;
        //std::cout<<"elementVector size : "<<elementVector_.size()<<std::endl;
        
        if (vectorID >= elementVector_.size())
        {
            //std::cout << "Warning: Setting an element vector that was not preallocated. If this is expected, please allocate more element vectors in the mesh generator" << std::endl;
            elementVector_.resize(vectorID + 1);
        }
        elementVector_[vectorID] = vector;
    }
    
    LinearAlgebra::NumericalVector ElementData::getElementVector(std::size_t vectorID) const
    {
        logger.assert(vectorID < elementVector_.size(), "insufficient element vectors stored");
        return elementVector_[vectorID];
    }

    void ElementData::setNumberOfBasisFunctions(std::size_t number)
    {
        nrOfBasisFunctions_ = number;
    }    
   
    std::size_t ElementData::getNrOfBasisFunctions() const
    {
        return nrOfBasisFunctions_;
    }
    
    /// \param[in] timeLevel Index corresponding to the time level.
    /// \param[in] unknown Index corresponding to the variable.
    /// \param[in] basisFunction Index corresponding to the basisFunction.
    /// \param[in] val Value to set the expansionCoeffient.
    void ElementData::setData(std::size_t timeLevel, std::size_t unknown, std::size_t basisFunction, double val)
    {
        if (timeLevel < timeLevels_ && unknown < nrOfUnknowns_ && basisFunction < nrOfBasisFunctions_)
        {
            if (expansionCoefficients_[timeLevel].size() != nrOfUnknowns_ * nrOfBasisFunctions_)
            {
                expansionCoefficients_[timeLevel].resize(nrOfUnknowns_ * nrOfBasisFunctions_);
            }
            expansionCoefficients_[timeLevel](convertToSingleIndex(basisFunction, unknown)) = val;
        }
        else
        {
            throw "Error: Asked for a time level, or unknown, greater than the amount of time levels";
        }
    }
    
    /// \param[in] timeLevel Index corresponding to the time level.
    /// \param[in] unknown Index corresponding to the variable.
    /// \param[in] basisFunction Index corresponding to the basisFunction.
    double ElementData::getData(std::size_t timeLevel, std::size_t unknown, std::size_t basisFunction) const
    {
        if (timeLevel < timeLevels_ && unknown < nrOfUnknowns_ && basisFunction < nrOfBasisFunctions_)
        {
            logger.assert(expansionCoefficients_[timeLevel].size() == nrOfUnknowns_
            * nrOfBasisFunctions_, "Wrong number of expansion coefficients.");
            return expansionCoefficients_[timeLevel](convertToSingleIndex(basisFunction, unknown));

        }
        else
        {
            throw "Error: Asked for a time level, or unknown, greater than the amount of time levels";
        }
    }
    
    /// \details Rewrite with swap!!! This method is slow and inefficient. It is therefore advised to use getTimeLevelDataVector instead.
    /// \param[in] timeLevel Index corresponding to the time level.
    /// \param[in] unknown Index corresponding to the variable.
    /// \param[in] val Vector of values to set the expansionCoeffient corresponding to the given unknown and time level.
    void ElementData::setTimeLevelData(std::size_t timeLevel, std::size_t unknown, const LinearAlgebra::NumericalVector& val)
    {
        if (timeLevel < timeLevels_ && unknown < nrOfUnknowns_)
        {
            if (expansionCoefficients_[timeLevel].size() != nrOfUnknowns_ * nrOfBasisFunctions_)
            {
                expansionCoefficients_[timeLevel].resize(nrOfUnknowns_ * nrOfBasisFunctions_);
            }

            for (std::size_t iB = 0; iB < val.size(); ++iB) // iB = iBasisFunction
            {
                expansionCoefficients_[timeLevel](convertToSingleIndex(iB,unknown)) = val[iB];
            }
        }
        else
        {
            throw "Error: Asked for a time level, or unknown, greater than the amount of time levels";
        }
    }
    
    void ElementData::setTimeLevelData(std::size_t timeLevel, const LinearAlgebra::NumericalVector& val)
    {
        setTimeLevelData(timeLevel, 0, val);
    }    

    /// \param[in] timeLevel Index corresponding to the time level.
    /// \param[in] unknown Index corresponding to the variable.
    const LinearAlgebra::NumericalVector ElementData::getTimeLevelData(std::size_t timeLevel, std::size_t unknown) const
    {

        if (timeLevel < timeLevels_)
        {
            LinearAlgebra::NumericalVector timeLevelData(nrOfBasisFunctions_);
            for(std::size_t iB = 0; iB < nrOfBasisFunctions_; iB++) // iB = iBasisFunction
            {
                timeLevelData(iB) = expansionCoefficients_[timeLevel](convertToSingleIndex(iB, unknown));
            }
            return timeLevelData;
        }
        else
        {
            throw "Error: Asked for a time level greater than the amount of time levels";
        }
    } 
    
    
    /**
      \details This method returns the TimeLevelData present in this Element for the given timeLevel in the form of a matrix. If the data does not exist yet (or better said, is of dimension 0), it will be initialised with the proper dimension. Tbis method is slow since it needs to reshape a vector to a matrix and returns a copy. Therefore it is advised to use getTimeLevelDataVector.
      
      \param[in] timeLevel Index corresponding to the time level.
      \return A matrix M such that M(iV,iB) is the expansion coefficient corresponding to variable iV and basisfunction iB at the given time level.
    */
    LinearAlgebra::Matrix ElementData::getTimeLevelDataMatrix(std::size_t timeLevel)
    {
        if (timeLevel < timeLevels_)
        {
            // The vector can be of dimension 0 if it hasn't been used before.
            // So lets resize it first!
            if (expansionCoefficients_[timeLevel].size() != nrOfUnknowns_ * nrOfBasisFunctions_)
            {
                expansionCoefficients_[timeLevel].resize(nrOfUnknowns_ * nrOfBasisFunctions_);
            }
            LinearAlgebra::Matrix M(nrOfUnknowns_, nrOfBasisFunctions_);
            for(std::size_t iV = 0; iV < nrOfUnknowns_; iV++)   // iV = iVariable
            {
                for(std::size_t iB = 0; iB < nrOfBasisFunctions_; iB++) // iB = iBasisFunction
                {
                    M(iV,iB) = expansionCoefficients_[timeLevel](convertToSingleIndex(iB,iV));
                }
            }
            return M;
        }
        else
        {
            throw "Error: Asked for a time level greater than the amount of time levels!";
        }
    }
    
    /// \param[in] timeLevel Index corresponding to the time level.
    /// \param[in] val Vector of values to set the expansionCoeffient corresponding to the given unknown and time level.
    void ElementData::setTimeLevelDataVector(std::size_t timeLevel, LinearAlgebra::NumericalVector &val)
    {
        if (timeLevel < timeLevels_)
        {
            // The vector can be of dimension 0 if it hasn't been used before.
            // So lets resize it first!
            if (expansionCoefficients_[timeLevel].size() != nrOfUnknowns_ * nrOfBasisFunctions_)
            {
                expansionCoefficients_[timeLevel].resize(nrOfUnknowns_ * nrOfBasisFunctions_);
            }
            expansionCoefficients_[timeLevel] = val;
        }
        else
        {
            throw "Error: Asked for a time level greater than the amount of time levels!";
        }
    }
    
    /// \param[in] timeLevel Index corresponding to the time level.
    /// \return The expansion coefficient corresponding to the given time level.
    const LinearAlgebra::NumericalVector & ElementData::getTimeLevelDataVector(std::size_t timeLevel) const
    {
        if (timeLevel < timeLevels_)
        {
            logger.assert(expansionCoefficients_[timeLevel].size() == nrOfUnknowns_ * nrOfBasisFunctions_,
                          "Wrong number of expansion coefficients.");
            return expansionCoefficients_[timeLevel];
        }
        else
        {
            throw "Error: Asked for a time level greater than the amount of time levels!";
        }
    }

    LinearAlgebra::NumericalVector & ElementData::getTimeLevelDataVector(std::size_t timeLevel)
    {
        if (timeLevel < timeLevels_)
        {
            // The vector can be of dimension 0 if it hasn't been used before.
            // So lets resize it first!
            if (expansionCoefficients_[timeLevel].size() != nrOfUnknowns_ * nrOfBasisFunctions_)
            {
                expansionCoefficients_[timeLevel].resize(nrOfUnknowns_ * nrOfBasisFunctions_);
            }
            return expansionCoefficients_[timeLevel];
        }
        else
        {
            throw "Error: Asked for a time level greater than the amount of time levels!";
        }
    }
    
    void ElementData::setCurrentData(const LinearAlgebra::NumericalVector& data)
    {
        currentData_ = data;
    }
    
    LinearAlgebra::NumericalVector& ElementData::getCurrentData()
    {
        // The vector can be of dimension 0 if it hasn't been used before.
        // So lets resize it first!
        if (currentData_.size() != nrOfBasisFunctions_)
        {
            currentData_.resize(nrOfBasisFunctions_);
        }
        return currentData_;
    }
    
    std::size_t ElementData::getNrOfUnknows() const
    {
        return nrOfUnknowns_;
    }    
    
    void ElementData::setResidue(LinearAlgebra::NumericalVector& residue)
    {
        residue_ = residue;
    }
    
    const typename LinearAlgebra::NumericalVector& ElementData::getResidue() const
    {
        return residue_;
    }
    
    void ElementData::setUserData(UserElementData* data)
    {
        userData_ = data;
    }
    
    UserElementData* ElementData::getUserData() const
    {
        return userData_;
    }
}
