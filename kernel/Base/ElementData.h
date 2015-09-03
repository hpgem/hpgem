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
#ifndef ElementData_h
#define ElementData_h
//----------------------------------------------------------------
#include <vector>
#include "LinearAlgebra/MiddleSizeMatrix.h"
#include "LinearAlgebra/MiddleSizeVector.h"

namespace LinearAlgebra
{
    class MiddleSizeVector;
}

struct UserElementData;

namespace Base
{
    
    class ElementData
    {
        
    public:
        using MatrixT = LinearAlgebra::MiddleSizeMatrix;
        using VectorOfMatrices = std::vector<LinearAlgebra::MiddleSizeMatrix>;

        ElementData(std::size_t timeLevels, std::size_t numberOfUnkowns, std::size_t numberOfBasisFunctions, std::size_t numberOfElementMatrixes = 0, std::size_t numberOfElementVectors = 0);
        
        ElementData(const ElementData& other);

        virtual ~ElementData() = default;
        
        /// \brief Set/update the element matrix.
        void setElementMatrix(const LinearAlgebra::MiddleSizeMatrix &, std::size_t matrixID = 0);

        /// \brief Get the element matrix corresponding to the given matrixiD.
        const LinearAlgebra::MiddleSizeMatrix &getElementMatrix(std::size_t matrixID = 0) const;

        /// \brief Get the element matrix corresponding to the given matrixiD.
        LinearAlgebra::MiddleSizeMatrix & getElementMatrix(std::size_t matrixID = 0);

        /// \brief Set the element vector corresponding to the given vectorID.
        void setElementVector(const LinearAlgebra::MiddleSizeVector &vector, std::size_t vectorID = 0);

        /// \brief Get the element vector corresponding to the given vectorID.
        LinearAlgebra::MiddleSizeVector getElementVector(std::size_t vectorID = 0) const;

        
        /// \brief Set the expansion coefficients corresponding to the given time level.
        void setTimeLevelDataVector(std::size_t timeLevel, LinearAlgebra::MiddleSizeVector &val);

        /// \brief Returns a reference to the expansion coefficients corresponding to the given time level.
        const LinearAlgebra::MiddleSizeVector& getTimeLevelDataVector(std::size_t timeLevel) const;
        LinearAlgebra::MiddleSizeVector& getTimeLevelDataVector(std::size_t timeLevel);

        /// \brief Specify a time level index and variable index, return a vector containing the corresponding expansion coefficients.
        const LinearAlgebra::MiddleSizeVector getTimeLevelData(std::size_t timeLevel, std::size_t unknown = 0) const;

        /// \brief Specify a time level index and a variable index (unknown), set the corresponding expansionCoefficients. Better use getTimeLevelDataVector if possible because that's faster!
        void setTimeLevelData(std::size_t timeLevel, std::size_t unknown, const LinearAlgebra::MiddleSizeVector &val);
        void setTimeLevelData(std::size_t timeLevel, const LinearAlgebra::MiddleSizeVector &val);

        /// \brief Specify a time level index, a variable index and a basis function index, return the corresponding expansionCoefficient (double).
        LinearAlgebra::MiddleSizeVector::type getData(std::size_t timeLevel, std::size_t unknown, std::size_t basisFunction) const;
        
        /// \brief Specify a time level index, a variable index and a basis function index, set the corresponding expansionCoefficient (double).
        void setData(std::size_t timeLevel, std::size_t unknown, std::size_t basisFunction, double val);
         
        
        /// \brief Set the number of time integration vectors.
        void setNumberOfTimeIntegrationVectors(const std::size_t numberOfTimeIntegrationVectors)
        {
            timeIntegrationVectors_.resize(numberOfTimeIntegrationVectors, LinearAlgebra::MiddleSizeVector(numberOfUnknowns_ * numberOfBasisFunctions_));
        }
        
        /// \brief Get the number of time integration vectors.
        std::size_t getNumberOfTimeIntegrationVectors()
        {
            return timeIntegrationVectors_.size();
        }
        
        /// \brief Return a reference to a time integration vector.
        const LinearAlgebra::MiddleSizeVector & getTimeIntegrationVector(std::size_t timeIntegrationVectorId) const;
        LinearAlgebra::MiddleSizeVector& getTimeIntegrationVector(std::size_t timeIntegrationVectorId);
        
        /// \brief Set a time integration vector corresponding to the given id.
        void setTimeIntegrationVector(std::size_t timeIntegrationVectorId, LinearAlgebra::MiddleSizeVector &val);

        /// \brief Return a subvector corresponding to the given unknown (variable id).
        const LinearAlgebra::MiddleSizeVector getTimeIntegrationSubvector(std::size_t timeIntegrationVectorId, std::size_t unknown) const;
        
        /// \brief Set a subvector corresponding to the given unknown (variable id).
        void setTimeIntegrationSubvector(std::size_t timeIntegrationVectorId, std::size_t unknown, LinearAlgebra::MiddleSizeVector val);
    
        /// \brief Return the value corresponding to the given time integration vector index, unknown (variable id) and basis function index
        LinearAlgebra::MiddleSizeVector::type getTimeIntegrationData(std::size_t timeIntegrationVectorId, std::size_t unknown, std::size_t basisFunction);
        
        /// \brief Set the value corresponding to the given time integration vector index, unknown (variable id) and basis function index
        void setTimeIntegrationData(std::size_t timeIntegrationVectorId, std::size_t unknown, std::size_t basisFunction, double val);
        
        
        ///\deprecated Spelling mistake, please use getNumberOfUnknowns instead.
        std::size_t getNrOfUnknows() const;

        ///\deprecated Does not follow naming guidelines, use getNumberOfUnknowns instead.
        std::size_t getNrOfUnknowns() const;

        ///\deprecated Does not follow naming guidelines, use getNumberOfBasisFunctions instead.
        std::size_t getNrOfBasisFunctions() const;

        std::size_t getNumberOfUnknowns() const;

        ///Returns the number of basis functions that are non-zero inside this element
        std::size_t getNumberOfBasisFunctions() const;

        void setUserData(UserElementData* data);

        UserElementData* getUserData() const;

        /// \brief Convert the index corresponding to the basis function (iBasisFunction) 
        /// and the index corresponding to the variable (iVar) to a single index.
        /// \param[in] iVar The index corresponding to the variable.
        /// \param[in] iBasisFunction The index corresponding to the basisfunction.
        std::size_t convertToSingleIndex(std::size_t iBasisFunction, std::size_t iVar = 0) const
        {
            return iVar * numberOfBasisFunctions_ + iBasisFunction;
        }
        
    protected:
        void setNumberOfBasisFunctions(std::size_t number);

    private:
        /// The number of time levels. Multiple time levels can be used to store intermediate results, help variables and residues.
        std::size_t timeLevels_;
        
        /// The number of variables (unknowns).
        std::size_t numberOfUnknowns_;
        
        /// The number of basis functions
        std::size_t numberOfBasisFunctions_;

        /// \brief Stores the expansion coefficients.
        /// \details The value expansionCoefficients_[iT](iVB) is the expansion
        /// coefficient corresponding to the solution at time level iT and vector-
        /// basisfunction iVB. Index iVB satisfies iVB = convertToSingleIndex(iB,iV), 
        /// where iB is the index corresponding to the scalar basis function and iV the
        /// index corresponding to the variable.
        std::vector<LinearAlgebra::MiddleSizeVector> expansionCoefficients_;
        
        /// \brief Vectors used for the time integration
        /// \details The value timeIntegrationVectors_[i](iVB) is the expansion
        /// coefficient corresponding to vector i and vector-basisfunction iVB.
        /// Index iVB satisfies iVB = convertToSingleIndex(iB,iV),
        /// where iB is the index corresponding to the scalar basis function and iV the
        /// index corresponding to the variable.
        std::vector<LinearAlgebra::MiddleSizeVector> timeIntegrationVectors_;

        ///Stores polymorphic pointer to UserDefined Data, internally not used.
        ///Used only outside of the Kernel.
        UserElementData* userData_;

        ///Stores element matrix(es) for this element
        VectorOfMatrices elementMatrix_;
        
        ///Stores element vector(s) for this element
        std::vector<LinearAlgebra::MiddleSizeVector> elementVector_;
    };
}
#endif
