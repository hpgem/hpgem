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
#ifndef ElementData_hpp
#define ElementData_hpp
//----------------------------------------------------------------
#include <vector>
#include "LinearAlgebra/Matrix.hpp"
#include "LinearAlgebra/NumericalVector.hpp"

namespace LinearAlgebra
{
    class NumericalVector;
}

struct UserElementData;

namespace Base
{

    class ElementData
    {
        /*!
         * ElementData contains the data for every element in the vector elementData_.
         * This is a two dimensional vector where the first dimensions is the time level,
         * and the second the number of unknowns times the number of basis functions.
         */
    public:
        using MatrixT = LinearAlgebra::Matrix;
        using VectorOfMatrices = std::vector<LinearAlgebra::Matrix>;

        ElementData(std::size_t timeLevels, std::size_t nrOfUnkowns, std::size_t nrOfBasisFunctions, std::size_t nrOfElementMatrixes = 0, std::size_t nrOfElementVectors = 0);

        virtual ~ ElementData() { }

        ///Set/update the element matrix.
        void setElementMatrix(const LinearAlgebra::Matrix&, int matrixID = 0);

        /// Get the element matrix corresponding to the given matrixiD.
        virtual const LinearAlgebra::Matrix & getElementMatrix(std::size_t matrixID = 0) const;
        
        // Dated version of getElementMatrix. Use the version that returns a reference instead!
        // virtual void getElementMatrix(LinearAlgebra::Matrix&, int matrixID = 0) const;

        /// \brief Returns (and creates if unavailable) the expansion coefficients correspdonding to a given timelevel in the form of a matrix. Better use getTimeLevelDataVector instead!
        LinearAlgebra::Matrix getTimeLevelDataMatrix(std::size_t timeLevel);
    
        virtual void setElementVector(const LinearAlgebra::NumericalVector& vector, int vectorID = 0);

        virtual void getElementVector(LinearAlgebra::NumericalVector&, int vectorID = 0) const;

        /// \brief Sets (and creates if unavailable) the expansion coefficients correspdonding to the given time level.
        void setTimeLevelDataVector(std::size_t timeLevel, LinearAlgebra::NumericalVector &val);

        /// \brief Returns a reference to the expansion coefficients correspdonding to the given time level.
        const LinearAlgebra::NumericalVector & getTimeLevelDataVector(std::size_t timeLevel) const;
        LinearAlgebra::NumericalVector & getTimeLevelDataVector(std::size_t timeLevel);
        
        /// \brief Specify a time level index and variable index, return a vector containing the corresponding expansion coefficients.
        virtual const LinearAlgebra::NumericalVector getTimeLevelData(std::size_t timeLevel, std::size_t unknown = 0) const;

        /// \brief Specify a time level index and a variable index (unknown), set the corresponding expansionCoefficients. Better use getTimeLevelDataVector if possible!
        void setTimeLevelData(std::size_t timeLevel, std::size_t unknown, const LinearAlgebra::NumericalVector &val);
        void setTimeLevelData(std::size_t timeLevel, const LinearAlgebra::NumericalVector &val);
        
        void setCurrentData(const LinearAlgebra::NumericalVector& data);
        
        LinearAlgebra::NumericalVector& getCurrentData();

        /// \brief Specify a time level index, a variabale index and a basis function index, return the corresponding expansionCoefficient (double).
        virtual double getData(std::size_t timeLevel, std::size_t unknown, std::size_t basisFunction) const;

        /// \brief Specify a time level index, a variable index and a basis function index, set the corresponding expansionCoefficient (double).
        void setData(std::size_t timeLevel, std::size_t unknown, std::size_t basisFunction, double val);

        virtual int getNrOfUnknows() const;

        virtual int getNrOfBasisFunctions() const;

        //this needs to store information about all variables, so it needs to be a matrix (?)
        virtual const LinearAlgebra::NumericalVector& getResidue() const;

        void setResidue(LinearAlgebra::NumericalVector& residue);

        void setUserData(UserElementData* data);

        virtual UserElementData* getUserData() const;
        
        /// \brief Convert the index corresponding to the basis function (iBasisFunction) and the index corresponding to the variable (iVar) to a single index.
        /// \param[in] iVar The index corresponding to the variable.
        /// \param[in] iBasisFunction The index corresponding to the basisfunction.
        virtual const std::size_t convertToSingleIndex(std::size_t iBasisFunction, std::size_t iVar = 0) const
        {
            return iVar * nrOfBasisFunctions_ + iBasisFunction;
        }

    protected:
        void setNumberOfBasisFunctions(std::size_t number);

    private:
        std::size_t timeLevels_;
        std::size_t nrOfUnknowns_;
        std::size_t nrOfBasisFunctions_;

        /// \brief Stores the expansion coefficients.
        /// \details The value expansionCoefficients_(iT)(iVB) is the expansion coefficient corresponding to the solution at time level iT and vector basisfunction iVB. Index iVB satisfies iVB = convertToSingleIndex(iB,iV), where iB is the index corresponding to the basis function and iVB the index corresponding to the variable.
        std::vector<LinearAlgebra::NumericalVector> expansionCoefficients_;

        ///Stores the result of an element integration
        LinearAlgebra::NumericalVector residue_;
        
        ///Store temporary data
        LinearAlgebra::NumericalVector currentData_;

        ///Stores polymorphic pointer to UserDefined Data, internally not used! Used only outside of the Kernel!!!
        UserElementData* userData_;

        ///Stores element matrix(es) for this element
        VectorOfMatrices elementMatrix_;

        std::vector<LinearAlgebra::NumericalVector> elementVector_;
    };
}
#endif
