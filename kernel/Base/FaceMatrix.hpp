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
#ifndef FaceMatrix_hpp
#define FaceMatrix_hpp
//----------------------------------------------------------------
#include "Base/Side.hpp"
#include "LinearAlgebra/Matrix.hpp"


namespace Base
{
    /// \class FaceMatrix
    /// \brief This class stores and can manipulate data for the face matrices.
    ///
    /// \details
    /// This class stores and can manipulate data for the face matrices. the entire face matrix can be accessed using the member function getEntireMatrix. The submatrix correspoding to two adjacent elements can be obtained using the member function getElementMatrix. Since getElementMatrix is significantly faster than getEntireMatrix, it is advised to use getElementMatrix instead of getEntireMatrix where possible.
    /// There are also several operators defined for addition and multiplication for this class, so for these operations you do not first need to convert it to a standard matrix.
    class FaceMatrix
    {
    public:
        // Constructors
        FaceMatrix(){}
        
        FaceMatrix(const std::size_t nDOFLeft, const std::size_t nDOFRight);
        
        FaceMatrix(const FaceMatrix &other);
        
        // Operators
        /// \brief Defines the operator (iSide, iVarBasisFunc, jSide, jVarBasisFunc) such that a reference to data (iVarBasisFunc, jVarBasisFunc) from the element matrix corresponding to (iSide, jSide) will be returned.
        double & operator() (Side iSide, Side jSide, std::size_t iVarBasisFunction, std::size_t jVarBasisFunction);
        
        /// \brief Defines the operator (i,j) such that a reference to data (i,j) from the face matrix will be returned.
        double & operator() (std::size_t i, std::size_t j);
        
        /// \Sets the face matrix equal to another face matrix.
        FaceMatrix & operator= (const FaceMatrix &other);
        
        /// \biref Adds another face matrix to this face matrix.
        FaceMatrix & operator+= (const FaceMatrix &other);
        
        /// \brief Multiplies the face matrix by a scalar.
        FaceMatrix & operator*= (const double &scalar);
        
        // Other member functions
        /// \brief Gets the number of degrees of freedom (usually the amount of (vector)-basis functions) corresponding to the element at side iSide.
        const std::size_t getNrOfDegreesOfFreedom(Side iSide) const
        {
            if(iSide == Side::LEFT)
                return M_LeftLeft_.getNRows();
            else
                return M_RightRight_.getNRows();
        }
        
        /// \brief Resizes the element matrices.
        void resize(const std::size_t nDOFLeft, const std::size_t nDOFRight);
        
        /// \brief Returns a reference of the submatrix corresponding to a combination of two elements connected to the face.
        const LinearAlgebra::Matrix & getElementMatrix(Side iSide, Side jSide) const;
        
        /// \brief Sets the submatrix corresponding to a combination of two elements connected to the face.
        void setElementMatrix(const LinearAlgebra::Matrix & elementMatrix, Side iSide, Side jSide);
        
        /// \brief Returns the complete face matrix as a standard matrix. It is advised to use getElementMatrix instead when possible.
        const LinearAlgebra::Matrix getEntireMatrix() const;
        
        /// \brief Sets the FaceMatrix using the complete face matrix as input. It is advised to use setElementMatrix instead when possible.
        void setEntireMatrix(const LinearAlgebra::Matrix & entireMatrix);
        
        /// \brief Applies the operation y=ax + y, where a is a scalar and x another FaceMatrix.
        void axpy(const double &a, const FaceMatrix &x);
        
    private:
        /// \brief Four matrices, corresponding to the submatrices of the face matrix, where each submatrix is an element matrix corresponding to (iSide, jSide).
        LinearAlgebra::Matrix M_LeftLeft_;
        LinearAlgebra::Matrix M_LeftRight_;
        LinearAlgebra::Matrix M_RightLeft_;
        LinearAlgebra::Matrix M_RightRight_;
    };
};
#endif

