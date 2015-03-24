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
//------------------------------------------------------------------------------
#ifndef MATRIX_HH
#define MATRIX_HH

// System includes
#include <iostream>
#ifdef LA_STL_VECTOR
#include <vector>
#else
#include <valarray>
#endif

#include "NumericalVector.h"
#include <complex>

namespace LinearAlgebra
{
    class NumericalVector;
    //We need the ostream for outputting and we encapsulate from valarray.
#ifdef LA_STL_VECTOR
    using std::vector;
#else
    using std::valarray;
#endif
    /// \class Matrix
    /// \brief Data type for small dense matrix.
    /// 
    /// \details Stores small dense matrix efficiently.
    /// Since this class is inherited from std::valarray, it also inherits members of std::valarray
    /// Note, valarray was speed tested and was shown to be quicker than stl vector and it very light weight
    /// It only store doubles as this is the main type linear algebra is done on in hpGEM
    /// It stores the matrix in fortran style (column-major) to give quicker access to extern BLAS libraries.
    /// For example,  the order they are stored in a 2x2 matrix is 
    ///     0   2 
    ///     1   3
    /// Examples for the implementation are given in the unit test (../tests/unit/LinearAlgebra/MatrixUnitTest).
    /// \bug valarray does not compile, since data() is not defined on it
    /// \todo number of rows/columns is saved as a std::size_t, but to BLAS is limited to 32 bits (not sure if signed or unsigned)
    /// implement fall-back routines for matrices that are too large. Note that in this situation dense data storage may not be the best option
    /// \todo Complete the set of operators. Missing: operator-(), operator-(const Matrix&)
    class Matrix
    {
    public:
        
        /// \brief Default Matrix constructor : Simply creates a zero size matrix
        Matrix();

        /// \brief Constructs a matrix of size n-rows by m-columns.
        Matrix(const std::size_t n, const std::size_t m);

        /// \brief Constructs a matrix of size n-rows by m-columns and initialises all entry to a constant 
        Matrix(const std::size_t n, const std::size_t m, const double& c);

        /// \brief Construct and copy Matrix from another Matrix i.e. B(A) where B and A are both matrices
        Matrix(const Matrix& other);
        
        /// \brief Move Matrix from another Matrix
        Matrix(Matrix&& other);

        /// \brief defines the operator (n,m) to access the element on row n and column m        
        double& operator()(std::size_t n, std::size_t m)
        {
            logger.assert(n < nRows_, "Requested row number % for a matrix with only % rows", n, nRows_);
            logger.assert(m < nCols_, "Requested column number % for a matrix with only % columns", m, nCols_);
            return data_[n + m * nRows_];
        }
        
        /// \brief defines the operator (n,m) to access the element on row n and column m        
        const double& operator()(std::size_t n, std::size_t m) const
        {
            logger.assert(n < nRows_, "Requested row number % for a matrix with only % rows", n, nRows_);
            logger.assert(m < nCols_, "Requested column number % for a matrix with only % columns", m, nCols_);
            return data_[n + m * nRows_];
        }
        
        /// \brief Access the n linear element in the matrix. 
        double& operator[](const std::size_t n);

        const double& operator[](const std::size_t n) const;

        /// \brief Defines Matrix A times vector B and return vector C i.e. C_,j= A_ij B_,j
        NumericalVector operator*(NumericalVector& right) const;

        /// \brief Does matrix A_ij = B_ik * C_kj
        Matrix operator*(const Matrix &other);
        Matrix operator*(const Matrix &other) const;

        Matrix& operator+=(const Matrix& other);

        /// \brief Does matrix A_ij=scalar*A_ij
        Matrix& operator*=(const double &scalar);

        /// \brief this does element by divided by a scalar 
        Matrix& operator/=(const double& scalar);

        /// \brief this does element by divided by a scalar
        Matrix operator/(const double& scalar);

        /// \brief Assigns the Matrix by a scalar
        Matrix& operator=(const double& c);

        /// \brief Assigns one matrix to another.
        Matrix& operator=(const Matrix& right);

        /// \brief Assigns one matrix to another.
        Matrix& operator=(Matrix&& right);

        /// \brief computeWedgeStuffVector.
        NumericalVector computeWedgeStuffVector() const;

        /// \brief Applies the matrix y=ax + y, where x is another matrix and a is a scalar
        void axpy(double a, const Matrix& x);

        /// \brief Resize the Matrix to be n-Rows by m-columns
        void resize(std::size_t n, std::size_t m);

        /// \brief Glues two matrices with the same number of columns together
        void concatenate(const Matrix& other);

        /// \brief Get total number of Matrix entries
        const std::size_t size() const;

        /// \brief Get the number of rows
        const std::size_t getNRows() const;

        /// \brief Get the number of columns
        const std::size_t getNCols() const;

        /// \brief get the j^th column
        /// If someone knows how to do this such that it returns a reference, please
        ///implement it.
        LinearAlgebra::NumericalVector getColumn(std::size_t j) const;

        /// \brief get the i^th row
        /// If someone knows how to do this such that it returns a reference, please
        ///implement it.
        LinearAlgebra::NumericalVector getRow(std::size_t i) const;

        /// \brief Return the LUfactorisation of the matrix
        Matrix LUfactorisation() const;

        /// \brief return the inverse in the vector result. The size of result matches the matrix.
        Matrix inverse() const;

        /// \brief solves Ax=B where A is the current matrix and B is passed in. The result is returned in B.
        void solve(Matrix& B) const;

        /// \brief solves Ax=b where A is the current matrix and NumericalVector b 
        /// is the input parameter. The result is returned in b.
        void solve(NumericalVector& b) const;

#ifdef HPGEM_USE_COMPLEX_PETSC
        std::complex<double>* data();
        const std::complex<double>* data() const;

#else
        double* data();
        const double* data() const;
#endif
        
    private:
        /// The actually data of the matrix class
#ifdef LA_STL_VECTOR
        vector<double> data_;
#else
#error "valarray is broken, please turn on hpGEM_USE_STL_VECTOR_FOR_LA"
        valarray<double> data_;
#endif
        
        /// Stores the number of rows of the matrix
        std::size_t nRows_;

        /// Store the number of columns of the matrix
        std::size_t nCols_;
        
    };
    
    /// Writes nicely formatted entries of the Matrix A to the stream os.
    std::ostream& operator<<(std::ostream& os, const Matrix& A);
    
    ///Adds two matrices
    Matrix operator+(const Matrix& mat1, const Matrix& mat2);
    
    ///Multiplies a matrix with a double
    Matrix operator*(const double d, const Matrix& mat);

}
#endif
