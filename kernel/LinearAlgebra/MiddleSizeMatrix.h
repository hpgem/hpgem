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
#ifndef MIDDLESIZEMATRIX_HH
#define MIDDLESIZEMATRIX_HH

// System includes
#include <iostream>
#include <vector>

#include "MiddleSizeVector.h"
#include <complex>

namespace LinearAlgebra
{
    template<std::size_t numberOfRows, std::size_t numberOfColumns>
    class SmallMatrix;
    class MiddleSizeVector;
    /// \class MiddleSizeMatrix
    /// \brief Data type for dense matrix.
    /// 
    /// \details Stores dense matrix efficiently.
    /// It is templated to store complex numbers when PETSc needs those and stores doubles otherwise
    /// It stores the matrix in fortran style (column-major) to give quicker access to extern BLAS libraries.
    /// For example,  the order they are stored in a 2x2 matrix is 
    ///     0   2 
    ///     1   3
    /// Examples for the implementation are given in the unit test (../tests/unit/LinearAlgebra/MatrixUnitTest).
    /// implement fall-back routines for matrices that are too large. Note that in this situation dense data storage may not be the best option
    class MiddleSizeMatrix
    {
    public:
#ifdef HPGEM_USE_COMPLEX_PETSC
        using type = std::complex<double>;
#else
        using type = double;
#endif
        
        /// \brief Default Matrix constructor : Simply creates a zero size matrix
        MiddleSizeMatrix();

        /// \brief Constructs a matrix of size n-rows by m-columns.
        MiddleSizeMatrix(const std::size_t n, const std::size_t m);

        /// \brief Constructs a matrix of size n-rows by m-columns and initialises all entry to a constant 
        MiddleSizeMatrix(const std::size_t n, const std::size_t m, const type& c);

        /// \brief Construct and copy Matrix from another Matrix i.e. B(A) where B and A are both matrices
        MiddleSizeMatrix(const MiddleSizeMatrix& other);
        
        //implemented with SmallMatrix for dependency reasons
        template<std::size_t numberOfRows, std::size_t nCols>
        MiddleSizeMatrix(const SmallMatrix<numberOfRows, nCols>& other);

        /// \brief construct a matrix by placing some vectors next to each other. Note that vectors in hpGEM are column vectors
        MiddleSizeMatrix(const MiddleSizeVector& other);

        /// \brief Glues one or more matrices with the same number of rows together
        MiddleSizeMatrix(std::initializer_list<MiddleSizeMatrix>);

        /// \brief Move Matrix from another Matrix
        MiddleSizeMatrix(MiddleSizeMatrix&& other);

        /// \brief defines the operator (n,m) to access the element on row n and column m        
        type& operator()(std::size_t n, std::size_t m)
        {
            logger.assert_debug(n < numberOfRows_, "Requested row number % for a matrix with only % rows", n, numberOfRows_);
            logger.assert_debug(m < numberOfColumns_, "Requested column number % for a matrix with only % columns", m, numberOfColumns_);
            return data_[n + m * numberOfRows_];
        }
        
        /// \brief defines the operator (n,m) to access the element on row n and column m        
        const type& operator()(std::size_t n, std::size_t m) const
        {
            logger.assert_debug(n < numberOfRows_, "Requested row number % for a matrix with only % rows", n, numberOfRows_);
            logger.assert_debug(m < numberOfColumns_, "Requested column number % for a matrix with only % columns", m, numberOfColumns_);
            return data_[n + m * numberOfRows_];
        }
        
        /// \brief Access the n linear element in the matrix. 
        type& operator[](const std::size_t n);

        const type& operator[](const std::size_t n) const;

        /// \brief Defines Matrix A times vector B and return vector C i.e. C_,j= A_ij B_,j
        MiddleSizeVector operator*(MiddleSizeVector& right);
        MiddleSizeVector operator*(MiddleSizeVector& right) const;

        /// \brief Does matrix A_ij=scalar*B_ij
        MiddleSizeMatrix operator*(const type& right) const;

        /// \brief Does matrix A_ij = B_ik * C_kj
        MiddleSizeMatrix operator*(const MiddleSizeMatrix &other);
        MiddleSizeMatrix operator*(const MiddleSizeMatrix &other) const;

        MiddleSizeMatrix& operator+=(const MiddleSizeMatrix& other);
        MiddleSizeMatrix& operator-=(const MiddleSizeMatrix& other);

        MiddleSizeMatrix operator+(const MiddleSizeMatrix &other) const;
        MiddleSizeMatrix operator-(const MiddleSizeMatrix &other) const;
        MiddleSizeMatrix operator-() const;

        /// \brief Does matrix A_ij=scalar*A_ij
        MiddleSizeMatrix& operator*=(const type &scalar);

        /// \brief Does matrix A_ij = A_ik * B_kj
        MiddleSizeMatrix& operator*=(const MiddleSizeMatrix &other);

        /// \brief Does matrix A_ij=scalar*A_ij
        MiddleSizeMatrix& operator/=(const type &scalar);

        /// \brief this does element by divided by a scalar
        MiddleSizeMatrix operator/(const type& scalar) const;

        /// \brief Assigns the Matrix by a scalar
        MiddleSizeMatrix& operator=(const type& c);

        /// \brief Assigns one matrix to another.
        MiddleSizeMatrix& operator=(const MiddleSizeMatrix& right);

        /// \brief Assigns one matrix to another.
        MiddleSizeMatrix& operator=(MiddleSizeMatrix&& right);

        /// \brief computeWedgeStuffVector.
        MiddleSizeVector computeWedgeStuffVector() const;

        /// \brief Applies the matrix y=ax + y, where x is another matrix and a is a scalar
        void axpy(type a, const MiddleSizeMatrix& x);

        /// \brief Resize the Matrix to be n-Rows by m-columns
        void resize(std::size_t n, std::size_t m);

        /// \brief Glues two matrices with the same number of columns together
        void concatenate(const MiddleSizeMatrix& other);

        /// \brief Get total number of Matrix entries
        std::size_t size() const;

        /// \brief Get the number of rows
        std::size_t getNumberOfRows() const;

        ///\deprecated Does not conform the naming conventions, please use getNumberOfRows instead.
        std::size_t getNRows() const
        {
          return getNumberOfRows();
        }

        /// \brief Get the number of columns
        std::size_t getNumberOfColumns() const;

        ///\deprecated Does not conform the naming conventions, please use getNumberOfColumns instead.
        std::size_t getNCols() const
        {
          return getNumberOfColumns();
        }

        /// \brief get the j^th column
        /// If someone knows how to do this such that it returns a reference, please
        ///implement it.
        LinearAlgebra::MiddleSizeVector getColumn(std::size_t j) const;

        /// \brief get the i^th row
        /// If someone knows how to do this such that it returns a reference, please
        ///implement it.
        //Not a good idea: Returning a reference means that the data that gets returned persists internally, but the data is stored in a column major ordering, so we should not keep row vectors -FB
        LinearAlgebra::MiddleSizeVector getRow(std::size_t i) const;

        /// \brief Return the LUfactorisation of the matrix
        MiddleSizeMatrix LUfactorisation() const;

        /// \brief return the inverse in the vector result. The size of result matches the matrix.
        MiddleSizeMatrix inverse() const;

        MiddleSizeMatrix transpose() const;

        /// \brief solves Ax=B where A is the current matrix and B is passed in. The result is returned in B.
        void solve(MiddleSizeMatrix& B) const;

        /// \brief solves Ax=b where A is the current matrix and NumericalVector b 
        /// is the input parameter. The result is returned in b.
        void solve(MiddleSizeVector& b) const;

        /// \brief Computes the Cholesky decomposition in place.
        ///
        /// Computes the Cholesky decomposition L*L^H = A of this matrix (A),
        /// assuming but not checking that this matrix is Hermitian positive
        /// definite. The contents of this matrix is overwritten with the
        /// lower triangular part of the Cholesky decomposition. The strictly
        /// upper triangular part is left untouched.
        void cholesky();

        /// \brief Solve a system (L*L^H)*X = B after Cholesky decomposition
        ///
        /// Solves the system (L^H *L)*X = B, where L is a lower triangular
        /// matrix such as obtained from cholesky(). This is not checked.
        /// \param B The rhs and result matrix.
        void solveCholesky(MiddleSizeMatrix &B) const;

        /// \brief Solve the system L*X = B or X*L^H = B
        ///
        /// Solve the system L*X = B or X*L^H = B, where L is the current matrix
        /// that is assumed to be lower triangular. The strictly upper triagonal
        /// part of the matrix is not referenced.
        ///
        /// \param B The right hand side, that will be overwritten by the
        /// solution.
        /// \param left whether to solve L*X = B (true) or X*L^H = B (false)
        void solveLowerTriangular(MiddleSizeMatrix &B, bool left = true) const;

        /// \brief Computes the minimum norm solution to a real linear least squares problem.
        void pseudoSolve(MiddleSizeVector& b, double rCond) const;
        
        type* data();
        const type* data() const;

        /// \brief Checks if the matrix is Hermitian up to a tolerance
        ///
        /// \param tol Difference of absolute magnitude larger than this is seen
        ///   as non hermitian.
        /// \return Whether the matrix is hermitian up to the specified tolerance
        bool isHermitian(double tol) const;

    private:
        /// The actually data of the matrix class
        std::vector<type> data_;

        
        /// Stores the number of rows of the matrix
        std::size_t numberOfRows_;

        /// Store the number of columns of the matrix
        std::size_t numberOfColumns_;
        
    };
    
    /// Writes nicely formatted entries of the Matrix A to the stream os.
    std::ostream& operator<<(std::ostream& os, const MiddleSizeMatrix& A);

    ///Multiplies a matrix with a double
    MiddleSizeMatrix operator*(const MiddleSizeMatrix::type d, const MiddleSizeMatrix& mat);
    
    ///Multiplies a matrix with a double
    MiddleSizeVector operator*(MiddleSizeVector& vec, MiddleSizeMatrix& mat);

}
#endif
