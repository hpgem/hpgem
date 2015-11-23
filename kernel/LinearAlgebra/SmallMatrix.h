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

#ifndef SMALLMATRIX_H_
#define SMALLMATRIX_H_

#include <array>
#include "SmallVector.h"
#include "MiddleSizeMatrix.h"

namespace LinearAlgebra
{
    /// \class SmallMatrix
    /// \brief Data type for small dense matrix.
    ///
    /// \details Stores small dense matrix efficiently.
    /// It only store doubles as this is the main type linear algebra is done on in hpGEM
    /// It stores the matrix in fortran style (column-major) to give quicker access to extern BLAS libraries.
    /// For example,  the order they are stored in a 2x2 matrix is
    ///     0   2
    ///     1   3
    /// Examples for the implementation are given in the unit test (../tests/unit/LinearAlgebra/MatrixUnitTest).
    template<std::size_t nRows, std::size_t nCols>
    class SmallMatrix
    {
    public:

        /// \brief Constructs a matrix of size n-rows by m-columns.
        SmallMatrix()
            : data_()
        {
        }

        SmallMatrix(const SmallVector<nRows>& other)
            : data_()
        {
            logger.assert(nCols == 1, "Trying to construct a matrix with more than 1 columns from a vector");
            std::copy(other.data(), other.data() + nRows, data_.begin());
        }

        /// \brief Constructs a matrix of size n-rows by m-columns and initialises all entry to a constant
        SmallMatrix(const double& c)
            : data_()
        {
            data_.fill(c);
        }

        /// \brief Construct and copy Matrix from another Matrix i.e. B(A) where B and A are both matrices
        SmallMatrix(const SmallMatrix& other)
            : data_()
        {
            std::copy(other.data_.begin(), other.data_.end(), data_.begin());
        }

        /// \brief Construct and copy Matrix from another Matrix i.e. B(A) where B and A are both matrices
        SmallMatrix(const MiddleSizeMatrix& other)
            : data_()
        {
            logger.assert(other.getNRows() == nRows, "expected a matrix with % rows, but got a matrix with % rows", nRows, other.getNRows());
            logger.assert(other.getNCols() == nCols, "expected a matrix with % columns, but got a matrix with % columns", nCols, other.getNCols());
            for(std::size_t i = 0; i < nRows * nCols; ++i)
            {
                data_[i] = std::real(other[i]);
            }
        }

        /// \brief Glues one or more vectors with the same number of rows together
        SmallMatrix(std::array<SmallVector<nRows>, nCols> entries)
            : data_()
        {
            for(std::size_t i = 0; i < nRows; ++i)
            {
                for(std::size_t j = 0; j < nCols; ++j)
                {
                    (*this)(i, j) = entries[j][i];
                }
            }
        }

        /// \brief Move Matrix from another Matrix
        SmallMatrix(SmallMatrix&& other)
            : data_(std::move(other.data_))
        {
        }

        /// \brief defines the operator (n,m) to access the element on row n and column m
        double& operator()(std::size_t n, std::size_t m)
        {
            logger.assert(n < nRows, "Requested row number % for a matrix with only % rows", n, nRows);
            logger.assert(m < nCols, "Requested column number % for a matrix with only % columns", m, nCols);
            return data_[n + m * nRows];
        }

        /// \brief defines the operator (n,m) to access the element on row n and column m
        const double& operator()(std::size_t n, std::size_t m) const
        {
            logger.assert(n < nRows, "Requested row number % for a matrix with only % rows", n, nRows);
            logger.assert(m < nCols, "Requested column number % for a matrix with only % columns", m, nCols);
            return data_[n + m * nRows];
        }

        /// \brief Access the n linear element in the matrix.
        double& operator[](const std::size_t n)
        {
            logger.assert(n < nRows * nCols, "Requested entry % for a matrix with only % entries", n, nRows * nCols);
            return data_[n];
        }

        const double& operator[](const std::size_t n) const
        {
            logger.assert(n < nRows * nCols, "Requested entry % for a matrix with only % entries", n, nRows * nCols);
            return data_[n];
        }

        /// \brief Defines Matrix A times vector B and return vector C i.e. C_,j= A_ij B_,j
        SmallVector<nRows> operator*(SmallVector<nCols>& right);
        SmallVector<nRows> operator*(SmallVector<nCols>& right) const;

        /// \brief Does matrix A_ij=scalar*B_ij
        SmallMatrix operator*(const double& right) const
        {
            SmallMatrix result;
            std::transform(data_.begin(), data_.end(), result.data_.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, right));
            return result;
        }

        /// \brief Does matrix A_ij = B_ik * C_kj
        template<std::size_t K>
        SmallMatrix<nRows, K> operator*(const SmallMatrix<nCols, K> &other);
        template<std::size_t K>
        SmallMatrix<nRows, K> operator*(const SmallMatrix<nCols, K> &other) const;

        SmallMatrix& operator+=(const SmallMatrix& other)
        {
            std::transform(data_.begin(), data_.end(), other.data_.begin(), data_.begin(), std::plus<double>());
            return *this;
        }

        SmallMatrix& operator-=(const SmallMatrix& other)
        {
            std::transform(data_.begin(), data_.end(), other.data_.begin(), data_.begin(), std::minus<double>());
            return *this;
        }

        SmallMatrix operator+(const SmallMatrix &other) const
        {
            SmallMatrix result;
            std::transform(data_.begin(), data_.end(), other.data_.begin(), result.data_.begin(), std::plus<double>());
            return result;
        }

        SmallMatrix operator-(const SmallMatrix &other) const
        {
            SmallMatrix result;
            std::transform(data_.begin(), data_.end(), other.data_.begin(), result.data_.begin(), std::minus<double>());
            return result;
        }

        SmallMatrix operator-() const
        {
            return *this * -1.;
        }

        /// \brief Does matrix A_ij=scalar*A_ij
        SmallMatrix& operator*=(const double &scalar)
        {
            std::transform(data_.begin(), data_.end(), data_.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, scalar));
            return *this;
        }

        /// \brief Does matrix A_ij = A_ik * B_kj
        /// note that other must be square because this is a fixed-size matrix
        SmallMatrix& operator*=(const SmallMatrix<nCols, nCols> &other);

        /// \brief Does matrix A_ij=scalar*A_ij
        SmallMatrix& operator/=(const double &scalar)
        {
            std::transform(data_.begin(), data_.end(), data_.begin(), std::bind(std::divides<double>(), std::placeholders::_1, scalar));
            return *this;
        }

        /// \brief this does element by divided by a scalar
        SmallMatrix operator/(const double& scalar) const
        {
            SmallMatrix result;
            std::transform(data_.begin(), data_.end(), result.data_.begin(), std::bind(std::divides<double>(), std::placeholders::_1, scalar));
            return result;
        }

        /// \brief Assigns one matrix to another.
        SmallMatrix& operator=(const SmallMatrix& right)
        {
            std::copy(right.data_.begin(), right.data_.end(), data_.begin());
            return *this;
        }

        /// \brief Assigns one matrix to another.
        SmallMatrix& operator=(SmallMatrix&& right)
        {
            std::move(right.data_.begin(), right.data_.end(), data_.begin());
            return *this;
        }

        /// \brief computeWedgeStuffVector.
        SmallVector<nRows> computeWedgeStuffVector() const;

        /// \brief Applies the matrix y=ax + y, where x is another matrix and a is a scalar
        void axpy(double a, const SmallMatrix& x)
        {
            for(std::size_t i = 0; i < nRows * nCols; ++i)
            {
                data_[i] += a * x[i];
            }
        }

        /// \brief Get total number of Matrix entries
        std::size_t size() const
        {
            return nRows * nCols;
        }

        /// \brief Get the number of rows
        std::size_t getNRows() const
        {
            return nRows;
        }

        /// \brief Get the number of columns
        std::size_t getNCols() const
        {
            return nCols;
        }

        /// \brief get the j^th column
        SmallVector<nRows> getColumn(std::size_t j) const
        {
            logger.assert(j < nCols, "Asked for column %, but there are only % columns", j, nCols);
            return SmallVector<nRows>(data() + j * nRows);
        }

        /// \brief get the i^th row
        SmallVector<nCols> getRow(std::size_t i) const
        {
            logger.assert(i < nRows, "Asked for row %, but there are only % rows", i, nRows);
            SmallVector<nCols> result;
            for(std::size_t j = 0; j < nCols; ++j)
            {
                result[j] = (*this)(i,j);
            }
            return result;
        }

        /// \brief Return the LUfactorisation of the matrix
        SmallMatrix LUfactorisation() const;

        double determinant() const;

        /// \brief return the inverse in the vector result. The size of result matches the matrix.
        SmallMatrix inverse() const;

        SmallMatrix<nCols, nRows> transpose() const
        {
            SmallMatrix<nCols, nRows> result;
            for(std::size_t i = 0; i < nRows; ++i)
            {
                for(std::size_t j = 0; j < nCols; ++j)
                {
                    result(j, i) = (*this)(i, j);
                }
            }
            return result;
        }

        /// \brief solves Ax=B where A is the current matrix and B is passed in. The result is returned in B.
        template<std::size_t nRHS>
        void solve(SmallMatrix<nRows, nRHS>& B) const;

        /// \brief solves Ax=b where A is the current matrix and NumericalVector b
        /// is the input parameter. The result is returned in b.
        void solve(SmallVector<nRows>& b) const;

        double* data()
        {
            return data_.data();
        }

        const double* data() const
        {
            return data_.data();
        }

    private:
        /// The actually data of the matrix class
        std::array<double, nRows * nCols> data_;
    };

    /// Writes nicely formatted entries of the Matrix A to the stream os.
    template<std::size_t nRows, std::size_t nCols>
    std::ostream& operator<<(std::ostream& os, const SmallMatrix<nRows, nCols>& A)
    {
        for(std::size_t i = 0; i < nRows; ++i)
        {
            os << A.getRow(i) << std::endl;
        }
        return os;
    }

    ///Multiplies a matrix with a double
    template<std::size_t nRows, std::size_t nCols>
    SmallMatrix<nRows, nCols> operator*(const double d, const SmallMatrix<nRows, nCols>& mat)
    {
        return mat * d;
    }

    ///Multiplies a matrix with a vector
    template<std::size_t nRows, std::size_t nCols>
    SmallVector<nCols> operator*(SmallVector<nRows>& vec, SmallMatrix<nRows, nCols>& mat);

    template<std::size_t nRows, std::size_t nCols>
    MiddleSizeMatrix::MiddleSizeMatrix(const SmallMatrix<nRows, nCols>& other)
        : data_(nRows * nCols), nRows_(nRows), nCols_(nCols)
    {
        std::copy(other.data(), other.data() + nRows * nCols, data_.begin());
    }

}

#include "SmallMatrix_impl.h"


#endif /* SMALLMATRIX_H_ */
