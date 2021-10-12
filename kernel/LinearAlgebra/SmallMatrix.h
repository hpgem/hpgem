/*
 This file forms part of hpGEM. This package has been developed over a number of
 years by various people at the University of Twente and a full list of
 contributors can be found at http://hpgem.org/about-the-code/team

 This code is distributed using BSD 3-Clause License. A copy of which can found
 below.


 Copyright (c) 2014, University of Twente
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.

 3. Neither the name of the copyright holder nor the names of its contributors
 may be used to endorse or promote products derived from this software without
 specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef HPGEM_KERNEL_SMALLMATRIX_H
#define HPGEM_KERNEL_SMALLMATRIX_H

#include <array>
#include "SmallVector.h"
#include "MiddleSizeMatrix.h"

#include "Eigen/Eigen"

namespace hpgem {
namespace LinearAlgebra {
/// \class SmallMatrix
/// \brief Data type for small dense matrix.
///
/// \details Stores small dense matrix efficiently.
/// It only store doubles as this is the main type linear algebra is done on in
/// hpGEM It stores the matrix in fortran style (column-major) to give quicker
/// access to extern BLAS libraries. For example,  the order they are stored in
/// a 2x2 matrix is
///     0   2
///     1   3
/// Examples for the implementation are given in the unit test
/// (../tests/unit/LinearAlgebra/MatrixUnitTest).
template <std::size_t numberOfRows, std::size_t numberOfColumns>
class SmallMatrix {
   public:
    using EigenType = Eigen::Matrix<double, numberOfRows, numberOfColumns>;

    /// \brief Constructs a matrix of size n-rows by m-columns.
    SmallMatrix() : data2_(EigenType::Zero()) {}

    explicit SmallMatrix(const SmallVector<numberOfRows>& other) : data2_() {
        logger.assert_debug(numberOfColumns == 1,
                            "Trying to construct a matrix with more than 1 "
                            "columns from a vector");
        std::copy(other.data(), other.data() + numberOfRows, data2_.begin());
    }

    /// \brief Constructs a matrix of size n-rows by m-columns and initialises
    /// all entry to a constant
    explicit SmallMatrix(const double& c) : data2_(EigenType::Constant(c)) {}

    SmallMatrix(const std::initializer_list<SmallVector<numberOfRows>>& entries)
        : data2_() {
        logger.assert_debug(
            entries.size() == numberOfColumns,
            "expected a matrix with % columns, but got a matrix with % columns",
            numberOfColumns, entries.size());
        std::size_t column = 0;
        for (const SmallVector<numberOfRows>& entrie : entries) {
            for (std::size_t i = 0; i < numberOfRows; ++i) {
                (*this)(i, column) = entrie[i];
            }
            ++column;
        }
    }

    /// \brief Construct and copy Matrix from another Matrix i.e. B(A) where B
    /// and A are both matrices
    SmallMatrix(const SmallMatrix& other) : data2_(other.data2_) {
    }

    /// \brief Construct and copy Matrix from another Matrix i.e. B(A) where B
    /// and A are both matrices
    SmallMatrix(const MiddleSizeMatrix& other) : data2_() {
        logger.assert_debug(
            other.getNumberOfRows() == numberOfRows,
            "expected a matrix with % rows, but got a matrix with % rows",
            numberOfRows, other.getNumberOfRows());
        logger.assert_debug(
            other.getNumberOfColumns() == numberOfColumns,
            "expected a matrix with % columns, but got a matrix with % columns",
            numberOfColumns, other.getNumberOfColumns());
        for (std::size_t i = 0; i < numberOfRows * numberOfColumns; ++i) {
            data2_(i) = std::real(other[i]);
        }
    }

    /// \brief Glues one or more vectors with the same number of rows together
    SmallMatrix(std::array<SmallVector<numberOfRows>, numberOfColumns> entries)
        : data2_() {
        for (std::size_t i = 0; i < numberOfRows; ++i) {
            for (std::size_t j = 0; j < numberOfColumns; ++j) {
                (*this)(i, j) = entries[j][i];
            }
        }
    }

    /// \brief Move Matrix from another Matrix
    SmallMatrix(SmallMatrix&& other) : data2_(std::move(other.data2_)) {}

    SmallMatrix(EigenType rawData) : data2_(rawData) {}

    /// \brief defines the operator (n,m) to access the element on row n and
    /// column m
    double& operator()(std::size_t n, std::size_t m) {
        return data2_(n,m);
    }

    /// \brief defines the operator (n,m) to access the element on row n and
    /// column m
    const double& operator()(std::size_t n, std::size_t m) const {
        return data2_(n,m);
    }

    /// \brief Access the n linear element in the matrix.
    double& operator[](const std::size_t n) {
        return data2_(n);
    }

    const double& operator[](const std::size_t n) const {
        return data2_(n);
    }

    /// \brief Defines Matrix A times vector B and return vector C i.e. C_,j=
    /// A_ij B_,j
    SmallVector<numberOfRows> operator*(SmallVector<numberOfColumns>& right);
    SmallVector<numberOfRows> operator*(
        SmallVector<numberOfColumns>& right) const;

    /// \brief Does matrix A_ij=scalar*B_ij
    SmallMatrix operator*(const double& right) const {
        EigenType result = data2_ * right;
        return result;
    }

    /// \brief Does matrix A_ij = B_ik * C_kj
    template <std::size_t K>
    SmallMatrix<numberOfRows, K> operator*(
        const SmallMatrix<numberOfColumns, K>& other);
    template <std::size_t K>
    SmallMatrix<numberOfRows, K> operator*(
        const SmallMatrix<numberOfColumns, K>& other) const;

    SmallMatrix& operator+=(const SmallMatrix& other) {
        data2_ += other.data2_;
        return *this;
    }

    SmallMatrix& operator-=(const SmallMatrix& other) {
        data2_ -= other.data2_;
        return *this;
    }

    SmallMatrix operator+(const SmallMatrix& other) const {
        EigenType res = data2_ + other.data2_;
        return res;
    }

    SmallMatrix operator-(const SmallMatrix& other) const {
        EigenType res = data2_ - other.data2_;
        return res;
    }

    SmallMatrix operator-() const { return *this * -1.; }

    /// \brief Does matrix A_ij=scalar*A_ij
    SmallMatrix& operator*=(const double& scalar) {
        data2_ *= scalar;
        return *this;
    }

    /// \brief Does matrix A_ij = A_ik * B_kj
    /// note that other must be square because this is a fixed-size matrix
    SmallMatrix& operator*=(
        const SmallMatrix<numberOfColumns, numberOfColumns>& other);

    /// \brief Does matrix A_ij=scalar*A_ij
    SmallMatrix& operator/=(const double& scalar) {
        data2_ /= scalar;
        return *this;
    }

    /// \brief this does element by divided by a scalar
    SmallMatrix operator/(const double& scalar) const {
        EigenType res = data2_ / scalar;
        return res;
    }

    /// \brief Assigns one matrix to another.
    SmallMatrix& operator=(const SmallMatrix& right) {
        data2_ = right.data2_;
        return *this;
    }

    /// \brief Assigns one matrix to another.
    SmallMatrix& operator=(SmallMatrix&& right) {
        data2_ = right.data2_;
        return *this;
    }

    /// \brief computeWedgeStuffVector.
    SmallVector<numberOfRows> computeWedgeStuffVector() const;

    /// \brief Applies the matrix y=ax + y, where x is another matrix and a is a
    /// scalar
    void axpy(double a, const SmallMatrix& x) {
        data2_ += a * x.data2_;
    }

    /// \brief Get total number of Matrix entries
    std::size_t size() const { return numberOfRows * numberOfColumns; }

    /// \brief Get the number of rows
    std::size_t getNumberOfRows() const { return numberOfRows; }

    ///\deprecated Does not conform naming convention, please use
    /// getNumberOfRows instead.
    std::size_t getNRows() const { return getNumberOfRows(); }

    /// \brief Get the number of columns
    std::size_t getNumberOfColumns() const { return numberOfColumns; }

    std::size_t getNCols() const { return getNumberOfColumns(); }

    /// \brief get the j^th column
    SmallVector<numberOfRows> getColumn(std::size_t j) const {
        typename SmallVector<numberOfRows>::EigenType res = data2_.col(j);
        return res;
    }

    /// \brief get the i^th row
    SmallVector<numberOfColumns> getRow(std::size_t i) const {
        typename SmallVector<numberOfColumns>::EigenType res = data2_.row(i);
        return res;
    }

    double determinant() const;

    /// \brief return the inverse in the vector result. The size of result
    /// matches the matrix.
    SmallMatrix inverse() const;

    SmallMatrix<numberOfColumns, numberOfRows> transpose() const {
        Eigen::Matrix<double, numberOfColumns, numberOfRows> res;
        res = data2_.transpose();
        return res;
    }

    /// \brief solves Ax=B where A is the current matrix and B is passed in. The
    /// result is returned in B.
    template <std::size_t numberOfRightHandSideColumns>
    void solve(
        SmallMatrix<numberOfRows, numberOfRightHandSideColumns>& B) const;

    /// \brief solves Ax=b where A is the current matrix and NumericalVector b
    /// is the input parameter. The result is returned in b.
    void solve(SmallVector<numberOfRows>& b) const;

    double* data() { return data2_.data(); }

    const double* data() const { return data2_.data(); }

    operator EigenType& () {
        return data2_;
    }

    operator const EigenType& () const {
        return data2_;
    }

   private:

    /// The actually data of the matrix class
    EigenType data2_;
};

/// Writes nicely formatted entries of the Matrix A to the stream os.
template <std::size_t numberOfRows, std::size_t numberOfColumns>
std::ostream& operator<<(std::ostream& os,
                         const SmallMatrix<numberOfRows, numberOfColumns>& A) {
    for (std::size_t i = 0; i < numberOfRows; ++i) {
        os << A.getRow(i) << std::endl;
    }
    return os;
}

/// Multiplies a matrix with a double
template <std::size_t numberOfRows, std::size_t numberOfColumns>
SmallMatrix<numberOfRows, numberOfColumns> operator*(
    const double d, const SmallMatrix<numberOfRows, numberOfColumns>& mat) {
    return mat * d;
}

/// Multiplies a matrix with a vector
template <std::size_t numberOfRows, std::size_t numberOfColumns>
SmallVector<numberOfColumns> operator*(
    SmallVector<numberOfRows>& vec,
    SmallMatrix<numberOfRows, numberOfColumns>& mat);

template <std::size_t numberOfRows, std::size_t numberOfColumns>
MiddleSizeMatrix::MiddleSizeMatrix(
    const SmallMatrix<numberOfRows, numberOfColumns>& other)
    : data_(numberOfRows * numberOfColumns),
      numberOfRows_(numberOfRows),
      numberOfColumns_(numberOfColumns) {
    std::copy(other.data(), other.data() + numberOfRows * numberOfColumns,
              data_.begin());
}

}  // namespace LinearAlgebra
}  // namespace hpgem
#include "SmallMatrix_impl.h"

#endif  // HPGEM_KERNEL_SMALLMATRIX_H
