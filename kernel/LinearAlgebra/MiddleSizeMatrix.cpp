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

#include "MiddleSizeMatrix.h"
#include "MiddleSizeVector.h"
#include "Logger.h"
#include <algorithm>
#include <complex>

namespace LinearAlgebra
{
    
    extern "C"
    {

        /// This does matrix times vector and is from blas level 2.
        void dgemv_(const char* trans, int* m, int* n, double* alpha, double* A, int* LDA, double *x, int *incx, double *beta, double *y, int *incy);
        /// This does matrix times vector and is from blas level 2.
        void zgemv_(const char* trans, int* m, int* n, std::complex<double>* alpha, std::complex<double>* A, int* LDA, std::complex<double> *x, int *incx, std::complex<double> *beta, std::complex<double> *y, int *incy);

        ///This is the gernal matrix multiplication from blas level 3
        int dgemm_(const char *transA, const char *transB, int *M, int *N, int *k, double *alpha, double *A, int *LDA, double *B, int *LDB, double *beta, double *C, int *LDC);
        ///This is the gernal matrix multiplication from blas level 3
        int zgemm_(const char *transA, const char *transB, int *M, int *N, int *k, std::complex<double> *alpha, std::complex<double> *A, int *LDA, std::complex<double> *B, int *LDB, std::complex<double> *beta, std::complex<double> *C, int *LDC);

        ///This is the gerneral scalar times vector + vector from blas, hence from blas level 1. Here we also use on a matrix by treating as a vector
        int daxpy_(int* N, double* DA, double* DX, int* INCX, double* DY, int* INCY);
        ///This is the gerneral scalar times vector + vector from blas, hence from blas level 1. Here we also use on a matrix by treating as a vector
        int zaxpy_(int* N, std::complex<double>* DA, std::complex<double>* DX, int* INCX, std::complex<double>* DY, int* INCY);

        /// This is LU factorisation of the matrix A. This has been taken from LAPACK 
        void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
        /// This is LU factorisation of the matrix A. This has been taken from LAPACK
        void zgetrf_(int* M, int *N, std::complex<double>* A, int* lda, int* IPIV, int* INFO);

        /// This is the inverse calulation also from LAPACK. Calculates inverse if you pass it the LU factorisation.
        void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
        /// This is the inverse calulation also from LAPACK. Calculates inverse if you pass it the LU factorisation.
        void zgetri_(int* N, std::complex<double>* A, int* lda, int* IPIV, std::complex<double>* WORK, int* lwork, int* INFO);

        /// This is used for solve Ax=B for x. Again this is from LAPACK.
        void dgesv_(int* N, int* NRHS, double* A, int* lda, int* IPIV, double* B, int* LDB, int* INFO);
        /// This is used for solve Ax=B for x. Again this is from LAPACK.
        void zgesv_(int* N, int* NRHS, std::complex<double>* A, int* lda, int* IPIV, std::complex<double>* B, int* LDB, int* INFO);
    }
    
    MiddleSizeMatrix::MiddleSizeMatrix()
            : nRows_(0), nCols_(0)
    {
    }
    
    /// \param[in]  n The number of rows the matrix will have
    /// \param[in]  m The number of columns the matrix will have
    MiddleSizeMatrix::MiddleSizeMatrix(const std::size_t n, const std::size_t m)
            : data_(n * m), nRows_(n), nCols_(m)
    {
        logger.assert(n <= std::numeric_limits<int>::max() && m <= std::numeric_limits<int>::max(), "Dense linear algebra is not supported on this system for matrices that are this large");
    }
    
    /// \param[in]  n The number of rows the matrix will have
    /// \param[in]  m The number of columns the matrix will have
    /// \param[in]  c The value all the entries of the matrix will contain
    ///
    /// \details Example usage : Matrix<double> A(4,3,2) with create a 4 by 3 
    /// matrix called A with entries equal to 2.
    MiddleSizeMatrix::MiddleSizeMatrix(const std::size_t n, const std::size_t m, const type& c)
            :
#ifdef LA_STL_VECTOR
                    data_(n * m, c),
#else
                    data_(c, n*m),
#endif
                    nRows_(n), nCols_(m)
    {
        logger.assert(n <= std::numeric_limits<int>::max() && m <= std::numeric_limits<int>::max(), "Dense linear algebra is not supported on this system for matrices that are this large");
    }
    
    /// \param[in] Matrix A i.e. the matrix to be copies.
    MiddleSizeMatrix::MiddleSizeMatrix(const MiddleSizeMatrix& other)
            : data_(other.data_), nRows_(other.nRows_), nCols_(other.nCols_)
    {
    }
    
    MiddleSizeMatrix::MiddleSizeMatrix(MiddleSizeMatrix&& other)
            : data_(std::move(other.data_)), nRows_(other.nRows_), nCols_(other.nCols_)
    {
    }
    
    MiddleSizeMatrix::MiddleSizeMatrix(const MiddleSizeVector& list)
            : data_(list.data(), list.data() + list.size()), nRows_(list.size()), nCols_(1)
    {
        logger.assert(list.size() <= std::numeric_limits<int>::max(), "Dense linear algebra is not supported on this system for matrices that are this large");
    }

    MiddleSizeMatrix::MiddleSizeMatrix(std::initializer_list<MiddleSizeMatrix> list)
            : data_(0), nRows_(list.begin()->getNRows()), nCols_(0)
    {
        for(MiddleSizeMatrix mat : list)
        {
            logger.assert(nRows_ == mat.getNRows(), "Can only construct a matrix from vectors of the same size");
            nCols_ += mat.getNCols();
        }
        logger.assert(nRows_ <= std::numeric_limits<int>::max() && nCols_ <= std::numeric_limits<int>::max(), "Dense linear algebra is not supported on this system for matrices that are this large");
        data_.resize(nRows_ * nCols_);
        auto inserter = data_.begin();
        for(MiddleSizeMatrix mat : list)
        {
            inserter = std::copy(mat.data(), mat.data() + mat.size(), inserter);
        }
    }

    /// \param[in] n The number of the row you want the element from
    /// \return double i.e. the value of the element you requested
    ///
    /// \details Recall that the matrix is stored in fortran style i.e. columns first and then rows
    MiddleSizeMatrix::type& MiddleSizeMatrix::operator[](const std::size_t n)
    {
        logger.assert(n < data_.size(), "Requested entry % for a matrix with only % entries", n, data_.size());
        return data_[n];
    }
    
    const MiddleSizeMatrix::type& MiddleSizeMatrix::operator[](const std::size_t n) const
    {
        logger.assert(n < data_.size(), "Requested entry % for a matrix with only % entries", n, data_.size());
        return data_[n];
    }
    
    /// \param[in] other : the Matrix that is added to this matrix
    /// \return Matrix
    MiddleSizeMatrix& MiddleSizeMatrix::operator+=(const MiddleSizeMatrix& other)
    {
        //Make sure the matrices are the same size
        logger.assert(size() == other.size() && nCols_ == other.nCols_, "Dimensions of matrices are not the same.");
        
        //add the matrices element-wise
        for (std::size_t i = 0; i < size(); ++i)
        {
            data_[i] += other[i];
        }
        logger(DEBUG, "Called for matrix addition");
        return (*this);
    }
    
    MiddleSizeMatrix& MiddleSizeMatrix::operator-=(const MiddleSizeMatrix& other)
    {
        //Make sure the matrices are the same size
        logger.assert(size() == other.size() && nCols_ == other.nCols_, "Dimensions of matrices are not the same.");

        //add the matrices element-wise
        for (std::size_t i = 0; i < size(); ++i)
        {
            data_[i] -= other[i];
        }
        return (*this);
    }

    MiddleSizeMatrix& MiddleSizeMatrix::operator *=(const MiddleSizeMatrix& other)
    {
        ///\todo BLAS has no in-place matrix-matrix multiplication
        return (*this) = (*this) * other;
    }

    /// \param[in] scalar : A double that each element of the matrix is multiplied by
    /// \return Matrix
    MiddleSizeMatrix& MiddleSizeMatrix::operator*=(const type &scalar)
    {
#ifdef LA_STL_VECTOR
        for (type& d : data_)
            d *= scalar;
#else
        data_ *= scalar;
#endif
        return *this;
    }

    /// \param[in] scalar : A double that each element of the matrix is divided by
    /// \return Matrix
    MiddleSizeMatrix& MiddleSizeMatrix::operator/=(const type& scalar)
    {
#ifdef LA_STL_VECTOR
        for (type& d : data_)
            d /= scalar;
#else
        data_/=scalar;
#endif
        return *this;
    }

    MiddleSizeMatrix MiddleSizeMatrix::operator +(const MiddleSizeMatrix& other) const
    {
        MiddleSizeMatrix result(*this);
        return result += other;
    }

    MiddleSizeMatrix MiddleSizeMatrix::operator -(const MiddleSizeMatrix& other) const
    {
        MiddleSizeMatrix result(*this);
        return result -= other;
    }

    MiddleSizeMatrix MiddleSizeMatrix::operator -() const
    {
        return (*this) * -1.;
    }


    /*! \details Computes Matrix * vector and return the vector
     This is done by calling the BLAS (level 2) routine dgemv.
     */
    MiddleSizeVector MiddleSizeMatrix::operator*(MiddleSizeVector& right) const
    {
        logger.assert(nCols_ == right.size(), "Matrix-vector multiplication with mismatching sizes");
        
        if (nRows_ == 0)
        {
            logger(WARN, "Trying to multiply a vector with a matrix without any rows.");
            return MiddleSizeVector(0);
        }
        if (nCols_ == 0)
        {
            logger(WARN, "Trying to multiply a vector with a matrix without any columns.");
            return MiddleSizeVector(nRows_);
        }
        int nr = nRows_;
        int nc = nCols_;
        
        int i_one = 1;
        type d_one = 1.0;
        type d_zero = 0.0;
        
        MiddleSizeVector result(nr);
        
        logger(DEBUG, "Matrix size: % x % \n Vector size: %", nr, nc, right.size());
#ifdef HPGEM_USE_COMPLEX_PETSC
        zgemv_("N", &nr, &nc, &d_one, ((*(const_cast<MiddleSizeMatrix *>(this))).data()), &nr, right.data(), &i_one, &d_zero, result.data(), &i_one);
#else
        dgemv_("N", &nr, &nc, &d_one, ((*(const_cast<MiddleSizeMatrix *>(this))).data()), &nr, right.data(), &i_one, &d_zero, result.data(), &i_one);
#endif
        return result;
    }


    /*! \details Computes Matrix * vector and return the vector
     This is done by calling the BLAS (level 2) routine dgemv.
     */
    MiddleSizeVector MiddleSizeMatrix::operator*(MiddleSizeVector& right)
    {
        logger.assert(nCols_ == right.size(), "Matrix-vector multiplication with mismatching sizes");

        if (nRows_ == 0)
        {
            logger(WARN, "Trying to multiply a vector with a matrix without any rows.");
            return MiddleSizeVector(0);
        }
        if (nCols_ == 0)
        {
            logger(WARN, "Trying to multiply a vector with a matrix without any columns.");
            return MiddleSizeVector(nRows_);
        }

        int nr = nRows_;
        int nc = nCols_;

        int i_one = 1;
        type d_one = 1.0;
        type d_zero = 0.0;

        MiddleSizeVector result(nr);

        logger(DEBUG, "Matrix size: % x % \n Vector size: %", nr, nc, right.size());

#ifdef HPGEM_USE_COMPLEX_PETSC
        zgemv_("N", &nr, &nc, &d_one, this->data(), &nr, right.data(), &i_one, &d_zero, result.data(), &i_one);
#else
        dgemv_("N", &nr, &nc, &d_one, this->data(), &nr, right.data(), &i_one, &d_zero, result.data(), &i_one);
#endif
        return result;
    }
    
    /// \param[in] other : Matrix on the right of the multiplication
    /// \return Matrix
    /// \details
    /*! This uses the BLAS level 3 libaray dgemm to undertake the calculated
     Note it create the matrix that is return, but this is required as the return matrix may be a different size.
     */
    MiddleSizeMatrix MiddleSizeMatrix::operator*(const MiddleSizeMatrix &other)
    {
        logger.assert(nCols_ == other.nRows_, "Inner dimensions not equal.");

        if (nCols_ == 0)
        {
            logger(WARN, "Trying to multiply a matrix with a matrix without any columns.");
            return MiddleSizeMatrix(nRows_, other.getNCols());
        }
        int i = nRows_;
        int j = nCols_;
        int k = other.getNCols();
        
        ///The result of the matrix is left.Nrows, right.NCols()
        MiddleSizeMatrix C(i, k);
        
        type d_one = 1.0;
        type d_zero = 0.0;
        
        //Let the actual multiplication be done by Fortran
#ifdef HPGEM_USE_COMPLEX_PETSC
        zgemm_("N", "N", &i, &k, &j, &d_one, ((*this).data()), &i, const_cast<type*>(other.data()), &j, &d_zero, C.data(), &i);
#else
        dgemm_("N", "N", &i, &k, &j, &d_one, ((*this).data()), &i, const_cast<type*>(other.data()), &j, &d_zero, C.data(), &i);
#endif
        
        return C;
    }
    
    MiddleSizeMatrix MiddleSizeMatrix::operator*(const MiddleSizeMatrix &other) const
    {
        
        logger.assert(nCols_ == other.nRows_, "Inner dimensions are not the same.");

        if (nCols_ == 0)
        {
            logger(WARN, "Trying to multiply a matrix with a matrix without any columns.");
            return MiddleSizeMatrix(nRows_, other.getNCols());
        }
        int i = nRows_;
        int j = nCols_;
        int k = other.getNCols();
        
        //The result of the matrix is left.Nrows, right.NCols()
        MiddleSizeMatrix C(i, k);
        
        type d_one = 1.0;
        type d_zero = 0.0;
        
        //Let the actual multiplication be done by Fortran
#ifdef HPGEM_USE_COMPLEX_PETSC
        zgemm_("N", "N", &i, &k, &j, &d_one, ((*(const_cast<MiddleSizeMatrix *>(this))).data()), &i, (((const_cast<MiddleSizeMatrix&>(other))).data()), &j, &d_zero, C.data(), &i);
#else
        dgemm_("N", "N", &i, &k, &j, &d_one, ((*(const_cast<MiddleSizeMatrix *>(this))).data()), &i, (((const_cast<MiddleSizeMatrix&>(other))).data()), &j, &d_zero, C.data(), &i);
#endif
        
        return C;
    }
    
    /// \param[in] scalar : A double that each element of the matrix is multiplied by
    /// \return Matrix
    MiddleSizeMatrix MiddleSizeMatrix::operator*(const type &scalar) const
    {
        MiddleSizeMatrix result(*this);
        return (result *= scalar);
    }

    /// \param[in] scalar : A double that each element of the matrix is divided by
    /// \return Matrix
    MiddleSizeMatrix MiddleSizeMatrix::operator/(const type &scalar) const
    {
        MiddleSizeMatrix result(*this);
        return (result /= scalar);
    }
    
    /// \param [in] double c the  value all the entries of the matrix are set to
    ///
    /// \details Sets all the entries in the matrix equal to the scalar c.
    MiddleSizeMatrix& MiddleSizeMatrix::operator=(const type& c)
    {
        if (size() != 1)
        {
            nRows_ = 1;
            nCols_ = 1;
            data_.resize(1);
        }
#ifdef LA_STL_VECTOR
        data_[0] = c;
#else
        data_=c;
#endif
        return *this;
    }
    
    /// \param[in] Matrix : this is the matrix of the right hand side of the assignment
    MiddleSizeMatrix& MiddleSizeMatrix::operator=(const MiddleSizeMatrix& right)
    {
        data_ = (right.data_);
        nRows_ = right.nRows_;
        nCols_ = right.nCols_;
        return *this;
    }
    
    MiddleSizeMatrix& MiddleSizeMatrix::operator=(MiddleSizeMatrix&& right)
    {
        data_ = std::move(right.data_);
        nRows_ = right.nRows_;
        nCols_ = right.nCols_;
        return *this;
    }

    ///\return NumericalVector : The answer is return in this vector which is created by this function call
    ///\details This is the computeWedgeStuffVector magic taken directly from hpGEM version 1.0
    ///Here we repeat the orginal commment from hpGEM 1.0
    ///Compute the wedge product to find a vector which completes the
    ///vectors in the columns of the Jacobian matrix to a basis.
     
    ///Wedge product means: the i-th component of the solution vector is found by
    ///augmenting the Jacobian matrix on the right with the i-th standard basis
    ///vector, and computing the determinant of this square matrix.  At least for
    ///dimension 2 and 3 I do not form the square matrices, since the
    ///evaluation of the determinant is easy and can be inserted directly.
    MiddleSizeVector MiddleSizeMatrix::computeWedgeStuffVector() const
    {
        logger.assert(nCols_ == nRows_ - 1, "Matrix has wrong dimensions to construct the wedge stuff vector");
        MiddleSizeVector result(nRows_);
        
        switch (nRows_)
        {
            case 2:
                result[0] = -(*this)(1, 0);
                result[1] = +(*this)(0, 0);
                break;
            case 3:
                result[0] = (*this)(1, 0) * (*this)(2, 1) - (*this)(2, 0) * (*this)(1, 1);
                result[1] = (*this)(0, 1) * (*this)(2, 0) - (*this)(0, 0) * (*this)(2, 1); // includes minus sign already!
                result[2] = (*this)(0, 0) * (*this)(1, 1) - (*this)(1, 0) * (*this)(0, 1);
                break;
            case 4:
                result[0] = (*this)(1, 0) * (-(*this)(2, 1) * (*this)(3, 2) + (*this)(3, 1) * (*this)(2, 2)) + (*this)(2, 0) * ((*this)(1, 1) * (*this)(3, 2) - (*this)(3, 1) * (*this)(1, 2)) + (*this)(3, 0) * (-(*this)(1, 1) * (*this)(2, 2) + (*this)(2, 1) * (*this)(1, 2));
                
                result[1] = (*this)(0, 0) * ((*this)(2, 1) * (*this)(3, 2) - (*this)(3, 1) * (*this)(2, 2)) + (*this)(2, 0) * (-(*this)(0, 1) * (*this)(3, 2) + (*this)(3, 1) * (*this)(0, 2)) + (*this)(3, 0) * ((*this)(0, 1) * (*this)(2, 2) - (*this)(2, 1) * (*this)(0, 2));
                result[2] = (*this)(0, 0) * (-(*this)(1, 1) * (*this)(3, 2) + (*this)(3, 1) * (*this)(1, 2)) + (*this)(1, 0) * ((*this)(0, 1) * (*this)(3, 2) - (*this)(3, 1) * (*this)(0, 2)) + (*this)(3, 0) * (-(*this)(0, 1) * (*this)(1, 2) + (*this)(1, 1) * (*this)(0, 2));
                result[3] = (*this)(0, 0) * ((*this)(1, 1) * (*this)(2, 2) - (*this)(2, 1) * (*this)(1, 2)) + (*this)(1, 0) * (-(*this)(0, 1) * (*this)(2, 2) + (*this)(2, 1) * (*this)(0, 2)) + (*this)(2, 0) * ((*this)(0, 1) * (*this)(1, 2) - (*this)(1, 1) * (*this)(0, 2));
                break;
            default:
                logger(ERROR, "Wedge product not implemented for this dimension");
        } //end switch
        
        return (result);
        
    }

    /// \param[in] a : double scalar that is multiple by the matrix x
    /// \param[in] x : matrix that is multiple 
    ///
    /// \details Adds to the matrix a*X_ij where a is scalar and X is a matrix
    void MiddleSizeMatrix::axpy(type a, const MiddleSizeMatrix& x)
    {
        
        int size = nRows_ * nCols_;
        logger.assert(nRows_ == x.nRows_, "Dimensions are not the same.");
        logger.assert(nCols_ == x.nCols_, "Dimensions are not the same.");
        int i_one = 1;
        
#ifdef HPGEM_USE_COMPLEX_PETSC
        zaxpy_(&size, &a, const_cast<type *>(x.data()), &i_one, data(), &i_one);
#else
        daxpy_(&size, &a, const_cast<type *>(x.data()), &i_one, data(), &i_one);
#endif
        
    }
    
    /// \param[in] n the number of row in the new matrix
    /// \param[in] m the number of columns in the new matrix
    void MiddleSizeMatrix::resize(std::size_t n, std::size_t m)
    {
        logger.assert(n <= std::numeric_limits<int>::max() && m <= std::numeric_limits<int>::max(), "Dense linear algebra is not supported on this system for matrices that are this large");
        nRows_ = n;
        nCols_ = m;
        if (n * m != data_.size())
        {
            data_.resize(nRows_ * nCols_);
        }
    }
    
    /// If two matrices have the same number of columns, glue them together.
    /// \todo Find a more elegant way to do this.
    void MiddleSizeMatrix::concatenate(const MiddleSizeMatrix& other)
    {
        logger.assert(nCols_ == other.nCols_, "Number of columns is not the same.");
        logger.assert(nRows_+other.nRows_ <= std::numeric_limits<int>::max(), "Dense linear algebra is not supported on this system for matrices that are this large");
        
#ifdef LA_STL_VECTOR
        std::vector<type> data_new(nCols_ * (nRows_ + other.nRows_));
#else
        std::valarray<type> data_new(nCols_ * (nRows_ + other.nRows_));
#endif
        
        for (std::size_t col = 0; col < nCols_; ++col)
        {
            //First insert the values of this matrix, then of the other matrix.
            //Index row stands for the row number in the new matrix.
            for (std::size_t row = 0; row < nRows_; ++row)
            {
                data_new[row + col * (nRows_ + other.nRows_)] = data_[row + col * nRows_];
            }
            for (std::size_t row = nRows_; row < nRows_ + other.nRows_; ++row)
            {
                data_new[row + col * (nRows_ + other.nRows_)] = other.data_[(row - nRows_) + col * other.nRows_];
            }
        }
        nRows_ += other.nRows_;
        data_ = data_new;
    }
    
    /// \return the total number of entries
    std::size_t MiddleSizeMatrix::size() const
    {
        return nRows_ * nCols_;
    }
    
    /// \return the number of rows
    std::size_t MiddleSizeMatrix::getNRows() const
    {
        return nRows_;
    }
    
    /// \brief Get the number of columns
    /// \return int : the number of columns
    std::size_t MiddleSizeMatrix::getNCols() const
    {
        return nCols_;
    }
    
    LinearAlgebra::MiddleSizeVector MiddleSizeMatrix::getColumn(std::size_t j) const
    {
        logger.assert(j < nCols_, "Requested column %, but there are only % columns", j, nCols_);
        LinearAlgebra::MiddleSizeVector ret(nRows_);
        for (std::size_t i = 0; i < nRows_; ++i)
        {
            ret[i] = data_[j * nRows_ + i];
        }
        return ret;
    }
    
    LinearAlgebra::MiddleSizeVector MiddleSizeMatrix::getRow(std::size_t i) const
    {
        logger.assert(i < nRows_, "Requested row %, but there are only % rows", i, nRows_);
        LinearAlgebra::MiddleSizeVector ret(nCols_);
        for (std::size_t j = 0; j < nCols_; ++j)
        {
            ret[j] = data_[j * nRows_ + i];
        }
        return ret;
    }
    
    /// return Matrix which is the LUfactorisation of the current matrix
    MiddleSizeMatrix MiddleSizeMatrix::LUfactorisation() const
    {
        
        int nr = nRows_;
        int nc = nCols_;
        int nPivot = std::min(nRows_, nCols_);
        std::vector<int> iPivot(nPivot);
        
        MiddleSizeMatrix result(*this);
        
        int info;

#ifdef HPGEM_USE_COMPLEX_PETSC
        zgetrf_(&nr, &nc, result.data(), &nr, iPivot.data(), &info);
#else
        dgetrf_(&nr, &nc, result.data(), &nr, iPivot.data(), &info);
#endif
        
        return result;
    }
    
    /// \param[out] result this is the inverse of the current matrix
    MiddleSizeMatrix MiddleSizeMatrix::inverse() const
    {
        logger.assert(nRows_ == nCols_, "Cannot invert a non-square matrix");
        MiddleSizeMatrix result = (*this);
        
        int nr = nRows_;
        int nc = nCols_;
        
        int nPivot = std::min(nRows_, nCols_);
        std::vector<int> iPivot(nPivot);
        
        int info = 0;
        
        int lwork = nRows_ * nCols_;
        
        MiddleSizeMatrix work(nRows_, nCols_);

#ifdef HPGEM_USE_COMPLEX_PETSC
        zgetrf_(&nr, &nc, result.data(), &nr, iPivot.data(), &info);
        zgetri_(&nc, result.data(), &nc, iPivot.data(), work.data(), &lwork, &info);
#else
        dgetrf_(&nr, &nc, result.data(), &nr, iPivot.data(), &info);
        dgetri_(&nc, result.data(), &nc, iPivot.data(), work.data(), &lwork, &info);
#endif

        
        return result;
    }
    
    MiddleSizeMatrix MiddleSizeMatrix::transpose() const
    {
        MiddleSizeMatrix result(nCols_, nRows_);
        for(std::size_t i = 0; i < nRows_; ++i)
        {
            for(std::size_t j = 0; j < nCols_; ++j)
            {
                result(j, i) = (*this)(i, j);
            }
        }
        return result;
    }

    /// \param[in,out] B. On enter is B in Ax=B and on exit is x.
    void MiddleSizeMatrix::solve(MiddleSizeMatrix& B) const
    {
        logger.assert(nRows_ == nCols_, "can only solve for square matrixes");
        logger.assert(nRows_ == B.nRows_, "size of the RHS does not match the size of the matrix");
        
        int n = nRows_;
        int nrhs = B.getNCols();
        int info;
        
        std::vector<int> IPIV(n);
        
        MiddleSizeMatrix matThis = *this;

#ifdef HPGEM_USE_COMPLEX_PETSC
        zgesv_(&n, &nrhs, matThis.data(), &n, IPIV.data(), B.data(), &n, &info);
#else
        dgesv_(&n, &nrhs, matThis.data(), &n, IPIV.data(), B.data(), &n, &info);
#endif
    }

    void MiddleSizeMatrix::solve(MiddleSizeVector& b) const
    {
        logger.assert(nRows_ == nCols_, "can only solve for square matrixes");
        logger.assert(nRows_ == b.size(), "size of the RHS does not match the size of the matrix");
        
        int n = nRows_;
        int nrhs = 1;
        int info;
        
        std::vector<int> IPIV(n);
        
        MiddleSizeMatrix matThis = *this;

#ifdef HPGEM_USE_COMPLEX_PETSC
        zgesv_(&n, &nrhs, matThis.data(), &n, IPIV.data(), b.data(), &n, &info);
#else
        dgesv_(&n, &nrhs, matThis.data(), &n, IPIV.data(), b.data(), &n, &info);
#endif
    }
    
    MiddleSizeMatrix::type* MiddleSizeMatrix::data()
    {
        return data_.data();
    }
    
    const MiddleSizeMatrix::type* MiddleSizeMatrix::data() const
    {
        return data_.data();
    }
    
    ///Print the matrix with () around each line and [] around the matrix.
    std::ostream& operator<<(std::ostream& os, const MiddleSizeMatrix& A)
    {
        std::size_t nRows = A.getNRows();
        std::size_t nCols = A.getNCols();
        os << "[" << std::endl;
        for (std::size_t i = 0; i < nRows; ++i)
        {
            os << "(";
            for (std::size_t j = 0; j < nCols; ++j)
            {
                os << A(i, j) << "\t ";
            }
            os << ")" << std::endl;
        }
        os << "]";
        return os;
    }
    
    MiddleSizeMatrix operator*(const MiddleSizeMatrix::type d, const MiddleSizeMatrix& mat)
    {
        MiddleSizeMatrix matNew = mat;
        matNew *= d;
        return matNew;
    }


    MiddleSizeVector operator*(MiddleSizeVector& left, MiddleSizeMatrix& right)
    {
        logger.assert(right.getNRows() == left.size(), "Matrix-vector multiplication with mismatching sizes");

        if (right.getNCols() == 0)
        {
            logger(WARN, "Trying to multiply a vector with a matrix without any columns.");
            return MiddleSizeVector(0);
        }
        int nr = right.getNRows();
        int nc = right.getNCols();

        int i_one = 1;
        MiddleSizeMatrix::type d_one = 1.0;
        MiddleSizeMatrix::type d_zero = 0.0;

        MiddleSizeVector result(nc);

        logger(DEBUG, "Matrix size: % x % \n Vector size: %", nr, nc, left.size());

#ifdef HPGEM_USE_COMPLEX_PETSC
        zgemv_("T", &nr, &nc, &d_one, right.data(), &nr, left.data(), &i_one, &d_zero, result.data(), &i_one);
#else
        dgemv_("T", &nr, &nc, &d_one, right.data(), &nr, left.data(), &i_one, &d_zero, result.data(), &i_one);
#endif
        return result;
    }
}
