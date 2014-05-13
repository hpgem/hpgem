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
#include "GlobalNamespaceLinearAlgebra.hpp"
#include "Matrix.hpp"




namespace LinearAlgebra
{
    
    
    extern "C"
    {
        
        /// This does matrix times vector and is from blas level 2.
        void dgemv_(const char* trans, int* m, int* n, double* alpha, double* A, int* LDA, double *x, int *incx, double *beta, double *y, int *incy);
        
        ///This is the gernal matrix multiplication from blas level 3
        int dgemm_(const char *transA, const char *transB, int *M, int *N, int *k, double *alpha, double *A, int *LDA, double *B, int *LDB, double *beta, double *C, int *LDC); 
        
        ///This is the gerneral scalar times vector + vector from blas, hence from blas level 1. Here we also use on a matrix by treating as a vector
        int daxpy_(unsigned int* N, double* DA, double* DX,unsigned int* INCX, double* DY, unsigned int* INCY);
        
        /// This is LU factorisation of the matrix A. This has been taken from LAPACK 
        void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
        
        /// This is the inverse calulation also from LAPACK. Calculates inverse if you pass it the LU factorisation.
        void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
        
        /// This is used for solve Ax=B for x. Again this is from LAPACK.
        void dgesv_(int* N, int* NRHS, double* A, int* lda,  int* IPIV, double* B, int* LDB, int* INFO);
        
        
        //Tito's def
        void sgetrf_(const int*, const int*, float*, const int*, int*, int*);
        //void dgetrf_(const int*, const int*, double*, const int*, int*, int*);
        
        void sgetrs_(const char*, const int*, const int*, const float*, const int*, const int*, float*, const int*, int*);
        void dgetrs_(const char*, const int*, const int*, const double*, const int*, const int*, double*, const int*, int*);
        
        void sgesv_(const int*, const int*, float*, const int*, int*, float*, const int*, int*);
       // void dgesv_(const int*, const int*, double*, const int*, int*, double*, const int*, int*);
        
        
        
        
        // call dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
        
        //int dgemm_(const char *transA, const char *transB, int *M, int *N, int *k, double *alpha, double *A, int *LDA, double *B, int *LDB, double *beta, double *C, int *LDC);
        
    }
    
    
    Matrix::Matrix():
    	nRows_(0),
    	nCols_(0)
    {
    }
    
    
    /// \param[in]  n The number of rows the matrix will have
    /// \param[in]  m The number of columns the matrix will have
    Matrix::Matrix(const int n, const int m):
    	data_(n*m),
    	nRows_(n),
    	nCols_(m)
    {
    }
    
    
    /// \param[in]  n The number of rows the matrix will have
    /// \param[in]  m The number of columns the matrix will have
    /// \param[in]  c The value all the entries of the matrix will contain
    ///
    /// \details
    /// Example usage : Matrix<double> A(4,3,2) with create a 4 by 3 matrix called A with entries equal to 2.
    Matrix::Matrix(const int n, const int m, const double& c):
    	data_(c, n*m),
    	nRows_(n),
    	nCols_(m)
    {
    
    }
    
    
    /// \param[in] Matrix A i.e. the matrix to be copies.
    Matrix::Matrix(const Matrix& other):
    	data_(other.data_),
    	nRows_(other.nRows_),
    	nCols_(other.nCols_)
    {
    }
    
    
//    /// \param[in] n The number of the row you want the element from
//    /// \param[in] m The numebr of the column you want the element from
//    /// \return double ref i.e. the value of the element you requested
//    /// \bug Range checking has not been added yet.
//*****S.N need to be brought to the hpp, due to compiling error in realease mode. Strange bug, need to be revistied

//    inline double& Matrix::operator()(int n, int m)
//    {
//    	return data_[n + m*nRows_];
//    }
    
    
    /// \param[in] n The number of the row you want the element from
    /// \param[in] m The numebr of the column you want the element from
    /// \return double ref i.e. the value of the element you requested
    /// \bug Range checking has not been added yet.
    //*****S.N need to be brought to the hpp, due to compiling error in realease mode. Strange bug, need to be revistied
    //inline const double& Matrix::operator() (int n, int m) const {return data_[n + m*nRows_];}
    
    
    /// \param[in] n The number of the row you want the element from
    /// \return double i.e. the value of the element you requested
    /// \bug Range checking has not been added yet.
    ///
    /// \details
    /// Recall that the matrix is stored in fortran style i.e. columns first and then rows
    double& Matrix::operator[](const int n){return data_[n];}
    
    const double& Matrix::operator[](const int n) const  {return data_[n];}
    
    /// \details
    /*! Computes Matrix * vector and return the vector
     This is done by calling the BLAS (level 2) routine dgemv.
     */
    NumericalVector Matrix::operator*(NumericalVector& right)
    {
        int nr=nRows_;
        int nc=nCols_;            
        
        int i_one=1;
        double d_one=1.0;
        double d_zero=0.0;
        
        //NumericalVector result(other);
        NumericalVector result(nc);
        
       dgemv_("N", &nr, &nc, &d_one, &((*this)[0]), &nr,&right[0],&i_one, &d_zero, &result[0], &i_one);
        
        return result;        
    }
    
    /// \details
    /*! Computes Matrix * vector and return the vector
     This is done by calling the BLAS (level 2) routine dgemv.
     */
    NumericalVector Matrix::operator*( NumericalVector& right) const
    {
        int nr=nRows_;
        int nc=nCols_;            
        
        int i_one=1;
        double d_one=1.0;
        double d_zero=0.0;
        
        NumericalVector result(nc);
        
        dgemv_("N", &nr, &nc, &d_one, &((*(const_cast<Matrix *> (this)))[0]), &nr,&right[0],&i_one, &d_zero, &result[0], &i_one);
        return result;
    }
    
    
    /// \param[in] other : Matrix on the right of the multiplication
    /// \return Matrix
    /// \details
    /// \bug Direct calculation may be quicker this needs testing.
    /*! This uses the BLAS level 3 libaray dgemm to undertake the calculated
     Note it create the matrix that is return, but this is required as the return matrix may be a different size.
     */
    Matrix Matrix::operator* (Matrix &other )
    {
        
        if (nCols_!=other.nRows_)
        {
            /// \bug this need fixed when we have real error handling.
            throw(10);
        }
        
        
        int i=nRows_;
        int j=nCols_;
        int k=other.getNCols();
        
        ///The result of the matrix is left.Nrows, right.NCols()
        Matrix C(i,k);
        
        int i_one=1;
        double d_one=1.0;
        double d_zero=0.0;
        
        
        
        
        dgemm_("N","N",&i,&k,&j,&d_one,&((*this)[0]),&i,&other[0],&j,&d_zero,&C[0],&i);
        
        
        
        return C;
    }
    
    Matrix Matrix::operator* (const Matrix &other )const
    {
        
        if (nCols_!=other.nRows_)
        {
                /// \bug this need fixed when we have real error handling.
            throw(10);
        }
        
        
        int i=nRows_;
        int j=nCols_;
        int k=other.getNCols();
        
            ///The result of the matrix is left.Nrows, right.NCols()
        Matrix C(i,k);
        
        int i_one=1;
        double d_one=1.0;
        double d_zero=0.0;
        
        
        
        
        dgemm_("N","N",&i,&k,&j,&d_one,&((*(const_cast<Matrix *> (this)))[0]),&i,&(((const_cast<Matrix&> (other)))[0]),&j,&d_zero,&C[0],&i);
        
        
        
        return C;
    }
    

    /// \param[in] scale : A double that each element of the matrix is multiplied bu
    /// \return Matrix
    Matrix Matrix::operator*= (const double &scalar)
    {
        data_*=scalar;
        return *this;
    }
    
    
    Matrix& Matrix::operator/= (const double& scalar)
    {
        data_/=scalar;
        return *this;
    }
    
    Matrix Matrix::operator/ (const double &scalar)
    {
        Matrix result(*this);
        return (result/=scalar);
    }
    
    
    /// \param [in] double c the  value all the entries of the matrix are set to
    ///
    /// \details
    /// Sets all the entries in the matrix equal to the scalar c.
    Matrix& Matrix::operator=(const double& c)
    {
        if (size()==1){nRows_=1; nCols_=1; data_.resize(1);}
        data_=c;
        return *this;
    }
    
    
    /// \param[in] Matrix : this is the matrix of the right hand side of the assigment
    /// \bug Error checking needs to be added
    Matrix& Matrix::operator=(const Matrix& right) 
    	{
    	   data_=(right.data_);
     	   return *this;
    	};
    
    
    
    
    /// \param[out] result : This returns the result in a passed in NumericalVector 
    /// \details
    /// \see computeWedgeStuffVector()
    /*! This is the computeWedgeStuffVector magic taken directly from hpGEM version 1.0
     Here we repeat the orginal commment from hpGEM 1.0
     Compute the wedge product to find a vector which completes the
     vectors in the columns of the Jacobian matrix to a basis.
     
     Wedge product means: the i-th component of the solution vector is found by
     augmenting the Jacobian matrix on the right with the i-th standard basis
     vector, and computing the determinant of this square matrix.  At least for
     dimension 2 and three I do not form the square matrices, since the
     evaluation of the determinant is easy and can be inserted directly.
     */
    void Matrix::computeWedgeStuffVector(NumericalVector& result)
    {
        
        //if (nRows_ != nCols_){throw("Wedge product only defined for square matrices");}///Wrong...
        
        if (nRows_ != result.size()){throw("Passed vector is the wrong size for a wedge product");}
        
        
        switch (nRows_)
        {
            case 2 :
                result[0] = - (*this)(1,0);
                result[1] = + (*this)(0,0);
                break;
            case 3:
                result[0] = (*this)(1,0) * (*this)(2,1) - (*this)(2,0) * (*this)(1,1);
                result[1] = (*this)(0,1) * (*this)(2,0) - (*this)(0,0) * (*this)(2,1); // includes minus sign already!
                result[2] = (*this)(0,0) * (*this)(1,1) - (*this)(1,0) * (*this)(0,1);
                break;
            case 4:
                result[0] = (*this)(1,0) * (-(*this)(2,1)*(*this)(3,2) + (*this)(3,1)*(*this)(2,2)) +
                (*this)(2,0) * ( (*this)(1,1)*(*this)(3,2) - (*this)(3,1)*(*this)(1,2)) +
                (*this)(3,0) * (-(*this)(1,1)*(*this)(2,2) + (*this)(2,1)*(*this)(1,2));
                
                result[1] = (*this)(0,0) * ( (*this)(2,1)*(*this)(3,2) - (*this)(3,1)*(*this)(2,2)) +
                (*this)(2,0) * (-(*this)(0,1)*(*this)(3,2) + (*this)(3,1)*(*this)(0,2)) +
                (*this)(3,0) * ( (*this)(0,1)*(*this)(2,2) - (*this)(2,1)*(*this)(0,2));
                result[2] = (*this)(0,0) * (-(*this)(1,1)*(*this)(3,2) + (*this)(3,1)*(*this)(1,2)) +
                (*this)(1,0) * ( (*this)(0,1)*(*this)(3,2) - (*this)(3,1)*(*this)(0,2)) +
                (*this)(3,0) * (-(*this)(0,1)*(*this)(1,2) + (*this)(1,1)*(*this)(0,2));
                result[3] = (*this)(0,0) * ( (*this)(1,1)*(*this)(2,2) - (*this)(2,1)*(*this)(1,2)) +
                (*this)(1,0) * (-(*this)(0,1)*(*this)(2,2) + (*this)(2,1)*(*this)(0,2)) +
                (*this)(2,0) * ( (*this)(0,1)*(*this)(1,2) - (*this)(1,1)*(*this)(0,2));
                break;
            default:
                throw("Wedge product not defined for this dimension");
        }//end switch
        
        
    }
    
    /// \details
    /// \param[out] result : This returns the result in a passed in NumericalVector 
    /// \see computeWedgeStuffVector()
    /*! This is the computeWedgeStuffVector magic taken directly from hpGEM version 1.0
     Here we repeat the orginal commment from hpGEM 1.0
     Compute the wedge product to find a vector which completes the
     vectors in the columns of the Jacobian matrix to a basis.
     
     Wedge product means: the i-th component of the solution vector is found by
     augmenting the Jacobian matrix on the right with the i-th standard basis
     vector, and computing the determinant of this square matrix.  At least for
     dimension 2 and three I do not form the square matrices, since the
     evaluation of the determinant is easy and can be inserted directly.
     */
    void Matrix::computeWedgeStuffVector(NumericalVector& result) const
    {
        
//         if (nRows_ != nCols_){std::cout<<"Wedge product only defined for square matrices"<<std::endl;}
        
//         if (nRows_ != result.size()){std::cout<<"Passed vector is the wrong size for a wedge product"<<std::endl;}
        
        
        switch (nRows_)
        {
            case 2 :
                result[0] = - (*this)(1,0);
                result[1] = + (*this)(0,0);
                break;
            case 3:
                result[0] = (*this)(1,0) * (*this)(2,1) - (*this)(2,0) * (*this)(1,1);
                result[1] = (*this)(0,1) * (*this)(2,0) - (*this)(0,0) * (*this)(2,1); // includes minus sign already!
                result[2] = (*this)(0,0) * (*this)(1,1) - (*this)(1,0) * (*this)(0,1);
                break;
            case 4:
                result[0] = (*this)(1,0) * (-(*this)(2,1)*(*this)(3,2) + (*this)(3,1)*(*this)(2,2)) +
                (*this)(2,0) * ( (*this)(1,1)*(*this)(3,2) - (*this)(3,1)*(*this)(1,2)) +
                (*this)(3,0) * (-(*this)(1,1)*(*this)(2,2) + (*this)(2,1)*(*this)(1,2));
                
                result[1] = (*this)(0,0) * ( (*this)(2,1)*(*this)(3,2) - (*this)(3,1)*(*this)(2,2)) +
                (*this)(2,0) * (-(*this)(0,1)*(*this)(3,2) + (*this)(3,1)*(*this)(0,2)) +
                (*this)(3,0) * ( (*this)(0,1)*(*this)(2,2) - (*this)(2,1)*(*this)(0,2));
                result[2] = (*this)(0,0) * (-(*this)(1,1)*(*this)(3,2) + (*this)(3,1)*(*this)(1,2)) +
                (*this)(1,0) * ( (*this)(0,1)*(*this)(3,2) - (*this)(3,1)*(*this)(0,2)) +
                (*this)(3,0) * (-(*this)(0,1)*(*this)(1,2) + (*this)(1,1)*(*this)(0,2));
                result[3] = (*this)(0,0) * ( (*this)(1,1)*(*this)(2,2) - (*this)(2,1)*(*this)(1,2)) +
                (*this)(1,0) * (-(*this)(0,1)*(*this)(2,2) + (*this)(2,1)*(*this)(0,2)) +
                (*this)(2,0) * ( (*this)(0,1)*(*this)(1,2) - (*this)(1,1)*(*this)(0,2));
                break;
            default:
                std::cout<<"Wedge product not defined for this dimension"<<std::endl;
        }//end switch
        
        
    }
    
    /// \return NumericalVector : The answer is return in this vector which is created by this function call
    /// \see computeWedgeStuffVector (NumericalVector)
    NumericalVector Matrix::computeWedgeStuffVector()
    {
        NumericalVector result(nRows_);
        computeWedgeStuffVector(result);
        return(result);
        
    }
    
    
    /// \param[in] a : double scalar that is multiple by the matrix x
    /// \param[in] x : matrix that is multiple 
    /// \bug Does not undertake range checking.
    ///
    /// \details
    /*!
     *  Adds to the matrix a*X_ij where a is scalar and X is a matrix
     !*/
    void Matrix::axpy(double a, const Matrix& x)
    {
     
        unsigned int size=nRows_*nCols_;
     
        unsigned int i_one=1;
     
        
        daxpy_(&size, &a, &((*(const_cast<Matrix *> (&x)))[0]), &i_one, &((*this)[0]) , &i_one);
        
        
    }
    
    
  
    
    
    
    /// \param[in] n the number of row in the new matrix
    /// \param[in] m the number of columns in the new matrix
    void Matrix::resize(int n, int m){nRows_=n; nCols_=m; data_.resize(nRows_*nCols_);}
    
    
    /// \return int : the total number of entries
    const int Matrix::size() {return nRows_*nCols_;}
    
    
    /// \return int : the number of rows
    const int Matrix::getNRows() {return nRows_;}
    
    
    /// \brief Get the number of columns
    /// \return int : the number of columns
    const int Matrix::getNCols() {return nCols_;}
    
    
    
    /// \return int : the number of rows
    const int Matrix::getNRows() const {return nRows_;}
    
    
    /// \brief Get the number of columns
    /// \return int : the number of columns
    const int Matrix::getNCols() const {return nCols_;}
    
    /// return Matrix which is the LUfactorisation of the current matrix
    Matrix Matrix::LUfactorisation() const
    {
        
        int nr=nRows_;
        int nc=nCols_; 
        int nPivot=std::min(nRows_,nCols_);
        int iPivot[nPivot];
        
        Matrix result(*this);
    
        int info;
        
        dgetrf_(&nr,&nc,&result[0],&nr,iPivot,&info);
        
        
        return result;
    }
    
    /// \param[out] result this is the inverse of the current matrix
    /// \bug if the dimensations of result are correct is not checked.
    void Matrix::inverse(Matrix& result) const
    {
     
        result=(*this);
        
        int nr=nRows_;
        int nc=nCols_;
        
        int nPivot=std::min(nRows_,nCols_);
        int iPivot[nPivot];
        
        int info=0;
        
        dgetrf_(&nr,&nc,&result[0],&nr,iPivot,&info);
        
        int lwork = nRows_*nCols_;
        
        double work[lwork];
        
        dgetri_(&nc,&result[0],&nc,iPivot,&work[0],&lwork,&info);

        
        
            //std::cout << "Mautrix="<<result<<std::endl;
        
//        std::cout << "ENDINVERSE"<<std::endl;
        
    }
    
    /// \param[in,out] B. On enter is B in Ax=B and on exit is x.
    void Matrix::solve(Matrix& B) const
    {
        
        
        int n=nRows_;
        int nrhs=B.getNCols();
        int info;
        
        int IPIV;
    
        
        dgesv_(&n,&nrhs,&((*(const_cast<Matrix *> (this)))[0]),&n,&IPIV,&B[0],&n,&info);
    
        
    }
    
    
    ostream& operator<<(ostream& os, const Matrix& A)
    {
        unsigned int nRows=A.getNRows();
        unsigned int nCols=A.getNCols();
        os << "["<<std::endl;
        for (int i=0; i<nRows; ++i)
        {
            os << "(";

             for(int j=0; j<nCols; ++j)
                os<< A(i,j) <<"\t ";
            os << ")"<<std::endl;
        }
         os << "]";
        
//        for (int i=0; i<nRows-1; ++i) 
//        {
//            os << "[(";
//            for(int j=0; j<nCols-1; ++j)
//                os<< A(i,j) <<",";
//            os << A(i, nCols-1) << "),\n ";
//        }
//    
//        os << "(";
//        for(int j=0; j<nCols-1; ++j)
//            os<< A(nRows-1,j) <<",";
//        os << A(nRows-1, nCols-1) << ")";
//        
//        os << "]";
        return os;
    }
    
    

    
    
    
    
}
/// \brief Defines vector B times matrix A and updates the values in A i.e. A_ij=B_,j A_ij
/// \param [in] NumericalVector
/// This should be in vector
// void operator*=(NumericalVector& left);
/// \details
/// Computes vector * matrix and updates the value in the matrix
/// This is done by calling the BLAS (level 2) routine gdemv
/// Note this operator does still require a copy of the Matrix
/// \bug This may be able to done quicker using directly i.e. not using hte BLAS libararies. 
//  void Matrix::operator*=(NumericalVector& left)
//    {
//        int nr=nRows_;
//        int nc=nCols_;
//        
//        
//        int i_one=1;
//        double d_one=1.0;
//        double d_zero=0.0;
//        
//        Matrix temp_copy(*this);
//        
//        dgemv_("N", &nr, &nc, &d_one,&temp_copy[0] , &nr,&left[0],&i_one, &d_zero, &((*this)[0]), &i_one);
//        
//        
//    }









    
