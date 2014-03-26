//------------------------------------------------------------------------------
// File: Matrix.hh
//
// March 2010,  M.T. Julianto
// July 2010: wrapper for lapack's gesv routines added
// Coverted and updated for hpGEM version 2 Feb 2013 (During the hpGEM intensive week)
//------------------------------------------------------------------------------
#ifndef MATRIX_HH
#define MATRIX_HH

// System includes
#include <iostream>
#include <valarray>
#include "NumericalVector.hpp"

namespace LinearAlgebra
{
    
    //We need the ostream for outputting and we encapulate from valarray.
    using std::ostream;
    using std::valarray;
    /// \class Matrix
    /// \brief Data type for small dense matrix.
    /// 
    /// \details
    /// Stores small dense matrix eficiently.
    /// Since this class is inherited from std::valarray, it also inherits members of std::valarray
    /// Note, valarry was speed tested and was shown to be quicker than stl vector and it very light weight
    /// It only store doubles as this is the main type linear algebra is done on in hpGEM
    /// It stores the matrix in fortran style to give quicker access to extern BLAS libaries.
    class Matrix
    {
    public:
        
        /// \brief Default Matrix constructor : Simply creates a zero size matrix
        Matrix();
        
        /// \brief Constructs a matrix of size n-rows by m-columns.
        Matrix(const int n, const int m);
        
        /// \brief Constructs a matrix of size n-rows by m-columns and initialises all entery to a constant 
        Matrix(const int n, const int m, const double& c);
    
        /// \brief Construct and copy Matrix from another Matrix i.e. B(A) where B and A are both matrices
        Matrix(const Matrix& other);
    
        /// \brief defines the operator (n,m) to access the element on row n and column m
        double& operator()(int n, int m)
            {
                return data_[n + m*nRows_];
            }
                
        /// \brief defines the operator (n,m) to access the element on row n and column m
        const double& operator() (int n, int m) const
        {
            return data_[n + m*nRows_];
        }

        /// WRONG!!
        /// \brief Access the n linear element in the matrix. 
        double& operator[](const int n);
        
            /// WRONG!!
        const double& operator[](const int n) const;
        
        /// \brief Defines Matrix A times vector B and return vector C i.e. C_,j= A_ij B_,j
        NumericalVector operator*(NumericalVector& right);
        
        /// \brief Defines Matrix A times vector B and return vector C i.e. C_,j= A_ij B_,j (constant version)
        NumericalVector operator*(NumericalVector& right) const;
        
        
            /// NEED ONE WITHOUT COPY!!!!
        /// \brief Does matrix A_ij = B_ik * C_kj
        Matrix operator* (Matrix &other);
        Matrix operator* (const Matrix &other)const;
        
        /// \brief Does matrix A_ij=scalar*A_ij
        //THIS SHOULD NOT RETURN ANYTHING, in PLACE S.N
        Matrix operator*= (const double &scalar);
        
        //THIS SHOULD NOT RETURN ANYTHING, in PLACE S.N
        /// \breif this does element by divived by a scalar 
        Matrix& operator/= (const double& scalar);
        
        /// \breif this does element by divived by a scalar
        Matrix operator/   (const double& scalar);
        
        /// \brief Assigns the Matrix by a scalar
        Matrix& operator=(const double& c);
        
        /// \brief Assigns one matrix to another.
        Matrix& operator=(const Matrix& right);
                    
        /// \brief computeWedgeStuffVector. The answer is return in result which you should precreate.
        void computeWedgeStuffVector(NumericalVector& result);
        
        /// \brief computeWedgeStuffVector. The answer is return in result which you should precreate.
        void computeWedgeStuffVector(NumericalVector& result) const;
        
        /// \brief computerWedgeStuffVector and create and returns it in a vector. 
        NumericalVector computeWedgeStuffVector();
        
        /// \brief Applies the matrix y=ax + y, where x is another matrix and a is a scalar
        void axpy(double a, const Matrix& x);
        
        /// \brief Resize the Matrix to be n-Rows by m-columns
        void resize(int n, int m);
        
        /// \brief Get total number of Matrix entries
        const int size();
        
        /// \brief Get the number of rows
        const int getNRows();
        
        /// \brief Get the number of columns
        const int getNCols();
        
        /// \brief Get the number of rows
        const int getNRows() const;
        
        /// \brief Get the number of columns
        const int getNCols() const;
        
        /// \brief Return the LUfactorisation of the matrix
        Matrix LUfactorisation() const;
        
        /// \brief return the inverse in the vector result. The size of result must match the matrix.
        void inverse(Matrix& result) const;
        
        /// \brief solves Ax=B where A is the current matrix and B is passed in. The result is returned in B.
        void solve(Matrix& B) const;
        
    private:
        /// The actually data of the matrix class
        valarray<double> data_;
        
        /// Stores the number of rows of the matrix
        unsigned int nRows_;
        
        /// Store the number of columns of the matrix
        unsigned int nCols_;
        
    };
    
    /// Writes nices format enteries of the Matrix A to the stream os.
    ostream& operator<<(ostream& os, const Matrix& A);
    
}
#endif
    


//
//	//----------------------------------------------------------------------
//	// Unary operators
//	//----------------------------------------------------------------------
//
//	/*! \brief Opposite sign of Matrix
//	 *
//	 *  Returns a matrix instance contains all entries of the matrix with the opposite sign.
//	 * 
//	 * <H3>Input Parameter(s):</H3>
//	 * <TABLE border="0" cellpadding="0" cellspacing="0">
//	 * <TR><TD WIDTH=50></TD> <TD><TD>none</TD></TR>
//	 * </TABLE>
//	 * <H3>Output Parameter(s):</H3>
//	 * <TABLE border="0" cellpadding="0" cellspacing="0">
//	 * <TR><TD WIDTH=50></TD> <TD><TD>returns an instance of Base::Matrix contains entries with the opposite sign</TD></TR>
//	 * </TABLE>
//	 * <H3>Note(s):</H3>
//	 * <TABLE border="0" cellpadding="0" cellspacing="0">
//	 * <TR><TD WIDTH=50></TD> <TD><TD>the entries of Base::Matrix <b>*this</b> itself are unchanged. </TD></TR>
//	 * </TABLE>
//	 * <H3>Example(s):</H3>
//	 * <PRE>
//	 * \#include <Base/Matrix.hh> 
//	 * 
//	 * Base::Matrix<double> A(5,3);  // matrix of dimension 5-rows by 3-columns
//	 * Base::Matrix<double> B(5,3);  // matrix of dimension 5-rows by 3-columns
//	 * A = 4.25;   // sets all entries of A to 4.25
//	 * B = -A;  // assigns B by entries of A with the opposite sign, the matrix A itself is unchanged
//	 * </PRE>
//	 **********************************************************************/
//	Matrix operator-()
//	{
//	   Matrix result(*this);
//	   result.valarray<EntryType>::operator*=(-1);
//	   return result;
//	};
//
//	// Other unary operators will be just inherited from the parent
//
//      private:
//        //! Private variable, the number of rows
//        int nrows_;
//
//        //! Private variable, the number of collumns
//        int ncols_;
//
//    };  // class Matrix
//
//
//    //==========================================================================
//    // Globally defined operations: (these might incomplete yet)
//    //--------------------------------------------------------------------------
//
//    /*! \brief Addition operation of two Matrices
//     *
//     *  Performs addition of two Base::Matrix.  Both matrices must have the same size.
//     *  An instance of Base::Matrix is returned.
//     * 
//     * <H3>Input Parameter(s):</H3>
//     * <TABLE border="0" cellpadding="0" cellspacing="0">
//     * <TR><TD WIDTH=50></TD> <TD><b>A</b></TD> <TD WIDTH=10></TD> <TD> - a Base::Matrix<EntryType></TD></TR>
//     * <TR><TD WIDTH=50></TD> <TD><b>B</b></TD> <TD WIDTH=10></TD> <TD> - a Base::Matrix<EntryType></TD></TR>
//     * </TABLE>
//     * <H3>Output Parameter(s):</H3>
//     * <TABLE border="0" cellpadding="0" cellspacing="0">
//     * <TR><TD WIDTH=50></TD><TD>returns an instance of Base::Matrix<EntryType> contains result of <b>A + B</b></TD></TR>
//     * </TABLE>
//     * <H3>Example(s):</H3>
//     * <PRE>
//     * \#include <Base/Matrix.hh> 
//     * 
//     * Base::Matrix<double> P(5,3);  // matrix of dimension 5-rows by 3-columns
//     * Base::Matrix<double> Q(5,3);  // matrix of dimension 5-rows by 3-columns
//     * Base::Matrix<double> R(5,3);  // matrix of dimension 5-rows by 3-columns
//     * P = 4.25;   // sets all entries of P to 4.25
//     * Q = -0.75;  // sets all entries of Q to -0.75
//     * R = P + Q;  // assigns R by the result of addition P + Q
//     * </PRE>
//     **************************************************************************/
//    template <class EntryType>
//    Matrix<EntryType> operator+(const Matrix<EntryType>& A, const Matrix<EntryType>& B) {
//	//TEST_ERROR_DEBUG((A.dim1() == B.dim1()) && (A.dim2() == B.dim2()), 
//	//	const char*, "Matrix: both matrices must have the same size!");
//		
//		/// \bug Error checking not yet implemented.
//	Matrix<EntryType> result(A);
//	result += B;
//	return result;
//    }
//
//
//    /*! \brief Addition operation of scallar and Matrix
//     *
//     *  Performs addition of a scallar and a Base::Matrix.  An instance of Base::Matrix is returned.
//     * 
//     * <H3>Input Parameter(s):</H3>
//     * <TABLE border="0" cellpadding="0" cellspacing="0">
//     * <TR><TD WIDTH=50></TD> <TD><b>c</b></TD> <TD WIDTH=10></TD> <TD> - a scallar</TD></TR>
//     * <TR><TD WIDTH=50></TD> <TD><b>A</b></TD> <TD WIDTH=10></TD> <TD> - a Base::Matrix<EntryType></TD></TR>
//     * </TABLE>
//     * <H3>Output Parameter(s):</H3>
//     * <TABLE border="0" cellpadding="0" cellspacing="0">
//     * <TR><TD WIDTH=50></TD><TD>returns an instance of Base::Matrix<EntryType> contains result of addition <b>c + A</b></TD></TR>
//     * </TABLE>
//     * <H3>Example(s):</H3>
//     * <PRE>
//     * \#include <Base/Matrix.hh> 
//     * 
//     * Base::Matrix<double> P(5,3);  // matrix of dimension 5-rows by 3-columns
//     * Base::Matrix<double> Q(5,3);  // matrix of dimension 5-rows by 3-columns
//     * double c = 0.91;  // a scallar
//     * P = 4.25;   // sets all entries of P to 4.25
//     * Q = c + P;  // assigns Q by the result of addition c + P
//     * </PRE>
//     **************************************************************************/
//    template <class EntryType>
//    Matrix<EntryType> operator+(const EntryType c, const Matrix<EntryType>& A) {
//	Matrix<EntryType> result(A);
//	result += c;
//	return result;
//    }
//
//    /*! \brief Addition operation of Matrix and scallar
//     *
//     *  Performs addition of a scallar and a Base::Matrix.  An instance of Base::Matrix is returned.
//     * 
//     * <H3>Input Parameter(s):</H3>
//     * <TABLE border="0" cellpadding="0" cellspacing="0">
//     * <TR><TD WIDTH=50></TD> <TD><b>A</b></TD> <TD WIDTH=10></TD> <TD> - a Base::Matrix<EntryType></TD></TR>
//     * <TR><TD WIDTH=50></TD> <TD><b>c</b></TD> <TD WIDTH=10></TD> <TD> - a scallar</TD></TR>
//     * </TABLE>
//     * <H3>Output Parameter(s):</H3>
//     * <TABLE border="0" cellpadding="0" cellspacing="0">
//     * <TR><TD WIDTH=50></TD><TD>returns an instance of Base::Matrix<EntryType> contains result of addition <b>A + c</b></TD></TR>
//     * </TABLE>
//     * <H3>Example(s):</H3>
//     * <PRE>
//     * \#include <Base/Matrix.hh> 
//     * 
//     * Base::Matrix<double> P(5,3);  // matrix of dimension 5-rows by 3-columns
//     * Base::Matrix<double> Q(5,3);  // matrix of dimension 5-rows by 3-columns
//     * double c = 0.91;  // a scallar
//     * P = 4.25;   // sets all entries of P to 4.25
//     * Q = P + c;  // assigns Q by the result of addition P + c
//     * </PRE>
//     **************************************************************************/
//    template <class EntryType>
//    Matrix<EntryType> operator+(const Matrix<EntryType>& A, const EntryType c) {
//	Matrix<EntryType> result(A);
//	result += c;
//	return result;
//    }
//
//    /*! \brief Substraction operation of two Matrices
//     *
//     *  Performs substraction of two Base::Matrix.    Both matrices must have the same size.
//     *  An instance of Base::Matrix is returned.
//     * 
//     * <H3>Input Parameter(s):</H3>
//     * <TABLE border="0" cellpadding="0" cellspacing="0">
//     * <TR><TD WIDTH=50></TD> <TD><b>A</b></TD> <TD WIDTH=10></TD> <TD> - a Base::Matrix<EntryType></TD></TR>
//     * <TR><TD WIDTH=50></TD> <TD><b>B</b></TD> <TD WIDTH=10></TD> <TD> - a Base::Matrix<EntryType></TD></TR>
//     * </TABLE>
//     * <H3>Output Parameter(s):</H3>
//     * <TABLE border="0" cellpadding="0" cellspacing="0">
//     * <TR><TD WIDTH=50></TD><TD>returns an instance of Base::Matrix<EntryType> contains result of <b>A - B</b></TD></TR>
//     * </TABLE>
//     * <H3>Example(s):</H3>
//     * <PRE>
//     * \#include <Base/Matrix.hh> 
//     * 
//     * Base::Matrix<double> P(5,3);  // matrix of dimension 5-rows by 3-columns
//     * Base::Matrix<double> Q(5,3);  // matrix of dimension 5-rows by 3-columns
//     * Base::Matrix<double> R(5,3);  // matrix of dimension 5-rows by 3-columns
//     * P = 4.25;   // sets all entries of P to 4.25
//     * Q = -0.75;  // sets all entries of Q to -0.75
//     * R = P - Q;  // assigns R by the result of substraction  P - Q
//     * </PRE>
//     **************************************************************************/
//    template <class EntryType>
//    Matrix<EntryType> operator-(const Matrix<EntryType>& A, const Matrix<EntryType>& B) {
//	//TEST_ERROR_DEBUG((A.dim1() == B.dim1()) && (A.dim2() == B.dim2()), 
//	//	const char*, "Matrix: both matrices must have the same size!");
//		
//		/// \bug Range checking not yet implemented
//	Matrix<EntryType> result(A);
//	result -= B;
//	return result;
//    }
//
//
//    /*! \brief Substraction operation of scallar and Matrix
//     *
//     *  Performs substraction of a scallar by a Base::Matrix.  An instance of Base::Matrix is returned.
//     * 
//     * <H3>Input Parameter(s):</H3>
//     * <TABLE border="0" cellpadding="0" cellspacing="0">
//     * <TR><TD WIDTH=50></TD> <TD><b>c</b></TD> <TD WIDTH=10></TD> <TD> - a scallar</TD></TR>
//     * <TR><TD WIDTH=50></TD> <TD><b>A</b></TD> <TD WIDTH=10></TD> <TD> - a Base::Matrix<EntryType></TD></TR>
//     * </TABLE>
//     * <H3>Output Parameter(s):</H3>
//     * <TABLE border="0" cellpadding="0" cellspacing="0">
//     * <TR><TD WIDTH=50></TD><TD>returns an instance of Base::Matrix<EntryType> contains result of <b>c - A</b></TD></TR>
//     * </TABLE>
//     * <H3>Example(s):</H3>
//     * <PRE>
//     * \#include <Base/Matrix.hh> 
//     * 
//     * Base::Matrix<double> P(5,3);  // matrix of dimension 5-rows by 3-columns
//     * Base::Matrix<double> Q(5,3);  // matrix of dimension 5-rows by 3-columns
//     * double c = 1.59;  // a scallar
//     * P = 4.25;   // sets all entries of P to 4.25
//     * Q = c - P;  // assigns Q by the result of substraction  c - P
//     * </PRE>
//     **************************************************************************/
//    template <class EntryType>
//    Matrix<EntryType> operator-(const EntryType c, const Matrix<EntryType>& A) {
//	Matrix<EntryType> result(A);
//	result.valarray<EntryType>::operator*=(-1);
//	result += c;
//	return result;
//    }
//
//    /*! \brief Substraction operation of Matrix and scallar
//     *
//     *  Performs substraction of a Base::Matrix by a scallar.  An instance of Base::Matrix is returned.
//     * 
//     * <H3>Input Parameter(s):</H3>
//     * <TABLE border="0" cellpadding="0" cellspacing="0">
//     * <TR><TD WIDTH=50></TD> <TD><b>A</b></TD> <TD WIDTH=10></TD> <TD> - a Base::Matrix<EntryType></TD></TR>
//     * <TR><TD WIDTH=50></TD> <TD><b>c</b></TD> <TD WIDTH=10></TD> <TD> - a scallar</TD></TR>
//     * </TABLE>
//     * <H3>Output Parameter(s):</H3>
//     * <TABLE border="0" cellpadding="0" cellspacing="0">
//     * <TR><TD WIDTH=50></TD><TD>returns an instance of Base::Matrix<EntryType> contains result of <b>A - c</b></TD></TR>
//     * </TABLE>
//     * <H3>Example(s):</H3>
//     * <PRE>matrix_prod
//     * \#include <Base/Matrix.hh> 
//     * 
//     * Base::Matrix<double> P(5,3);  // matrix of dimension 5-rows by 3-columns
//     * Base::Matrix<double> Q(5,3);  // matrix of dimension 5-rows by 3-columns
//     * double c = 1.59;  // a scallar
//     * P = 4.25;   // sets all entries of P to 4.25
//     * Q = P - c;  // assigns Q by the result of substraction  P - c
//     * </PRE>
//     **************************************************************************/
//    template <class EntryType>
//    Matrix<EntryType> operator-(const Matrix<EntryType>& A, const EntryType c) {
//	Matrix<EntryType> result(A);
//	result -= c;
//	return result;
//    }
//
//    /*! \brief Multiplication operation of two Matrices
//     *
//     *  Performs multiplication of two Base::Matrix.  The size of matrices must be appropriate, i.e.
//     *  C_{p x r} = A_{p x q} * B_{q x r}
//     *
//     * <H3>Input Parameter(s):</H3>
//     * <TABLE border="0" cellpadding="0" cellspacing="0">
//     * <TR><TD WIDTH=50></TD> <TD><b>A</b></TD> <TD WIDTH=10></TD> <TD> - a Base::Matrix<EntryType></TD></TR>
//     * <TR><TD WIDTH=50></TD> <TD><b>B</b></TD> <TD WIDTH=10></TD> <TD> - a Base::Matrix<EntryType></TD></TR>
//     * </TABLE>
//     * <H3>Output Parameter(s):</H3>
//     * <TABLE border="0" cellpadding="0" cellspacing="0">
//     * <TR><TD WIDTH=50></TD> <TD><b>C</b></TD> <TD WIDTH=10></TD> <TD> - a Base::Matrix<EntryType> contains result of <b>A * B</b></TD></TR>
//     * </TABLE>
//     * <H3>Example(s):</H3>
//     * <PRE>
//     * \#include <Base/Matrix.hh> 
//     * 
//     * Base::Matrix<double> P(5,3);  // matrix of dimension 5-rows by 3-columns
//     * Base::Matrix<double> Q(3,4);  // matrix of dimension 3-rows by 4-columns
//     * Base::Matrix<double> R(5,4);  // matrix of dimension 5-rows by 4-columns
//     * P = 4.25;   // sets all entries of P to 4.25
//     * Q = 2.0;    // sets all entries of Q to 2.0
//     * matrix_prod(P, Q, R);  // assigns R by the result of multiplication P * Q
//     * </PRE>
//     **************************************************************************/
//    template <typename EntryType>
//    void matrix_prod(const Matrix<EntryType>& A, const Matrix<EntryType>& B, Matrix<EntryType>& C) {
////	TEST_ERROR_DEBUG(((A.dim1() == C.dim1()) && (A.dim2() == B.dim1()) && (B.dim2() == C.dim2())), 
////		const char*, "Matrix: incompatible matrix size!");
//		
//		/// \bug Range checking not implemented
//		
//	int nrows_=A.dim1(), ncols_=A.dim2(), d3=B.dim2();
//	std::slice slc_i(0,nrows_,1);
//	for (int j=0; j<d3;++j) {
//	   std::slice slc_j(j*nrows_,nrows_,1);
//	   C[slc_j] = (B[ncols_*j] * A[slc_i]);
//	}
//	for (int i=1; i<ncols_;++i) {
//	   std::slice slc_i(i*nrows_,nrows_,1);
//	   for (int j=0; j<d3;++j) {
//	      std::slice slc_j(j*nrows_,nrows_,1);
//	      C[slc_j] += (B[i+ncols_*j] * A[slc_i]);
//	   }
//	}
//    }
//	
//
//    /*! \brief Solve linear system A*X = B
//     *
//     *  Solve system of linear equation \b A*X \b = \b B using lapack's ?gesv 
//     *  routines, which use Gauss elimination method with partial pivoting.  
//     *  On return, the coeficient matrix A is replaced by its LU factor: 
//     *  lower tridiagonal L and upper tridiagonal U.  Note that diagonal elements 
//     *  of L are not stored.
//     *
//     *  This function is overloaded for \c float and \b double \b EntryType .
//     *
//     * <H3>Input Parameter(s):</H3>
//     * <TABLE border="0" cellpadding="0" cellspacing="0">
//     * <TR><TD WIDTH=50></TD> <TD><b>A</b></TD> <TD WIDTH=10></TD> <TD> - coeficient matrix<EntryType></TD></TR>
//     * <TR><TD WIDTH=50></TD> <TD><b>B</b></TD> <TD WIDTH=10></TD> <TD> - right hand side(s)<EntryType></TD></TR>
//     * </TABLE>
//     * <H3>Output Parameter(s):</H3>
//     * <TABLE border="0" cellpadding="0" cellspacing="0">
//     * <TR><TD WIDTH=50></TD> <TD><b>A</b></TD> <TD WIDTH=10></TD> <TD> - L and U factors of the coeficient matrix<EntryType></TD></TR>
//     * <TR><TD WIDTH=50></TD> <TD><b>B</b></TD> <TD WIDTH=10></TD> <TD> - solution of the system<EntryType></TD></TR>
//     * <TR><TD WIDTH=50></TD> <TD>returns flag of computation status: <br>
//     *                          0  ==> succeeded, <br>
//     *                          <0 ==> negative \b i means illegal value on i-th argument, <br>
//     *                          >0 ==> positive \b i means U(i,i) is zero, i.e. it is a singular system.</TD></TR>
//     * </TABLE>
//     * <H3>Example(s):</H3>
//     * <PRE>
//     * \#include <Base/Matrix.hh> 
//     * 
//     * Base::Matrix<double> A(5,5);  // coeficient matrix of the linear system
//     * Base::Matrix<double> B(5,3);  // right hand sides
//     * // then fill in the entries of matrices A & B
//     *
//     * // Now, solve the system of linear equations
//     * int info = Base::gesv(A, B);
//     * if (info == 0) {
//     *    cout << "LU matrix = " << A << endl;
//     *    cout << "Solution = " << B << endl;
//     * }
//     * else cout << "Failed to solve the system\n"
//     * </PRE>
//     **************************************************************************/
//    template<typename EntryType>
//    inline int gesv(Matrix<EntryType>&, Matrix<EntryType>&);
//
//    //! Base::gesv for double precission matrix entries
//    template<>
//    inline int gesv<double>(Matrix<double>& A, Matrix<double>& B){
//	int ipiv[A.nrows_], info;
//	dgesv_ (&A.nrows_, &B.ncols_, 
//		&static_cast<valarray<double>&>(A)[0], &A.nrows_, ipiv, 
//		&static_cast<valarray<double>&>(B)[0], &B.nrows_, &info);
//	return info;
//    }
//
//    //! Base::gesv for single precission matrix entries
//    template<>
//    inline int gesv<float>(Matrix<float>& A, Matrix<float>& B){
//	int ipiv[A.nrows_], info;
//	sgesv_ (&A.nrows_, &B.ncols_, 
//		&static_cast<valarray<float>&>(A)[0], &A.nrows_, ipiv, 
//		&static_cast<valarray<float>&>(B)[0], &B.nrows_, &info);
//	return info;
//    }
//
//    /*! \brief Prints content of the Matrix
//     *
//     *  Prints matrix contents to an ostream.
//     * 
//     * <H3>Input Parameter(s):</H3>
//     * <TABLE border="0" cellpadding="0" cellspacing="0">
//     * <TR><TD WIDTH=50></TD> <TD><b>A</b></TD> <TD WIDTH=10></TD> <TD> - a Base::Matrix<EntryType></TD></TR>
//     * </TABLE>
//     * <H3>Output Parameter(s):</H3>
//     * <TABLE border="0" cellpadding="0" cellspacing="0">
//     * <TR><TD WIDTH=50></TD> <TD>returns reference to an ostream instance</TD></TR>
//     * </TABLE>
//     * <H3>Example(s):</H3>
//     * <PRE>
//     * \#include <Base/Matrix.hh> 
//     * 
//     * Base::Matrix<double> P(5,3);  // matrix of dimension 5-rows by 3-columns
//     * P = 3.87;  // assigns P by a scallar
//     * std::cout << P << std::endl;  // display matrix P on stdout
//     * </PRE>
//     **************************************************************************/
//    template <class EntryType>
//    ostream& operator<<(ostream& os, const Matrix<EntryType>& A)
//    {
//	os << "[" << A.size1() << "," << A.size2() << "](";
//	for (int i=0; i<A.dim1()-1; ++i) {
//	   os << "(";
//	   for(int j=0; j<A.dim2()-1; ++j)
//		os<< A(i,j) <<",";
//	   os << A(i, A.dim2()-1) << "),";
//	}
//
//	os << "(";
//	for(int j=0; j<A.dim2()-1; ++j)
//	   os<< A(A.dim1()-1,j) <<",";
//	os << A(A.dim1()-1, A.dim2()-1) << ")";
//
//	os << ")";
//	return os;
//    }
//    
//    
//    
//
//}
//#endif
////------------------------------------------------------------------------------
//// Local variables:
//// mode:c++
//// comment-column: 48
//// End:
