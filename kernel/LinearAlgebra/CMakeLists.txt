include_directories(${hpGEM_SOURCE_DIR}/kernel/ .)
add_library(LinearAlgebra SHARED 
		daxpy.f lsame.f xerbla.f dgemm.f dgemv.f 
		#Next line for LAPACK stuff inverse and LU factorisation
		#First the BLAS required parts
		dger.f dtrsm.f dscal.f dswap.f idamax.f dtrmm.f dtrmv.f
		#Now the actually LAPACK parts
		ieeeck.f iparmq.f dlamch.f dlaswp.f ilaenv.f dgetf2.f dgetrf.f dtrtri.f  dgetri.f dtrti2.f
		
		#These are the extra from LAPACK for solve Ax=B
		dgesv.f dgetrs.f
		NumericalVector.cpp
		Matrix.cpp)
target_link_libraries(LinearAlgebra Logger)