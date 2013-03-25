add_library(Matrix SHARED 
		daxpy.f lsame.f xerbla.f dgemm.f dgemv.f Matrix.cpp)
