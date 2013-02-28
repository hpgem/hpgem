add_library(Matrix SHARED 
		daxpy.f lsame.f xerbla.f dgemm.f dgemv.f Matrix.cpp)
#find_library(gfortran)
target_link_libraries(Matrix)
#SET_TARGET_PROPERTIES(Matrix PROPERTIES LINKER_LANGUAGE CXX)
