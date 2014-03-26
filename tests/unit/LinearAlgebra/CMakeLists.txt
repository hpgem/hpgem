include_directories(${hpGEM_SOURCE_DIR}/kernel/ .)


#Test 1 - Simple Matrix Algebra Test : Check it both complies and runs
add_executable(MatrixTest.out MatrixTest.cpp)
target_link_libraries(MatrixTest.out Matrix)
add_test(MatrixTest MatrixTest.out)



#Test 2 - Numerical Vector test
add_executable(VectorTest.out 
	VectorTest.cpp 
	${hpGEM_SOURCE_DIR}/kernel/Geometry/Point.hpp 
	${hpGEM_SOURCE_DIR}/kernel/LinearAlgebra/daxpy.f
	${hpGEM_SOURCE_DIR}/kernel/LinearAlgebra/lsame.f 
	${hpGEM_SOURCE_DIR}/kernel/LinearAlgebra/xerbla.f 
	${hpGEM_SOURCE_DIR}/kernel/LinearAlgebra/dgemv.f 
	${hpGEM_SOURCE_DIR}/kernel/LinearAlgebra/NumericalVector.hpp)
target_link_libraries(VectorTest.out gfortran)
add_test(VectorTest VectorTest.out)  
