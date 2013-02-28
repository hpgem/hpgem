include_directories(${hpGEM_SOURCE_DIR}/src/ .)


#First declare the rules that will form the executable tests#
#############################################################

#Define MatrixTest.out dependencies
add_executable(MatrixTest.out MatrixTest.cpp)
target_link_libraries(MatrixTest.out Matrix)

#add_executable(MatrixSpeedTest.out
#	${hpGEM_SOURCE_DIR}/src/LinearAlgebra/daxpy.f
#	${hpGEM_SOURCE_DIR}/src/LinearAlgebra/lsame.f
#	${hpGEM_SOURCE_DIR}/src/LinearAlgebra/xerbla.f
#	${hpGEM_SOURCE_DIR}/src/LinearAlgebra/dgemm.f MatrixSpeedTest.cpp)
#target_link_libraries(MatrixSpeedTest.out gfortran)

add_executable(VectorTest.out 
	VectorTest.cpp 
	${hpGEM_SOURCE_DIR}/src/Geometry/Point.hpp 
	${hpGEM_SOURCE_DIR}/src/LinearAlgebra/daxpy.f
	${hpGEM_SOURCE_DIR}/src/LinearAlgebra/lsame.f 
	${hpGEM_SOURCE_DIR}/src/LinearAlgebra/xerbla.f 
	${hpGEM_SOURCE_DIR}/src/LinearAlgebra/dgemv.f 
	${hpGEM_SOURCE_DIR}/src/LinearAlgebra/NumericalVector.hpp)
target_link_libraries(VectorTest.out gfortran)

#Now create the actually tests#
###############################

#Test 1 - Simple Matrix Algebra Test : Check it both complies and runs
add_test(BuildTest "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target MatrixTest.out)
add_test(MatrixTest MatrixTest.out)
set_tests_properties(MatrixTest PROPERTIES DEPENDS BuildTest)

#Test 2 - To be cleaned up
add_test(VectorTest VectorTest.out)  
