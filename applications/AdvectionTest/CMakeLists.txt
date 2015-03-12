include_directories(. ${hpGEM_SOURCE_DIR}/kernel/)

#First declare the rules that will form the executable tests#
#############################################################
   
add_executable(AdvectionTest.out 
		main.cpp
      		)
target_link_libraries(AdvectionTest.out hpGEM_Base Output ${PETSC_LIBRARIES})

