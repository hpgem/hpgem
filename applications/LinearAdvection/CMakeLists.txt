include_directories(. ${hpGEM_SOURCE_DIR}/kernel/)

#First declare the rules that will form the executable tests#
#############################################################
   
add_executable(LinearAdvection.out 
		main.cpp
      		)
target_link_libraries(LinearAdvection.out hpGEM_Base Output ${PETSC_LIBRARIES})

