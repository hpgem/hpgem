include_directories(. ${hpGEM_SOURCE_DIR}/kernel/)

#First declare the rules that will form the executable tests#
#############################################################

if(hpGEM_USE_PETSC)

    include_directories(${PETSC_INCLUDES})
    
add_executable(Laplace.out 
		main.cpp
      		)
target_link_libraries(Laplace.out hpGEM_Base Output ${PETSC_LIBRARIES})

else() 
	add_custom_target(Laplace.out 
		COMMENT "\t ERROR:: You need PETSC to compile Laplace.out. \n\t ERROR:: Please rerun cmake and enable PETSC support")
endif()
