include_directories(. ${hpGEM_SOURCE_DIR}/kernel/)

#First declare the rules that will form the executable tests#
#############################################################
if(hpGEM_USE_PETSC)

    include_directories(${PETSC_INCLUDES})
    
add_executable(TutorialLaplace.out 
		main.cpp
      		)
target_link_libraries(TutorialLaplace.out hpGEM_Base Output ${PETSC_LIBRARIES})

else() 
	add_custom_target(TutorialLaplace.out 
		COMMENT "\t ERROR:: You need PETSC to compile LaplaceEasier.out. \n\t ERROR:: Please rerun cmake and enable PETSC support")
endif()
