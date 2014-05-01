include_directories(. ${hpGEM_SOURCE_DIR}/kernel/)

#First declare the rules that will form the executable tests#
#############################################################


if(hpGEM_INCLUDE_PETSC_SUPPORT)

	configure_file(bCoeff.txt
			${PROJECT_BINARY_DIR}/bCoeff.txt @ONLY IMMEDIATE)
    include_directories(${PETSC_INCLUDES})

	add_executable(hEulerApp.out
        main.cpp
        HEuler.cpp
        InitialConditions.cpp
        )
	
		target_link_libraries(hEulerApp.out Output ${PETSC_LIBRARIES})
else() 
	add_custom_target(hEulerApp.out 
		COMMENT "\t ERROR:: You need PETSC to compile hEuler.out. \n\t ERROR:: Please rerun cmake and enable PETSC support \n\t ERROR:: If using raw cmake add the flag -DhpGEM_INCLUDE_PETSC_SUPPORT='on'")
endif()

