#First declare the rules that will form the executable tests#
#############################################################


if(hpGEM_USE_PETSC)

	configure_file(bCoeff.txt
			${PROJECT_BINARY_DIR}/bCoeff.txt @ONLY IMMEDIATE)

	add_executable(hEulerApp.out
        main.cpp
        HEuler.cpp
        InitialConditions.cpp
        )
	
		target_link_libraries(hEulerApp.out HPGEM::HPGEM PETSc::PETSc)
else() 
	add_custom_target(hEulerApp.out 
		COMMENT "\t ERROR:: You need PETSC to compile hEuler.out. \n\t ERROR:: Please rerun cmake and enable PETSC support \n\t ERROR:: If using raw cmake add the flag -DhpGEM_INCLUDE_PETSC_SUPPORT='on'")
endif()

