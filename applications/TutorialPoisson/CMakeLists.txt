
#First declare the rules that will form the executable tests#
#############################################################
if(hpGEM_USE_PETSC)

add_executable(TutorialPoisson.out 
		main.cpp
      		)
target_link_libraries(TutorialPoisson.out HPGEM::HPGEM PETSc::PETSc)

else() 
	add_custom_target(TutorialPoisson.out 
		COMMENT "\t ERROR:: You need PETSC to compile TutorialPoisson.out. \n\t ERROR:: Please rerun cmake and enable PETSC support")
endif()
