
#First declare the rules that will form the executable tests#
#############################################################

if(hpGEM_USE_PETSC)
add_executable(LinearPotentialFlow.out 
		main.cpp
      		)
target_link_libraries(LinearPotentialFlow.out HPGEM::HPGEM PETSc::PETSc)

else() 
	add_custom_target(LinearPotentialFlow.out 
		COMMENT "\t ERROR:: You need PETSC to compile Laplace.out. \n\t ERROR:: Please rerun cmake and enable PETSC support")
endif()
