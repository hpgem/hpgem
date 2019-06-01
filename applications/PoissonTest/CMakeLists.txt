#First declare the rules that will form the executable tests#
#############################################################

if(hpGEM_USE_PETSC)

    
add_executable(PoissonTest.out 
		main.cpp
      		)
target_link_libraries(PoissonTest.out HPGEM::HPGEM PETSc::PETSc)

else() 
	add_custom_target(PoissonTest.out 
		COMMENT "\t ERROR:: You need PETSC to compile PoissonTest.out. \n\t ERROR:: Please rerun cmake and enable PETSC support")
endif()
