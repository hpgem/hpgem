
#First declare the rules that will form the executable tests#
#############################################################
    if(NOT hpGEM_USE_COMPLEX_PETSC)
    
add_executable(ExampleMultipleVariableProblem.out
		main.cpp
      		)
target_link_libraries(ExampleMultipleVariableProblem.out HPGEM::HPGEM)
else()
    add_custom_target(ExampleMultipleVariableProblem.out 
        COMMENT "\t ERROR:: You need real numbers to compile ExampleMultipleVariableProblem.out. \n\t ERROR:: Please rerun cmake and disable complex PETSC support")
endif()

