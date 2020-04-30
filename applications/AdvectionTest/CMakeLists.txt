

#First declare the rules that will form the executable tests#
#############################################################
if(NOT hpGEM_USE_COMPLEX_PETSC)
   
add_executable(AdvectionTest.out 
		main.cpp
      		)
target_link_libraries(AdvectionTest.out HPGEM::HPGEM)
else()
    add_custom_target(AdvectionTest.out 
        COMMENT "\t ERROR:: You need real numbers to compile AdvectionTest.out. \n\t ERROR:: Please rerun cmake and disable complex PETSC support")
endif()

