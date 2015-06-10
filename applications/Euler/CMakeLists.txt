include_directories(. ${hpGEM_SOURCE_DIR}/kernel/)

#First declare the rules that will form the executable tests#
#############################################################
 if(NOT hpGEM_USE_COMPLEX_PETSC)
     
add_executable(Euler.out 
		main.cpp
      		)
target_link_libraries(Euler.out hpGEM_Base Output)
else()
    add_custom_target(Euler.out 
        COMMENT "\t ERROR:: You need real numbers to compile Euler.out. \n\t ERROR:: Please rerun cmake and disable complex PETSC support")
endif()

