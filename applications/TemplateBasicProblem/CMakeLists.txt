include_directories(. ${hpGEM_SOURCE_DIR}/kernel/)

#First declare the rules that will form the executable tests#
#############################################################
if(NOT hpGEM_USE_COMPLEX_PETSC)
   
#change the name of the |
#executable here        V
    add_executable(ExampleProblem.out 
    		main.cpp
          		)
#and here                    |
#                            V
    target_link_libraries(ExampleProblem.out hpGEM_Base Output)
endif()
