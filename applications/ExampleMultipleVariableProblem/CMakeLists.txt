include_directories(. ${hpGEM_SOURCE_DIR}/kernel/)

#First declare the rules that will form the executable tests#
#############################################################
   
add_executable(ExampleMultipleVariableProblem.out
		main.cpp
      		)
target_link_libraries(ExampleMultipleVariableProblem.out hpGEM_Base Output)

