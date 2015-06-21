include_directories(. ${hpGEM_SOURCE_DIR}/kernel/)

#First declare the rules that will form the executable tests#
#############################################################
   
add_executable(ExampleProblem.out 
		main.cpp
      		)
target_link_libraries(ExampleProblem.out hpGEM_Base Output)

