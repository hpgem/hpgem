include_directories(. ${hpGEM_SOURCE_DIR}/kernel/)

#First declare the rules that will form the executable tests#
#############################################################
   
add_executable(Euler.out 
		Euler.cpp
		main.cpp
      		)
target_link_libraries(Euler.out hpGEM_Base Output)

