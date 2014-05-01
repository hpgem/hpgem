include_directories(. ${hpGEM_SOURCE_DIR}/kernel/)

#First declare the rules that will form the executable tests#
#############################################################

add_executable(app1.out 
		app1.cpp
      		)
target_link_libraries(app1.out hpGEM_Base QuadratureRules)
