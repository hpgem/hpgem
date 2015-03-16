include_directories(. ${hpGEM_SOURCE_DIR}/kernel/)

#First declare the rules that will form the executable tests#
#############################################################
   
add_executable(TutorialAdvection.out 
		main.cpp
      		)
target_link_libraries(TutorialAdvection.out hpGEM_Base Output)

