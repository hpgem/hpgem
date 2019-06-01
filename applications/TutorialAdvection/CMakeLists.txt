
#First declare the rules that will form the executable tests#
#############################################################
   
add_executable(TutorialAdvection.out 
		main.cpp
      		)
target_link_libraries(TutorialAdvection.out HPGEM::HPGEM)

