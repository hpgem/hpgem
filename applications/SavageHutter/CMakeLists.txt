include_directories(. ${hpGEM_SOURCE_DIR}/kernel/)

#First declare the rules that will form the executable tests#
#############################################################
   
add_executable(SavageHutter.out
    SavageHutterRightHandSideComputer.cpp
		SavageHutter.cpp
		main.cpp
      		)
target_link_libraries(SavageHutter.out hpGEM_Base Output)

