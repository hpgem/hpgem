include_directories(. ${hpGEM_SOURCE_DIR}/kernel/)

#First declare the rules that will form the executable tests#
#############################################################
   
add_executable(CompressibleNavierStokes.out
		Inviscid.cpp
		Viscous.cpp 
		CompressibleNavierStokes.cpp
		main.cpp
      		)
target_link_libraries(CompressibleNavierStokes.out hpGEM_Base Output)

