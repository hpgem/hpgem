include_directories(. ${hpGEM_SOURCE_DIR}/kernel/)

#First declare the rules that will form the executable tests#
#############################################################
if(NOT hpGEM_USE_COMPLEX_PETSC)
    
    add_executable(CompressibleNavierStokes.out
		Inviscid.cpp
		Viscous.cpp 
		CompressibleNavierStokes.cpp
		main.cpp
      		)
    target_link_libraries(CompressibleNavierStokes.out hpGEM_Base Output)

endif()