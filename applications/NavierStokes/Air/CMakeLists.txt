include_directories(. ${hpGEM_SOURCE_DIR}/kernel/)

#First declare the rules that will form the executable tests#
#############################################################
if(NOT hpGEM_USE_COMPLEX_PETSC)
    
    add_executable(Air.out
		StateCoefficientsStructAir.cpp
		Air.cpp
		main.cpp
      		)
    target_link_libraries(Air.out hpGEM_Base Output)

endif()
