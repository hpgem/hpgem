include_directories(. ${hpGEM_SOURCE_DIR}/kernel/)

#First declare the rules that will form the executable tests#
#############################################################
    if(NOT hpGEM_USE_COMPLEX_PETSC)
    
add_executable(SavageHutter.out
    HelperFunctions.cpp
    TvbLimiterWithDetector1D.cpp
    PositiveLayerLimiter.cpp
    SavageHutterRightHandSideComputer.cpp
		SavageHutter.cpp
		main.cpp
      		)
target_link_libraries(SavageHutter.out hpGEM_Base Output)
else()
    add_custom_target(SavageHutter.out 
        COMMENT "\t ERROR:: You need real numbers to compile SavageHutter.out. \n\t ERROR:: Please rerun cmake and disable complex PETSC support")
endif()

