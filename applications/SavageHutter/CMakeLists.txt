include_directories(. ${hpGEM_SOURCE_DIR}/kernel/)

#First declare the rules that will form the executable tests#
#############################################################
add_executable(SavageHutter.out
    HelperFunctions.cpp
    TvbLimiterWithDetector1D.cpp
    PositiveLayerLimiter.cpp
    SavageHutterRightHandSideComputer.cpp
    SavageHutterRHS2D.cpp
		SavageHutter.cpp
    SavageHutterBase.cpp
		main.cpp
      		)
target_link_libraries(SavageHutter.out hpGEM_Base Output)

