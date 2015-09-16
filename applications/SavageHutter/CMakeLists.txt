include_directories(. ${hpGEM_SOURCE_DIR}/kernel/)

if(NOT hpGEM_USE_COMPLEX_PETSC)

add_subdirectory(SlopeLimiters)
add_subdirectory(HeightLimiters)

add_executable(SavageHutter.out
    HelperFunctions.cpp
    SavageHutterRightHandSideComputer.cpp
    BidisperseRHS1D.cpp
    SavageHutterRHS2D.cpp
		SavageHutter.cpp
    SavageHutterBase.cpp
		main.cpp
      		)
target_link_libraries(SavageHutter.out hpGEM_Base Output SlopeLimiters HeightLimiters)

endif()