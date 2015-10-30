#include_directories(. ${hpGEM_SOURCE_DIR}/kernel/)

#When hpGEM is configured to interact with a complex build of PETSc some of its internal data
#structures will also use complex numbers. This is why applications that don't need PETSc
#are disabled for this case. If you need to interact with complex PETSc, make sure the code
#compiles both with PETSc turned on and PETSc turned off.
#if(NOT hpGEM_USE_COMPLEX_PETSC)

#add_subdirectory(SlopeLimiters)
#add_subdirectory(HeightLimiters)

#add_executable(SavageHutter.out
#    HelperFunctions.cpp
#    SavageHutter1DBasic.cpp
#    SavageHutter1DBidispersed.cpp
#    SavageHutter1DWidthAveraged.cpp
#    SavageHutter1DWidthHAndU.cpp
#    SavageHutter2DBasic.cpp
#    SavageHutter1DBase.cpp
#    SavageHutter2DBase.cpp
#		main.cpp
#      		)
#target_link_libraries(SavageHutter.out hpGEM_Base Output SlopeLimiters HeightLimiters)

#endif()