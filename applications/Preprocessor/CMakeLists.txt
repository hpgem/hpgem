include_directories(. ${hpGEM_SOURCE_DIR}/kernel/)

if(hpGEM_USE_PETSC OR hpGEM_USE_COMPLEX_PETSC)
    include_directories(${PETSC_INCLUDES})
endif()

if(hpGEM_USE_SLEPC)
    include_directories(${SLEPC_INCLUDE_DIRS})
endif()

#the preprocessor does not perform floating poing arithmetic, so it doesn't need special Petsc/not petsc guards
add_executable(Preprocessor.out
               main.cpp
               centaur.cpp
               hpgem.cpp
               ../../kernel/Logger.cpp
               ../../kernel/Base/CommandLineOptions.cpp
               ../../kernel/Base/MpiContainer.cpp)

set_property(TARGET Preprocessor.out PROPERTY CXX_STANDARD 14)
set(CURRENT_TARGET Preprocessor.out)
include(${CMAKE_SOURCE_DIR}/conf/cmake/CLionFix.cmake)
set(CURRENT_TARGET)

target_link_libraries(Preprocessor.out ${METIS_LIBRARIES})

if(hpGEM_USE_MPI)
    set( DEP_COMPILE_FLAGS "${DEP_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}" )
    set( DEP_LINKER_FLAGS "${DEP_LINKER_FLAGS} ${MPI_LINK_FLAGS}")
    target_link_libraries(Preprocessor.out "${MPI_LIBRARIES}")
endif()

if(hpGEM_USE_METIS)
    target_link_libraries(Preprocessor.out "${METIS_LIBRARIES}")
endif()

if(hpGEM_USE_PETSC OR hpGEM_USE_COMPLEX_PETSC)
    target_link_libraries(Preprocessor.out "${PETSC_LIBRARIES}")
endif()

if(hpGEM_USE_SLEPC)
    target_link_libraries(Preprocessor.out "${SLEPC_LIBRARIES}")
endif()
