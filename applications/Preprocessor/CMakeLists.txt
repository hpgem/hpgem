include_directories(. ${hpGEM_SOURCE_DIR}/kernel/)

#the preprocessor does not perform floating poing arithmetic, so it doesn't need special Petsc/not petsc guards
add_executable(Preprocessor.out
               main.cpp
               centaur.cpp
               FortranUnformattedFile.cpp
               hpgem.cpp
        SparseUnionFind.cpp
               ../../kernel/Logger.cpp
               ../../kernel/Base/CommandLineOptions.cpp
               ../../kernel/Base/MpiContainer.cpp)

add_executable(ExtractSurface.out
               extractSurface.cpp
               centaur.cpp
               FortranUnformattedFile.cpp
               hpgem.cpp
        SparseUnionFind.cpp
               ../../kernel/Logger.cpp
               ../../kernel/Base/CommandLineOptions.cpp
               ../../kernel/Base/MpiContainer.cpp)

set(CURRENT_TARGET Preprocessor.out)
include(${CMAKE_SOURCE_DIR}/conf/cmake/CLionFix.cmake)
set(CURRENT_TARGET ExtractSurface.out)
include(${CMAKE_SOURCE_DIR}/conf/cmake/CLionFix.cmake)
set(CURRENT_TARGET)

if(hpGEM_USE_MPI)
    target_link_libraries(Preprocessor.out MPI::MPI_CXX)
    target_link_libraries(ExtractSurface.out MPI::MPI_CXX)
endif()

if(hpGEM_USE_METIS)
    target_link_libraries(Preprocessor.out METIS::METIS)
    target_link_libraries(ExtractSurface.out METIS::METIS)
endif()

if(hpGEM_USE_ANY_PETSC)
    target_link_libraries(Preprocessor.out PETSc::PETSc)
    target_link_libraries(ExtractSurface.out PETSc::PETSc)
endif()

if(hpGEM_USE_SLEPC)
    target_link_libraries(Preprocessor.out SLEPc::SLEPc)
    target_link_libraries(ExtractSurface.out SLEPc::SLEPc)
endif()
