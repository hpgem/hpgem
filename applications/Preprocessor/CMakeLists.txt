
add_library(Preprocessor-lib STATIC
lib/centaur.cpp
lib/hpgem.cpp
lib/gmsh.cpp
	)
	# TODO: Lacking source/include separation
	target_include_directories(Preprocessor-lib
		PUBLIC
			${CMAKE_CURRENT_SOURCE_DIR}/lib
	)
	target_link_libraries(Preprocessor-lib
		PUBLIC
            HPGEM::HPGEM)
            

if(hpGEM_USE_MPI)
    target_link_libraries(Preprocessor-lib PUBLIC MPI::MPI_CXX)
endif()

if(hpGEM_USE_METIS)
    target_link_libraries(Preprocessor-lib PUBLIC METIS::METIS)
endif()

if(hpGEM_USE_ANY_PETSC)
    target_link_libraries(Preprocessor-lib PUBLIC PETSc::PETSc)
endif()

if(hpGEM_USE_SLEPC)
    target_link_libraries(Preprocessor-lib PUBLIC SLEPc::SLEPc)
endif()

#the preprocessor does not perform floating poing arithmetic, so it doesn't need special Petsc/not petsc guards
add_executable(Preprocessor.out
               main.cpp
               ../../kernel/Logger.cpp
               ../../kernel/Base/CommandLineOptions.cpp
               ../../kernel/Base/MpiContainer.cpp)
               target_link_libraries(Preprocessor.out Preprocessor-lib)
add_executable(ExtractSurface.out
               extractSurface.cpp
               ../../kernel/Logger.cpp
               ../../kernel/Base/CommandLineOptions.cpp
               ../../kernel/Base/MpiContainer.cpp)

               target_link_libraries(ExtractSurface.out Preprocessor-lib)
