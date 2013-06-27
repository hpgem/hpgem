cmake_minimum_required (VERSION 2.8)
project (hpGEM CXX Fortran)

set(hpGEM_VERSION_MAJOR 2)
set(hpGEM_VERSION_MINOR 0)

#Directory were the external cmake modules are stored
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")


enable_testing()

add_subdirectory(src)
add_subdirectory(tests)


#This section generates a configuration file that contains the CMake varibles#
#############################################

configure_file(conf/CMakeDefinitions.conf
		${PROJECT_BINARY_DIR}/conf/CMakeDefinitions.hpp)


#This section is about the build of documentation#
##############################################

option(hpGEM_BUILD_DOCUMENTATION "Use Doxygen to create the HTML based API documentation" OFF)
option(hpGEM_INCLUDE_PETSC_SUPPORT "Include PETSC this is required for GlobalAssemble" OFF)
if(hpGEM_BUILD_DOCUMENTATION)
    FIND_PACKAGE(Doxygen)
    if (NOT DOXYGEN_FOUND)
        message(FATAL_ERROR 
            "Doxygen is needed to build the documentation. Please install it correctly or turn the option")
    endif()

  if (DOXYGEN_FOUND)

    configure_file(conf/doxygen.conf 
        	 	${PROJECT_BINARY_DIR}/conf/doxygen.conf  @ONLY IMMEDIATE)
    #-- Add custom targets to both make and clean (delete) the documentation
        add_custom_target 	(doc
                            COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/conf/doxygen.conf
                            SOURCES ${PROJECT_BINARY_DIR}/conf/doxygen.conf)

	add_custom_target	(docClean
                        COMMAND rm -r ${PROJECT_BINARY_DIR}/docs/*
		       	COMMENT "Cleaning (deleting) the documentation"	) 
  endif()
endif()

if(hpGEM_INCLUDE_PETSC_SUPPORT)
	FIND_PACKAGE(PETSc REQUIRED)
	if (NOT PETSC_FOUND)
		message(FATAL_ERROR
			"The option you have chosen requires PETSc and you do not have this installed. Please install")
	endif()


#endif PETSc required
endif()
