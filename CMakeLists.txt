cmake_minimum_required (VERSION 2.8)
project (hpGEM CXX Fortran)

enable_testing()

add_subdirectory(src)
add_subdirectory(tests)


#This section generates a configuration file that contains the CMake varibles#
#############################################

configure_file(conf/CMakeDefinitions.conf
		${PROJECT_BINARY_DIR}/conf/CMakeDefinitions.hpp)


#Extra gfortran checks
#######################


	########
		#Check for propery gfortran lib

		#find_library(GFORTRAN_LIBRARY gfortran)
		#
		#find_program(_GFORTRAN_EXECUTABLES NAMES gfortran)
		#
		#if (_GFORTRAN_EXECUTABLE)
		#EXECUTE_PROCESS(COMMAND ${_GFORTRAN_EXECUTABLE} -print-file-name=libgfortran.so
		#OUTPUT_VARIABLE _libgfortran_path
		#OUTPUT_STRIP_TRAILING_WHITESPACE																				)
		#endif()

		#if(EXISTS ${_libgfortran_path})
		#	set(GFORTRAN_LIBRARY ${_libgfortran_path})
			#else()
			#find_library(GFORTRAN_LIBRARY gfortran)
			#endif()

			#if (NOT GFORTRAN_LIBRARY)
			#message(FATAL_ERROR "gfortran is need but some of library could not be found, please install correctly")
			#endif()



#This section is about the build of documentation#
##############################################
option(hpGEM_BUILD_DOCUMENTATION "Use Doxygen to create the HTML based API documentation" OFF)
option(hpGEM_DOCUMENTATION_INCLUDE_UML_DIAGRAMS "Generate UML diagrams in the doxygen output, require DOT to be installed",OFF)
if(hpGEM_BUILD_DOCUMENTATION)
    FIND_PACKAGE(Doxygen)
    if (NOT DOXYGEN_FOUND)
        message(FATAL_ERROR 
            "Doxygen is needed to build the documentation. Please install it correctly or turn the option")
    endif()

  if (DOXYGEN_FOUND)

	  #	  FIND_PACKAGE(GraphViz)

	  # THis package does not exist but there is help on the web to find it  
	#if (GRAPHVIZ_FOUND)
	#	message(FATAL_ERROR
	#		"Graph Viz is required to generate diagrams for the documentation")
	#	endif()
  
    #Copy the doxygen configuration file in the build directory
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
