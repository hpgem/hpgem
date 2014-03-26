cmake_minimum_required (VERSION 2.8)
project (hpGEM CXX Fortran)

set(hpGEM_VERSION_MAJOR 2)
set(hpGEM_VERSION_MINOR 0)

#Directory were the external cmake modules are stored
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/conf/cmake/Modules/")


enable_testing()

add_subdirectory(kernel)
add_subdirectory(tests)


#This section generates a configuration file that contains the CMake varibles#
#############################################

configure_file(conf/CMakeDefinitions.conf
		${PROJECT_BINARY_DIR}/conf/CMakeDefinitions.hpp)


#This section is about the build of documentation#
##############################################

option(hpGEM_BUILD_DOCUMENTATION "Use Doxygen to create the HTML based API documentation" OFF)
option(hpGEM_INCLUDE_PETSC_SUPPORT "Include PETSC this is required for GlobalAssemble" OFF)
option(hpGEM_INCLUDE_SLEPC_SUPPORT "Include SLEPC this is needed for some applications" OFF)
if(hpGEM_BUILD_DOCUMENTATION)
    FIND_PACKAGE(Doxygen)
    if (NOT DOXYGEN_FOUND)
        message(FATAL_ERROR 
            "Doxygen is needed to build the documentation. Please install it correctly or turn the option")
    endif()

  if (DOXYGEN_FOUND)
    
    #This is the configure file for normal doxygen builds	
    configure_file(conf/doxygen.conf 
        	 	${PROJECT_BINARY_DIR}/conf/doxygen.conf  @ONLY IMMEDIATE)

    #The next four and for website doxygen builds. The should be hinded from the public cmake at some point
    configure_file(conf/web_doxygen.conf
	    		${PROJECT_BINARY_DIR}/conf/web_doxygen.conf @ONLY IMMEDIATE)

    configure_file(conf/hpg.css
	    		${PROJECT_BINARY_DIR}/conf/hpg.css @ONLY IMMEDIATE)

    configure_file(conf/new_footer.html
	    		${PROJECT_BINARY_DIR}/conf/new_footer.html @ONLY IMMEDIATE)

    configure_file(conf/new_header.html
	    		${PROJECT_BINARY_DIR}/conf/new_header.html @ONLY IMMEDIATE)





    #-- Add custom targets to both make and clean (delete) the documentation
        add_custom_target 	(doc
                            COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/conf/doxygen.conf
                            SOURCES ${PROJECT_BINARY_DIR}/conf/doxygen.conf)

	add_custom_target	(docClean
                        COMMAND rm -r ${PROJECT_BINARY_DIR}/docs/*
			COMMENT "Cleaning (deleting) the documentation"	) 

	add_custom_target	(docWeb
			COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/conf/web_doxygen.conf
			SOURCES ${PROJECT_BINARY_DIR}/conf/web_doxygen.conf)

	add_custom_target	(docPublishAlpha
			COMMAND rsync -vr docs/html/* WebAdmin@einder.ewi.utwente.nl:/home/www/html/hpGEM/assets/doxygen/alpha/
			COMMENT "Publushing (uploading) the curent documentation as ALPHA realese" )
  endif()
endif()

if(hpGEM_INCLUDE_PETSC_SUPPORT)
        enable_language(C)
	FIND_PACKAGE(PETSc REQUIRED)
	#FIND_PACKAGE(SLEPc REQUIRED)
        if (NOT PETSC_FOUND)
		message(FATAL_ERROR
			"The option you have choosen requires PETSc and you do not have this installed. Please install")
	endif()

#?	message(FATAL_ERROR 
#?		"You should not ahve choose this option")

#endif PETSc required
endif()

if(hpGEM_INCLUDE_SLEPC_SUPPORT)
        enable_language(C)
        if(NOT hpGEM_INCLUDE_PETSC_SUPPORT)
               message(FATAL_ERROR
                       "SLEPc depends on PETSc, please also enable PETSc support")
        endif()
        FIND_PACKAGE(SLEPc REQUIRED)
        if(NOT SLEPC_FOUND)
               message(FATAL_ERROR
                       "The option you have chosen requires SLEPc and you do not have this installed. Please install")
        endif()
endif()


