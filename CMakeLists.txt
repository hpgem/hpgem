cmake_minimum_required (VERSION 2.8)
project (hpGEM C CXX Fortran)

set(hpGEM_VERSION_MAJOR 2)
set(hpGEM_VERSION_MINOR 0)

#Directory were the external cmake modules are stored
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/conf/cmake/Modules/")

#allow dynamically locating dynamic libraries if cmake and the OS support it
if(${CMAKE_VERSION} VERSION_GREATER 2.8.12)
    set(CMAKE_MACOSX_RPATH 1)
endif()

enable_testing()

#Here is the check for CX11 support : We now use some features of this, so turn it on if possible
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)

#Add the flag to compiler options if the current compiler is found
if (NOT CMAKE_CXX_FLAGS)
	if (COMPILER_SUPPORTS_CXX11)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
	elseif(COMPILER_SUPPORTS_CXX0X)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
	else()	
	message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. We require a compiler that has support for CXX0X, for example GCC 4.3 or later")
	endif()
endif(NOT CMAKE_CXX_FLAGS)




#This section generates a configuration file that contains the CMake varibles#
#############################################

configure_file(${hpGEM_SOURCE_DIR}/conf/cmake/CMakeDefinitions.conf ${hpGEM_SOURCE_DIR}/conf/cmake/CMakeDefinitions.h ESCAPE_QUOTES)
include_directories(${hpGEM_SOURCE_DIR}/conf/cmake/)


#Make the default build release is this in what a normal user wants
###################################################################
if(NOT CMAKE_BUILD_TYPE)
	        set(CMAKE_BUILD_TYPE "Release" CACHE 
		STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
        set_property( CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Release Debug RelWithDebInfo MinSizeRel )
endif(NOT CMAKE_BUILD_TYPE)

#Add the custom target make fullTest, that clean and check the whole build system
#################################################################################
add_custom_target(fullTest
	        COMMAND make rebuild_cache      
		COMMAND make clean
		COMMAND make
		COMMAND make test
		COMMENT "This will clean the project, rebuild and test." )


################################################################################
#This section is about the logger, debugging and other stuff
################################################################################

#==========[ Assert Config ]============
set(hpGEM_DISABLE_ASSERTS OFF CACHE BOOL "Disables asserts - use at your own risk!")
mark_as_advanced( FORCE hpGEM_DISABLE_ASSERTS )
if (hpGEM_DISABLE_ASSERTS)
  add_definitions( -DHPGEM_NO_ASSERTS )
endif()

set(hpGEM_ENFORCE_ASSERTS OFF CACHE BOOL "Enables asserts even in release mode - may be slow!")
mark_as_advanced( FORCE hpGEM_ENFORCE_ASSERTS )
if (hpGEM_ENFORCE_ASSERTS)
  if(hpGEM_DISABLE_ASSERTS)
    message(FATAL_ERROR, "Cannot disable the assertions when they are enforced!!!")
  endif()
  add_definitions( -DHPGEM_FORCE_ASSERTS )
endif()

#=======[ The default loglevel ]========
set(hpGEM_LOGLEVEL "DEFAULT" CACHE STRING "Verbosity of hpGEM. Default is recommended.")
set_property( CACHE hpGEM_LOGLEVEL PROPERTY STRINGS FATAL ERROR WARN INFO DEFAULT VERBOSE DEBUG )
mark_as_advanced( FORCE hpGEM_LOGLEVEL )
add_definitions( -DHPGEM_LOGLEVEL=Log::${hpGEM_LOGLEVEL} )

#===========[ Stacktraces ]=============
#we need this to check if we can
#even make stacktraces
include(CheckIncludeFiles)
include(CheckIncludeFileCXX)
include(CheckLibraryExists)

#Now, we need execinfo() to do it a bit
#in a sane way (nobody likes assembly)
CHECK_INCLUDE_FILES(execinfo.h  HAS_EXECINFO_H)

#Now, we would also like to have some
#name demangling going on...
#dynlinker
CHECK_INCLUDE_FILES(dlfcn.h     HAS_DLFCN_H)
#CXX demangler api
CHECK_INCLUDE_FILE_CXX(cxxabi.h HAS_CXXABI_H)
CHECK_LIBRARY_EXISTS(dl dladdr "" HAS_DLADDR_LIB)

#So, now see what we've end up with.
if (${HAS_EXECINFO_H})
  #the bare minimum is there...
  if (CMAKE_BUILD_TYPE)
    if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
      set( hpGEM_STACKTRACE_ENABLE ON CACHE BOOL "Enable stacktraces for log messages")
    else()
      set( hpGEM_STACKTRACE_ENABLE OFF CHACHE BOOL "Enable stacktraces for log messages")
    endif()
  endif()
  #now, check whether or not we can do demangling
  set( USE_DEMANGLE ON )
  if (NOT HAS_DLFCN_H)
    set( USE_DEMANGLE OFF )
  endif()
  if (NOT HAS_DLADDR_LIB)
    set( USE_DEMANGLE OFF )
  endif()
  if (NOT HAS_CXXABI_H)
    set( USE_DEMANGLE OFF )
  endif()
  set( hpGEM_STACKTRACE_DEMANGLE ${USE_DEMANGLE} CACHE BOOL "Enable name demangling for stacktraces")
else()
  set( hpGEM_STACKTRACE_ENABLE   OFF CACHE INTERNAL "execinfo.h not found - no stacktraces possible")
  set( hpGEM_STACKTRACE_DEMANGLE OFF CACHE INTERNAL "demangling is not going to happen, no bt!")
endif()
mark_as_advanced( FORCE hpGEM_STACKTRACE_ENABLE hpGEM_STACKTRACE_DEMANGLE )
#Well, that should work now...

if (hpGEM_STACKTRACE_ENABLE)
  add_definitions(-DHPGEM_STACKTRACE_ENABLE)
endif()
if (hpGEM_STACKTRACE_DEMANGLE)
  add_definitions(-DHPGEM_STACKTRACE_DEMANGLE)
endif()

#This section is about the build of documentation#
##############################################

option(hpGEM_BUILD_DOCUMENTATION "Use Doxygen to create the HTML based API documentation" OFF)
option(hpGEM_USE_METIS "Include metis this will enable parallel computation" OFF)
option(hpGEM_USE_MPI "Include mpi this will enable parallel computation" OFF)
option(hpGEM_USE_QHULL "Include QHull this will allow unstructured mesh generation and r-refinement" OFF)
option(hpGEM_USE_PETSC "Include PETSC this is required for GlobalAssemble" OFF)
option(hpGEM_USE_COMPLEX_PETSC "Include a more complex version of PETSc instead" OFF)
option(hpGEM_USE_SLEPC "Include SLEPC this is needed for some applications" OFF)
option(hpGEM_USE_EXTERNAL_BLAS "Try to use a version of BLAS optimised for your system" ON)
option(hpGEM_USE_SUNDIALS "Include SUNDIALS containing robust time integrators and nonlinear solvers" OFF)
mark_as_advanced(hpGEM_USE_EXTERNAL_BLAS)
option(hpGEM_USE_EXTERNAL_LAPACK "Try to use a version of LAPACK optimised for your system" ON)
mark_as_advanced(hpGEM_USE_EXTERNAL_LAPACK)
option(hpGEM_USE_STL_VECTOR_FOR_LA "Will swap the internal storage of the linear algebra from being based on valarray to STL vector" ON)
mark_as_advanced(hpGEM_USE_STL_VECTOR_FOR_LA)
if(hpGEM_USE_STL_VECTOR_FOR_LA)
	add_definitions(-DLA_STL_VECTOR)
endif(hpGEM_USE_STL_VECTOR_FOR_LA)
if (NOT hpGEM_USE_STL_VECTOR_FOR_LA)
	message (FATAL_ERROR
		"Val array is no longer support; please select use hpGEM_USE_STK_VECTOR_FOR_LA")
endif()
if(hpGEM_BUILD_DOCUMENTATION)
    FIND_PACKAGE(Doxygen)
    if (NOT DOXYGEN_FOUND)
        message(FATAL_ERROR 
            "Doxygen is needed to build the documentation. Please install it correctly or turn the option")
    endif()

  if (DOXYGEN_FOUND)
    
    #This is the configure file for normal doxygen builds	
    configure_file(conf/doxygen/doxygen.conf 
        	 	${PROJECT_BINARY_DIR}/conf/doxygen/doxygen.conf  @ONLY IMMEDIATE)

    #The next four and for website doxygen builds. The should be hinded from the public cmake at some point
    configure_file(conf/doxygen/web_doxygen.conf
	    		${PROJECT_BINARY_DIR}/conf/doxygen/web_doxygen.conf @ONLY IMMEDIATE)

    configure_file(conf/doxygen/hpg.css
	    		${PROJECT_BINARY_DIR}/conf/doxygen/hpg.css @ONLY IMMEDIATE)

    configure_file(conf/doxygen/new_footer.html
	    		${PROJECT_BINARY_DIR}/conf/doxygen/new_footer.html @ONLY IMMEDIATE)

    configure_file(conf/doxygen/new_header.html
	    		${PROJECT_BINARY_DIR}/conf/doxygen/new_header.html @ONLY IMMEDIATE)





    #-- Add custom targets to both make and clean (delete) the documentation
        add_custom_target 	(doc
                            COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/conf/doxygen/doxygen.conf
                            SOURCES ${PROJECT_BINARY_DIR}/conf/doxygen/doxygen.conf)

	add_custom_target	(docClean
                        COMMAND rm -r ${PROJECT_BINARY_DIR}/docs/*
			COMMENT "Cleaning (deleting) the documentation"	) 

	add_custom_target	(docWeb
			COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/conf/doxygen/web_doxygen.conf
			SOURCES ${PROJECT_BINARY_DIR}/conf/doxygen/web_doxygen.conf)

	add_custom_target	(docPublishAlpha
			COMMAND rsync -vr docs/html/* WebAdmin@einder.ewi.utwente.nl:/home/www/html/hpGEM/assets/doxygen/alpha/
			COMMENT "Publushing (uploading) the curent documentation as ALPHA realese" )
  endif()
endif()

if(hpGEM_USE_MPI)
    if(NOT hpGEM_USE_METIS)
        message(FATAL_ERROR "The option you have chosen requires some form of mesh decomposition, please enable metis")
    endif()
    FIND_PACKAGE(MPI REQUIRED)
    add_definitions(-DHPGEM_USE_MPI)
    include_directories(${MPI_CXX_INCLUDE_PATH})
    if(NOT MPI_FOUND)
        message(FATAL_ERROR "The option you have chosen requires mpi and you do not have this installed. Please install")
    endif()
endif()

if(hpGEM_USE_METIS)
    if(NOT hpGEM_USE_MPI)
        message(FATAL_ERROR "metis is useless without mpi, please use MPI or do not use metis")
    endif()
    FIND_PACKAGE(METIS REQUIRED)
    add_definitions(-DHPGEM_USE_METIS)
    include_directories(${METIS_INCLUDE_DIR})
    if(NOT METIS_FOUND)
        message(FATAL_ERROR "The option you have chosen requires metis and you do not have this installed. Please install")
    endif()
endif()

if(hpGEM_USE_QHULL)
    FIND_PACKAGE(QHULL REQUIRED)
    add_definitions(-DHPGEM_USE_QHULL)
    if(NOT QHULL_FOUND)
        message(FATAL_ERROR "The option you have chosen requires QHull and you do not have this installed. Please install")
    endif()
endif()

if (hpGEM_USE_PETSC OR hpGEM_USE_COMPLEX_PETSC)
	set(PETSC_ARCH "NOTFOUND" CACHE STRING "PETSc arch string should be set be the user; he will know this if he has installed PETSc")
	set(hpGEM_LOGLEVEL "DEFAULT" CACHE STRING "Verbosity of hpGEM. Default is recommended.")
endif()

if(hpGEM_USE_PETSC) 
    if(hpGEM_USE_COMPLEX_PETSC)
        message(FATAL_ERROR "The options 'hpGEM_USE_PETSC' and 'hpGEM_USE_COMPLEX_PETSC' are mutually exclusive, please disable one of them")
    endif()
        enable_language(C)
	FIND_PACKAGE(PETSc REQUIRED COMPONENTS CXX)
        add_definitions(-DHPGEM_USE_PETSC)
        if (NOT PETSC_FOUND)
		message(FATAL_ERROR
			"The option you have choosen requires PETSc and you do not have this installed. Please install")
	endif()

#endif PETSc required
endif()

if(hpGEM_USE_COMPLEX_PETSC)
    if(hpGEM_USE_PETSC)
        message(FATAL_ERROR "The options 'hpGEM_USE_PETSC' and 'hpGEM_USE_COMPLEX_PETSC' are mutually exclusive, please disable one of them")
    endif()
        enable_language(C)
	FIND_PACKAGE(PETSc REQUIRED COMPONENTS CXX)
        add_definitions(-DHPGEM_USE_COMPLEX_PETSC)
        if (NOT PETSC_FOUND)
		message(FATAL_ERROR
			"The option you have choosen requires PETSc and you do not have this installed. Please install")
	endif()

#endif PETSc required
endif()

if(hpGEM_USE_SLEPC)
        enable_language(C)
        if(NOT hpGEM_USE_PETSC AND NOT hpGEM_USE_COMPLEX_PETSC)
               message(FATAL_ERROR
                       "SLEPc depends on PETSc, please also enable PETSc support")
        endif()
        FIND_PACKAGE(SLEPc REQUIRED)
        add_definitions(-DHPGEM_USE_SLEPC)
        if(NOT SLEPC_FOUND)
               message(FATAL_ERROR
                       "The option you have chosen requires SLEPc and you do not have this installed. Please install")
        endif()
endif()

if(hpGEM_USE_EXTERNAL_BLAS)
    FIND_PACKAGE(BLAS QUIET)
    if(NOT BLAS_FOUND)
        message(SYSTEM "Could not find external BLAS, using reference implementation instead")
        set(hpGEM_USE_EXTERNAL_BLAS OFF CACHE BOOL "Try to use a version of BLAS optimised for your system" FORCE)
    endif()
endif()

if(hpGEM_USE_EXTERNAL_LAPACK)
    FIND_PACKAGE(LAPACK QUIET)
    if(NOT LAPACK_FOUND)
        message(SYSTEM "Could not find external LAPACK, using reference implementation instead")
        set(hpGEM_USE_EXTERNAL_LAPACK OFF CACHE BOOL "Try to use a version of LAPACK optimised for your system" FORCE)
    endif()
endif()

if(hpGEM_USE_SUNDIALS)
    if(hpGEM_USE_MPI)
        message(FATAL_ERROR "SUNDIALS is only serial (for now)")
    endif()
    #enable_language(C)
    FIND_PACKAGE(SUNDIALS QUIET)
    add_definitions(-DHPGEM_USE_SUNDIALS)
    if(NOT SUNDIALS_FOUND)
           message(FATAL_ERROR
                   "The option you have choosen requires SUNDIALS and you do not have this installed. Please install")
    endif()
endif()



add_subdirectory(kernel)
add_subdirectory(tests)
add_subdirectory(applications)

