cmake_minimum_required (VERSION 3.3)
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

#Here is the check for CXX14 support : We now use some features of this, so turn it on if possible
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 14)




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

#==========[ assert_debug Config ]============
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
  link_libraries(dl)
endif()

#This section is about the build of documentation#
##############################################

# Optional dependencies
option(hpGEM_BUILD_DOCUMENTATION "Use Doxygen to create the HTML based API documentation" OFF)
option(hpGEM_USE_METIS "Include metis this will enable parallel computation" OFF)
option(hpGEM_USE_MPI "Include mpi this will enable parallel computation" OFF)
option(hpGEM_USE_QHULL "Include QHull this will allow unstructured mesh generation and r-refinement" OFF)
option(hpGEM_USE_PETSC "Include PETSC with real numbers, which is required for assembling a global system" OFF)
option(hpGEM_USE_COMPLEX_PETSC "Include PETSc with complex numbers" OFF)
option(hpGEM_USE_SLEPC "Include SLEPC this is needed for some applications" OFF)
option(hpGEM_USE_SUNDIALS "Include SUNDIALS containing robust time integrators and nonlinear solvers" OFF)

# By default use system BLAS, LAPACK
option(hpGEM_USE_EXTERNAL_BLAS "Try to use a version of BLAS optimised for your system" ON)
option(hpGEM_USE_EXTERNAL_LAPACK "Try to use a version of LAPACK optimised for your system" ON)
mark_as_advanced(hpGEM_USE_EXTERNAL_BLAS) # These should only be used in special cases
mark_as_advanced(hpGEM_USE_EXTERNAL_LAPACK) # These should only be used in special cases

# Old option (partially broken) to use valarray in the linear algebra classes
option(hpGEM_USE_STL_VECTOR_FOR_LA "Will swap the internal storage of the linear algebra from being based on valarray to STL vector" ON)
mark_as_advanced(hpGEM_USE_STL_VECTOR_FOR_LA)
if(hpGEM_USE_STL_VECTOR_FOR_LA)
	add_definitions(-DLA_STL_VECTOR)
endif(hpGEM_USE_STL_VECTOR_FOR_LA)
if (NOT hpGEM_USE_STL_VECTOR_FOR_LA)
	message (FATAL_ERROR
		"Val array is no longer support; please select use hpGEM_USE_STK_VECTOR_FOR_LA")
endif()





################################################################################
####### CONFIGURE DOCUMENTATION & DEPENDENCIES
################################################################################

# Build documentation if requested, including HTML documentation
if(hpGEM_BUILD_DOCUMENTATION)
    include(conf/cmake/ConfigureDocumentation.cmake)
endif()

# Find all the optional dependencies
include(conf/cmake/ConfigureOptionalDependencies.cmake)




add_subdirectory(kernel)
add_subdirectory(tests)
add_subdirectory(applications)

