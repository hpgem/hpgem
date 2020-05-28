cmake_minimum_required (VERSION 3.3)
project (hpGEM C CXX Fortran)

set(hpGEM_VERSION_MAJOR 2)
set(hpGEM_VERSION_MINOR 0)

## Set CMAKE options

#Directory were the external cmake modules are stored
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/conf/cmake/Modules/")

#Make the default build release is this in what a normal user wants
###################################################################
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE "Release" CACHE
			STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
	set_property( CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Release Debug RelWithDebInfo MinSizeRel )
endif(NOT CMAKE_BUILD_TYPE)


#allow dynamically locating dynamic libraries if cmake and the OS support it
if(${CMAKE_VERSION} VERSION_GREATER 2.8.12)
    set(CMAKE_MACOSX_RPATH 1)
endif()

enable_testing()




#Here is the check for CXX14 support : We now use some features of this, so turn it on if possible
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 14)



################################################################################
####### CUSTOM TARGETS
#################################################################################
# Add the custom target make fullTest, that clean and check the whole build system
add_custom_target(fullTest
		COMMAND make rebuild_cache
		COMMAND make clean
		COMMAND make
		COMMAND make test
		COMMENT "This will clean the project, rebuild and test." )

################################################################################
####### OPTIONAL FLAGS FOR hpGEM
################################################################################

# Optional dependencies
option(hpGEM_BUILD_DOCUMENTATION "Use Doxygen to create the HTML based API documentation" OFF)
option(hpGEM_USE_METIS "Include metis this will enable parallel computation" OFF)
option(hpGEM_USE_MPI "Include mpi this will enable parallel computation" OFF)
option(hpGEM_USE_QHULL "Include QHull this will allow unstructured mesh generation and r-refinement" OFF)
option(hpGEM_USE_PETSC "Include PETSC with real numbers, which is required for assembling a global system" OFF)
option(hpGEM_USE_COMPLEX_PETSC "Include PETSc with complex numbers" OFF)
option(hpGEM_USE_SLEPC "Include SLEPC this is needed for some applications" OFF)
option(hpGEM_USE_SUNDIALS "Include SUNDIALS containing robust time integrators and nonlinear solvers" OFF)


################################################################################
####### LOGGER CONFIGURATION
################################################################################

# The logger configuration configures the
# - Whether to use assersions
# - What log messages to print
# - Possible printing of stacktraces and demangling of the stacktrace

include(conf/cmake/ConfigureLogger.cmake)

################################################################################
####### CONFIGURE DOCUMENTATION & DEPENDENCIES
################################################################################

# Build documentation if requested, including HTML documentation
if(hpGEM_BUILD_DOCUMENTATION)
    include(conf/cmake/ConfigureDocumentation.cmake)
endif()

# Find all the optional dependencies
include(conf/cmake/ConfigureDependencies.cmake)


################################################################################
####### EXPOSE CMAKE VARIABLES
################################################################################

# Make a header file which can contain relevant cmake variables
configure_file(
		${hpGEM_SOURCE_DIR}/conf/cmake/CMakeDefinitions.conf
		${hpGEM_SOURCE_DIR}/conf/cmake/CMakeDefinitions.h ESCAPE_QUOTES)
# And ensure it is included by the targets
include_directories(${hpGEM_SOURCE_DIR}/conf/cmake/)

################################################################################
####### KERNEL, TESTS, APPLICATIONS
################################################################################
add_subdirectory(kernel)
add_subdirectory(tests)
add_subdirectory(applications)

