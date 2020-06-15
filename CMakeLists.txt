cmake_minimum_required (VERSION 3.10)
project (hpGEM C CXX)

set(hpGEM_VERSION_MAJOR 2)
set(hpGEM_VERSION_MINOR 0)

## Set CMAKE options

#Directory were the external cmake modules are stored
include(FeatureSummary)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/conf/cmake/Modules/")
#Make the default build release is this in what a normal user wants
###################################################################
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE "Release" CACHE
			STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
	set_property( CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Release Debug RelWithDebInfo MinSizeRel )
endif(NOT CMAKE_BUILD_TYPE)


#allow dynamically locating dynamic libraries if cmake and the OS support it
set(CMAKE_MACOSX_RPATH 1)


# Enable testing via flag
option(ENABLE_TESTING "Build and copy testing stuff" OFF)
add_feature_info(ENABLE_TESTING ENABLE_TESTING "Build and copy testing stuff")
if(ENABLE_TESTING)
  enable_testing()
endif()

#Here is the check for CXX14 support : We now use some features of this, so turn it on if possible
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 14)

option(ENABLE_WARNING_FLAGS "Inject more warning flags" OFF)
if(ENABLE_WARNING_FLAGS)
  include(CheckCXXCompilerFlag)
  foreach(_FLAG -Wall -Wextra -Wpedantic -Wshadow -Wconversion) 
    check_cxx_compiler_flag("${_FLAG}" COMPILER_SUPPORTS${_FLAG})
    if(COMPILER_SUPPORTS${_FLAG})
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${_FLAG}")
    endif()
  endforeach()
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-sign-conversion")
  endif()
  if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wcheck")
  endif()
endif()

option(ENABLE_WERROR "Inject -Werror" OFF)
if(ENABLE_WERROR)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
endif()

option(ENABLE_COVERAGE_BUILD "Do a coverage build" OFF)
if(ENABLE_COVERAGE_BUILD)
    message(STATUS "Enabling coverage build")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --coverage")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} --coverage")
endif()



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
add_feature_info(hpGEM_BUILD_DOCUMENTATION hpGEM_BUILD_DOCUMENTATION "Use Doxygen to create the HTML based API documentation")
option(hpGEM_USE_METIS "Include metis this will enable parallel computation" OFF)
add_feature_info(hpGEM_USE_METIS hpGEM_USE_METIS "Include metis this will enable parallel computation")
option(hpGEM_USE_MPI "Include mpi this will enable parallel computation" OFF)
add_feature_info(hpGEM_USE_MPI hpGEM_USE_MPI "Include mpi this will enable parallel computation")
option(hpGEM_USE_QHULL "Include QHull this will allow unstructured mesh generation and r-refinement" OFF)
add_feature_info(hpGEM_USE_QHULL hpGEM_USE_QHULL "Include QHull this will allow unstructured mesh generation and r-refinement")
option(hpGEM_USE_PETSC "Include PETSC with real numbers, which is required for assembling a global system" OFF)
add_feature_info(hpGEM_USE_PETSC hpGEM_USE_PETSC "Include PETSC with real numbers, which is required for assembling a global system")
option(hpGEM_USE_COMPLEX_PETSC "Include PETSc with complex numbers" OFF)
add_feature_info(hpGEM_USE_COMPLEX_PETSC hpGEM_USE_COMPLEX_PETSC "Include PETSc with complex numbers")
option(hpGEM_USE_SLEPC "Include SLEPC this is needed for some applications" OFF)
add_feature_info(hpGEM_USE_SLEPC hpGEM_USE_SLEPC "Include SLEPC this is needed for some applications")


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
if(ENABLE_TESTING)
	add_subdirectory(tests)
endif()
add_subdirectory(applications)





feature_summary(INCLUDE_QUIET_PACKAGES WHAT ALL)
