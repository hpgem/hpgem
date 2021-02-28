# - Try to find SLEPc
# Once done this will define
#
#  SLEPc_FOUND        - system has SLEPc
#  SLEPc_INCLUDE_DIRS - the SLEPc include directories
#  SLEPc_LIBRARIES    - Link these to use SLEPc
#  SLEPc_ARCH         - Compiler switches for using SLEPc

find_package(PkgConfig)
list(APPEND CMAKE_PREFIX_PATH "${SLEPC_DIR}/${PETSC_ARCH}")

# SLEPC decided to change the casing of the slepc module from SLEPc to slepc,
# probably with version 3.14. Unfortunately there does not seem to be a case
# insensitive version of pkg_check_modules, so we have to do it this way.
pkg_check_modules(SLEPc_PKG slepc)
if (NOT SLEPc_PKG_FOUND)
    pkg_check_modules(SLEPc_PKG SLEPc)
    if (NOT SLEPc_PKG_FOUND)
        message(FATAL_ERROR "SLEPc pkgconfig file not found")
    endif()
endif()

find_path(SLEPc_INCLUDE_DIR slepc.h HINTS ${SLEPc_PKG_INCLUDE_DIRS})
find_path(SLEPc_INCLUDE_DIR2 slepcconf.h HINTS ${SLEPc_PKG_INCLUDE_DIRS})
list(APPEND SLEPc_INCLUDE_DIR "${SLEPc_INCLUDE_DIR2}")
find_library(SLEPc_LIBRARY NAMES slepc HINTS ${SLEPc_PKG_LIBRARY_DIRS} )
include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set SLEPc_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(SLEPc DEFAULT_MSG SLEPc_LIBRARY SLEPc_INCLUDE_DIR )

if (SLEPc_FOUND)
  include(CheckSymbolExists)
  foreach(SLEPc_dir ${SLEPc_INCLUDE_DIR})
      set(SLEPc_SLEPcconf_path "${SLEPc_dir}/slepcconf.h")
      if(EXISTS "${SLEPc_SLEPcconf_path}")
        file(STRINGS ${SLEPc_SLEPcconf_path} _contents REGEX "#define SLEPC_PETSC_ARCH[ \\t]+")
          if(_contents)
            string(REGEX REPLACE ".*#define SLEPC_PETSC_ARCH[ \\t]+\"?([^\\t \"]*)\"?" "\\1" SLEPc_ARCH "${_contents}")
          endif()
      endif()
  endforeach()
endif ()

# Copy the results to the output variables and target.
if(SLEPc_FOUND)
  set(SLEPc_LIBRARIES ${SLEPc_LIBRARY} )
  set(SLEPc_INCLUDE_DIRS ${SLEPc_INCLUDE_DIR} )

  if(SLEPc_PKG_FOUND)
        add_library(SLEPc::SLEPc UNKNOWN IMPORTED)
        set_target_properties(SLEPc::SLEPc PROPERTIES 
            IMPORTED_LOCATION "${SLEPc_LIBRARY}" 
            INTERFACE_INCLUDE_DIRECTORIES "${SLEPc_INCLUDE_DIR}")
  endif()

endif()
mark_as_advanced(SLEPc_INCLUDE_DIRS SLEPc_LIBRARIES SLEPc_ARCH)