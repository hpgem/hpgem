# - Try to find SLEPc
# Once done this will define
#
#  SLEPc_FOUND        - system has SLEPc
#  SLEPc_INCLUDE_DIRS - the SLEPc include directories
#  SLEPc_LIBRARIES    - Link these to use SLEPc
#  SLEPc_ARCH         - Compiler switches for using SLEPc

find_package(PkgConfig)

pkg_check_modules(SLEPc_PKG SLEPc)

find_path(SLEPc_INCLUDE_DIR slepc.h HINTS ${SLEPc_PKG_INCLUDE_DIRS})

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
          string(REGEX REPLACE ".*#define SLEPC_PETSC_ARCH[ \\t]+\"?([^\\t \"]+)\"?" "\\1" SLEPc_ARCH "${_contents}")
          message(STATUS "${SLEPc_ARCH}")
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