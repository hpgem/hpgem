# - Try to find PETSc
# Once done this will define
#
#  PETSc_FOUND        - system has PETSc
#  PETSc_INCLUDE_DIRS - the PETSc include directories
#  PETSc_LIBRARIES    - Link these to use PETSc
#  PETSc_IS_COMPLEX   - supports complex numbers
#  PETSc_ARCH  - Compiler switches for using PETSc
#
#

find_package(PkgConfig)
list(APPEND CMAKE_PREFIX_PATH "${PETSC_DIR}/${PETSC_ARCH}")
pkg_check_modules(PETSC_PKG PETSc)

find_path(PETSc_INCLUDE_DIR petsc.h HINTS ${PETSC_PKG_INCLUDE_DIRS})
find_path(PETSc_INCLUDE_DIR2 petscconf.h HINTS ${PETSC_PKG_INCLUDE_DIRS})
list(APPEND PETSc_INCLUDE_DIR "${PETSc_INCLUDE_DIR2}")
find_library(PETSc_LIBRARY NAMES petsc HINTS ${PETSC_PKG_LIBRARY_DIRS} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set PETSC_FOUND to TRUE
# if all listed variables are TRUE

find_package_handle_standard_args(PETSc DEFAULT_MSG PETSc_LIBRARY PETSc_INCLUDE_DIR )

message(STATUS "${PETSc_INCLUDE_DIR}")

if (PETSc_FOUND)
  include(CheckSymbolExists)
  foreach(petsc_dir ${PETSc_INCLUDE_DIR})
      
      set(petsc_petscconf_path "${petsc_dir}/petscconf.h")
      if(EXISTS "${petsc_petscconf_path}")
          set(CMAKE_REQUIRED_QUIET 1)
          check_symbol_exists(PETSC_USE_COMPLEX ${petsc_petscconf_path} PETSc_IS_COMPLEX)
          file(STRINGS ${petsc_petscconf_path} _contents REGEX "#define PETSC_ARCH[ \\t]+")
          if(_contents)
              string(REGEX REPLACE ".*#define PETSC_ARCH[ \\t]+\"?([^\\t \"]*)\"?" "\\1" PETSc_ARCH "${_contents}")
          endif()

      endif()
  endforeach()
endif ()

# Copy the results to the output variables and target.
if(PETSc_FOUND)
  set(PETSc_LIBRARIES ${PETSc_LIBRARY} )
  set(PETSc_INCLUDE_DIRS ${PETSc_INCLUDE_DIR} )

  if(PETSC_PKG_FOUND)
        add_library(PETSc::PETSc UNKNOWN IMPORTED)
        set_target_properties(PETSc::PETSc PROPERTIES 
            IMPORTED_LOCATION "${PETSc_LIBRARY}" 
            INTERFACE_INCLUDE_DIRECTORIES "${PETSc_INCLUDE_DIR}")
  endif()

endif()
mark_as_advanced(PETSc_INCLUDE_DIRS PETSc_LIBRARIES PETSc_IS_COMPLEX)

