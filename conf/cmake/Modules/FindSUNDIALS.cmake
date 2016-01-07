# - Find KINSOL
# This module looks for KINSOL (Message Passing Interface) support
# it will define the following values
#  KINSOL_INCLUDE_PATH = cached location of kinsol.h
#  KINSOL_LIBRARIES    = cached list of libraries to link in (mpi mpich etc)

#FIND_PATH(KINSOL_INCLUDE_PATH kinsol.h
#  PATHS /usr/local/include
#  /usr/local/include/kinsol 
#  /usr/include 
#  /usr/include/kinsol
#  )
#
#FIND_LIBRARY(KINSOL_LIBRARY
#  NAMES sundials_kinsol 
#  PATHS /usr/lib /usr/local/lib
#  /usr/lib/kinsol /usr/local/lib/kinsol
#  )

#if ( KINSOL_INCLUDE_PATH AND KINSOL_LIBRARY)
#    set(KINSOL_FOUND TRUE)
#else()
#    set(KINSOL_FOUND FALSE)
#endif()

#FIND_PATH(NVECTOR_INCLUDE_PATH nvector_serial.h
#  PATHS /usr/local/include
#  /usr/local/include/nvector
#  /usr/include
#  /usr/include/nvector
#  )

message(STATUS "Looking for sundials")

find_path(SUNDIALS_INCLUDE_DIR
    nvector/nvector_serial.h
    sundials/sundials_dense.h
    sundials/sundials_types.h
    kinsol/kinsol.h
    kinsol/kinsol_dense.h
    kinsol/kinsol_band.h
    kinsol/kinsol_spgmr.h
    kinsol/kinsol_spbcgs.h
    kinsol/kinsol_sptfqmr.h
    kinsol/kinsol_impl.h
    idas/idas.h
    idas/idas_dense.h
    idas/idas_band.h
    idas/idas_spgmr.h
    idas/idas_spbcgs.h
    idas/idas_sptfqmr.h
    idas/idas_impl.h
    cvodes/cvodes.h
    cvodes/cvodes_dense.h
    cvodes/cvodes_band.h
    cvodes/cvodes_spgmr.h
    cvodes/cvodes_spbcgs.h
    cvodes/cvodes_sptfqmr.h
    cvodes/cvodes_impl.h
  HINTS $ENV{SUNDIALS}/include
  )

if(SUNDIALS_INCLUDE_DIR)
  message(STATUS "Found sundials include directory: ${SUNDIALS_INCLUDE_DIR}")
else()
  message(STATUS "Could not fund sundials include dir")
endif()

# libraries
set(SUNDIALS_LIBS_LIST
#  sundials_cvode
  sundials_cvodes
#  sundials_fnvecserial
#  sundials_ida
  sundials_idas
  sundials_kinsol
  sundials_nvecserial
  )

set(SUNDIALS_LIBRARIES)
foreach(LIB ${SUNDIALS_LIBS_LIST})
  find_library(SUNDIALS_LIB_${LIB}
    NAMES ${LIB}
    HINTS $ENV{SUNDIALS}/lib)
  if(SUNDIALS_LIB_${LIB})
#    message(STATUS "Found ${LIB}: ${SUNDIALS_LIB_${LIB}}")
    set(SUNDIALS_LIBRARIES ${SUNDIALS_LIBRARIES} ${SUNDIALS_LIB_${LIB}})
#  else()
#    message(STATUS "Could not find lib${LIB}")
  endif()
endforeach()

if(SUNDIALS_LIBRARIES)
  set(SUNDIALS_LIBRARIES ${SUNDIALS_LIBRARIES})
  message(STATUS "Found sundials libraries ${SUNDIALS_LIBRARIES}")
  set(SUNDIALS_FOUND_LIBS TRUE)
else()
  message(STATUS "Could not find sundials libraries")
endif()

if(SUNDIALS_INCLUDE_DIR AND SUNDIALS_FOUND_LIBS)
  set(SUNDIALS_FOUND TRUE)
else()
  message(STATUS "SUNDIALS: Cound not find sundials. Try setting SUNDIALS env var.")
endif()



