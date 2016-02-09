include_directories(${hpGEM_SOURCE_DIR}/kernel/ .)

if(hpGEM_USE_METIS)
  include_directories(${METIS_INCLUDE_DIR})
endif()

if(hpGEM_USE_PETSC OR hpGEM_USE_COMPLEX_PETSC)
  include_directories(${PETSC_INCLUDES})
endif()

if(hpGEM_USE_SLEPC)
  include_directories(${SLEPC_INCLUDE_DIRS})
endif()

if(hpGEM_USE_QHULL)
  include_directories(${QHULL_INCLUDE_DIR})
endif()

if(hpGEM_USE_SUNDIALS) 
  include_directories(${SUNDIALS_INCLUDE_DIR})
endif()

#add_library(hpGEM_Serialization
# 		#H1ConformingTransformation.cpp
#		)

#target_link_libraries(hpGEM_Serialization hpGEM_Base Output Geometry TimeIntegration #boost_serialization)

#set_target_properties(hpGEM_Serialization PROPERTIES POSITION_INDEPENDENT_CODE true)

set( DEP_COMPILE_FLAGS "" )
set( DEP_LINKER_FLAGS "" )
set( DEP_LIBRARIES "")

#if(hpGEM_USE_MPI)
#set( DEP_COMPILE_FLAGS "${DEP_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}" )
#set( DEP_LINKER_FLAGS "${DEP_LINKER_FLAGS} ${MPI_LINK_FLAGS}")
#target_link_libraries(hpGEM_Serialization "${MPI_LIBRARIES}")
#endif()

#if(hpGEM_USE_METIS)
#target_link_libraries(hpGEM_Serialization "${METIS_LIBRARIES}")
#endif()

#if(hpGEM_USE_PETSC OR hpGEM_USE_COMPLEX_PETSC)
#target_link_libraries(hpGEM_Serialization "${PETSC_LIBRARIES}")
#endif()

#if(hpGEM_USE_SLEPC)
#target_link_libraries(hpGEM_Serialization "${SLEPC_LIBRARIES}")
#endif()

#if(hpGEM_USE_QHULL)
#    if(CMAKE_BUILD_TYPE MATCHES "Debug")
#        #if it is possible to link to the debug libraries, do so
#        #if there are no debug libraries, this will try to link against the normal #libraries instead
#        target_link_libraries(hpGEM_Serialization "${QHULL_DEBUG_LIBRARIES}")
#    else()
#        target_link_libraries(hpGEM_Serialization "${QHULL_LIBRARIES}")
#    endif()
#endif()

#if(hpGEM_USE_SUNDIALS)
#target_link_libraries(hpGEM_Serialization "${SUNDIALS_LIB_sundials_cvodes}")
#target_link_libraries(hpGEM_Serialization "${SUNDIALS_LIB_sundials_idas}")
#target_link_libraries(hpGEM_Serialization "${SUNDIALS_LIB_sundials_kinsol}")
#target_link_libraries(hpGEM_Serialization "${SUNDIALS_LIB_sundials_nvecserial}")
#endif()

#set_target_properties(hpGEM_Serialization PROPERTIES COMPILE_FLAGS “${DEP_COMPILE_FLAGS}")
#set_target_properties(hpGEM_Serialization PROPERTIES LINK_FLAGS "${DEP_LINKER_FLAGS}")
