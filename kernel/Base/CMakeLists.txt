include_directories(${hpGEM_SOURCE_DIR}/kernel/ .)
add_subdirectory(TimeIntegration)

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

add_library(hpGEM_Base STATIC
 		${hpGEM_SOURCE_DIR}/kernel/Base/BasisFunctionSet.cpp	
 		BaseBasisFunction.cpp	  
		CommandLineOptions.cpp
 		Element.cpp
 		ElementData.cpp
        ElementFactory.cpp
 		Face.cpp
 		FaceData.cpp
        FaceFactory.cpp
		FaceMatrix.cpp
 		Edge.cpp
 		Node.cpp
		HpgemAPISimplified.cpp
		#miscellaneous required stuff
		${hpGEM_SOURCE_DIR}/kernel/Base/MeshManipulatorBase.cpp
        Submesh.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Base/AssembleBasisFunctionSet.cpp
		${hpGEM_SOURCE_DIR}/kernel/Base/L2Norm.cpp
        MpiContainer.cpp
		        
        ${hpGEM_SOURCE_DIR}/kernel/Utilities/BasisFunctions1DH1ConformingLine.cpp
		${hpGEM_SOURCE_DIR}/kernel/Utilities/BasisFunctions2DH1ConformingSquare.cpp
		${hpGEM_SOURCE_DIR}/kernel/Utilities/BasisFunctions2DH1ConformingTriangle.cpp
		${hpGEM_SOURCE_DIR}/kernel/Utilities/BasisFunctions3DH1ConformingCube.cpp
		${hpGEM_SOURCE_DIR}/kernel/Utilities/BasisFunctions3DH1ConformingPrism.cpp
		${hpGEM_SOURCE_DIR}/kernel/Utilities/BasisFunctions3DH1ConformingTetrahedron.cpp
		${hpGEM_SOURCE_DIR}/kernel/Utilities/BasisFunctions3DNedelec.cpp
		${hpGEM_SOURCE_DIR}/kernel/Utilities/BasisFunctions3DAinsworthCoyle.cpp
		${hpGEM_SOURCE_DIR}/kernel/Utilities/helperFunctions.cpp
		${hpGEM_SOURCE_DIR}/kernel/Utilities/GlobalMatrix.cpp
		${hpGEM_SOURCE_DIR}/kernel/Utilities/GlobalVector.cpp
		${hpGEM_SOURCE_DIR}/kernel/Serialization/H1ConformingTransformation.cpp)

target_link_libraries(hpGEM_Base Output Geometry TimeIntegration boost_serialization)

set_target_properties(hpGEM_Base PROPERTIES POSITION_INDEPENDENT_CODE true)

set( DEP_COMPILE_FLAGS "" )
set( DEP_LINKER_FLAGS "" )
set( DEP_LIBRARIES "")

if(hpGEM_USE_MPI)
set( DEP_COMPILE_FLAGS "${DEP_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}" )
set( DEP_LINKER_FLAGS "${DEP_LINKER_FLAGS} ${MPI_LINK_FLAGS}")
target_link_libraries(hpGEM_Base "${MPI_LIBRARIES}")
endif()

if(hpGEM_USE_METIS)
target_link_libraries(hpGEM_Base "${METIS_LIBRARIES}")
endif()

if(hpGEM_USE_PETSC OR hpGEM_USE_COMPLEX_PETSC)
target_link_libraries(hpGEM_Base "${PETSC_LIBRARIES}")
endif()

if(hpGEM_USE_SLEPC)
target_link_libraries(hpGEM_Base "${SLEPC_LIBRARIES}")
endif()

if(hpGEM_USE_QHULL)
    if(CMAKE_BUILD_TYPE MATCHES "Debug")
        #if it is possible to link to the debug libraries, do so
        #if there are no debug libraries, this will try to link against the normal libraries instead
        target_link_libraries(hpGEM_Base "${QHULL_DEBUG_LIBRARIES}")
    else()
        target_link_libraries(hpGEM_Base "${QHULL_LIBRARIES}")
    endif()
endif()

if(hpGEM_USE_SUNDIALS)
target_link_libraries(hpGEM_Base "${SUNDIALS_LIB_sundials_cvodes}")
target_link_libraries(hpGEM_Base "${SUNDIALS_LIB_sundials_idas}")
target_link_libraries(hpGEM_Base "${SUNDIALS_LIB_sundials_kinsol}")
target_link_libraries(hpGEM_Base "${SUNDIALS_LIB_sundials_nvecserial}")
endif()

set_target_properties(hpGEM_Base PROPERTIES COMPILE_FLAGS "${DEP_COMPILE_FLAGS}")
set_target_properties(hpGEM_Base PROPERTIES LINK_FLAGS "${DEP_LINKER_FLAGS}")
