include_directories(${hpGEM_SOURCE_DIR}/kernel/ .)
add_subdirectory(TimeIntegration)

if(hpGEM_USE_METIS)
  include_directories(${METIS_INCLUDE_DIR})
endif()

if(hpGEM_USE_PETSC)
  include_directories(${PETSC_INCLUDES})
endif()

add_library(hpGEM_Base SHARED
		${hpGEM_SOURCE_DIR}/kernel/Output/TecplotPhysicalGeometryIterator.cpp	
                ${hpGEM_SOURCE_DIR}/kernel/Output/TecplotDiscontinuousSolutionWriter.cpp
                ${hpGEM_SOURCE_DIR}/kernel/Output/VTKTimeDependentWriter.cpp
                ${hpGEM_SOURCE_DIR}/kernel/Output/VTKSpecificTimeWriter.cpp
                ${hpGEM_SOURCE_DIR}/kernel/Output/base64.cpp
 		${hpGEM_SOURCE_DIR}/kernel/Base/BasisFunctionSet.cpp	
 		BaseBasisFunction.cpp	  
		CommandLineOptions.cpp
 		PhysGradientOfBasisFunction.cpp
 		Element.cpp
 		ElementData.cpp
 		ElementCacheData.cpp
                ElementFactory.cpp
 		${hpGEM_SOURCE_DIR}/kernel/Integration/ElementIntegral.cpp
 		Face.cpp
 		FaceData.cpp
 		FaceCacheData.cpp
                FaceFactory.cpp
 		${hpGEM_SOURCE_DIR}/kernel/Integration/FaceIntegral.cpp
 		Edge.cpp
 		Node.cpp
 		HpgemUI.cpp
                HpgemUISimplified.cpp
		#miscellaneous required stuff
		${hpGEM_SOURCE_DIR}/kernel/Base/MeshManipulator.cpp
                Mesh.cpp
                ProblemDescriptorTimeDependent.cpp
                Submesh.cpp
        	${hpGEM_SOURCE_DIR}/kernel/Base/AssembleBasisFunctionSet.cpp
		${hpGEM_SOURCE_DIR}/kernel/Base/L2Norm.cpp
 		${hpGEM_SOURCE_DIR}/kernel/Base/Norm2.cpp
		${hpGEM_SOURCE_DIR}/kernel/Base/Logger.cc
                MpiContainer.cpp
		ShortTermStorageElementBase.cpp
		ShortTermStorageElementH1.cpp
		ShortTermStorageFaceBase.cpp
		ShortTermStorageFaceH1.cpp
	) 

target_link_libraries(hpGEM_Base Geometry Reference_geometries TimeIntegration)

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

if(hpGEM_USE_PETSC)
target_link_libraries(hpGEM_Base "${PETSC_LIBRARIES}")
endif()

set_target_properties(hpGEM_Base PROPERTIES COMPILE_FLAGS "${DEP_COMPILE_FLAGS}")
set_target_properties(hpGEM_Base PROPERTIES LINK_FLAGS "${DEP_LINKER_FLAGS}")
