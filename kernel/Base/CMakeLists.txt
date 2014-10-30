include_directories(${hpGEM_SOURCE_DIR}/kernel/ .)

if(hpGEM_INCLUDE_METIS_SUPPORT)
  include_directories(${METIS_INCLUDE_DIR})
endif()

add_library(hpGEM_Base SHARED
		${hpGEM_SOURCE_DIR}/kernel/Output/TecplotPhysicalGeometryIterator.cpp	
 		${hpGEM_SOURCE_DIR}/kernel/Base/BasisFunctionSet.cpp	
 		BaseBasisFunction.cpp	  
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
 		HpgemUI.cpp
		#miscellaneous required stuff
		${hpGEM_SOURCE_DIR}/kernel/Base/MeshManipulator.cpp
                Mesh.cpp
                Submesh.cpp
        	${hpGEM_SOURCE_DIR}/kernel/Base/AssembleBasisFunctionSet.cpp
		${hpGEM_SOURCE_DIR}/kernel/Base/L2Norm.cpp
 		${hpGEM_SOURCE_DIR}/kernel/Base/Norm2.cpp
		${hpGEM_SOURCE_DIR}/kernel/Base/Logger.cc

		ShortTermStorageElementBase.cpp
		ShortTermStorageElementH1.cpp
		ShortTermStorageFaceBase.cpp
		ShortTermStorageFaceH1.cpp
	) 

if(hpGEM_INCLUDE_METIS_SUPPORT)
target_link_libraries(hpGEM_Base Geometry Reference_geometries ${METIS_LIBRARY})
else()
target_link_libraries(hpGEM_Base Geometry Reference_geometries)
endif()
