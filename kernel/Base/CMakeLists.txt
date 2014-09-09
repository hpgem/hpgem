include_directories(${hpGEM_SOURCE_DIR}/kernel/ .)

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
                mesh.cpp
                submesh.cpp
        	${hpGEM_SOURCE_DIR}/kernel/Base/AssembleBasisFunctionSet.cpp
		${hpGEM_SOURCE_DIR}/kernel/Base/L2Norm.cpp
 		${hpGEM_SOURCE_DIR}/kernel/Base/Norm2.cpp
		${hpGEM_SOURCE_DIR}/kernel/Base/Logger.cc

		ShortTermStorageElementBase.cpp
		ShortTermStorageElementH1.cpp
		ShortTermStorageFaceBase.cpp
		ShortTermStorageFaceH1.cpp
	) 

target_link_libraries(hpGEM_Base Geometry Reference_geometries)
