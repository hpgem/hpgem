include_directories(${hpGEM_SOURCE_DIR}/kernel/ .)

add_library(hpGEM_Base SHARED
		${hpGEM_SOURCE_DIR}/kernel/Output/TecplotPhysicalGeometryIterator.cpp	
 		${hpGEM_SOURCE_DIR}/kernel/Base/BasisFunctionSet.cpp		  
 		Element.cpp
 		ElementData.cpp
 		ElementCacheData.cpp
 		${hpGEM_SOURCE_DIR}/kernel/Integration/ElementIntegral.cpp
 		Face.cpp
 		FaceData.cpp
 		FaceCacheData.cpp
 		${hpGEM_SOURCE_DIR}/kernel/Integration/FaceIntegral.cpp
 		HpgemUI.cpp
		#miscellaneous required stuff
		${hpGEM_SOURCE_DIR}/kernel/Base/MeshManipulator.cpp
        	${hpGEM_SOURCE_DIR}/kernel/Base/AssembleBasisFunctionSet.cpp
		${hpGEM_SOURCE_DIR}/kernel/Base/L2Norm.cpp
 		${hpGEM_SOURCE_DIR}/kernel/Base/Norm2.cpp


	) 

target_link_libraries(hpGEM_Base Geometry Reference_geometries)
