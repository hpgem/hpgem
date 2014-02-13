include_directories(${hpGEM_SOURCE_DIR}/src/ .)

add_library(hpGEM_Base SHARED
		${hpGEM_SOURCE_DIR}/src/Output/TecplotPhysicalGeometryIterator.cpp	
 		${hpGEM_SOURCE_DIR}/src/Base/BasisFunctionSet.cpp		  
 		Element.cpp
 		ElementData.cpp
 		ElementCacheData.cpp
 		${hpGEM_SOURCE_DIR}/src/Integration/ElementIntegral.cpp
 		Face.cpp
 		FaceData.cpp
 		FaceCacheData.cpp
 		${hpGEM_SOURCE_DIR}/src/Integration/FaceIntegral.cpp
 		HpgemUI.cpp
		#miscellaneous required stuff
		${hpGEM_SOURCE_DIR}/src/Base/MeshManipulator.cpp
        	${hpGEM_SOURCE_DIR}/src/Base/AssembleBasisFunctionSet.cpp
		${hpGEM_SOURCE_DIR}/src/Base/L2Norm.cpp
 		${hpGEM_SOURCE_DIR}/src/Base/Norm2.cpp
	) 

target_link_libraries(hpGEM_Base Geometry Reference_geometries Matrix)