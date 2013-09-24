include_directories(${hpGEM_SOURCE_DIR}/src/ .)

add_library(hpGEM_Base SHARED
		${hpGEM_SOURCE_DIR}/src/Output/TecplotPhysicalGeometryIterator.cpp				  
		#miscellaneous required stuff
		${hpGEM_SOURCE_DIR}/src/Base/MeshManipulator.cpp
        	${hpGEM_SOURCE_DIR}/src/Base/AssembleBasisFunctionSet.cpp
		${hpGEM_SOURCE_DIR}/src/Base/L2Norm.cpp
 		${hpGEM_SOURCE_DIR}/src/Base/Norm2.cpp
	)

