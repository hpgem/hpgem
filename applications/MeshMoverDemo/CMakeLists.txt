include_directories(. ${hpGEM_SOURCE_DIR}/kernel/)

add_executable(MeshMoverDemo.out 
            MeshMoverDemo.cpp
            ${hpGEM_SOURCE_DIR}/kernel/Base/MeshMoverBase.hpp
            ${hpGEM_SOURCE_DIR}/kernel/Base/AssembleBasisFunctionSet.cpp
            ${hpGEM_SOURCE_DIR}/kernel/Output/TecplotPhysicalGeometryIterator.cpp
            ${hpGEM_SOURCE_DIR}/kernel/Output/TecplotDiscontinuousSolutionWriter.cpp
  
	    #miscellaneous required stuff
	    ${hpGEM_SOURCE_DIR}/kernel/Base/AssembleBasisFunctionSet.cpp
	    ${hpGEM_SOURCE_DIR}/kernel/Base/L2Norm.cpp
	    )

target_link_libraries(MeshMoverDemo.out hpGEM_Base Geometry Matrix QuadratureRules)                
