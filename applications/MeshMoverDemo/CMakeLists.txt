include_directories(. ${hpGEM_SOURCE_DIR}/src/)

add_executable(MeshMoverDemo.out 
            MeshMoverDemo.cpp
            ${hpGEM_SOURCE_DIR}/src/Base/MeshMoverBase.hpp
            ${hpGEM_SOURCE_DIR}/src/Base/AssembleBasisFunctionSet.cpp
            ${hpGEM_SOURCE_DIR}/src/Output/TecplotPhysicalGeometryIterator.cpp
            ${hpGEM_SOURCE_DIR}/src/Output/TecplotDiscontinuousSolutionWriter.cpp
  
	    #miscellaneous required stuff
	    ${hpGEM_SOURCE_DIR}/src/Base/AssembleBasisFunctionSet.cpp
	    ${hpGEM_SOURCE_DIR}/src/Base/L2Norm.cpp
	    )

target_link_libraries(MeshMoverDemo.out hpGEM_Base Geometry Matrix QuadratureRules)                
