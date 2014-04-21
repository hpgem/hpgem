include_directories(${hpGEM_SOURCE_DIR}/kernel/ .)

add_library(Output SHARED
	${hpGEM_SOURCE_DIR}/kernel/Output/TecplotPhysicalGeometryIterator.cpp		${hpGEM_SOURCE_DIR}/kernel/Output/TecplotDiscontinuousSolutionWriter.cpp
	)

target_link_libraries(Output hpGEM_Base QuadratureRules)
