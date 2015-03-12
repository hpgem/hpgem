include_directories(${hpGEM_SOURCE_DIR}/kernel/ .)

add_library(Output STATIC
	${hpGEM_SOURCE_DIR}/kernel/Output/TecplotPhysicalGeometryIterator.cpp		
        ${hpGEM_SOURCE_DIR}/kernel/Output/TecplotDiscontinuousSolutionWriter.cpp
	${hpGEM_SOURCE_DIR}/kernel/Output/GNUPlotDiscontinuousSolutionWriter.cpp
    ${hpGEM_SOURCE_DIR}/kernel/Output/VTKSpecificTimeWriter.cpp
    ${hpGEM_SOURCE_DIR}/kernel/Output/VTKTimeDependentWriter.cpp
    ${hpGEM_SOURCE_DIR}/kernel/Output/base64.cpp
	)

target_link_libraries(Output hpGEM_Base)
set_target_properties(Output PROPERTIES POSITION_INDEPENDENT_CODE true)
