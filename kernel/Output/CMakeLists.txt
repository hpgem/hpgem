include_directories(${hpGEM_SOURCE_DIR}/kernel/ .)

add_library(Output OBJECT
	${hpGEM_SOURCE_DIR}/kernel/Output/TecplotPhysicalGeometryIterator.cpp
    ${hpGEM_SOURCE_DIR}/kernel/Output/base64.cpp
	VTKLagrangeCurve.cpp
	VTKLagrangeTriangle.cpp
	VTKStandardElements.cpp
)

target_link_libraries(Output)
set_target_properties(Output PROPERTIES POSITION_INDEPENDENT_CODE true)
