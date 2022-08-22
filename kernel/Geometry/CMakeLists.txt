add_subdirectory(Mappings)

include_directories(. ${hpGEM_SOURCE_DIR}/kernel/)

add_library(Reference_geometries OBJECT
	${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferencePoint.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceLine.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceTriangle.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceSquare.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceTetrahedron.cpp 
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferencePyramid.cpp 
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceTriangularPrism.cpp 
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceCube.cpp 
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceHypercube.cpp
		ReferenceCurvilinearLine.cpp
		ReferenceCurvilinearTriangle.cpp
		ReferenceCurvilinearTetrahedron.cpp

		ReferenceGeometryFactory.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceGeometry.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/PointPhysicalBase.cpp
	# Mappings

        ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefCubeToHypercube.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefFaceToPyramid.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefFaceToTriangularPrism.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefTriangleToTetrahedron.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefSquareToSquare.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefSquareToCube.cpp
	${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefCubeToCube.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefTriangleToTriangle.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefLineToSquare.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefLineToTriangle.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefLineToLine.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefPointToLine.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/ConcatenatedMapping.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefPointToPoint.cpp
		Mappings/MappingRefLineToTetrahedron.cpp
	)

set_target_properties(Reference_geometries PROPERTIES POSITION_INDEPENDENT_CODE true)


add_library(Geometry OBJECT

        ${hpGEM_SOURCE_DIR}/kernel/Geometry/ElementGeometry.cpp
 		${hpGEM_SOURCE_DIR}/kernel/Geometry/FaceGeometry.cpp
        # Physical geometries
		${hpGEM_SOURCE_DIR}/kernel/Geometry/PhysicalGeometryBase.cpp

        # Mappings
		Mappings/MappingToPhysTetrahedronQuadratic.cpp
        Mappings/MappingToPhysTriangleQuadratic.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToPhysHypercubeLinear.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToPhysPyramid.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToPhysTriangularPrism.cpp
      )

target_link_libraries(Reference_geometries
	PUBLIC
		LinearAlgebra)
target_link_libraries(Geometry PUBLIC Reference_geometries)
