add_subdirectory(Mappings)

include_directories(. ${hpGEM_SOURCE_DIR}/kernel/)

add_library(Reference_geometries STATIC
	${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferencePoint.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceLine.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceTriangle.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceSquare.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceTetrahedron.cpp 
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferencePyramid.cpp 
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceTriangularPrism.cpp 
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceCube.cpp 
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceHypercube.cpp 
	
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceGeometry.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/PointReferenceBase.cpp
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
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefPointToPoint.cpp)

set_target_properties(Reference_geometries PROPERTIES POSITION_INDEPENDENT_CODE true)


add_library(Geometry SHARED

        ${hpGEM_SOURCE_DIR}/kernel/Geometry/ElementGeometry.cpp
 		${hpGEM_SOURCE_DIR}/kernel/Geometry/FaceGeometry.cpp
        # Physical geometries
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/PhysicalLine.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/PhysicalTriangle.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/PhysicalQuadrilateral.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/PhysicalPyramid.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/PhysicalTriangularPrism.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/PhysicalTetrahedron.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/PhysicalHexahedron.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/PhysicalOctachoron.cpp

        # Mappings
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToPhysHypercubeLinear.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToPhysPyramid.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToPhysTriangularPrism.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/OutwardNormalVectorSign.cpp

        # Refinements
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/RefinementLine.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/RefinementTriangle.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/RefinementQuadrilateral.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/RefinementTriangularPrism.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Geometry/RefinementHexahedron.cpp

      )

target_link_libraries(Reference_geometries LinearAlgebra QuadratureRules)
target_link_libraries(Geometry Reference_geometries)
