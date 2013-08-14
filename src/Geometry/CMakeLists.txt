add_subdirectory(Mappings)

include_directories(. ${hpGEM_SOURCE_DIR}/src/)

add_library(Reference_geometries SHARED
	${hpGEM_SOURCE_DIR}/src/Geometry/ReferencePoint.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceLine.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceTriangle.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceSquare.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceTetrahedron.cpp 
        ${hpGEM_SOURCE_DIR}/src/Geometry/ReferencePyramid.cpp 
        ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceTriangularPrism.cpp 
        ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceCube.cpp 
        ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceHypercube.cpp 
	
	# Mappings

        ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefFaceToPyramid.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefFaceToTriangularPrism.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefTriangleToTetrahedron.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefSquareToSquare.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefSquareToCube.cpp
	${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefCubeToCube.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefTriangleToTriangle.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefLineToSquare.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefLineToTriangle.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefLineToLine.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefPointToLine.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefCubeToHypercube.cpp)



add_library(Geometry SHARED

        ${hpGEM_SOURCE_DIR}/src/Geometry/ElementGeometry.cpp
        # Physical geometries
        ${hpGEM_SOURCE_DIR}/src/Geometry/PhysicalLine.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/PhysicalTriangle.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/PhysicalQuadrilateral.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/PhysicalPyramid.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/PhysicalTriangularPrism.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/PhysicalTetrahedron.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/PhysicalHexahedron.cpp
#         ${hpGEM_SOURCE_DIR}/src/Geometry/PhysicalOctachoron.cpp

        # Mappings
        ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToPhysHypercubeLinear.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToPhysPyramid.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefCubeToHypercube.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToPhysTriangularPrism.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/OutwardNormalVectorSign.cpp

        # Refinements
        ${hpGEM_SOURCE_DIR}/src/Geometry/RefinementLine.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/RefinementTriangle.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/RefinementQuadrilateral.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/RefinementTriangularPrism.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/RefinementHexahedron.cpp

      )

target_link_libraries(Reference_geometries Matrix)
target_link_libraries(Geometry Reference_geometries)

