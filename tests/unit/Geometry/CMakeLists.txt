include_directories(${hpGEM_SOURCE_DIR}/kernel/ .)


#First declare the rules that will form the executable tests#
#############################################################

#Define Polytope_UnitTest.out dependencies

#-----REFERNECE-GEOMETRIES-----REFERNECE-GEOMETRIES-----REFERNECE-GEOMETRIES-----REFERNECE-GEOMETRIES
#0-d#
add_executable(ReferencePointUnitTest.out
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferencePoint.cpp
     		   ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefPointToPoint.cpp
               ReferencePointUnitTest.cpp)
target_link_libraries(ReferencePointUnitTest.out Reference_geometries)
add_test(ReferencePointUnitTest ReferencePointUnitTest.out)

#1-d#
add_executable(ReferenceLineUnitTest.out
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferencePoint.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceLine.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefPointToLine.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefLineToLine.cpp
      		   ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefPointToPoint.cpp
               ReferenceLineUnitTest.cpp)
target_link_libraries(ReferenceLineUnitTest.out Reference_geometries Matrix)
add_test(ReferenceLineUnitTest ReferenceLineUnitTest.out)

#2-d#
add_executable(ReferenceTriangleUnitTest.out
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferencePoint.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceLine.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceTriangle.cpp
    	       ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefPointToPoint.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefPointToLine.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefLineToLine.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefLineToTriangle.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefTriangleToTriangle.cpp
               ReferenceTriangleUnitTest.cpp)
target_link_libraries(ReferenceTriangleUnitTest.out Reference_geometries Matrix)
add_test(ReferenceTriangleUnitTest ReferenceTriangleUnitTest.out)

add_executable(ReferenceSquareUnitTest.out
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferencePoint.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceLine.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceSquare.cpp
     		   ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefPointToPoint.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefPointToLine.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefLineToLine.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefLineToSquare.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefSquareToSquare.cpp
               ReferenceSquareUnitTest.cpp)
target_link_libraries(ReferenceSquareUnitTest.out Reference_geometries Matrix)
add_test(ReferenceSquareUnitTest ReferenceSquareUnitTest.out)

#3-d#
add_executable(ReferenceTetrahedronUnitTest.out
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferencePoint.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceLine.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceTriangle.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceTetrahedron.cpp
     		   ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefPointToPoint.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefPointToLine.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefLineToLine.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefLineToTriangle.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefTriangleToTriangle.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefTriangleToTetrahedron.cpp
               ReferenceTetrahedronUnitTest.cpp)
target_link_libraries(ReferenceTetrahedronUnitTest.out Reference_geometries Matrix)
add_test(ReferenceTetrahedronUnitTest ReferenceTetrahedronUnitTest.out)

add_executable(ReferencePyramidUnitTest.out
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferencePoint.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceLine.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceTriangle.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceSquare.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferencePyramid.cpp
     		   ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefPointToPoint.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefPointToLine.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefLineToLine.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefLineToTriangle.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefLineToSquare.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefTriangleToTriangle.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefSquareToSquare.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefFaceToPyramid.cpp
               ReferencePyramidUnitTest.cpp)
target_link_libraries(ReferencePyramidUnitTest.out Reference_geometries Matrix)
add_test(ReferencePyramidUnitTest ReferencePyramidUnitTest.out)

add_executable(ReferenceTriangularPrismUnitTest.out
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferencePoint.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceLine.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceTriangle.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceSquare.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceTriangularPrism.cpp
     		   ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefPointToPoint.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefPointToLine.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefLineToLine.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefLineToTriangle.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefLineToSquare.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefTriangleToTriangle.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefSquareToSquare.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefFaceToTriangularPrism.cpp
               ReferenceTriangularPrismUnitTest.cpp)
target_link_libraries(ReferenceTriangularPrismUnitTest.out Reference_geometries Matrix)        
add_test(ReferenceTriangularPrismUnitTest ReferenceTriangularPrismUnitTest.out)

add_executable(ReferenceCubeUnitTest.out
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferencePoint.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceLine.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceSquare.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceCube.cpp
     		   ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefPointToPoint.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefPointToLine.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefLineToLine.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefLineToSquare.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefSquareToSquare.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefSquareToCube.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefCubeToCube.cpp
               ReferenceCubeUnitTest.cpp)
target_link_libraries(ReferenceCubeUnitTest.out Reference_geometries Matrix)
add_test(ReferenceCubeUnitTest ReferenceCubeUnitTest.out)

#4-d#
add_executable(ReferenceHypercubeUnitTest.out 
               ReferenceHypercubeUnitTest.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferencePoint.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceLine.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceSquare.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceCube.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceHypercube.cpp
     		   ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefPointToPoint.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefPointToLine.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefLineToLine.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefLineToSquare.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefSquareToSquare.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefSquareToCube.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefCubeToCube.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefCubeToHypercube.cpp
               ReferenceHypercubeUnitTest.cpp)
target_link_libraries(ReferenceHypercubeUnitTest.out Reference_geometries Matrix)
add_test(ReferenceHypercubeUnitTest ReferenceHypercubeUnitTest.out)


#---------POINTS---------POINTS---------POINTS---------POINTS---------POINTS---------POINTS---------

add_executable(PointUnitTest.out
               PointUnitTest.cpp
               #${hpGEM_SOURCE_DIR}/kernel/Geometry/Point.hpp //headers go in source files
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Point.cpp
               )
target_link_libraries(PointUnitTest.out Matrix)
add_test(PointUnitTest PointUnitTest.out)

add_executable(PointPhysicalUnitTest.out
               PointPhysicalUnitTest.cpp
               #${hpGEM_SOURCE_DIR}/kernel/Geometry/Point.hpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Point.cpp
               #${hpGEM_SOURCE_DIR}/kernel/Geometry/PointPhysical.hpp
               )
       target_link_libraries(PointPhysicalUnitTest.out Matrix)
       
add_test(PointPhysicalUnitTest PointPhysicalUnitTest.out)


add_executable(PointReferenceUnitTest.out
               PointReferenceUnitTest.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/PointReference.hpp
               #${hpGEM_SOURCE_DIR}/kernel/Geometry/Point.hpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Point.cpp
               )
       target_link_libraries(PointReferenceUnitTest.out Matrix)

add_test(PointReferenceUnitTest PointReferenceUnitTest.out)

#--------ELEMENT--------ELEMENTS--------ELEMENTS--------ELEMENTS--------ELEMENTS--------ELEMENTS----

#add_executable(ElementGeometry_UnitTest.out
#               ElementGeometry_UnitTest.cpp
#               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingHypercubeLinear.cpp
#               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingReferenceToPhysical.hpp
#               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingInterface.hpp
#               ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceSquare.cpp
#               ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceSquare.hpp
#               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Jacobian.hpp)

#--------MAPPINGS--------MAPPINGS--------MAPPINGS--------MAPPINGS--------MAPPINGS--------MAPPINGS--
               
add_executable(MappingHypercubeLinear_UnitTest.out
               MappingHypercubeLinear_UnitTest.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferencePoint.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceLine.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/PhysicalLine.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceSquare.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/PhysicalQuadrilateral.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceCube.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/PhysicalHexahedron.cpp
     		   ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefPointToPoint.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefPointToLine.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefLineToLine.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefLineToSquare.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefSquareToSquare.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefSquareToCube.cpp
	       ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefCubeToCube.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToPhysHypercubeLinear.cpp
               )
target_link_libraries(MappingHypercubeLinear_UnitTest.out Reference_geometries Matrix)
               
add_executable(MappingSimplexLinear_UnitTest.out
               MappingSimplexLinear_UnitTest.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferencePoint.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceLine.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceTriangle.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceSquare.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/PhysicalTriangle.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/ReferenceCube.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/PhysicalTetrahedron.cpp
     		   ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefPointToPoint.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefPointToLine.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefLineToLine.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefLineToTriangle.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefLineToSquare.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefTriangleToTriangle.cpp
	       ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefCubeToCube.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefSquareToSquare.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefSquareToCube.cpp
               )
target_link_libraries(MappingSimplexLinear_UnitTest.out Reference_geometries Matrix)

add_executable(MappingPointToLine_UnitTest.out
               MappingPointToLine_UnitTest.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Point.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefPointToLine.cpp
               )
target_link_libraries(MappingPointToLine_UnitTest.out Matrix)
               
add_executable(MappingLineToSquare_UnitTest.out
               MappingLineToSquare_UnitTest.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Point.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/PointReference.hpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefLineToSquare.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefLineToSquare.hpp
               )
target_link_libraries(MappingLineToSquare_UnitTest.out Matrix)

add_executable(MappingSquareToCube_UnitTest.out
               MappingSquareToCube_UnitTest.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/PointReference.hpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Point.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefSquareToCube.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Geometry/Mappings/MappingToRefSquareToCube.hpp)
target_link_libraries(MappingSquareToCube_UnitTest.out Matrix)
                              
add_executable(FaceGeometry.out
               FaceGeometryUnitTest.cpp)                
target_link_libraries(FaceGeometry.out Geometry)



