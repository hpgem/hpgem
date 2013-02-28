include_directories(${hpGEM_SOURCE_DIR}/src/ .)


#First declare the rules that will form the executable tests#
#############################################################

#Define Polytope_UnitTest.out dependencies

#-----REFERNECE-GEOMETRIES-----REFERNECE-GEOMETRIES-----REFERNECE-GEOMETRIES-----REFERNECE-GEOMETRIE

add_executable(ReferencePointUnitTest.out
               ${hpGEM_SOURCE_DIR}/src/Geometry/ReferencePoint.cpp
               ReferencePointUnitTest.cpp)

add_executable(ReferenceLineUnitTest.out
               ${hpGEM_SOURCE_DIR}/src/Geometry/ReferencePoint.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceLine.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefPointToLine.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefLineToLine.cpp
               ReferenceLineUnitTest.cpp)
target_link_libraries(ReferenceLineUnitTest.out Matrix)


add_executable(ReferenceTriangleUnitTest.out
               ${hpGEM_SOURCE_DIR}/src/Geometry/ReferencePoint.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceLine.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceTriangle.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefPointToLine.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefLineToLine.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefLineToTriangle.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefTriangleToTriangle.cpp
               ReferenceTriangleUnitTest.cpp)
target_link_libraries(ReferenceTriangleUnitTest.out Matrix)

add_executable(ReferenceSquareUnitTest.out
               ${hpGEM_SOURCE_DIR}/src/Geometry/ReferencePoint.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceLine.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceSquare.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefPointToLine.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefLineToLine.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefLineToSquare.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefSquareToSquare.cpp
               ReferenceSquareUnitTest.cpp)
target_link_libraries(ReferenceSquareUnitTest.out Matrix)

add_executable(ReferenceTetrahedronUnitTest.out
               ${hpGEM_SOURCE_DIR}/src/Geometry/ReferencePoint.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceLine.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceTriangle.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceTetrahedron.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefPointToLine.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefLineToLine.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefLineToTriangle.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefTriangleToTriangle.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefTriangleToTetrahedron.cpp
               ReferenceTetrahedronUnitTest.cpp)
target_link_libraries(ReferenceTetrahedronUnitTest.out Matrix)

add_executable(ReferenceCubeUnitTest.out
               ${hpGEM_SOURCE_DIR}/src/Geometry/ReferencePoint.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceLine.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceSquare.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceCube.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefPointToLine.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefLineToLine.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefLineToSquare.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefSquareToSquare.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefSquareToCube.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefCubeToCube.cpp
               ReferenceCubeUnitTest.cpp)
target_link_libraries(ReferenceCubeUnitTest.out Matrix)
               
add_executable(ReferencePyramidUnitTest.out
               ${hpGEM_SOURCE_DIR}/src/Geometry/ReferencePoint.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceLine.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceTriangle.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceSquare.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/ReferencePyramid.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefPointToLine.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefLineToLine.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefLineToTriangle.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefLineToSquare.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefTriangleToTriangle.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefSquareToSquare.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefFaceToPyramid.cpp
               ReferencePyramidUnitTest.cpp)
target_link_libraries(ReferencePyramidUnitTest.out Matrix)

add_executable(ReferenceTriangularPrismUnitTest.out
               ${hpGEM_SOURCE_DIR}/src/Geometry/ReferencePoint.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceLine.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceTriangle.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceSquare.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceTriangularPrism.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefPointToLine.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefLineToLine.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefLineToTriangle.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefLineToSquare.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefTriangleToTriangle.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefSquareToSquare.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefFaceToTriangularPrism.cpp
               ReferenceTriangularPrismUnitTest.cpp)
target_link_libraries(ReferenceTriangularPrismUnitTest.out Matrix)        

add_executable(ReferenceHypercubeUnitTest.out 
               ReferenceHypercubeUnitTest.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/ReferencePoint.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceLine.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceSquare.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceCube.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceHypercube.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefPointToLine.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefLineToLine.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefLineToSquare.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefSquareToSquare.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefSquareToCube.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefCubeToCube.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefCubeToHypercube.cpp
               ReferenceHypercubeUnitTest.cpp)
target_link_libraries(ReferenceHypercubeUnitTest.out Matrix)

#add_executable(TriangularPrismTest.out
#                ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceTriangularPrism.cpp
#                ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceTriangularPrism.hpp
#                ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceGeometry.hpp
#                ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceGeometry_Impl.hpp
#                ${hpGEM_SOURCE_DIR}/src/Geometry/Point.hpp
#                ${hpGEM_SOURCE_DIR}/src/Geometry/Point_Impl.hpp
#                ReferenceTriangularPrismUnitTest.cpp)

#---------POINTS---------POINTS---------POINTS---------POINTS---------POINTS---------POINTS---------

add_executable(Point_UnitTest.out
               Point_UnitTest.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Point.hpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Point_Impl.hpp)

add_executable(PointPhysical_UnitTest.out
               PointPhysical_UnitTest.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Point.hpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/PointPhysical.hpp)

add_executable(PointReference_UnitTest.out
               PointReference_UnitTest.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/PointReference.hpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Point.hpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Point_Impl.hpp)

#--------ELEMENT--------ELEMENTS--------ELEMENTS--------ELEMENTS--------ELEMENTS--------ELEMENTS----

#add_executable(ElementGeometry_UnitTest.out
#               ElementGeometry_UnitTest.cpp
#               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingHypercubeLinear.cpp
#               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingReferenceToPhysical.hpp
#               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingInterface.hpp
#               ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceSquare.cpp
#               ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceSquare.hpp
#               ${hpGEM_SOURCE_DIR}/src/Geometry/Jacobian.hpp)

#--------MAPPINGS--------MAPPINGS--------MAPPINGS--------MAPPINGS--------MAPPINGS--------MAPPINGS--
               
add_executable(MappingHypercubeLinear_UnitTest.out
               MappingHypercubeLinear_UnitTest.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/ReferencePoint.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceLine.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/PhysicalLine.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceSquare.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/PhysicalQuadrilateral.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceCube.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/PhysicalHexahedron.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefPointToLine.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefLineToLine.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefLineToSquare.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefSquareToSquare.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefSquareToCube.cpp
	       ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefCubeToCube.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToPhysHypercubeLinear.cpp
               )
target_link_libraries(MappingHypercubeLinear_UnitTest.out Matrix)
               
add_executable(MappingSimplexLinear_UnitTest.out
               MappingSimplexLinear_UnitTest.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/ReferencePoint.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceLine.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceTriangle.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceSquare.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/PhysicalTriangle.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceCube.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/PhysicalTetrahedron.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefPointToLine.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefLineToLine.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefLineToTriangle.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefLineToSquare.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefTriangleToTriangle.cpp
	       ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefCubeToCube.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefSquareToSquare.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefSquareToCube.cpp
               )
target_link_libraries(MappingSimplexLinear_UnitTest.out Matrix)

add_executable(MappingPointToLine_UnitTest.out
               MappingPointToLine_UnitTest.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefPointToLine.cpp
               )
target_link_libraries(MappingPointToLine_UnitTest.out Matrix)
               
add_executable(MappingLineToSquare_UnitTest.out
               MappingLineToSquare_UnitTest.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/PointReference.hpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefLineToSquare.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefLineToSquare.hpp
               )
target_link_libraries(MappingLineToSquare_UnitTest.out Matrix)

add_executable(MappingSquareToCube_UnitTest.out
               MappingSquareToCube_UnitTest.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/PointReference.hpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefSquareToCube.cpp
               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefSquareToCube.hpp)
target_link_libraries(MappingSquareToCube_UnitTest.out Matrix)
                              
add_executable(FaceGeometry.out
               FaceGeometryUnitTest.cpp)                
target_link_libraries(FaceGeometry.out Geometry)



#Now create the actually tests#
###############################
#Test -2 - Simple FaceGeometry Test : Check it both complies and runs
add_test(BuildTest "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target FaceGeometry.out)
add_test(FaceGeometry FaceGeometry.out)
set_tests_properties(FaceGeometry PROPERTIES DEPENDS BuildTest)


#Test 0 - Simple Hypercube Test : Check it both complies and runs
add_test(BuildTest "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target TetrahedronSampleTest.out)
add_test(TetrahedronSampleTest TetrahedronSampleTest.out)
set_tests_properties(TetrahedronSampleTest PROPERTIES DEPENDS BuildTest)

#Test -1 - Simple Polytope Test : Check it both complies and runs
add_test(BuildTest "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target PolytopeTest.out)
add_test(PolytopeTest PolytopeTest.out)
set_tests_properties(PolytopeTest PROPERTIES DEPENDS BuildTest)

#Test 1 - Simple Polytope Test : Check it both complies and runs
add_test(BuildTest "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target PolytopeTest.out)
add_test(PolytopeTest PolytopeTest.out)
set_tests_properties(PolytopeTest PROPERTIES DEPENDS BuildTest)

#Test 2 - Simple Line Test : Check it both complies and runs
add_test(BuildTest "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target LineTest.out)
add_test(LineTest LineTest.out)
set_tests_properties(LineTest PROPERTIES DEPENDS BuildTest)

#Test 3 - Simple Triangle Test : Check it both complies and runs
add_test(BuildTest "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target TriangleTest.out)
add_test(TriangleTest TriangleTest.out)
set_tests_properties(TriangleTest PROPERTIES DEPENDS BuildTest)

#Test 4 - Simple Square Test : Check it both complies and runs
add_test(BuildTest "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target SquareTest.out)
add_test(SquareTest SquareTest.out)
set_tests_properties(SquareTest PROPERTIES DEPENDS BuildTest)

#Test 5 - Simple Pyramid Test : Check it both complies and runs
add_test(BuildTest "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target PyramidTest.out)
add_test(PyramidTest PyramidTest.out)
set_tests_properties(PyramidTest PROPERTIES DEPENDS BuildTest)

#Test 6 - Simple TriangularPrism Test : Check it both complies and runs
add_test(BuildTest "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target TriangularPrismTest.out)
add_test(TriangularPrismTest TriangularPrismTest.out)
set_tests_properties(TriangularPrismTest PROPERTIES DEPENDS BuildTest)

#Test 7 - Simple Hypercube Test : Check it both complies and runs
add_test(BuildTest "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target HypercubeTest.out)
add_test(HypercubeTest HypercubeTest.out)
set_tests_properties(HypercubeTest PROPERTIES DEPENDS BuildTest)

#Test 8 - Simple Hypercube Test : Check it both complies and runs
add_test(BuildTest "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target TetrahedronTest.out)
add_test(TetrahedronTest TetrahedronTest.out)
set_tests_properties(TetrahedronTest PROPERTIES DEPENDS BuildTest)




