include_directories(${hpGEM_SOURCE_DIR}/src/ ${hpGEM_BINARY_DIR}/conf .)


#First declare the rules that will form the executable tests#
#############################################################

#Define MatrixTest.out dependencies
add_executable(RectangularMeshGeneratorTest.out 
        RectangularMeshGeneratorTest.cpp
        # Miscellaneous
        ${hpGEM_SOURCE_DIR}/src/Base/AssembleBasisFunctionSet.cpp
        ${hpGEM_SOURCE_DIR}/src/Base/L2Norm.cpp
        )
                
add_executable(CentaurMeshReaderTest.out
        CentaurMeshReaderTest.cpp

        # Reference geometries
        ${hpGEM_SOURCE_DIR}/src/Geometry/ReferencePoint.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceLine.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceTriangle.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceSquare.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceTetrahedron.cpp 
        ${hpGEM_SOURCE_DIR}/src/Geometry/ReferencePyramid.cpp 
        ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceTriangularPrism.cpp 
        ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceCube.cpp 
        ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceHypercube.cpp 

        # Physical geometries
        ${hpGEM_SOURCE_DIR}/src/Geometry/PhysicalLine.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/PhysicalTriangle.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/PhysicalQuadrilateral.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/PhysicalPyramid.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/PhysicalTriangularPrism.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/PhysicalTetrahedron.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/PhysicalHexahedron.cpp
#         ${hpGEM_SOURCE_DIR}/src/Geometry/PhysicalOctachoron.cpp

        # Gauss quadrature rules
#         ${hpGEM_SOURCE_DIR}/src/Integration/QuadratureRules/GaussQuadratureRulesForLine.cpp
#         ${hpGEM_SOURCE_DIR}/src/Integration/QuadratureRules/GaussQuadratureRulesForTriangle.cpp
#         ${hpGEM_SOURCE_DIR}/src/Integration/QuadratureRules/GaussQuadratureRulesForSquare.cpp
#         ${hpGEM_SOURCE_DIR}/src/Integration/QuadratureRules/GaussQuadratureRulesForTetrahedron.cpp
#         ${hpGEM_SOURCE_DIR}/src/Integration/QuadratureRules/GaussQuadratureRulesForPyramid.cpp
#         ${hpGEM_SOURCE_DIR}/src/Integration/QuadratureRules/GaussQuadratureRulesForTriangularPrism.cpp
#         ${hpGEM_SOURCE_DIR}/src/Integration/QuadratureRules/GaussQuadratureRulesForCube.cpp
#         ${hpGEM_SOURCE_DIR}/src/Integration/QuadratureRules/GaussQuadratureRulesForHypercube.cpp

        # Mappings
        ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefFaceToPyramid.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefFaceToTriangularPrism.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefTriangleToTetrahedron.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefSquareToSquare.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefSquareToCube.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefTriangleToTriangle.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefLineToSquare.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefLineToTriangle.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefLineToLine.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefPointToLine.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToPhysHypercubeLinear.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToPhysPyramid.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToRefCubeToHypercube.cpp
        ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToPhysTriangularPrism.cpp

        # Miscellaneous
        ${hpGEM_SOURCE_DIR}/src/Base/AssembleBasisFunctionSet.cpp
        ${hpGEM_SOURCE_DIR}/src/Base/L2Norm.cpp
        ${hpGEM_SOURCE_DIR}/src/Output/TecplotPhysicalGeometryIterator.cpp
        )


target_link_libraries(RectangularMeshGeneratorTest.out QuadratureRules Geometry)
target_link_libraries(CentaurMeshReaderTest.out QuadratureRules Geometry)

add_executable(LevelTreeTest.out LevelTreeUnitTest.cpp)

add_executable(BasisFunctionsTest.out 
               BasisFunctions_UnitTest.cpp 
               ${hpGEM_SOURCE_DIR}/src/Base/AssembleBasisFunctionSet.cpp)

#add_executable(Element_UnitTest.out
#               Element_UnitTest.cpp
#               ${hpGEM_SOURCE_DIR}/src/Base/ElementData.hpp
#               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToPhysHypercubeLinear.cpp
#               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingToPhysHypercubeLinear.hpp
#               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingReferenceToPhysical.hpp
#               ${hpGEM_SOURCE_DIR}/src/Geometry/Mappings/MappingInterface.hpp
#               ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceSquare.cpp
#               ${hpGEM_SOURCE_DIR}/src/Geometry/ReferenceSquare.hpp
#               ${hpGEM_SOURCE_DIR}/src/Geometry/Jacobian.hpp)

#Now create the actually tests#
###############################

#Test 1 - Simple Matrix Algebra Test : Check it both complies and runs
add_test(BuildTest "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target Element_UnitTest.out)
add_test(Element_UnitTest Element_UnitTest.out)
set_tests_properties(Element_UnitTest PROPERTIES DEPENDS BuildTest)
