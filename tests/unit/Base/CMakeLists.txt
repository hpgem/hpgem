include_directories(${hpGEM_SOURCE_DIR}/src/ ${hpGEM_BINARY_DIR}/conf .)


#Mesh Generator and readers tests
#################################


add_executable(RectangularMeshGeneratorUnitTest.out 
        RectangularMeshGeneratorUnitTest.cpp
        # Miscellaneous
        ${hpGEM_SOURCE_DIR}/src/Base/AssembleBasisFunctionSet.cpp
        ${hpGEM_SOURCE_DIR}/src/Base/MeshManipulator.cpp
        ${hpGEM_SOURCE_DIR}/src/Base/L2Norm.cpp
        )
target_link_libraries(RectangularMeshGeneratorUnitTest.out Geometry QuadratureRules)
add_test(RectangularMeshGeneratorUnitTest RectangularMeshGeneratorUnitTest.out)

add_executable(CentaurMeshReaderUnitTest.out
        CentaurMeshReaderUnitTest.cpp

        # Miscellaneous
        ${hpGEM_SOURCE_DIR}/src/Base/AssembleBasisFunctionSet.cpp
        ${hpGEM_SOURCE_DIR}/src/Base/L2Norm.cpp
        ${hpGEM_SOURCE_DIR}/src/Base/MeshManipulator.cpp
        ${hpGEM_SOURCE_DIR}/src/Output/TecplotPhysicalGeometryIterator.cpp
        )
target_link_libraries(CentaurMeshReaderUnitTest.out QuadratureRules Geometry)
add_test(CentaurMeshReaderUnitTest CentaurMeshReaderUnitTest.out)

#Basis functions and data structors
###################################

add_executable(LevelTreeUnitTest.out LevelTreeUnitTest.cpp)
add_test(LevelTreeUnitTest LevelTreeUnitTest.out)

add_executable(BasisFunctionsUnitTest.out 
               BasisFunctionsUnitTest.cpp 
	       ${hpGEM_SOURCE_DIR}/src/Base/AssembleBasisFunctionSet.cpp
           ${hpGEM_SOURCE_DIR}/src/Base/MeshManipulator.cpp)
add_test(BasisFunctionsUnitTest BasisFunctionsUnitTest.out)

add_executable(ElementUnitTest.out
               ElementUnitTest.cpp
               ${hpGEM_SOURCE_DIR}/src/Base/ElementData.hpp
               ${hpGEM_SOURCE_DIR}/src/Base/AssembleBasisFunctionSet.cpp
               ${hpGEM_SOURCE_DIR}/src/Base/MeshManipulator.cpp)
target_link_libraries(ElementUnitTest.out QuadratureRules Geometry)
add_test(ElementUnitTest ElementUnitTest.out)

