include_directories(${hpGEM_SOURCE_DIR}/kernel/ ${hpGEM_BINARY_DIR}/conf .)


#Mesh Generator and readers tests
#################################


add_executable(RectangularMeshGeneratorUnitTest.out 
        RectangularMeshGeneratorUnitTest.cpp
        # Miscellaneous
        ${hpGEM_SOURCE_DIR}/kernel/Base/AssembleBasisFunctionSet.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Base/L2Norm.cpp
        )
target_link_libraries(RectangularMeshGeneratorUnitTest.out hpGEM_Base Geometry QuadratureRules)
add_test(RectangularMeshGeneratorUnitTest RectangularMeshGeneratorUnitTest.out)

add_executable(CentaurMeshReaderUnitTest.out
        CentaurMeshReaderUnitTest.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Base/MeshManipulator.cpp
        # Miscellaneous
        ${hpGEM_SOURCE_DIR}/kernel/Base/AssembleBasisFunctionSet.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Base/L2Norm.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Output/TecplotPhysicalGeometryIterator.cpp
        )
target_link_libraries(CentaurMeshReaderUnitTest.out hpGEM_Base QuadratureRules Geometry)
add_test(CentaurMeshReaderUnitTest CentaurMeshReaderUnitTest.out)

#Basis functions and data structors
###################################

add_executable(LevelTreeUnitTest.out LevelTreeUnitTest.cpp)
add_test(LevelTreeUnitTest LevelTreeUnitTest.out)

add_executable(BasisFunctionsUnitTest.out 
               BasisFunctionsUnitTest.cpp 
	       ${hpGEM_SOURCE_DIR}/kernel/Base/AssembleBasisFunctionSet.cpp
	       ${hpGEM_SOURCE_DIR}/kernel/Base/BasisFunctionSet.cpp
	       ${hpGEM_SOURCE_DIR}/kernel/Geometry/Point.cpp
           #${hpGEM_SOURCE_DIR}/kernel/Base/MeshManipulator.cpp
           )
target_link_libraries(BasisFunctionsUnitTest.out Matrix)
add_test(BasisFunctionsUnitTest BasisFunctionsUnitTest.out)

add_executable(ElementUnitTest.out
               ElementUnitTest.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Base/ElementData.hpp
               ${hpGEM_SOURCE_DIR}/kernel/Base/AssembleBasisFunctionSet.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Base/MeshManipulator.cpp)
target_link_libraries(ElementUnitTest.out hpGEM_Base QuadratureRules Geometry)
add_test(ElementUnitTest ElementUnitTest.out)

