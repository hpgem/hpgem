include_directories(${hpGEM_SOURCE_DIR}/src/ .)

add_executable(TecplotOutputUnitTest.out
               TecplotOutputUnitTest.cpp
               ${hpGEM_SOURCE_DIR}/src/Output/TecplotPhysicalGeometryIterator.cpp
               ${hpGEM_SOURCE_DIR}/src/Output/TecplotDiscontinuousSolutionWriter.cpp
               ${hpGEM_SOURCE_DIR}/src/Base/AssembleBasisFunctionSet.cpp)
target_link_libraries(TecplotOutputUnitTest.out hpGEM_Base Geometry QuadratureRules)
add_test(TecplotOutputUnitTest TecplotOutputUnitTest.out)
