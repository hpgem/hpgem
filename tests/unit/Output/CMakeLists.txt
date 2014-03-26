include_directories(${hpGEM_SOURCE_DIR}/kernel/ .)

add_executable(TecplotOutputUnitTest.out
               TecplotOutputUnitTest.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Output/TecplotPhysicalGeometryIterator.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Output/TecplotDiscontinuousSolutionWriter.cpp
               ${hpGEM_SOURCE_DIR}/kernel/Base/AssembleBasisFunctionSet.cpp)
target_link_libraries(TecplotOutputUnitTest.out hpGEM_Base Geometry QuadratureRules)
add_test(TecplotOutputUnitTest TecplotOutputUnitTest.out)
