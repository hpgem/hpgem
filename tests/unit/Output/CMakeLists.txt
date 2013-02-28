include_directories(${hpGEM_SOURCE_DIR}/src/ .)

add_executable(TecplotOutputTest.out
               TecplotOutputTest.cpp
               ${hpGEM_SOURCE_DIR}/src/Output/TecplotPhysicalGeometryIterator.cpp
               ${hpGEM_SOURCE_DIR}/src/Base/AssembleBasisFunctionSet.cpp)
       target_link_libraries(TecplotOutputTest.out Geometry)

