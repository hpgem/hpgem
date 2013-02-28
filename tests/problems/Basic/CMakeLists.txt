include_directories(. ${hpGEM_SOURCE_DIR}/src/)

add_executable(Problem.out 
            Problem.cpp
            ${hpGEM_SOURCE_DIR}/src/Base/MeshMoverBase.hpp
            ${hpGEM_SOURCE_DIR}/src/Base/AssembleBasisFunctionSet.cpp
            ${hpGEM_SOURCE_DIR}/src/Output/TecplotPhysicalGeometryIterator.cpp
          )
target_link_libraries(Problem.out Geometry Matrix)                