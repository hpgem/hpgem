
add_library(Preprocessor-lib STATIC
lib/centaur.cpp
lib/hpgem.cpp
lib/gmsh.cpp
lib/ElementReorder.cpp
lib/SparseUnionFind.cpp
lib/FortranUnformattedFile.cpp
lib/ElementShapes.cpp
)

# TODO: Lacking source/include separation
target_include_directories(Preprocessor-lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/lib)
target_link_libraries(Preprocessor-lib PUBLIC HPGEM::HPGEM)

add_executable(Preprocessor.out main.cpp)
target_link_libraries(Preprocessor.out Preprocessor-lib)
add_executable(ExtractSurface.out extractSurface.cpp)
target_link_libraries(ExtractSurface.out Preprocessor-lib)
