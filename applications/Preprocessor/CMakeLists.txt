
add_library(Preprocessor-lib STATIC
lib/centaur.cpp
lib/hpgem.cpp
lib/gmsh.cpp
lib/ElementReorder.cpp
)

# TODO: Lacking source/include separation
target_include_directories(Preprocessor-lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/lib)
target_link_libraries(Preprocessor-lib PUBLIC HPGEM::HPGEM)


#the preprocessor does not perform floating poing arithmetic, so it doesn't need special Petsc/not petsc guards
add_executable(Preprocessor.out main.cpp)
target_link_libraries(Preprocessor.out Preprocessor-lib)
add_executable(ExtractSurface.out extractSurface.cpp)
target_link_libraries(ExtractSurface.out Preprocessor-lib)
