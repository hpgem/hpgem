
add_subdirectory(lib)

add_executable(Preprocessor.out main.cpp)
target_link_libraries(Preprocessor.out Preprocessor-lib)
add_executable(ExtractSurface.out extractSurface.cpp)
target_link_libraries(ExtractSurface.out Preprocessor-lib)
