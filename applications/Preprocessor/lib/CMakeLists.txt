add_library(Preprocessor-lib STATIC
        centaur.cpp
        hpgem.cpp
        gmsh.cpp
        ElementReorder.cpp
        FortranUnformattedFile.cpp

        mesh/ElementShapes.cpp
        utils/SparseUnionFind.cpp
        )
#target_include_directories(Preprocessor-lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/lib)
target_link_libraries(Preprocessor-lib PUBLIC HPGEM::HPGEM)
target_include_directories(Preprocessor-lib
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        )