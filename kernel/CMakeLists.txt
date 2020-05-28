add_library(Logger SHARED Logger.cpp)
set_target_properties(Logger PROPERTIES POSITION_INDEPENDENT_CODE true)

if( hpGEM_USE_MPI )
  target_link_libraries(Logger PRIVATE MPI::MPI_CXX)
endif()

add_subdirectory(LinearAlgebra)
add_subdirectory(Geometry)
add_subdirectory(Integration)
add_subdirectory(Base)
add_subdirectory(Output)


# Target which assembles the full kernel
# dummy.cpp needed as a library always needs at least one source file
add_library(HPGEM STATIC dummy.cpp)

# Include all parts, this unfortunately includes all source files, not just the headers
target_include_directories(HPGEM PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
# Link to all components
target_link_libraries(HPGEM
        Geometry
        hpGEM_Base
        LinearAlgebra
        Logger
        Output
        QuadratureRules
        Reference_geometries
        TimeIntegration
)

add_library(HPGEM::HPGEM ALIAS HPGEM)
