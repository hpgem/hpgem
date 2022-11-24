
################
#### LOGGER ####
################

# Shared logger and other base stuff
add_library(Logger STATIC Logger.cpp)

set_target_properties(Logger PROPERTIES POSITION_INDEPENDENT_CODE true)

if( hpGEM_USE_MPI )
  target_link_libraries(Logger PRIVATE MPI::MPI_CXX)
endif()

### Exposing CMAKE and build info ###

# Extract GIT hash
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
git_local_changes(GIT_DIRTY)

# Template cmake variables, e.g. the SHA
configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/hpgem-cmake.cpp.in
        ${CMAKE_CURRENT_BINARY_DIR}/hpgem-cmake.cpp ESCAPE_QUOTES)
target_sources(Logger PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/hpgem-cmake.cpp)
#Needed so that the hpgem-cmake.cpp can find the header files
target_include_directories(Logger PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})


##################
#### SUB LIBS ####
##################

add_subdirectory(LinearAlgebra)
add_subdirectory(EigenSolvers)
add_subdirectory(Geometry)
add_subdirectory(FE)
add_subdirectory(Integration)
add_subdirectory(Base)
add_subdirectory(Output)

################################
#### COMBINE INTO HPGEM LIB ####
################################


# Target which assembles the full kernel
# dummy.cpp needed as a library always needs at least one source file
add_library(HPGEM STATIC dummy.cpp)

# Include all parts, this unfortunately includes all source files, not just the headers
target_include_directories(HPGEM PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
# Link to all components
target_link_libraries(HPGEM
        FE
        Geometry
        hpGEM_Base
        LinearAlgebra
        Logger
        Output
        QuadratureRules
        Reference_geometries
        TimeIntegration
)

if(hpGEM_USE_SLEPC AND hpGEM_USE_MPI)
        target_link_libraries(HPGEM
                EigenSolvers)
endif()

add_library(HPGEM::HPGEM ALIAS HPGEM)
