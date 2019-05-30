add_library(Logger SHARED Logger.cpp)
set_target_properties(Logger PROPERTIES POSITION_INDEPENDENT_CODE true)

if( hpGEM_USE_MPI )
  target_link_libraries(Logger "${MPI_LIBRARIES}")
endif()

add_subdirectory(LinearAlgebra)
add_subdirectory(Geometry)
add_subdirectory(Integration)
add_subdirectory(Base)
add_subdirectory(Output)
