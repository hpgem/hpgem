add_library(Logger SHARED Logger.cc)
set_target_properties(Logger PROPERTIES POSITION_INDEPENDENT_CODE true)

if( hpGEM_STACKTRACE_DEMANGLE )
  target_link_libraries(Logger dl)
endif()

add_subdirectory(LinearAlgebra)
add_subdirectory(Geometry)
add_subdirectory(Integration)
add_subdirectory(Base)
add_subdirectory(Output)
