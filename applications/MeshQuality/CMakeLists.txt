if(hpGEM_USE_ANY_PETSC)
   
    add_executable(MeshQuality.out main.cpp)
    target_link_libraries(MeshQuality.out HPGEM::HPGEM)
endif()
