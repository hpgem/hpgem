add_subdirectory(LinearAlgebra)
add_subdirectory(Geometry)
add_subdirectory(Base)
add_subdirectory(Integration)
add_subdirectory(Output)
add_subdirectory(Utilities)
add_subdirectory(Preprocessor)

if (hpGEM_USE_COMPLEX_PETSC AND hpGEM_USE_SLEPC)
    add_subdirectory(DGMax)
endif(hpGEM_USE_COMPLEX_PETSC AND hpGEM_USE_SLEPC)

add_library (tests-main tests-main.cpp)
