include_directories(${hpGEM_SOURCE_DIR}/kernel/ .)

add_library(LinearAlgebra OBJECT
        MiddleSizeVector.cpp
        MiddleSizeMatrix.cpp
        SmallVector.cpp
        JacobiDavidsonMaxwell.cpp)

target_link_libraries(LinearAlgebra
    PUBLIC
        Logger
        MPI::MPI_CXX
        BLAS::BLAS
        LAPACK::LAPACK
        PETSc::PETSc SLEPc::SLEPc
)
