include_directories(${hpGEM_SOURCE_DIR}/kernel/ .)

add_library(EigenSolvers OBJECT
            JacobiDavidsonMaxwell.cpp)

target_link_libraries(EigenSolvers
    PUBLIC
        Logger
        MPI::MPI_CXX
        BLAS::BLAS
        LAPACK::LAPACK
        PETSc::PETSc SLEPc::SLEPc
)
