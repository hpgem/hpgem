include_directories(${hpGEM_SOURCE_DIR}/kernel/ .)

if(hpGEM_USE_ANY_PETSC AND hpGEM_USE_SLEPC AND hpGEM_USE_MPI)

        add_library(EigenSolvers OBJECT
                JacobiDavidsonMaxwell.cpp
                DoehlerMaxwell.cpp)

        target_link_libraries(EigenSolvers
        PUBLIC
                Logger
                MPI::MPI_CXX
                BLAS::BLAS
                LAPACK::LAPACK
                PETSc::PETSc SLEPc::SLEPc
        )
endif()
