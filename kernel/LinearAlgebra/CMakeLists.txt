include_directories(${hpGEM_SOURCE_DIR}/kernel/ .)

add_library(LinearAlgebra SHARED
        MiddleSizeVector.cpp
        MiddleSizeMatrix.cpp)
target_link_libraries(LinearAlgebra
    PUBLIC
        Logger
        BLAS::BLAS
        LAPACK::LAPACK
)
