include_directories(${hpGEM_SOURCE_DIR}/kernel/)

add_library(FE
    OBJECT
        BaseBasisFunction.cpp
        BasisFunctionSet.cpp
        BasisFunctions1DH1ConformingLine.cpp
        BasisFunctions2DH1ConformingSquare.cpp
        BasisFunctions2DH1ConformingTriangle.cpp
        BasisFunctions2DNedelec.cpp
        BasisFunctions3DH1ConformingCube.cpp
        BasisFunctions3DH1ConformingPrism.cpp
        BasisFunctions3DH1ConformingTetrahedron.cpp
        BasisFunctions3DAinsworthCoyle.cpp
        BasisFunctions3DNedelec.cpp
        BasisFunctions3DH1ConformingPyramid.cpp
        BasisFunctionsMonomials.cpp
        BasisFunctionsPiecewiseConstant.cpp
        helperFunctions.cpp
)

target_link_libraries(FE
    PUBLIC
        LinearAlgebra
        Reference_geometries)
