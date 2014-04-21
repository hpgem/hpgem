include_directories(. ${hpGEM_SOURCE_DIR}/kernel/)

add_library(QuadratureRules SHARED
        # Gauss quadrature rules
        ${hpGEM_SOURCE_DIR}/kernel/Integration/QuadratureRules/AllGaussQuadratureRules.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Integration/QuadratureRules/GaussQuadratureRulesForLine.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Integration/QuadratureRules/GaussQuadratureRulesForTriangle.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Integration/QuadratureRules/GaussQuadratureRulesForSquare.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Integration/QuadratureRules/GaussQuadratureRulesForTetrahedron.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Integration/QuadratureRules/GaussQuadratureRulesForPyramid.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Integration/QuadratureRules/GaussQuadratureRulesForTriangularPrism.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Integration/QuadratureRules/GaussQuadratureRulesForCube.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Integration/QuadratureRules/GaussQuadratureRulesForHypercube.cpp
        
        )
target_link_libraries(QuadratureRules Reference_geometries)

