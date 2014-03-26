include_directories(. ${hpGEM_SOURCE_DIR}/src/)

add_library(QuadratureRules SHARED
        # Gauss quadrature rules
        ${hpGEM_SOURCE_DIR}/src/Integration/QuadratureRules/AllGaussQuadratureRules.cpp
        ${hpGEM_SOURCE_DIR}/src/Integration/QuadratureRules/GaussQuadratureRulesForLine.cpp
        ${hpGEM_SOURCE_DIR}/src/Integration/QuadratureRules/GaussQuadratureRulesForTriangle.cpp
        ${hpGEM_SOURCE_DIR}/src/Integration/QuadratureRules/GaussQuadratureRulesForSquare.cpp
        ${hpGEM_SOURCE_DIR}/src/Integration/QuadratureRules/GaussQuadratureRulesForTetrahedron.cpp
        ${hpGEM_SOURCE_DIR}/src/Integration/QuadratureRules/GaussQuadratureRulesForPyramid.cpp
        ${hpGEM_SOURCE_DIR}/src/Integration/QuadratureRules/GaussQuadratureRulesForTriangularPrism.cpp
        ${hpGEM_SOURCE_DIR}/src/Integration/QuadratureRules/GaussQuadratureRulesForCube.cpp
        ${hpGEM_SOURCE_DIR}/src/Integration/QuadratureRules/GaussQuadratureRulesForHypercube.cpp
        
        )
target_link_libraries(QuadratureRules Reference_geometries)

