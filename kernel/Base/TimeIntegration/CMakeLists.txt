include_directories(. ${hpGEM_SOURCE_DIR}/kernel/)

add_library(TimeIntegration STATIC
        # Gauss quadrature rules
         ${hpGEM_SOURCE_DIR}/kernel/Base/TimeIntegration/ButcherTableau.cpp
         ${hpGEM_SOURCE_DIR}/kernel/Base/TimeIntegration/RK4Methods.cpp
         ${hpGEM_SOURCE_DIR}/kernel/Base/TimeIntegration/AllTimeIntegrators.cpp        
        )

set_target_properties(TimeIntegration PROPERTIES POSITION_INDEPENDENT_CODE true)