include_directories(. ${hpGEM_SOURCE_DIR}/kernel/)

add_library(TimeIntegration STATIC
        #Fully explicit time integrators that depend on butcher tableaus
         ${hpGEM_SOURCE_DIR}/kernel/Base/TimeIntegration/RK4Methods.cpp
         ${hpGEM_SOURCE_DIR}/kernel/Base/TimeIntegration/ForwardEuler.cpp
         ${hpGEM_SOURCE_DIR}/kernel/Base/TimeIntegration/ForwardEulerTVD.cpp
    	 ${hpGEM_SOURCE_DIR}/kernel/Base/TimeIntegration/RK2TVD.cpp
         ${hpGEM_SOURCE_DIR}/kernel/Base/TimeIntegration/RK3TVD.cpp     
         ${hpGEM_SOURCE_DIR}/kernel/Base/TimeIntegration/MidPoint.cpp
         ${hpGEM_SOURCE_DIR}/kernel/Base/TimeIntegration/AllTimeIntegrators.cpp   
        )

set_target_properties(TimeIntegration PROPERTIES POSITION_INDEPENDENT_CODE true)
