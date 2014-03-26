include_directories(. ${hpGEM_SOURCE_DIR}/kernel/)

#First declare the rules that will form the executable tests#
#############################################################

add_executable (QuadratureRules_UnitTest.out 
        QuadratureRules_UnitTest.cpp        
        ${hpGEM_SOURCE_DIR}/kernel/Integration/QuadratureRules/QuadratureRuleSet.cpp 
      )

target_link_libraries(QuadratureRules_UnitTest.out Matrix Geometry QuadratureRules)
