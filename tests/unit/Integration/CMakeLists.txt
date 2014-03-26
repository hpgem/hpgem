add_subdirectory(QuadratureRules)
include_directories(${hpGEM_SOURCE_DIR}/kernel/ .)


#add_executable(IntegratioinOfPysicalFunctionUnitTest.out
#		IntegrationOfPhysicalFunctionUnitTest.cpp
#       # Miscellaneous
#        ${hpGEM_SOURCE_DIR}/kernel/Base/AssembleBasisFunctionSet.cpp
#        ${hpGEM_SOURCE_DIR}/kernel/Base/L2Norm.cpp
#        )
#add_test(IntegratioinOfPysicalFunctionUnitTest IntegratioinOfPysicalFunctionUnitTest.out)

add_executable(IntegrationUnitTest.out 
        IntegrationUnitTest.cpp

        # Miscellaneous
        ${hpGEM_SOURCE_DIR}/kernel/Base/AssembleBasisFunctionSet.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Base/FaceCacheData.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Base/ElementCacheData.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Base/L2Norm.cpp
        )
target_link_libraries(IntegrationUnitTest.out Geometry hpGEM_Base QuadratureRules Matrix)
add_test(IntegrationUnitTest IntegrationUnitTest.out)


