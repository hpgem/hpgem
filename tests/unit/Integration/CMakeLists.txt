add_subdirectory(QuadratureRules)
include_directories(${hpGEM_SOURCE_DIR}/src/ .)


#add_executable(IntegratioinOfPysicalFunctionUnitTest.out
#		IntegrationOfPhysicalFunctionUnitTest.cpp
#       # Miscellaneous
#        ${hpGEM_SOURCE_DIR}/src/Base/AssembleBasisFunctionSet.cpp
#        ${hpGEM_SOURCE_DIR}/src/Base/L2Norm.cpp
#        )
#add_test(IntegratioinOfPysicalFunctionUnitTest IntegratioinOfPysicalFunctionUnitTest.out)

add_executable(IntegrationUnitTest.out 
        IntegrationUnitTest.cpp

        # Miscellaneous
        ${hpGEM_SOURCE_DIR}/src/Base/AssembleBasisFunctionSet.cpp
        ${hpGEM_SOURCE_DIR}/src/Base/L2Norm.cpp
        )
target_link_libraries(IntegrationUnitTest.out Geometry QuadratureRules Matrix)
add_test(IntegrationUnitTest IntegrationUnitTest.out)


