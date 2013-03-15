add_subdirectory(QuadratureRules)
include_directories(${hpGEM_SOURCE_DIR}/src/ .)

#First declare the rules that will form the executable tests#
#############################################################

add_executable(IntegratioinOfPysicalFunctionUnitTest.out
		IntegrationOfPhysicalFunctionUnitTest.cpp
        # Miscellaneous
        ${hpGEM_SOURCE_DIR}/src/Base/AssembleBasisFunctionSet.cpp
        ${hpGEM_SOURCE_DIR}/src/Base/L2Norm.cpp
        )

add_executable(IntegrationUnitTest.out 
        IntegrationUnitTest.cpp

        # Miscellaneous
        ${hpGEM_SOURCE_DIR}/src/Base/AssembleBasisFunctionSet.cpp
        ${hpGEM_SOURCE_DIR}/src/Base/L2Norm.cpp
        )
target_link_libraries(IntegrationUnitTest.out Geometry QuadratureRules Matrix)
target_link_libraries(IntegratioinOfPysicalFunctionUnitTest.out Geometry QuadratureRules Matrix)

