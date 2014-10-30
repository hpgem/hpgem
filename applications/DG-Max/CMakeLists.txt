include_directories(. ${hpGEM_SOURCE_DIR}/kernel/ ${PETSC_INCLUDES} ${SLEPC_INCLUDE_DIRS} ${hpGEM_SOURCE_DIR})


#First declare the rules that will form the executable tests#
#############################################################

if(hpGEM_USE_COMPLEX_PETSC AND hpGEM_USE_SLEPC)

    include_directories(${PETSC_INCLUDES})

add_executable( DG-Max.out
                DG-Max.cpp
                fillMatrices.cpp
                ElementInfos.cpp
                kspaceData.cpp
                BasisFunctionCollection_Curl.cpp
                BaseExtended.cpp
		${hpGEM_SOURCE_DIR}/kernel/Output/TecplotPhysicalGeometryIterator.cpp		
		${hpGEM_SOURCE_DIR}/kernel/Output/TecplotDiscontinuousSolutionWriter.cpp		  
		#miscellaneous required stuff
		${hpGEM_SOURCE_DIR}/kernel/Base/AssembleBasisFunctionSet.cpp
                ${hpGEM_SOURCE_DIR}/kernel/Integration/QuadratureRules/AllGaussQuadratureRules.cpp
                ${hpGEM_SOURCE_DIR}/kernel/Base/L2Norm.cpp
      		)
target_link_libraries(DG-Max.out Output ${PETSC_LIBRARIES} ${SLEPC_LIBRARIES} )

else() 
	add_custom_target(DG-Max.out 
		COMMENT "\t ERROR:: You need a very complex edition of PETSC to compile DG-Max.out. \n\t ERROR:: Please rerun cmake and enable complex PETSC support, also SLEPc")
endif()