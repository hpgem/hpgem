include_directories(. ${hpGEM_SOURCE_DIR}/src/ ${PETSC_INCLUDES} ${SLEPC_INCLUDE_DIRS} ${hpGEM_SOURCE_DIR})


#First declare the rules that will form the executable tests#
#############################################################


add_executable( DG-Max.out
                DG-Max.cpp
                fillMatrices.cpp
                ElementInfos.cpp
                kspaceData.cpp
                BasisFunctionCollection_Curl.cpp
                BaseExtended.cpp
		${hpGEM_SOURCE_DIR}/src/Output/TecplotPhysicalGeometryIterator.cpp		
		${hpGEM_SOURCE_DIR}/src/Output/TecplotDiscontinuousSolutionWriter.cpp		  
		#miscellaneous required stuff
		${hpGEM_SOURCE_DIR}/src/Base/AssembleBasisFunctionSet.cpp
                ${hpGEM_SOURCE_DIR}/src/Integration/QuadratureRules/AllGaussQuadratureRules.cpp
                ${hpGEM_SOURCE_DIR}/src/Base/L2Norm.cpp
      		)
target_link_libraries(DG-Max.out hpGEM_Base Geometry QuadratureRules Matrix ${PETSC_LIBRARIES} ${SLEPC_LIBRARIES} )
