include_directories(. ${hpGEM_SOURCE_DIR}/kernel/ ${PETSC_INCLUDES} ${SLEPC_INCLUDE_DIRS} ${hpGEM_SOURCE_DIR})


#First declare the rules that will form the executable tests#
#############################################################


add_executable( DG-Max.out
                DG-Max.cpp
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
${hpGEM_SOURCE_DIR}/kernel/Utilities/GlobalMatrix.cpp
${hpGEM_SOURCE_DIR}/kernel/Utilities/GlobalVector.cpp
      		)
target_link_libraries(DG-Max.out Output ${PETSC_LIBRARIES} ${SLEPC_LIBRARIES} )
