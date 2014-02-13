include_directories(. ${hpGEM_SOURCE_DIR}/src/)

#First declare the rules that will form the executable tests#
#############################################################

if(hpGEM_INCLUDE_PETSC_SUPPORT)

    include_directories(${PETSC_INCLUDES})
    
add_executable(Laplace.out 
		main.cpp
		GlobalAssembly.cpp
		
        	${hpGEM_SOURCE_DIR}/src/Base/MeshManipulator.cpp
		${hpGEM_SOURCE_DIR}/src/Output/TecplotPhysicalGeometryIterator.cpp
		${hpGEM_SOURCE_DIR}/src/Output/TecplotDiscontinuousSolutionWriter.cpp
        ${hpGEM_SOURCE_DIR}/src/Integration/QuadratureRules/AllGaussQuadratureRules.cpp		
							  
		#miscellaneous required stuff
		${hpGEM_SOURCE_DIR}/src/Base/AssembleBasisFunctionSet.cpp
		${hpGEM_SOURCE_DIR}/src/Base/Norm2.cpp
		${hpGEM_SOURCE_DIR}/src/Base/L2Norm.cpp
      		)
target_link_libraries(Laplace.out hpGEM_Base Geometry QuadratureRules Matrix ${PETSC_LIBRARIES})

else() 
	add_custom_target(Laplace.out 
		COMMENT "\t ERROR:: You need PETSC to compile Laplace.out. \n\t ERROR:: Please rerun cmake and enable PETSC support")
endif()