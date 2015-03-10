include_directories(. ${hpGEM_SOURCE_DIR}/kernel/)

#First declare the rules that will form the executable tests#
#############################################################

if(hpGEM_USE_PETSC)

    include_directories(${PETSC_INCLUDES})
    
add_executable(Laplace.out 
		main.cpp
		
        	${hpGEM_SOURCE_DIR}/kernel/Base/MeshManipulator.cpp
		${hpGEM_SOURCE_DIR}/kernel/Output/TecplotPhysicalGeometryIterator.cpp
		${hpGEM_SOURCE_DIR}/kernel/Output/TecplotDiscontinuousSolutionWriter.cpp
                ${hpGEM_SOURCE_DIR}/kernel/Output/VTKSpecificTimeWriter.cpp
        #	${hpGEM_SOURCE_DIR}/kernel/Integration/QuadratureRules/AllGaussQuadratureRules.cpp		
							  
		#miscellaneous required stuff
		${hpGEM_SOURCE_DIR}/kernel/Base/AssembleBasisFunctionSet.cpp
		${hpGEM_SOURCE_DIR}/kernel/Base/Norm2.cpp
		${hpGEM_SOURCE_DIR}/kernel/Base/L2Norm.cpp
      		)
target_link_libraries(Laplace.out hpGEM_Base Geometry QuadratureRules LinearAlgebra ${PETSC_LIBRARIES})

else() 
	add_custom_target(Laplace.out 
		COMMENT "\t ERROR:: You need PETSC to compile Laplace.out. \n\t ERROR:: Please rerun cmake and enable PETSC support")
endif()
