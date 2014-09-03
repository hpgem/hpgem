include_directories(. ${hpGEM_SOURCE_DIR}/kernel/)

#First declare the rules that will form the executable tests#
#############################################################

if(hpGEM_INCLUDE_PETSC_SUPPORT)

    include_directories(${PETSC_INCLUDES})
    
add_executable(PotentialFlow.out 
		main.cpp
		${hpGEM_SOURCE_DIR}/kernel/Utilities/GlobalMatrix.cpp
		${hpGEM_SOURCE_DIR}/kernel/Utilities/GlobalVector.cpp
		${hpGEM_SOURCE_DIR}/kernel/Utilities/BasisFunctions1DH1ConformingLine.cpp
		${hpGEM_SOURCE_DIR}/kernel/Utilities/BasisFunctions2DH1ConformingSquare.cpp
		${hpGEM_SOURCE_DIR}/kernel/Utilities/BasisFunctions2DH1ConformingTriangle.cpp
		${hpGEM_SOURCE_DIR}/kernel/Utilities/BasisFunctions3DH1ConformingCube.cpp
		${hpGEM_SOURCE_DIR}/kernel/Utilities/BasisFunctions3DH1ConformingTetrahedron.cpp
		${hpGEM_SOURCE_DIR}/kernel/Utilities/helperFunctions.cpp
		
        #	${hpGEM_SOURCE_DIR}/kernel/Integration/QuadratureRules/AllGaussQuadratureRules.cpp		
							  
		#miscellaneous required stuff
      		)
target_link_libraries(PotentialFlow.out hpGEM_Base Output Geometry QuadratureRules LinearAlgebra ${PETSC_LIBRARIES})

else() 
	add_custom_target(PotentialFlow.out 
		COMMENT "\t ERROR:: You need PETSC to compile Laplace.out. \n\t ERROR:: Please rerun cmake and enable PETSC support")
endif()
