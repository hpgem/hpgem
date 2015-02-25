include_directories(. ${hpGEM_SOURCE_DIR}/kernel/)

#First declare the rules that will form the executable tests#
#############################################################
if(hpGEM_USE_PETSC)

    include_directories(${PETSC_INCLUDES})
    
add_executable(TutorialLaplace.out 
		main.cpp
		${hpGEM_SOURCE_DIR}/kernel/Utilities/GlobalMatrix.cpp
		${hpGEM_SOURCE_DIR}/kernel/Utilities/GlobalVector.cpp
		${hpGEM_SOURCE_DIR}/kernel/Utilities/BasisFunctions1DH1ConformingLine.cpp
		${hpGEM_SOURCE_DIR}/kernel/Utilities/BasisFunctions2DH1ConformingSquare.cpp
		${hpGEM_SOURCE_DIR}/kernel/Utilities/BasisFunctions2DH1ConformingTriangle.cpp
		${hpGEM_SOURCE_DIR}/kernel/Utilities/BasisFunctions3DH1ConformingCube.cpp
		${hpGEM_SOURCE_DIR}/kernel/Utilities/BasisFunctions3DH1ConformingTetrahedron.cpp
		${hpGEM_SOURCE_DIR}/kernel/Utilities/helperFunctions.cpp
		
                ${hpGEM_SOURCE_DIR}/kernel/Base/MeshManipulator.cpp
		${hpGEM_SOURCE_DIR}/kernel/Output/TecplotPhysicalGeometryIterator.cpp
		${hpGEM_SOURCE_DIR}/kernel/Output/TecplotDiscontinuousSolutionWriter.cpp
		${hpGEM_SOURCE_DIR}/kernel/Output/GNUPlotDiscontinuousSolutionWriter.cpp
#       	${hpGEM_SOURCE_DIR}/kernel/Integration/QuadratureRules/AllGaussQuadratureRules.cpp		
							  
		#miscellaneous required stuff
		${hpGEM_SOURCE_DIR}/kernel/Base/AssembleBasisFunctionSet.cpp
		${hpGEM_SOURCE_DIR}/kernel/Base/Norm2.cpp
		${hpGEM_SOURCE_DIR}/kernel/Base/L2Norm.cpp
      		)
target_link_libraries(TutorialLaplace.out hpGEM_Base Geometry QuadratureRules LinearAlgebra ${PETSC_LIBRARIES})

else() 
        #message(WARNING "TutorialLaplace.out\nTutorialLaplace.out cannot build without PETSC support!\n\tPlease rerun cmake and enable PETSC support")
	add_custom_target(TutorialLaplace.out 
		COMMENT "\t ERROR:: You need PETSC to compile LaplaceEasier.out. \n\t ERROR:: Please rerun cmake and enable PETSC support")
endif()
