include_directories(. ${hpGEM_SOURCE_DIR}/kernel/)

#First declare the rules that will form the executable tests#
#############################################################
   
add_executable(LinearAdvection.out 
		main.cpp
                ${hpGEM_SOURCE_DIR}/kernel/Utilities/BasisFunctions2DH1ConformingTriangle.cpp
                ${hpGEM_SOURCE_DIR}/kernel/Utilities/helperFunctions.cpp

        	${hpGEM_SOURCE_DIR}/kernel/Base/MeshManipulator.cpp
		${hpGEM_SOURCE_DIR}/kernel/Output/TecplotPhysicalGeometryIterator.cpp
		${hpGEM_SOURCE_DIR}/kernel/Output/TecplotDiscontinuousSolutionWriter.cpp
        #	${hpGEM_SOURCE_DIR}/kernel/Integration/QuadratureRules/AllGaussQuadratureRules.cpp		
							  
		#miscellaneous required stuff
		${hpGEM_SOURCE_DIR}/kernel/Base/AssembleBasisFunctionSet.cpp
		${hpGEM_SOURCE_DIR}/kernel/Base/Norm2.cpp
		${hpGEM_SOURCE_DIR}/kernel/Base/L2Norm.cpp
      		)
target_link_libraries(LinearAdvection.out hpGEM_Base Geometry QuadratureRules LinearAlgebra ${PETSC_LIBRARIES})

