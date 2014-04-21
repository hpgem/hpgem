include_directories(. ${hpGEM_SOURCE_DIR}/kernel/)

#First declare the rules that will form the executable tests#
#############################################################


if(hpGEM_INCLUDE_PETSC_SUPPORT)

	configure_file(bCoeff.txt
			${PROJECT_BINARY_DIR}/bCoeff.txt @ONLY IMMEDIATE)


	#include_directories(${PETSC_DIR}/include)
	#include_directories(${PETSC_DIR}/${PETSC_ARCH}/include)
	#include_directories(${PETSC_DIR}/include/mpiuni)
    include_directories(${PETSC_INCLUDES})
   # add_library(heulerlib HEuler.cpp)

	add_executable(hEulerApp.out
        main.cpp
        HEuler.cpp
        InitialConditions.cpp
		${hpGEM_SOURCE_DIR}/kernel/Output/TecplotPhysicalGeometryIterator.cpp	
		${hpGEM_SOURCE_DIR}/kernel/Output/TecplotDiscontinuousSolutionWriter.cpp
		${hpGEM_SOURCE_DIR}/kernel/Integration/QuadratureRules/AllGaussQuadratureRules.cpp
		#miscellaneous required stuff
		${hpGEM_SOURCE_DIR}/kernel/Base/MeshManipulator.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Base/AssembleBasisFunctionSet.cpp
		${hpGEM_SOURCE_DIR}/kernel/Base/L2Norm.cpp
        ${hpGEM_SOURCE_DIR}/kernel/Base/Norm2.cpp
        )
	
		target_link_libraries(hEulerApp.out hpGEM_Base Geometry QuadratureRules LinearAlgebra ${PETSC_LIBRARIES})

else() 
	add_custom_target(hEulerApp.out 
		COMMENT "\t ERROR:: You need PETSC to compile hEuler.out. \n\t ERROR:: Please rerun cmake and enable PETSC support \n\t ERROR:: If using raw cmake add the flag -DhpGEM_INCLUDE_PETSC_SUPPORT='on'")
endif()

