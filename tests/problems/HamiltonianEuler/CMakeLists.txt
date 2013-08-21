include_directories(. ${hpGEM_SOURCE_DIR}/src/)

#First declare the rules that will form the executable tests#
#############################################################


if(hpGEM_INCLUDE_PETSC_SUPPORT)

	include_directories(${PETSC_DIR}/include)
	include_directories(${PETSC_DIR}/${PETSC_ARCH}/include)
	include_directories(${PETSC_DIR}/include/mpiuni)

	add_executable(hEuler.out 
		main.cpp
        	InitialConditions.cpp
		${hpGEM_SOURCE_DIR}/src/Output/TecplotPhysicalGeometryIterator.cpp				  
		#miscellaneous required stuff
		${hpGEM_SOURCE_DIR}/src/Base/AssembleBasisFunctionSet.cpp
		${hpGEM_SOURCE_DIR}/src/Base/L2Norm.cpp
        	${hpGEM_SOURCE_DIR}/src/Base/Norm2.cpp
      		)
	
		target_link_libraries(hEuler.out Geometry QuadratureRules Matrix ${PETSC_LIBRARIES})

else() 
	add_custom_target(hEuler.out 
		COMMENT "\t ERROR:: You need PETSC to compile hEuler.out. \n\t ERROR:: Please rerun cmake and enable PETSC support")
endif()

