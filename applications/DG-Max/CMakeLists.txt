include_directories(. ${hpGEM_SOURCE_DIR}/kernel/)


#First declare the rules that will form the executable tests#
if(hpGEM_USE_COMPLEX_PETSC AND hpGEM_USE_SLEPC)
 
    include_directories(${PETSC_INCLUDES})
    include_directories(${SLEPC_INCLUDES})

#############################################################


add_executable( DG-Max.out
                DG-Max.cpp
                ElementInfos.cpp
                kspaceData.cpp
                BasisFunctionCollection_Curl.cpp
                BaseExtended.cpp
      		)
target_link_libraries(DG-Max.out hpGEM_Base ${PETSC_LIBRARIES} ${SLEPC_LIBRARIES} )

else() 
	add_custom_target(DG-Max.out 
		COMMENT "\t ERROR:: You need a very complex edition of PETSC to compile DG-Max.out. \n\t ERROR:: Please rerun cmake and enable complex PETSC support, also SLEPc")
endif()
