include_directories(. ${hpGEM_SOURCE_DIR}/kernel/ ${PETSC_INCLUDES} ${SLEPC_INCLUDE_DIRS} ${hpGEM_SOURCE_DIR})


#First declare the rules that will form the executable tests#
if(hpGEM_USE_COMPLEX_PETSC AND hpGEM_USE_SLEPC)
 
    include_directories(${PETSC_INCLUDES})

	# DG-Max library #
	##################

	add_library(DGMax-lib STATIC
		ElementInfos.cpp
		DGMaxLogger.cpp
		kspaceData.cpp
        DGMaxProgramUtils.cpp

		Algorithms/DGMaxEigenValue.cpp
		Algorithms/DGMaxDiscretization.cpp
		Algorithms/DGMaxTimeIntegration.cpp
		Algorithms/DGMaxHarmonic.cpp
		Algorithms/DivDGMaxDiscretization.cpp
		Algorithms/DivDGMaxEigenValue.cpp
		Algorithms/DivDGMaxHarmonic.cpp

		ProblemTypes/Harmonic/SampleHarmonicProblems.cpp
		ProblemTypes/Time/DummyTestProblem.cpp
		ProblemTypes/Time/SampleTestProblems.cpp
		Utils/BandstructureGNUPlot.cpp
		Utils/BandstructureUtils.cpp
		Utils/HomogeneousBandStructure.cpp
		Utils/BraggStackBandstructure.cpp
		Utils/KSpacePath.cpp
	)

	target_link_libraries(DGMax-lib Output ${PETSC_LIBRARIES} ${SLEPC_LIBRARIES} )


	# Actual DG-Max program #
	#########################
	add_executable(DG-Max.out
		DG-Max.cpp
	)

	target_link_libraries(DG-Max.out DGMax-lib)

	# Add additional programs #
	###########################
	add_subdirectory(Programs)

else() 
	add_custom_target(DG-Max.out 
		COMMENT "\t ERROR:: You need a very complex edition of PETSC to compile DG-Max.out. \n\t ERROR:: Please rerun cmake and enable complex PETSC support, also SLEPc")
endif()
