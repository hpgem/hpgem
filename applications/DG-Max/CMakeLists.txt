
#First declare the rules that will form the executable tests#
if(hpGEM_USE_COMPLEX_PETSC AND hpGEM_USE_SLEPC)
	# DG-Max library #
	##################

	add_library(DGMax-lib STATIC
		ElementInfos.cpp
		DGMaxLogger.cpp
		kspaceData.cpp
        DGMaxProgramUtils.cpp

		Algorithms/DGMaxEigenvalue.cpp
		Algorithms/DGMaxDiscretization.cpp
		Algorithms/DGMaxTimeIntegration.cpp
		Algorithms/DGMaxHarmonic.cpp
		Algorithms/DivDGMaxDiscretization.cpp
		Algorithms/DivDGMaxEigenvalue.cpp
		Algorithms/DivDGMaxHarmonic.cpp

		ProblemTypes/Harmonic/SampleHarmonicProblems.cpp
		ProblemTypes/Time/DummyTestProblem.cpp
		ProblemTypes/Time/SampleTestProblems.cpp
		Utils/BandstructureGNUPlot.cpp
		Utils/BandstructureUtils.cpp
		Utils/HomogeneousBandStructure.cpp
		Utils/BraggStackBandstructure.cpp
		Utils/MatrixBlocks.cpp
		Utils/KSpacePath.cpp
		Utils/KPhaseShift.cpp
		Utils/CGDGMatrixKPhaseShiftBuilder.cpp
		Utils/FaceKPhaseShiftBuilder.cpp
		Utils/Verification/DGMaxEVConvergenceTest.cpp
		Utils/Verification/DivDGMaxEVConvergenceTest.cpp
		Utils/Verification/EVConvergenceResult.cpp
		Utils/Verification/AbstractEVConvergenceTest.cpp
	)
	# TODO: Lacking source/include separation
	target_include_directories(DGMax-lib
		PUBLIC
			${CMAKE_CURRENT_SOURCE_DIR}
	)
	target_link_libraries(DGMax-lib
		PUBLIC
			HPGEM::HPGEM PETSc::PETSc SLEPc::SLEPc
	)

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
