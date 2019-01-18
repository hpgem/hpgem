include_directories(. ${hpGEM_SOURCE_DIR}/kernel/ ${PETSC_INCLUDES} ${SLEPC_INCLUDE_DIRS} ${hpGEM_SOURCE_DIR})


#First declare the rules that will form the executable tests#
if(hpGEM_USE_COMPLEX_PETSC AND hpGEM_USE_SLEPC)
 
    include_directories(${PETSC_INCLUDES})

#############################################################


add_executable( DG-Max.out
                DG-Max.cpp
				DGMaxLogger.cpp
                ElementInfos.cpp
                kspaceData.cpp
                BaseExtended.cpp
		${hpGEM_SOURCE_DIR}/kernel/Output/TecplotPhysicalGeometryIterator.cpp		
		${hpGEM_SOURCE_DIR}/kernel/Output/TecplotDiscontinuousSolutionWriter_Impl.h		  
		#miscellaneous required stuff
		${hpGEM_SOURCE_DIR}/kernel/Base/AssembleBasisFunctionSet.cpp
                ${hpGEM_SOURCE_DIR}/kernel/Integration/QuadratureRules/AllGaussQuadratureRules.cpp
                ${hpGEM_SOURCE_DIR}/kernel/Base/L2Norm.cpp
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
                Utils/KSpacePath.cpp
                Utils/HomogeneousBandStructure.cpp
      		)
target_link_libraries(DG-Max.out Output ${PETSC_LIBRARIES} ${SLEPC_LIBRARIES} )

else() 
	add_custom_target(DG-Max.out 
		COMMENT "\t ERROR:: You need a very complex edition of PETSC to compile DG-Max.out. \n\t ERROR:: Please rerun cmake and enable complex PETSC support, also SLEPc")
endif()
