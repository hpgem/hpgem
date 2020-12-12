include_directories(${hpGEM_SOURCE_DIR}/kernel/ .)
add_subdirectory(TimeIntegration)

add_library(hpGEM_Base OBJECT
		CommandLineOptions.cpp
 		Element.cpp
		ElementBasisFunctions.cpp
 		ElementData.cpp
                ElementFactory.cpp
 		Face.cpp
 		FaceData.cpp
                FaceFactory.cpp
		FaceMatrix.cpp
 		Edge.cpp
 		Node.cpp
		HpgemAPISimplified.cpp
		#miscellaneous required stuff
		${hpGEM_SOURCE_DIR}/kernel/Base/MeshManipulatorBase.cpp
                Submesh.cpp
		${hpGEM_SOURCE_DIR}/kernel/Base/L2Norm.cpp
                MpiContainer.cpp

		${hpGEM_SOURCE_DIR}/kernel/Utilities/GlobalMatrix.cpp
		${hpGEM_SOURCE_DIR}/kernel/Utilities/GlobalVector.cpp
		${hpGEM_SOURCE_DIR}/kernel/Utilities/GlobalIndexing.cpp
		${hpGEM_SOURCE_DIR}/kernel/Utilities/ElementLocalIndexing.cpp
		${hpGEM_SOURCE_DIR}/kernel/Utilities/FaceLocalIndexing.cpp
		${hpGEM_SOURCE_DIR}/kernel/Utilities/SparsityEstimator.cpp
        GlobalUniqueIndex.cpp
)


set(CURRENT_TARGET hpGEM_Base)
include(${CMAKE_SOURCE_DIR}/conf/cmake/CLionFix.cmake)
set(CURRENT_TARGET)

target_link_libraries(hpGEM_Base FE Output Geometry TimeIntegration)

set_target_properties(hpGEM_Base PROPERTIES POSITION_INDEPENDENT_CODE true)

# All optional dependencies
if(hpGEM_USE_MPI)
	target_link_libraries(hpGEM_Base MPI::MPI_CXX)
endif()

if(hpGEM_USE_METIS)
	target_link_libraries(hpGEM_Base METIS::METIS)
endif()

if(hpGEM_USE_ANY_PETSC)
	target_link_libraries(hpGEM_Base PETSc::PETSc)
endif()

if(hpGEM_USE_SLEPC)
	target_link_libraries(hpGEM_Base SLEPc::SLEPc)
endif()

if(hpGEM_USE_QHULL)
	target_link_libraries(hpGEM_Base QHULL::QHULL)
endif()

