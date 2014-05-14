include(${CMAKE_SOURCE_DIR}/conf/cmake/MasterMake.cmake)

#temporary hack because I want to also test PETSc and GlobalAssembly
include_directories(${PETSC_INCLUDES})

#Part 1 : All cpp files will be made to exec files
#########################################
#Collect all the names of the cpp, note at the moment limited to demos, but will be fixed later
file(GLOB CPPFILES  "*.cpp")
#for every cpp found
foreach(CPPFILE ${CPPFILES})
	#extract the actually file name
	get_filename_component(FILENAME ${CPPFILE} NAME)
	#extract the filename minus the cpp. This will be the name of exe file
	get_filename_component(EXECNAME ${CPPFILE} NAME_WE)
	#Make the exe
	add_executable(${EXECNAME} ${FILENAME} 
	
	#temporary hack because I want to also test these basis functions and the Global assembly
	${hpGEM_SOURCE_DIR}/kernel/Utilities/BasisFunctions1DH1ConformingLine.cpp
	${hpGEM_SOURCE_DIR}/kernel/Utilities/BasisFunctions2DH1ConformingSquare.cpp
	${hpGEM_SOURCE_DIR}/kernel/Utilities/BasisFunctions2DH1ConformingTriangle.cpp
	${hpGEM_SOURCE_DIR}/kernel/Utilities/BasisFunctions3DH1ConformingCube.cpp
	${hpGEM_SOURCE_DIR}/kernel/Utilities/BasisFunctions3DH1ConformingTetrahedron.cpp
	${hpGEM_SOURCE_DIR}/kernel/Utilities/helperFunctions.cpp
	${hpGEM_SOURCE_DIR}/kernel/Utilities/GlobalMatrix.cpp
	${hpGEM_SOURCE_DIR}/kernel/Utilities/GlobalVector.cpp)

	target_link_libraries(${EXECNAME} hpGEM_Base Output ${PETSC_LIBRARIES})
endforeach()

