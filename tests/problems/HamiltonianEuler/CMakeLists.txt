include_directories(. ${hpGEM_SOURCE_DIR}/src/)

#First declare the rules that will form the executable tests#
#############################################################

add_executable(hEuler.out 
		main.cpp
		${hpGEM_SOURCE_DIR}/src/Output/TecplotPhysicalGeometryIterator.cpp				  
		#miscellaneous required stuff
		${hpGEM_SOURCE_DIR}/src/Base/AssembleBasisFunctionSet.cpp
		${hpGEM_SOURCE_DIR}/src/Base/L2Norm.cpp
      		)
target_link_libraries(hEuler.out Geometry QuadratureRules Matrix)
