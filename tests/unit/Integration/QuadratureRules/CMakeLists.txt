#include_directories(. ${hpGEM_SOURCE_DIR}/kernel/)
#
##First declare the rules that will form the executable tests#
##############################################################
#
#add_executable (QuadratureRules_UnitTest.out 
#        QuadratureRules_UnitTest.cpp        
#        ${hpGEM_SOURCE_DIR}/kernel/Integration/QuadratureRules/QuadratureRuleSet.cpp 
#      )
#
#target_link_libraries(QuadratureRules_UnitTest.out Geometry QuadratureRules)

include(${CMAKE_SOURCE_DIR}/conf/cmake/MasterMake.cmake)

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
	add_executable(${EXECNAME} ${FILENAME})

	target_link_libraries(${EXECNAME}  HPGEM::HPGEM tests-main)
endforeach()
