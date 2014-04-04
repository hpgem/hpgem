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
	target_link_libraries(${EXECNAME}  Geometry)
endforeach()



#include_directories(${hpGEM_SOURCE_DIR}/kernel/ .)
#
#
##First declare the rules that will form the executable tests#
##############################################################
#
##-----REFERNECE-GEOMETRIES-----REFERNECE-GEOMETRIES-----REFERNECE-GEOMETRIES-----REFERNECE-GEOMETRIES
##0-d#
#add_executable(ReferencePointUnitTest.out
#               ReferencePointUnitTest.cpp)
#target_link_libraries(ReferencePointUnitTest.out Geometry)
#add_test(ReferencePointUnitTest ReferencePointUnitTest.out)
#
##1-d#
#add_executable(ReferenceLineUnitTest.out
#               ReferenceLineUnitTest.cpp)
#target_link_libraries(ReferenceLineUnitTest.out Geometry)
#add_test(ReferenceLineUnitTest ReferenceLineUnitTest.out)
#
##2-d#
#add_executable(ReferenceTriangleUnitTest.out
#               ReferenceTriangleUnitTest.cpp)
#target_link_libraries(ReferenceTriangleUnitTest.out Geometry)
#add_test(ReferenceTriangleUnitTest ReferenceTriangleUnitTest.out)
#
#add_executable(ReferenceSquareUnitTest.out
#               ReferenceSquareUnitTest.cpp)
#target_link_libraries(ReferenceSquareUnitTest.out Geometry)
#add_test(ReferenceSquareUnitTest ReferenceSquareUnitTest.out)
#
##3-d#
#add_executable(ReferenceTetrahedronUnitTest.out
#               ReferenceTetrahedronUnitTest.cpp)
#target_link_libraries(ReferenceTetrahedronUnitTest.out Geometry)
#add_test(ReferenceTetrahedronUnitTest ReferenceTetrahedronUnitTest.out)
#
#add_executable(ReferencePyramidUnitTest.out
#               ReferencePyramidUnitTest.cpp)
#target_link_libraries(ReferencePyramidUnitTest.out Geometry)
#add_test(ReferencePyramidUnitTest ReferencePyramidUnitTest.out)
#
#add_executable(ReferenceTriangularPrismUnitTest.out
#               ReferenceTriangularPrismUnitTest.cpp)
#target_link_libraries(ReferenceTriangularPrismUnitTest.out Geometry)        
#add_test(ReferenceTriangularPrismUnitTest ReferenceTriangularPrismUnitTest.out)
#
#add_executable(ReferenceCubeUnitTest.out
#               ReferenceCubeUnitTest.cpp)
#target_link_libraries(ReferenceCubeUnitTest.out Geometry)
#add_test(ReferenceCubeUnitTest ReferenceCubeUnitTest.out)
#
#4-d#
#add_executable(ReferenceHypercubeUnitTest.out 
#               ReferenceHypercubeUnitTest.cpp
#	       )
#target_link_libraries(ReferenceHypercubeUnitTest.out Geometry)
#add_test(ReferenceHypercubeUnitTest ReferenceHypercubeUnitTest.out)
#
#
##---------POINTS---------POINTS---------POINTS---------POINTS---------POINTS---------POINTS---------
#
#add_executable(PointUnitTest.out
#               PointUnitTest.cpp
#               )
#target_link_libraries(PointUnitTest.out Geometry)
#add_test(PointUnitTest PointUnitTest.out)
#
#add_executable(PointPhysicalUnitTest.out
#               PointPhysicalUnitTest.cpp
#               )
#       target_link_libraries(PointPhysicalUnitTest.out Geometry)
#       
#add_test(PointPhysicalUnitTest PointPhysicalUnitTest.out)
#
#
#add_executable(PointReferenceUnitTest.out
#               PointReferenceUnitTest.cpp
#               )
#       target_link_libraries(PointReferenceUnitTest.out Geometry)
#
#add_test(PointReferenceUnitTest PointReferenceUnitTest.out)
#
##--------ELEMENT--------ELEMENTS--------ELEMENTS--------ELEMENTS--------ELEMENTS--------ELEMENTS----
#
#add_executable(ElementGeometry_UnitTest.out
#               ElementGeometry_UnitTest.cpp
#		)
#	add_test(ElementGeomertryUnitTest ElementGeometry_UnitTest.out)
#
##--------MAPPINGS--------MAPPINGS--------MAPPINGS--------MAPPINGS--------MAPPINGS--------MAPPINGS--
#               
#add_executable(MappingHypercubeLinear_UnitTest.out
#               MappingHypercubeLinear_UnitTest.cpp
#               )
#target_link_libraries(MappingHypercubeLinear_UnitTest.out Geometry LinearAlgebra)
#               
#add_executable(MappingSimplexLinear_UnitTest.out
#               MappingSimplexLinear_UnitTest.cpp
#	       )
#target_link_libraries(MappingSimplexLinear_UnitTest.out Geometry LinearAlgebra)
#
#add_executable(MappingPointToLine_UnitTest.out
#               MappingPointToLine_UnitTest.cpp
#               )
#target_link_libraries(MappingPointToLine_UnitTest.out Geometry)
#               
#add_executable(MappingLineToSquare_UnitTest.out
#               MappingLineToSquare_UnitTest.cpp
#               )
#       target_link_libraries(MappingLineToSquare_UnitTest.out Geometry)
#
#add_executable(MappingSquareToCube_UnitTest.out
#               MappingSquareToCube_UnitTest.cpp
#	       )
#target_link_libraries(MappingSquareToCube_UnitTest.out Geometry)
#                              
#add_executable(FaceGeometry.out
#               FaceGeometryUnitTest.cpp)                
#target_link_libraries(FaceGeometry.out Geometry)
#
#

