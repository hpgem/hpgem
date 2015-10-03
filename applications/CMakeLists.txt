#list all the directories
FILE( GLOB subdirs RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
		${CMAKE_CURRENT_SOURCE_DIR}/*)
SET( applications "" )
FOREACH( subdir ${subdirs} )
#well.. these could be files as well
#or empty directories. So yeah
  IF(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${subdir} AND
     EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/CMakeLists.txt )
    LIST( APPEND applications ${subdir} )
  ENDIF()
ENDFOREACH()

#now add all the application
FOREACH( application ${applications} )
  ADD_SUBDIRECTORY( ${application} )
