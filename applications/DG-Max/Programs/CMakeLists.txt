
# Definition of programs that are part of DGMax
# The current file only supports programs consisting of a single cpp file that
# is then linked to the DGMax library. If in the future more complex programs
# are needed, see the CMakeLists file in the 'applications' directory for an
# idea how to generically add subdirectories

file(GLOB CPPFILES  "*.cpp")
#for every cpp found
foreach(CPPFILE ${CPPFILES})
    #extract the actually file name
    get_filename_component(FILENAME ${CPPFILE} NAME)
    #extract the filename minus the cpp. This will be the name of exe file
    get_filename_component(EXECNAME ${CPPFILE} NAME_WE)
    #Make the exe
    add_executable(${EXECNAME} ${FILENAME})

    target_link_libraries(${EXECNAME}  DGMax-lib)
endforeach()


