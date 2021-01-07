# NOTE: Not including the MasterMake, instead registering the test ourselves

add_executable(vtk-test-output vtk-output.cpp)
target_link_libraries(vtk-test-output HPGEM::HPGEM)

add_test(NAME VTKOutputSelfTest
        COMMAND ${CMAKE_COMMAND}
            -Dgenerator_cmd=$<TARGET_FILE:vtk-test-output>
            -Dreference_path=${CMAKE_CURRENT_SOURCE_DIR}/vtk-references
            -P ${CMAKE_CURRENT_SOURCE_DIR}/vtk-test.cmake)
