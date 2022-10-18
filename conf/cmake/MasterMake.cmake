include_directories(	${hpGEM_SOURCE_DIR}/kernel
			${hpGEM_BINARY_DIR}/kernel
			.
			)


#Part 2 : Make run test for each of the demo files
##################################################

file(GLOB SELFTESTS "*SelfTest.cpp")
file(GLOB UNITTESTS "*UnitTest.cpp")
file(GLOB NEGATIVETESTS "*NegativeTest.cpp")
#for each demo add a test with the same name
foreach (TEST ${UNITTESTS})
        get_filename_component(EXECNAME ${TEST} NAME_WE)
        add_test(${EXECNAME} ${EXECNAME})
endforeach()

foreach (TEST ${SELFTESTS})
        get_filename_component(EXECNAME ${TEST} NAME_WE)
        add_test(${EXECNAME} ${EXECNAME})

        # For DivDGMax self tests, add parallel version of test
        if ((${TEST} MATCHES "DivDGMax") AND hpGEM_USE_MPI AND  hpGEM_USE_METIS)
                add_test(NAME "${EXECNAME}_parallel" COMMAND mpiexec -n 2 ${EXECNAME})
        endif()
endforeach()

foreach (TEST ${NEGATIVETESTS})
        get_filename_component(EXECNAME ${TEST} NAME_WE)
        add_test(${EXECNAME} ${EXECNAME})
        set_tests_properties(${EXECNAME} PROPERTIES WILL_FAIL true)
endforeach()

#Part 3 : Make tests for each of the selftest_data files
########################################################


file(GLOB TESTDATAFILES "${CMAKE_CURRENT_SOURCE_DIR}/SelfTestData/*")
#for each file in the selftest_data folder create a test. Which checks the data against this old data. The actually testing is done my the script self_test.
foreach(TESTFILE ${TESTDATAFILES})
        get_filename_component(TESTNAME ${TESTFILE} NAME)
        add_test(${TESTNAME} ${CMAKE_SOURCE_DIR}/Scripts/self_test ${TESTFILE} ${TESTNAME})
        #Add the newly created files to the clean target
        set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES "${TESTNAME}")
endforeach()
