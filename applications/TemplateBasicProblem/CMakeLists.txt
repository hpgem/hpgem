include_directories(. ${hpGEM_SOURCE_DIR}/kernel/)

#First declare the rules that will form the executable tests#
#############################################################

#choose if your application needs PETSc. Available options are
#'NOT hpGEM_USE_COMPLEX_PETSC', 'hpGEM_USE_PETSC'          and 'hpGEM_USE_COMPLEX_PETSC'
# (no PETSC) (default)        ,  (PETSc with real numbers) and (PETSC with complex numbers)
#When hpGEM is configured to interact with a complex build of PETSc some of its internal data
#structrues will also use complex numbers. This is why applications that don't need PETSc
#are disabled for this case.
if(NOT hpGEM_USE_COMPLEX_PETSC)
   
#change the name of the |
#executable here        V
    add_executable(ExampleProblem.out 
    		main.cpp
          		)
#and here                    |
#                            V
    target_link_libraries(ExampleProblem.out hpGEM_Base Output)
endif()
