add_subdirectory(Basic)
if(hpGEM_USE_PETSC)
	add_subdirectory(PETSc)
endif(hpGEM_USE_PETSC)

