add_subdirectory(Basic)
if(hpGEM_INCLUDE_PETSC_SUPPORT)
	add_subdirectory(PETSc)
endif(hpGEM_INCLUDE_PETSC_SUPPORT)

