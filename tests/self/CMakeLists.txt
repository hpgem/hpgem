add_subdirectory(Basic)
if(hpGEM_USE_PETSC OR hpGEM_USE_COMPLEX_PETSC)
	add_subdirectory(PETSc)
endif(hpGEM_USE_PETSC OR hpGEM_USE_COMPLEX_PETSC)

if (hpGEM_USE_COMPLEX_PETSC AND hpGEM_USE_SLEPC)
	add_subdirectory(DGMax)
endif(hpGEM_USE_COMPLEX_PETSC AND hpGEM_USE_SLEPC)

add_subdirectory(Output)