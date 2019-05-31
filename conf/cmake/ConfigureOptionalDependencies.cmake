
if(hpGEM_USE_MPI)
    FIND_PACKAGE(MPI REQUIRED)
    # Add OMPI_SKIP_MPICXX to exclude the deprecated c++ bindings of MPI
    add_definitions(-DHPGEM_USE_MPI -DOMPI_SKIP_MPICXX)
    include_directories(${MPI_CXX_INCLUDE_PATH})
    if(NOT MPI_FOUND)
        message(FATAL_ERROR "The option you have chosen requires mpi and you do not have this installed. Please install")
    endif()
endif()

if(hpGEM_USE_METIS)
    FIND_PACKAGE(METIS REQUIRED)
    add_definitions(-DHPGEM_USE_METIS)
    include_directories(${METIS_INCLUDE_DIR})
    if(NOT METIS_FOUND)
        message(FATAL_ERROR "The option you have chosen requires metis and you do not have this installed. Please install")
    endif()
endif()

if(hpGEM_USE_QHULL)
    FIND_PACKAGE(QHULL REQUIRED)
    add_definitions(-DHPGEM_USE_QHULL)
    if(NOT QHULL_FOUND)
        message(FATAL_ERROR "The option you have chosen requires QHull and you do not have this installed. Please install")
    endif()
endif()

if (hpGEM_USE_PETSC OR hpGEM_USE_COMPLEX_PETSC)
    set(PETSC_ARCH "" CACHE STRING "PETSc arch string should be set be the user; he will know this if he has installed PETSc")
    set(hpGEM_LOGLEVEL "DEFAULT" CACHE STRING "Verbosity of hpGEM. Default is recommended.")

    if (hpGEM_USE_PETSC AND hpGEM_USE_COMPLEX_PETSC)
        message(FATAL_ERROR "The options 'hpGEM_USE_PETSC' and 'hpGEM_USE_COMPLEX_PETSC' are mutually exclusive, please disable one of them")
    endif()

    enable_language(C)
    find_package(PETSc REQUIRED COMPONENTS CXX)

    # Check whether we found a petsc with complex or real numbers
    include(CheckSymbolExists)
    check_symbol_exists(PETSC_USE_COMPLEX petscconf.h PETSC_IS_COMPLEX)

    set(hpGEM_USE_ANY_PETSC ON "A version of PETSc has been configured")

    # Check if this matches with what the user expects
    if (PETSC_IS_COMPLEX)
        if (hpGEM_USE_PETSC)
            message(FATAL_ERROR "Requested PETSc with real numbers, found one with complex numbers")
        endif()
        add_definitions(-DHPGEM_USE_COMPLEX_PETSC -DHPGEM_USE_ANY_PETSC)

    else()
        if (hpGEM_USE_COMPLEX_PETSC)
            message(FATAL_ERROR "Requested PETSc with complex numbers, found one with real numbers")
        endif()
        add_definitions(-DHPGEM_USE_PETSC -DHPGEM_USE_ANY_PETSC)
    endif()
else()
    set(hpGEM_USE_ANY_PETSC OFF "No version of PETSc has been configured")
endif() # Petsc

if(hpGEM_USE_SLEPC)
    enable_language(C)
    if(NOT hpGEM_USE_PETSC AND NOT hpGEM_USE_COMPLEX_PETSC)
        message(FATAL_ERROR
                "SLEPc depends on PETSc, please also enable PETSc support")
    endif()
    FIND_PACKAGE(SLEPc REQUIRED)
    add_definitions(-DHPGEM_USE_SLEPC)
    if(NOT SLEPC_FOUND)
        message(FATAL_ERROR
                "The option you have chosen requires SLEPc and you do not have this installed. Please install")
    endif()
endif()

if(hpGEM_USE_EXTERNAL_BLAS)
    FIND_PACKAGE(BLAS QUIET)
    if(NOT BLAS_FOUND)
        message(SYSTEM "Could not find external BLAS, using reference implementation instead")
        set(hpGEM_USE_EXTERNAL_BLAS OFF CACHE BOOL "Try to use a version of BLAS optimised for your system" FORCE)
    endif()
endif()

if(hpGEM_USE_EXTERNAL_LAPACK)
    FIND_PACKAGE(LAPACK QUIET)
    if(NOT LAPACK_FOUND)
        message(SYSTEM "Could not find external LAPACK, using reference implementation instead")
        set(hpGEM_USE_EXTERNAL_LAPACK OFF CACHE BOOL "Try to use a version of LAPACK optimised for your system" FORCE)
    endif()
endif()

if(hpGEM_USE_SUNDIALS)
    if(hpGEM_USE_MPI)
        message(FATAL_ERROR "SUNDIALS is only serial (for now)")
    endif()
    #enable_language(C)
    FIND_PACKAGE(SUNDIALS QUIET)
    add_definitions(-DHPGEM_USE_SUNDIALS)
    if(NOT SUNDIALS_FOUND)
        message(FATAL_ERROR
                "The option you have choosen requires SUNDIALS and you do not have this installed. Please install")
    endif()
endif()