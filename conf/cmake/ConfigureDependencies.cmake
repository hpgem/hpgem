
### REQUIRED DEPENDENCIES ###
#############################

FIND_PACKAGE(BLAS REQUIRED)
if(BLAS_FOUND AND NOT TARGET BLAS::BLAS)
    # CREATING A BLAS::BLAS TARGET as FINDBLAS does not create one
    add_library(BLAS::BLAS INTERFACE IMPORTED)
    set_target_properties(BLAS::BLAS PROPERTIES INTERFACE_LINK_LIBRARIES "${BLAS_LIBRARIES}" INTERFACE_LINK_FLAGS "${BLAS_LINKER_FLAGS}" )
endif()

FIND_PACKAGE(LAPACK REQUIRED)
if(LAPACK_FOUND AND NOT TARGET LAPACK::LAPACK)
    # CREATING A LAPACK::LAPACK TARGET as FINDLAPACK does not create one
    add_library(LAPACK::LAPACK INTERFACE IMPORTED)
    set_target_properties(BLAS::BLAS PROPERTIES INTERFACE_LINK_LIBRARIES  "${LAPACK_LIBRARIES}" INTERFACE_LINK_FLAGS "${LAPACK_LINKER_FLAGS}" )
endif()



### OPTIONAL DEPENDENCIES ###
#############################

if(hpGEM_USE_MPI)
    FIND_PACKAGE(MPI REQUIRED)
    # Add OMPI_SKIP_MPICXX to exclude the deprecated c++ bindings of MPI
    add_definitions(-DHPGEM_USE_MPI -DOMPI_SKIP_MPICXX)
    # Create target for easy linking
    if(NOT TARGET MPI::MPI_CXX)
        # Note, from CMAKE 3.9 there is a proper target defined
        add_library(MPI::MPI_CXX INTERFACE IMPORTED)
        set_target_properties(MPI::MPI_CXX PROPERTIES 
        INTERFACE_LINK_LIBRARIES "${MPI_CXX_LIBRARIES}" 
        INTERFACE_COMPILE_OPTIONS "${MPI_CXX_COMPILE_FLAGS}"
        INTERFACE_LINK_OPTIONS "${MPI_CXX_LINK_FLAGS}" 
        INTERFACE_INCLUDE_DIRECTORIES "${MPI_CXX_INCLUDE_PATH}")
    endif()
endif()

if(hpGEM_USE_METIS)
    FIND_PACKAGE(METIS REQUIRED)
    add_definitions(-DHPGEM_USE_METIS)
    include_directories(${METIS_INCLUDE_DIR})
    # Create target for easy linking
    add_library(METIS::METIS INTERFACE IMPORTED)
    target_include_directories(METIS::METIS INTERFACE "${METIS_INCLUDE_DIR}")
    target_link_libraries(METIS::METIS INTERFACE "${METIS_LIBRARIES}")
endif()

if(hpGEM_USE_QHULL)
    FIND_PACKAGE(QHULL REQUIRED)
    add_definitions(-DHPGEM_USE_QHULL)
    # Create target for easy linking
    add_library(QHULL::QHULL INTERFACE IMPORTED)
    if(CMAKE_BUILD_TYPE MATCHES "Debug")
        set_target_properties(QHULL::QHULL PROPERTIES 
        INTERFACE_LINK_LIBRARIES "${QHULL_DEBUG_LIBRARIES}" 
        INTERFACE_INCLUDE_DIRECTORIES "${QHULL_INCLUDE_DIR}")
    else()
        set_target_properties(QHULL::QHULL PROPERTIES 
        INTERFACE_LINK_LIBRARIES "${QHULL_LIBRARIES}" 
        INTERFACE_INCLUDE_DIRECTORIES "${QHULL_INCLUDE_DIR}")
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

    set(hpGEM_USE_ANY_PETSC ON)

    # Check if this matches with what the user expects
    if (PETSC_IS_COMPLEX)
        if (hpGEM_USE_PETSC)
            message(FATAL_ERROR "Requested PETSc with real numbers, found one with complex numbers")
        endif()
        # TODO: This should be moved
        add_definitions(-DHPGEM_USE_COMPLEX_PETSC -DHPGEM_USE_ANY_PETSC)

    else()
        if (hpGEM_USE_COMPLEX_PETSC)
            message(FATAL_ERROR "Requested PETSc with complex numbers, found one with real numbers")
        endif()
        add_definitions(-DHPGEM_USE_PETSC -DHPGEM_USE_ANY_PETSC)
    endif()
    # Create target for easy linking
    if(PETSC_FOUND)
        add_library(PETSc::PETSc INTERFACE IMPORTED)
        target_include_directories(PETSc::PETSc INTERFACE "${PETSC_INCLUDES}")
        target_link_libraries(PETSc::PETSc INTERFACE "${PETSC_LIBRARIES}")
        if(PETSC_IS_COMPLEX)
            target_compile_definitions(PETSc::PETSc INTERFACE -DHPGEM_USE_COMPLEX_PETSC -DHPGEM_USE_ANY_PETSC)
        else()
            target_compile_definitions(PETSc::PETSc INTERFACE -DHPGEM_USE_PETSC -DHPGEM_USE_ANY_PETSC)
        endif()

    endif()

else()
    set(hpGEM_USE_ANY_PETSC OFF)
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
    # Create target for easy linking
    if(SLEPC_FOUND)
        add_library(SLEPc::SLEPc INTERFACE IMPORTED)
        target_include_directories(SLEPc::SLEPc INTERFACE "${SLEPC_INCLUDE_DIRS}")
        target_link_libraries(SLEPc::SLEPc INTERFACE "${SLEPC_LIBRARIES}")
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

    # Create target for easy linking
    add_library(SUNDIALS::SUNDIALS INTERFACE IMPORTED)
    target_include_directories(SUNDIALS::SUNDIALS ${SUNDIALS_INCLUDE_DIR})
    target_link_libraries(SUNDIALS::SUNDIALS
        INTERFACE
            "${SUNDIALS_LIB_sundials_cvodes}"
            "${SUNDIALS_LIB_sundials_idas}"
            "${SUNDIALS_LIB_sundials_kinsol}"
            "${SUNDIALS_LIB_sundials_nvecserial}"
    )
endif()
