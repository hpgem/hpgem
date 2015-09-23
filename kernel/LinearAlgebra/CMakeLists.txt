include_directories(${hpGEM_SOURCE_DIR}/kernel/ .)
if(NOT hpGEM_USE_EXTERNAL_BLAS)
    add_library(LinearAlgebra SHARED 
        daxpy.f lsame.f xerbla.f dgemm.f dgemv.f 
        #Next line for LAPACK stuff inverse and LU factorisation
        #First the BLAS required parts
        dger.f dtrsm.f dscal.f dswap.f idamax.f dtrmm.f dtrmv.f
        #Now the actually LAPACK parts
        ieeeck.f iparmq.f dlamch.f dlaswp.f ilaenv.f dgetf2.f dgetrf.f dtrtri.f  dgetri.f dtrti2.f
        zaxpy.f zgemm.f zgemv.f 
        #Next line for LAPACK stuff inverse and LU factorisation
        #First the BLAS required parts
        ztrsm.f zscal.f zswap.f ztrmm.f ztrmv.f
        #Now the actually LAPACK parts
        zlaswp.f zgetf2.f zgetrf.f ztrtri.f  zgetri.f ztrti2.f
		dcabs1.f izamax.f zgeru.f
		
        #These are the extra from LAPACK for solve Ax=B
        dgesv.f dgetrs.f
        #These are the extra from LAPACK for solve Ax=B
        zgesv.f zgetrs.f
		MiddleSizeVector.cpp
		MiddleSizeMatrix.cpp)
elseif(NOT hpGEM_USE_EXTERNAL_LAPACK)
    add_library(LinearAlgebra SHARED 
        #Next line for LAPACK stuff inverse and LU factorisation
        #Now the actually LAPACK parts
        ieeeck.f iparmq.f dlamch.f dlaswp.f ilaenv.f dgetf2.f dgetrf.f dtrtri.f  dgetri.f dtrti2.f
        #Next line for LAPACK stuff inverse and LU factorisation
        #Now the actually LAPACK parts
        zlaswp.f zgetf2.f zgetrf.f ztrtri.f  zgetri.f ztrti2.f
        
        #These are the extra from LAPACK for solve Ax=B
        dgesv.f dgetrs.f
        #These are the extra from LAPACK for solve Ax=B
        zgesv.f zgetrs.f
        MiddleSizeVector.cpp
        MiddleSizeMatrix.cpp)
    target_link_libraries(LinearAlgebra BLAS)
else()
    add_library(LinearAlgebra SHARED 
        MiddleSizeVector.cpp
        MiddleSizeMatrix.cpp)
    target_link_libraries(LinearAlgebra BLAS LAPACK)
endif()
target_link_libraries(LinearAlgebra Logger)
