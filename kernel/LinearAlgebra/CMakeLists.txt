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
		
	#Next line for LAPACK LDL factorisation for symmetric real matrices
	disnan.f dlaisnan.f dlasyf.f dsytf2.f dsytrf.f

	#Next line for LAPACK least square solver DGELSS
	dbdsqr.f dgebd2.f dgebrd.f dgelq2.f dgelqf.f dgelss.f dgeqr2.f dgeqrf.f
	dlabad.f dlabrd.f dlacpy.f dlange.f dlapy2.f dlarf.f dlarfb.f dlarfg.f
	dlarft.f dlartg.f dlas2.f dlascl.f dlaset.f dlasq1.f dlasq2.f dlasq3.f
	dlasq4.f dlasq5.f dlasq6.f dlasr.f dlasrt.f dlassq.f dlasv2.f dorg2r.f
	dorgbr.f dorgl2.f dorglq.f dorgqr.f dorm2r.f dormbr.f dorml2.f dormlq.f
	dormqr.f drscl.f iladlc.f iladlr.f 

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
		
	#Next line for LAPACK LDL factorisation for symmetric real matrices
	disnan.f dlaisnan.f dlasyf.f dsytf2.f dsytrf.f
        
	#Next line for LAPACK least square solver DGELSS
	dbdsqr.f dgebd2.f dgebrd.f dgelq2.f dgelqf.f dgelss.f dgeqr2.f dgeqrf.f
	dlabad.f dlabrd.f dlacpy.f dlange.f dlapy2.f dlarf.f dlarfb.f dlarfg.f
	dlarft.f dlartg.f dlas2.f dlascl.f dlaset.f dlasq1.f dlasq2.f dlasq3.f
	dlasq4.f dlasq5.f dlasq6.f dlasr.f dlasrt.f dlassq.f dlasv2.f dorg2r.f
	dorgbr.f dorgl2.f dorglq.f dorgqr.f dorm2r.f dormbr.f dorml2.f dormlq.f
	dormqr.f drscl.f iladlc.f iladlr.f 

        #These are the extra from LAPACK for solve Ax=B
        dgesv.f dgetrs.f
        #These are the extra from LAPACK for solve Ax=B
        zgesv.f zgetrs.f
        MiddleSizeVector.cpp
        MiddleSizeMatrix.cpp)
    target_link_libraries(LinearAlgebra ${BLAS_LIBRARIES})
else()
    add_library(LinearAlgebra SHARED 
        MiddleSizeVector.cpp
        MiddleSizeMatrix.cpp)
    target_link_libraries(LinearAlgebra ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
endif()
target_link_libraries(LinearAlgebra Logger)
