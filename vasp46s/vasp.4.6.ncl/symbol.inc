!-------- to be costumized by user (usually done in the makefile)-------
!#define vector              compile for vector machine
!#define essl                use ESSL instead of LAPACK
!#define single_BLAS         use single prec. BLAS

!#define wNGXhalf            gamma only wavefunctions (X-red)
!#define wNGZhalf            gamma only wavefunctions (Z-red)

!#define NGXhalf             charge stored in REAL array (X-red)
!#define NGZhalf             charge stored in REAL array (Z-red)
!#define NOZTRMM             replace ZTRMM by ZGEMM
!#define REAL_to_DBLE        convert REAL() to DBLE()
!#define MPI                 compile for parallel machine with MPI
!------------- end of user part         --------------------------------
#ifdef single_BLAS
!
!   use single precission BLAS and LAPACK calls (usually on Cray)
!
#define ZGEMM  CGEMM
#define ZGEMV  CGEMV
#define ZDSCAL CSSCAL
#define ZDOTC  CDOTC
#define ZAXPY  CAXPY
#define ZHEEV  CHEEV
#define ZHEGV  CHEGV
#define ZTRMM  CTRMM
#define ZTRTRI CTRTRI
#define ZPOTRF CPOTRF
#define PZPOTRF PCPOTRF
#define PZTRTRI PCTRTRI
#define PZHEEVX PCHEEVX
#define ZCOPY   CCOPY
#define ZGETRF  CGETRF
#define ZGETRS  CGETRS
#define ZGEFA   CGEFA
#define ZGECO   CGECO
#define ZGEDI   CGEDI
#define ZHEEVX  CHEEVX

#define DGEMM  SGEMM
#define DGEMV  SGEMV
#define DGEGV  SGEGV
#define DDOT   SDOT
#define DSCAL  SSCAL
#define DAXPY  SAXPY
#define DSYEV  SSYEV
#define DSYGV  SSYGV
#define DTRMM  STRMM
#define DTRTRI STRTRI
#define DPOTRF SPOTRF
#define DGETRF SGETRF
#define DGETRS SGETRS
#define PDPOTRF PSPOTRF
#define PDTRTRI PSTRTRI
#define PDSYEVX PSSYEVX
#define DCOPY   SCOPY
#define DSYEVX  SSYEVX

#endif
#ifdef NGZhalf
!
!   charge density: half grid mode Z direction
!
#define realmode
#define  zm(n)             (n/2+1)
#define  xm(n)              n
#define  xmin(n)          (-(n/2-1))
#define  SETN1fromN1P     N1=MOD(N1P+xmin(GRIDC%NGX)+GRIDC%NGX,GRIDC%NGX)
#define  SETFACT1         FACTM=1
#define  SETFACT          IF (N3 /= 1) FACTM=2
#define  SETSQRT          IF (N3 /= 1) FACTM=SQRT(2.)
#define  SETISQRT         IF (N3 /= 1) FACTM=1./SQRT(2.)
#define  MULFACT          FACTM*
#elif defined(NGXhalf)
!
!   charge density: half grid mode X direction
!
#define realmode
#define  xm(n)             (n/2+1)
#define  xmin(n)            0
#define  zm(n)              n
#define  SETN1fromN1P     N1=N1P
#define  SETFACT1         FACTM=1
#define  SETFACT          IF (N1 /= 1) FACTM=2
#define  SETSQRT          IF (N1 /= 1) FACTM=SQRT(2.)
#define  SETISQRT         IF (N1 /= 1) FACTM=1./SQRT(2.)
#define  MULFACT          FACTM*
#else
!
!   charge density: full grid mode
!
#define  zm(n)              n
#define  xm(n)              n
#define  xmin(n)          (-(n/2-1))
#define  SETN1fromN1P      N1=MOD(N1P+xmin(GRIDC%NGX)+GRIDC%NGX,GRIDC%NGX)
#define  SETFACT1          FACTM=1
#define  SETFACT
#define  SETSQRT
#define  SETISQRT
#define  MULFACT
#endif
#ifdef realmode
!
!   charge density real
!
#define RGRID REAL(q)
#define DIMREAL(n)  n*2
#else
!
!   charge density complex
!
#define RGRID COMPLEX(q)
#define DIMREAL(n)  n
#endif
#ifdef wNGZhalf
!
!   wavefunctions: half grid mode for Z direction
!
#define gammareal
#define  zw(n)            (n/2+1)
#define  yw(n)             n
#define  xw(n)             n
#define  WSETFACT1        FACTM=1
#define  WSETFACT(x)      IF (N1 /= 1.OR.N2 /= 1.OR.N3 /= 1) FACTM=x
#define  WMULFACT         FACTM*
#elif defined(wNGXhalf)
!
!   wavefunctions: half grid mode for X direction
!
#define gammareal
#define  xw(n)             (n/2+1)
#define  yw(n)              n
#define  zw(n)              n
#define  WSETFACT1        FACTM=1
#define  WSETFACT(x)      IF (N1 /= 1.OR.N2 /= 1.OR.N3 /= 1) FACTM=x
#define  WMULFACT         FACTM*
#else
!
!   wavefunctions: full grid mode
!
#define  xw(n)              n
#define  yw(n)              n
#define  zw(n)              n
#define  WSETFACT1         FACTM=1
#define  WSETFACT(x)
#define  WMULFACT
#endif
#ifdef  gammareal
!
!   wavefunctions real (gamma only)
!
#define GREAL(n)  REAL(n,KIND=q)
#define GDEF      REAL(q)
#define M_sum_g   M_sum_d
#define M_sumf_g M_sumf_d
#define M_bcast_g M_bcast_d
#define CPROJ     GPROJ
#define GCONJG(n) (n)
#define GGEMM     DGEMM
#define GGEMV     DGEMV
#define GGETRF    DGETRF
#define GGETRS    DGETRS
#define GDOTC     DDOT
#define GDSCAL    DSCAL
#define GAXPY     DAXPY
#define zero      0._q
#define one       1._q
#define m_        2*
#define ndata     1
#define trans     'T'
#define gammaarg   NPLWKC,NINDPC,WVSCAL,WVSCAI
#define gammaargNK NPLWKC(NK),NINDPC(1,NK),WVSCAL(1,NK),WVSCAI(1,NK)
#else
!
!   wavefunctions complex
!
#define GREAL(n)  (n)
#define GDEF      COMPLEX(q)
#define M_sumf_g  M_sumf_z
#define M_sum_g   M_sum_z
#define M_bcast_g M_bcast_z
#define GCONJG(n) CONJG(n)
#if defined (SCSL)
#define GGEMM     ZGEMM3M
#else
#define GGEMM     ZGEMM
#endif
#define GGEMV     ZGEMV
#define GGETRF    ZGETRF
#define GGETRS    ZGETRS
#define GDOTC     ZDOTC
#define GDSCAL    ZDSCAL
#define GAXPY     ZAXPY
#define zero      (0._q,0._q)
#define one       (1._q,0._q)
#define trans     'C'
#define m_
#define ndata     2
#define gammaarg
#define gammaargNK
#endif
!
!   common definitions
!
#ifdef realmode
#define  OVERLAP          REAL(q)
#else
#define  OVERLAP          COMPLEX(q)
#endif
#ifdef  REAL_to_DBLE
#define REAL(x) DBLE(x)
#endif
#if defined(MPI)
!
!   mpi parallel macros
!
#define CALLMPI(x) CALL x
#define CALLMPI_C(x) CALL x
#define STOP       CALL M_exit(); stop

#define io_begin IF (NODE_ME==IONODE) THEN
#define do_io    IF (NODE_ME==IONODE)
#define io_end   ENDIF
#elif defined(MPI_CHAIN)
!
!   mpi nudged chain version
!
#define CALLMPI(x)
#define CALLMPI_C(x) CALL x
#define STOP       CALL M_exit(); STOP

#define io_begin IF (NODE_ME==IONODE) THEN
#define do_io    IF (NODE_ME==IONODE)
#define io_end   ENDIF
#else
#define CALLMPI(x)
#define CALLMPI_C(x)
#define nodes
#define do_io
#define io_begin
#define io_end
#endif
#ifdef F90_T3D
!
!   T3D macros
!
!#define IVDEP
#endif
#if defined(T3D_SMA)
!
!   T3D SHMEM communication
!
#endif
#ifdef CRAY_MPP
!
!   CRAY_MPP macros
!
#define MPI_double_precision MPI_real
#define MPI_double_complex   MPI_complex
#endif
#ifdef USE_ERF
!
!  use instrinsic ERF
!
#define ERRF(x) ERF(x)
#define ERRFC(x) ERFC(x)
#endif

#define dwrite DWRITE
#ifdef debug
!
!  debugging primitives
!
#define DWRITE WRITE(0,*)
#define DWRITE0 do_io WRITE(0,*)
#else
#define DWRITE  !
#define DWRITE0 !
#endif

#define CW CPTWFP

#ifdef nonlr_single
!
! single precession non local projection operators  
!
#define RPROJ   RPROJS
#define qn      qs
#define DGEMVn  SGEMV
#define DGEMMn  SGEMM
#else
#define  qn    q
#define DGEMVn  DGEMV
#define DGEMMn  DGEMM
#endif

#if  ! (defined(NGXhalf)  || defined(NGYhalf) ||  defined(NGZhalf))
#define RHOLM_complex
#endif

#ifdef dotiming
!
!  some timing primitives
!
#define usetmg REAL(q) :: TV,TC,TV0,TC0
#define statmg CALLMPI( MPI_barrier( WDES%COMM%MPI_COMM, ierror )) ;  CALL VTIME(TV0,TC0)
#define stotmg CALLMPI( MPI_barrier( WDES%COMM%MPI_COMM, ierror )) ;  CALL VTIME(TV,TC) ; IF (NODE_ME==IONODE) WRITE(*,'(A,2F10.4)') 'timing ',TV-TV0,TC-TC0 ; CALL VTIME(TV0,TC0)
#define stotmgt(X) CALLMPI( MPI_barrier( WDES%COMM%MPI_COMM, ierror )) ;  CALL VTIME(TV,TC) ; IF (NODE_ME==IONODE) WRITE(*,'(A,2F10.4)') X,TV-TV0,TC-TC0 ; CALL VTIME(TV0,TC0)
#else
#define usetmg 
#define statmg 
#define stotmg 
#define stotmgt(X)
#endif
