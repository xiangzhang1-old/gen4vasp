!#define dotiming
!#define debug
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
!
!   charge density: full grid mode
!
!
!   charge density complex
!
!
!   wavefunctions: full grid mode
!
!
!   wavefunctions complex
!
!
!   common definitions
!






      MODULE subrot
      USE prec
      CONTAINS
!************************ SUBROUTINE EDDIAG ****************************
! RCS:  $Id: subrot.F,v 1.10 2003/06/27 13:22:23 kresse Exp kresse $
!
! this subroutine calculates the electronic eigenvalues and
! optionally performes a  sub-space diagonalisation
! i.e. unitary transforms the wavefunctions so that the Hamiltonian
!  becomes diagonal in the subspace spanned by the wavefunctions
! IFLAG:
!  0 only eigenvalues (without  diagonalisation no sub-space matrix)
!  1 only eigenvalues and sub-space matrix (no diagonalisation)
!  2 eigenvalues using diagonalisation of sub-space matrix
!    do not rotate wavefunctions
! 13+
!  3 eigenvalues and sub-space diagonalisation rotate wavefunctions
!    (for 13 no Jacoby algorithm is allowed)
!  4 eigenvalues and sub-space diagonalisation rotate wavefunctions
!    using Loewdin partubation theory (conserves ordering)
!  5 eigenvalues and sub-space diagonalisation + orthogonalization
!    unfortunately this option turns out to slow down the 
!    convergence of IALGO=48
!
!  for 2) CHAM is set to
!  CHAM(n2,n1) = < C(n2) | H | C(n1) >
!
!  28 .Feb 96 gK: changed to an more efficient blocked algorithm
!***********************************************************************

      SUBROUTINE EDDIAG(GRID,GRIDHF,LATT_CUR,NONLR_S,NONL_S,W,WDES, &
     &    LMDIM,CDIJ,CQIJ,IFLAG,LOVERL,LREAL,NBLOCK,SV,IU0,EXHF,LFOCK)
      USE prec
      USE wave_mpi
      USE wave
      USE dfast
      USE lattice
      USE mpimy
      USE mgrid
      USE nonl
      USE nonlr
      USE hamil
      USE constant
      USE jacobi
      USE scala
      USE main_mpi
!      USE fock

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRID
      TYPE (grid_3d)     GRIDHF
      TYPE (latt)        LATT_CUR
      TYPE (nonlr_struct) NONLR_S
      TYPE (nonl_struct) NONL_S
      TYPE (wavespin)    W
      TYPE (wavedes)     WDES
      LOGICAL LGEN, LFOCK

      COMPLEX(q)   SV(GRID%MPLWV,WDES%NCDIJ) ! local potential
      COMPLEX(q) CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ),CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)

      COMPLEX(q)      CSUM
      LOGICAL LREAL,LOVERL, DO_REDIS
! work arrays for ZHEEV (blocksize times number of bands)
      PARAMETER  (LWORK=32)
      COMPLEX(q)       CWRK(LWORK*WDES%NB_TOT)
      REAL(q)    R(WDES%NB_TOT)
      REAL(q)    RWORK(7*WDES%NB_TOT) ; INTEGER IWORK(5*WDES%NB_TOT), INFO(WDES%NB_TOT)
! work arrays (do max of 16 strips simultaneously)
      PARAMETER (NSTRIPD=16)

      TYPE (wavedes1)    WDES1          ! descriptor for (1._q,0._q) k-point
      TYPE (wavefun1)    W1             ! current wavefunction

      COMPLEX(q),ALLOCATABLE,TARGET::  CPROW(:,:),CPROW_S(:,:),CHAM(:,:),COVL(:,:)
      COMPLEX(q), ALLOCATABLE::  CFOCK(:,:)
      COMPLEX(q),ALLOCATABLE,TARGET::   CR(:),CH(:,:),CVR(:)
! redistributed plane wave coefficients
      COMPLEX(q), POINTER :: CW_RED(:,:),CH_RED(:,:)
      COMPLEX(q)   , POINTER :: CPROW_RED(:,:),CPROJ_RED(:,:),CPROW_RED_S(:,:)
      TYPE (REDIS_PW_CTR),POINTER :: H_PW1, H_PW2
      INTEGER :: icall=0
      icall=icall+1

      NODE_ME=0
      IONODE =0
      NCPU=1
!-----------------------------------------------------------------------
! determine whether redistribution is required
!-----------------------------------------------------------------------
      IF (NCPU /= 1) THEN

        DO_REDIS=.TRUE.
        NRPLWV_RED=WDES%NRPLWV/NCPU
        NPROD_RED =WDES%NPROD /NCPU

      ELSE

        DO_REDIS=.FALSE.
        NRPLWV_RED=WDES%NRPLWV
        NPROD_RED =WDES%NPROD

      ENDIF
      NB_TOT=WDES%NB_TOT
      NBANDS=WDES%NBANDS

      ! set NSTRIP between [1 and 32]
      NSTRIP=MAX(MIN(NSTRIPD,32/NCPU,NBANDS),1)

! allocate work space
      ALLOCATE(CPROW(WDES%NPROD,NBANDS), &
     &     CHAM(NB_TOT,NB_TOT), &
     &     CR(GRID%MPLWV*WDES%NRSPINORS),CVR(GRID%MPLWV*WDES%NRSPINORS), CH(WDES%NRPLWV,NSTRIP))

      IF (IFLAG==5) THEN
         ALLOCATE(CPROW_S(WDES%NPROD,NBANDS),COVL(NB_TOT,NB_TOT))
      ENDIF
      
      IF (LFOCK) THEN
         EXHF=0._q
         ALLOCATE(CFOCK(NB_TOT,NB_TOT))
      ENDIF

      W1%CR=>CR

!=======================================================================
      spin:  DO ISP=1,WDES%ISPIN
      kpoint: DO NK=1,WDES%NKPTS
!=======================================================================
! at this point calculate the Fock matrix if required 
      IF (LFOCK) THEN
!         CALL FOCK_EXOP(GRID, GRIDHF, LATT_CUR, W, WDES, NK, ISP, &
!              LOVERL, CFOCK, EXHF)
      ENDIF
      CALL SETWDES(WDES,WDES1,NK); CALL SETWGRID(WDES1,GRID)

!   get pointers for redistributed wavefunctions
!   I can not guarantee that this works with all f90 compilers
!   please see comments in wave_mpi.F
      IF (DO_REDIS) THEN
        CALL SET_WPOINTER(CW_RED,   NRPLWV_RED, NB_TOT, W%CPTWFP(1,1,NK,ISP))
        CALL SET_WPOINTER(CH_RED,   NRPLWV_RED, NCPU*NSTRIP, CH(1,1))
        CALL SET_GPOINTER(CPROJ_RED, NPROD_RED, NB_TOT, W%CPROJ(1,1,NK,ISP))
        CALL SET_GPOINTER(CPROW_RED, NPROD_RED, NB_TOT, CPROW(1,1))
        IF (IFLAG==5) &
          CALL SET_GPOINTER(CPROW_RED_S, NPROD_RED, NB_TOT, CPROW_S(1,1))
      ELSE
        CW_RED    => W%CPTWFP(:,:,NK,ISP)
        CH_RED    => CH(:,:)
        CPROJ_RED => W%CPROJ(:,:,NK,ISP)
        CPROW_RED => CPROW(:,:)
        IF (IFLAG==5) CPROW_RED_S => CPROW_S(:,:)
      ENDIF

!   set number of wavefunctions after redistribution
      NPL = WDES1%NPL     ! number of plane waves/node after data redistribution
      NPRO= WDES1%NPRO    ! number of projected wavef. after data redistribution
      NPRO_O=NPRO
      IF (.NOT. LOVERL) NPRO_O=0

      

      NGVECTOR=WDES1%NGVECTOR
!=======================================================================
!  IFLAG=0 calculate eigenvalues  only
!=======================================================================
      IF (IFLAG==0) THEN
         DO N=1,NBANDS
                                ! transform wavefunction to real space
                                ! and calculate eigenvalues calling ECCP, no redistribution !
            CALL SETWAV_(W,W1,N,NK,ISP) ! allocation for W1%CR done above
            DO ISPINOR=0,WDES%NRSPINORS-1
               CALL FFTWAV(NGVECTOR, WDES1%NINDPW(1), W1%CR(1+ISPINOR*GRID%MPLWV), W1%CPTWFP(1+ISPINOR*NGVECTOR), GRID)
            ENDDO
            CALL ECCP(WDES1,W1,W1,LMDIM,CDIJ(1,1,1,ISP),GRID,SV(1,ISP), W1%CELEN)
            W%CELEN(N,NK,ISP)=W1%CELEN
         ENDDO
         
         CYCLE kpoint
         
      ENDIF
!=======================================================================
!  IFLAG /= 0 calculate Hamiltonian CHAM
!=======================================================================
!  caclulate D |cfin_n> (D = non local strength of PP)
      IF (DO_REDIS .AND. LASYNC) THEN
        CALL REDIS_PW_ALLOC(WDES, NSTRIP, H_PW1)
        CALL REDIS_PW_ALLOC(WDES, NSTRIP, H_PW2)
        DO NPOS=1,NSTRIP
           CALL REDIS_PW_START(WDES, W%CPTWFP(1,NPOS,NK,ISP), NPOS, H_PW1)
        ENDDO
      ENDIF

      CALL OVERL(WDES1, .TRUE.,LMDIM,CDIJ(1,1,1,ISP), W%CPROJ(1,1,NK,ISP),CPROW)
      IF (IFLAG==5) THEN
         CALL OVERL(WDES1, LOVERL,LMDIM,CQIJ(1,1,1,ISP), W%CPROJ(1,1,NK,ISP),CPROW_S(1,1))
      ENDIF

    ! redistribute the projected wavefunctions
    ! wavefunctions are still required at this point
      IF (DO_REDIS) THEN
        CALL REDIS_PROJ(WDES1, NBANDS, CPROW(1,1))
        CALL REDIS_PROJ(WDES1, NBANDS, W%CPROJ(1,1,NK,ISP))
        IF (IFLAG==5) CALL REDIS_PROJ(WDES1, NBANDS, CPROW_S(1,1))
      ENDIF

      CHAM=0
      NDONE=0

  strip: DO NPOS=1,NBANDS,NSTRIP
        NSTRIP_ACT=MIN(NBANDS+1-NPOS,NSTRIP)

       !  calculate V_{local} |phi> + T | phi >
       !  for a block containing NSTRIP wavefunctions

        NGVECTOR=WDES1%NGVECTOR
        DO N=NPOS,NPOS+NSTRIP_ACT-1
          NDONE=NDONE+1
          NP=N-NPOS+1
          ! fft to real space
          DO ISPINOR=0,WDES%NRSPINORS-1
             CALL FFTWAV(NGVECTOR,WDES%NINDPW(1,NK),CR(1+ISPINOR*GRID%MPLWV),W%CPTWFP(1+ISPINOR*NGVECTOR,N,NK,ISP),GRID)
          ENDDO

          CALL VHAMIL(WDES1, GRID, SV(1,ISP), CR(1), CVR(1))
          
          ! fft of V(r) * phi(r) and  add the kinetic energy term
          spinors: DO ISPINOR=0,WDES%NRSPINORS-1
          CALL FFTEXT(NGVECTOR,WDES%NINDPW(1,NK),CVR(1+ISPINOR*GRID%MPLWV),CH(1+ISPINOR*NGVECTOR,NP),GRID,.FALSE.)
          DO M=1,NGVECTOR
             MM=M+ISPINOR*NGVECTOR
!-MM- changes to accommodate spin spirals
! original statement
!            CH(MM,NP)=CH(MM,NP)+W%CPTWFP(MM,N,NK,ISP)* WDES%DATAKE(M,NK)
             CH(MM,NP)=CH(MM,NP)+W%CPTWFP(MM,N,NK,ISP)* WDES%DATAKE(M,NK,ISPINOR+1)
!-MM- end of alterations
          ENDDO
          ENDDO spinors
          IF (DO_REDIS.AND. LASYNC) CALL REDIS_PW_START(WDES, CH(1,NP), N, H_PW2)
        ENDDO
    ! redistribute wavefunctions
    ! W%CPTWFP is then redistributed up to and including 1...NPOS+NSTRIP_ACT
      IF (DO_REDIS) THEN
        IF (LASYNC) THEN
          DO N=NPOS,NPOS+NSTRIP_ACT-1
             NP=N-NPOS+1
             CALL REDIS_PW_STOP (WDES, W%CPTWFP(1,N,NK,ISP), N, H_PW1)
             IF (N+NSTRIP<=NBANDS) &
             CALL REDIS_PW_START(WDES, W%CPTWFP(1,N+NSTRIP,NK,ISP), N+NSTRIP, H_PW1)
             CALL REDIS_PW_STOP (WDES, CH(1,NP), N, H_PW2)
          ENDDO
        ELSE
           CALL REDIS_PW(WDES1, NSTRIP_ACT, W%CPTWFP(1,NPOS,NK,ISP))
           CALL REDIS_PW(WDES1, NSTRIP_ACT, CH(1,1))
        ENDIF
      ENDIF

      NPOS_RED  =(NPOS-1)*NCPU+1
      NSTRIP_RED=NSTRIP_ACT*NCPU

      CALL ORTH1('U', &
        CW_RED(1,1),CH_RED(1,1),CPROJ_RED(1,1), &
        CPROW_RED(1,NPOS_RED),NB_TOT,NBLOCK, &
        NPOS_RED, NSTRIP_RED, NPL,NPRO,NRPLWV_RED,NPROD_RED,CHAM(1,1))

  ENDDO strip

      IF (DO_REDIS .AND. LASYNC) THEN
         CALL REDIS_PW_DEALLOC(H_PW1)
         CALL REDIS_PW_DEALLOC(H_PW2)
      ENDIF

!
!  just for safety everything ok ?
!  Mind that all wavefunctions are redistributed at this point
      IF (NDONE/=NBANDS .OR. NPOS_RED+NSTRIP_RED/=NB_TOT+1 ) THEN
       WRITE(*,*)'EDDIAG: fatal internal error (2):',NDONE,NBANDS,NPOS_RED+NSTRIP_RED,NB_TOT+1
       STOP
    ENDIF
      
      
      
! now CHAM and CFOCK are known, add them to CHAM (negative sign(!))
! every proc does this job because od the different distribution of
! CFOCK and CHAM before distribution to all 

      IF (LFOCK) THEN
         CHAM(:,:)=CHAM(:,:)-CFOCK(:,:)
!         CHAM(:,:)=CFOCK(:,:)
      ENDIF

!-----------------------------------------------------------------------
! calculate the overlap matrix
!-----------------------------------------------------------------------
      IF (IFLAG==5) THEN
        DO N=1,NB_TOT
        DO I=1,NB_TOT
          COVL(I,N)=(0._q,0._q)
        ENDDO; ENDDO

        ! set NSTRIP between [4 and 32] (approx NB_TOT/20)
        NSTRIP=2*(NB_TOT/20)
        NSTRIP=MAX(NSTRIP,4)
        NSTRIP=MIN(NSTRIP,32,NB_TOT)

        DO NPOS=1,NB_TOT-NSTRIP,NSTRIP
        CALL ORTH1('U',CW_RED(1,1),CW_RED(1,NPOS),CPROJ_RED(1,1), &
          CPROW_RED_S(1,NPOS),NB_TOT,NBLOCK, &
          NPOS,NSTRIP,NPL,NPRO_O,NRPLWV_RED,NPROD_RED,COVL(1,1))
        ENDDO

        CALL ORTH1('U',CW_RED(1,1),CW_RED(1,NPOS),CPROJ_RED(1,1), &
          CPROW_RED_S(1,NPOS),NB_TOT,NBLOCK, &
          NPOS,NB_TOT-NPOS+1,NPL,NPRO_O,NRPLWV_RED,NPROD_RED,COVL(1,1))
        
      ENDIF
!-----------------------------------------------------------------------
! dump some arrays
!-----------------------------------------------------------------------
!=======================================================================
! IFLAG =4 use Loewdin perturbation to get rotation matrix
! this preserves the ordering of the eigenvalues
! MIND: does not work for real matrices
!=======================================================================
      IF (IFLAG==4) THEN
      DIFMAX=0.001_q

      DO N2=1,NB_TOT
      DO N1=1,N2-1
         DIFCEL= REAL( CHAM(N2,N2)-CHAM(N1,N1) ,KIND=q)
         IF (ABS(DIFCEL)<DIFMAX) THEN
           CROT  =0
         ELSE
           CROT  =CONJG(CHAM(N1,N2))/DIFCEL
           IF (ABS(CROT)>0.1_q) THEN
             FAKT= 0.1_q/ABS(CROT)
             CROT  = CROT*FAKT
           ENDIF
         ENDIF
        CHAM(N2,N1) =-CROT
        CHAM(N1,N2) =-CONJG(CROT)
      ENDDO
      ENDDO
      DO N1=1,NB_TOT
         CHAM(N1,N1)=1
         W%CELTOT(N1,NK,ISP)=CHAM(N1,N1)
      ENDDO

      CALL  ORSP(NB_TOT,NB_TOT,NB_TOT,CHAM)
      ELSE
!=======================================================================
! IFLAG > 1 and IFLAG <4
! diagonalization of CHAM
! we have lots of choices for the parallel version that makes things
! rather complicated
! to allow for reasonable simple programming once the diagonalisation
! has been done I jump to line 100
!=======================================================================
      IF (IFLAG==1) GOTO 1000

      DO N1=1,NB_TOT
        IF (ABS(AIMAG(CHAM(N1,N1)))>1E-2_q .AND. IU0>=0) THEN
          WRITE(IU0,*)'WARNING: Sub-Space-Matrix is not hermitian subr', &
     &              AIMAG(CHAM(N1,N1)),N1
        ENDIF
        CHAM(N1,N1)= REAL( CHAM(N1,N1) ,KIND=q)
      ENDDO

!
! parallel versions
! if fast Jacobi method exists use it (T3D, T3E only)
! use the first line in that case
      

      IFAIL=0

!
!  seriell codes
!
      IF (IFLAG == 5) THEN
         CALL ZHEGV &
              (1,'V','U',NB_TOT,CHAM(1,1),NB_TOT,COVL(1,1),NB_TOT, &
              R,CWRK,LWORK*NB_TOT, RWORK,  IFAIL)
      ELSE
         ABSTOL=1E-10_q
         VL=0 ; VU=0 ; IL=0 ; IU=0
         ALLOCATE(COVL(NB_TOT,NB_TOT))
         CALL ZHEEVX( 'V', 'A', 'U', NB_TOT, CHAM(1,1) , NB_TOT, VL, VU, IL, IU, &
                            ABSTOL , NB_CALC , R, COVL(1,1), NB_TOT, CWRK, &
                            LWORK*NB_TOT, RWORK, IWORK, INFO, IFAIL )         
         CHAM=COVL
         DEALLOCATE(COVL)
      ENDIF
      ! T3D uses a global sum which does not guarantee to give the same results on all 
      ! the following line is required to make the code waterproof (we had problems)
      ! since we now use a propritary sum (see mpi.F) we should not require
      ! this broadcast anymore
      ! 
      ! 

 100  CONTINUE
      

      IF (IFAIL/=0) THEN
         WRITE(*,*) 'ERROR EDDIAG: Call to routine ZHEEV failed! '// &
     &              'Error code was ',IFAIL
         STOP
      ENDIF
      
      DO N1=1,NB_TOT
        W%CELTOT(N1,NK,ISP)=R(N1)
      ENDDO
      ENDIF


!=======================================================================
! IFLAG > 2
! rotate wavefunctions
!=======================================================================

      IF (IFLAG==2) GOTO 1000

      CALL LINCOM('F',CW_RED(1,1),CPROJ_RED(1,1),CHAM(1,1), &
             NB_TOT,NB_TOT,NPL,NPRO,NRPLWV_RED,NPROD_RED,NB_TOT, &
             NBLOCK,CW_RED(1,1),CPROJ_RED(1,1))
      ! "lincom ok"

 1000 CONTINUE

     !  back redistribution

      IF (DO_REDIS) THEN
        CALL REDIS_PROJ(WDES1, NBANDS, W%CPROJ(1,1,NK,ISP))
        IF (LASYNC) THEN
           W%OVER_BAND=.TRUE.
        ELSE
           CALL REDIS_PW  (WDES1, NBANDS, W%CPTWFP   (1,1,NK,ISP))
        ENDIF
        ! "redis ok"
      ENDIF

!=======================================================================
      ENDDO kpoint
      ENDDO spin
!=======================================================================
      DEALLOCATE(CPROW,CHAM,CR,CVR,CH)
      IF (IFLAG==5) DEALLOCATE(CPROW_S,COVL)
      IF (LFOCK) THEN
         DEALLOCATE(CFOCK)
         WRITE(*,*)
      ENDIF
      RETURN
      END SUBROUTINE
      END MODULE
