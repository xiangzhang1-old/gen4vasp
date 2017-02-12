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





      MODULE rmm_diis
      USE prec
      CONTAINS
!************************ SUBROUTINE EDDRMM *****************************
! RCS:  $Id: rmm-diis.F,v 1.7 2002/08/14 13:59:42 kresse Exp $
!
! this subroutine performes an optimization of the trial wavefunctions
! minimizing the expactation value of the Hamiltonian  i.e.
!     < phi | H |  phi >
! or the norm of the residual vector  i.e.
!    r^2 = < phi | H -e S | H - e S | phi >
! or using an inverse iteration method. In the last case
!    || ( H - e_initial S)  | phi > - | phi_initial > ||
! is optimized.
! The full name of the residual vector  minimization method
! is residual vector minimiziation method-
! direct inversion of the iterative subspace (RMM-DIIS)
!    see: D. M. Wood and A. Zunger, J. Phys. A, 1343 (1985)
!    and  P. Pulay,  Chem. Phys. Lett. 73, 393 (1980).
!
!
!  INFO%IALGO   determine type of preconditioning and the algorithm
!    6    rms-minimization          +  TAP preconditioning
!    7    rms-minimization          +  no preconditioning
!    8    precond rms-minimization  +  TAP preconditioning
!    9    precond rms-minimization  +  Jacobi like preconditioning
!    (TAP Teter Alan Payne)
!   LDELAY=.TRUE.
!          steepest descent eigenvalue minimization
!          maximum number of steps is 2
!  WEIMIN  treshhold for total energy minimisation
!    is the fermiweight of a band < WEIMIN,
!    minimisation will break after a maximum of two iterations
!  EBREAK  absolut break condition
!    intra-band minimisation is stopped if DE is < EBREAK
!  DEPER   intra-band break condition (see below)
!  ICOUEV  number of intraband evalue minimisations
!  DESUM   total change in eigenvalues
!  RMS     norm of residual vector
!
!***********************************************************************

      SUBROUTINE EDDRMM(GRID,INFO,LATT_CUR,NONLR_S,NONL_S,W,WDES, &
        LMDIM,CDIJ,CQIJ, RMS,DESUM,ICOUEV, SV,IU6,IU0, LDELAY)
      USE prec

      USE wave
      USE dfast
      USE base
      USE lattice
      USE mpimy
      USE mgrid
      USE nonl
      USE nonlr
      USE hamil
      USE constant
      USE wave_mpi
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRID
      TYPE (info_struct) INFO
      TYPE (latt)        LATT_CUR
      TYPE (nonlr_struct) NONLR_S
      TYPE (nonl_struct) NONL_S
      TYPE (wavespin)    W
      TYPE (wavedes)     WDES

      LOGICAL LDELAY
      COMPLEX(q)   SV(GRID%MPLWV,WDES%NCDIJ) ! local potential
      COMPLEX(q) CDIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ), CQIJ(LMDIM,LMDIM,WDES%NIONS,WDES%NCDIJ)

!----- local work arrays
      TYPE (wavedes1)    WDES1          ! descriptor for (1._q,0._q) k-point
      TYPE (wavefun1)    W1(WDES%NSIM)  ! current wavefunction
      TYPE (wavefun1)    WTMP(WDES%NSIM)! temporary trial wavefunction
      REAL(q) R(INFO%NDAV)
      PARAMETER  (LWORK=20)
      DIMENSION CWORK(LWORK*INFO%NDAV)
      DIMENSION RWORK(3*INFO%NDAV)
! work arrays
      COMPLEX(q),ALLOCATABLE:: CF(:,:,:),CF_INI(:,:)
      COMPLEX(q),ALLOCATABLE::    CPROF(:,:,:)
      REAL(q),ALLOCATABLE:: PRECON(:,:)
      COMPLEX(q),ALLOCATABLE:: CHAM(:,:),CTMP(:,:),CWORK1(:)
      COMPLEX(q),ALLOCATABLE:: AM(:,:),B(:),AM_(:,:),B_(:)
      INTEGER,ALLOCATABLE :: IPIV(:)
      INTEGER :: NB(WDES%NSIM)         ! contains a list of bands currently optimized
      REAL(q) :: EVALUE_INI(WDES%NSIM) ! eigenvalue of that band at the beginning
      REAL(q) :: EVALUE(WDES%NSIM)     ! eigenvalue during optimization
      REAL(q) :: DEIT(WDES%NSIM)       ! relative break criterion for that band
      REAL(q) :: IT(WDES%NSIM)         ! current iteration for this band
      LOGICAL :: LDO(WDES%NSIM)        ! band finished
      REAL(q) :: FPRE(WDES%NSIM)       ! norm of residual vector for each band
      REAL(q) :: TRIAL(WDES%NSIM)      ! trial step for each band
      LOGICAL :: LSTOP
      TYPE (REDIS_PW_CTR),POINTER :: H_PW


      NSIM=WDES%NSIM

      NODE_ME=0
      IONODE =0
!=======================================================================
!  INITIALISATION:
! maximum  number of iterations
! NRES position where H - E S| trial vector > is stored
!=======================================================================
      NITER=INFO%NDAV
      IF (LDELAY) NITER=MIN(NITER,1)
      IF (LDELAY .AND. INFO%IALGO ==0) NITER=1
      NRES =INFO%NDAV

      DESUM =0
      RMS   =0
      ICOUEV=0

      SLOCAL=0
      DO I=1,GRID%RL%NP
        SLOCAL=SLOCAL+SV(I,1)
      ENDDO

      
      SLOCAL=SLOCAL/GRID%NPLWV

      ALLOCATE(PRECON(WDES%NRPLWV,NSIM),CF_INI(WDES%NRPLWV,NSIM), &
     &        CF(WDES%NRPLWV,NRES*2,NSIM),CPROF(WDES%NPROD,NRES*2,NSIM), &
     &        CHAM(NRES,NRES),CTMP(NRES,NRES),CWORK1(NRES), &
     &        AM(NRES,NRES),B(NRES),AM_(NRES,NRES),B_(NRES),IPIV(NRES))
      LD=WDES%NRPLWV*NRES*2

      DO NP=1,NSIM
         CALL NEWWAV(WTMP(NP), WDES, GRID%MPLWV*WDES%NRSPINORS, .TRUE.)
         ALLOCATE(W1(NP)%CR(GRID%MPLWV*WDES%NRSPINORS))
      ENDDO
      CTMP=0
      CHAM=0

!=======================================================================
! do we have to distribute the wavefunctions back ?
!=======================================================================
      IF (W%OVER_BAND) THEN
         NCPU=1
         NSTRIP=MIN(NSIM,WDES%NBANDS)*2
         CALL REDIS_PW_ALLOC(WDES, NSTRIP, H_PW)
      ENDIF

!=======================================================================
      spin:    DO ISP=1,WDES%ISPIN
      kpoints: DO NK=1,WDES%NKPTS
!=======================================================================
!  first initiate communication between bands
      IF (W%OVER_BAND) THEN
         DO N=1,NSTRIP
            CALL REDIS_PW_START(WDES, W%CPTWFP(1,N,NK,ISP), N, H_PW)
         ENDDO
      ENDIF
!=======================================================================
      DE_ATT=ABS(W%CELEN(WDES%NBANDS,NK,ISP)-W%CELEN(1,NK,ISP))/4

      CALL SETWDES(WDES,WDES1,NK); CALL SETWGRID(WDES1,GRID)

      NPL=WDES1%NPL
      NGVECTOR=WDES1%NGVECTOR

      IF (INFO%LREAL) THEN
        CALL PHASER(GRID,LATT_CUR,NONLR_S,NK,WDES,0.0_q,0.0_q,0.0_q)
      ELSE
        CALL PHASE(WDES,NONL_S,NK)
      ENDIF
      NB=0          ! empty the list of bands, which are optimized currently
      NB_DONE=0     ! index the bands allready optimised
!=======================================================================
      bands: DO
!
!  check the NB list, whether there is any empty slot
!  fill in a not yet optimized wavefunction into the slot
!
      newband: DO NP=1,NSIM
      IF (NB(NP)==0 .AND.  NB_DONE < WDES%NBANDS ) THEN
        NB_DONE=NB_DONE+1
        N     =NB_DONE
        NB(NP)=NB_DONE

        IF (W%OVER_BAND) THEN
           CALL REDIS_PW_STOP (WDES, W%CPTWFP(1,N,NK,ISP), N, H_PW)
           IF (N+NSTRIP<=WDES%NBANDS) &
           CALL REDIS_PW_START(WDES, W%CPTWFP(1,N+NSTRIP,NK,ISP), N+NSTRIP, H_PW)
        ENDIF

        CALL SETWAV_(W,W1(NP),N,NK,ISP)  ! fill band N into W1(NP)


        DO M=1,NPL
           CF_INI(M,NP)=W1(NP)%CPTWFP(M)
        ENDDO

        IDUMP=0

        IF (IDUMP==2) WRITE(*,'(I3,1X)',ADVANCE='NO') N

        !===============================================================
        ! start with the exact evaluation of the eigenenergy
        !===============================================================
        DO ISPINOR=0,WDES%NRSPINORS-1
           CALL FFTWAV(NGVECTOR, WDES%NINDPW(1,NK),W1(NP)%CR(1+ISPINOR*GRID%MPLWV),W1(NP)%CPTWFP(1+ISPINOR*NGVECTOR),GRID)
        ENDDO
        CALL ECCP(WDES1,W1(NP),W1(NP),LMDIM,CDIJ(1,1,1,ISP),GRID,SV(1,ISP), W%CELEN(N,NK,ISP))
        EVALUE_INI(NP)=W%CELEN(N,NK,ISP)

        IF (IDUMP==2) WRITE(*,'(F9.4)',ADVANCE='NO') REAL( W%CELEN(N,NK,ISP) ,KIND=q)

        !===============================================================
        ! calculate the preconditioning matrix
        !===============================================================
        IF (INFO%IALGO==0 .OR. INFO%IALGO==8 .OR. INFO%IALGO==6 .OR. &
            (INFO%IALGO==9 .AND. LDELAY)) THEN
          EKIN=0
!DIR$ IVDEP
!OCL NOVREL

          DO ISPINOR=0,WDES%NRSPINORS-1
          DO M=1,NGVECTOR
             MM=M+ISPINOR*NGVECTOR
!-MM- changes to accommodate spin spirals
! original statements
!            IF (LDELAY .AND. WDES%DATAKE(M,NK)>INFO%ENINI) W1(NP)%CPTWFP(MM)=0
!            CPT=W1(NP)%CPTWFP(MM)
!            EKIN =EKIN+ REAL( CPT*CONJG(CPT) ,KIND=q) * WDES%DATAKE(M,NK)
             IF (LDELAY .AND. WDES%DATAKE(M,NK,ISPINOR+1)>INFO%ENINI) W1(NP)%CPTWFP(MM)=0
             CPT=W1(NP)%CPTWFP(MM)
             EKIN =EKIN+ REAL( CPT*CONJG(CPT) ,KIND=q) * WDES%DATAKE(M,NK,ISPINOR+1)
!-MM- end of alterations
          ENDDO
          ENDDO
          

          IF (EKIN<2.0_q) EKIN=2.0_q
          EKIN=EKIN*1.5_q
          IF (IDUMP==2)  WRITE(*,'(E9.2,"E")',ADVANCE='NO') EKIN

          FAKT=2._q/EKIN

          DO ISPINOR=0,WDES%NRSPINORS-1
          DO M=1,NGVECTOR
             MM=M+ISPINOR*NGVECTOR
!-MM- changes to accommodate spin spirals
! original statement
!            X=WDES%DATAKE(M,NK)/EKIN
             X=WDES%DATAKE(M,NK,ISPINOR+1)/EKIN
!-MM- end of alterations
             X2= 27+X*(18+X*(12+8*X))
             PRECON(MM,NP)=X2/(X2+16*X*X*X*X)*FAKT
          ENDDO
          ENDDO
        ELSE IF (INFO%IALGO==9) THEN
          DO ISPINOR=0,WDES%NRSPINORS-1
          DO M=1,NGVECTOR
            MM=M+ISPINOR*NGVECTOR
!-MM- changes to accommodate spin spirals
! original statement
!           X=MAX(WDES%DATAKE(M,NK)+SLOCAL-EVALUE_INI(NP),0._q)
            X=MAX(WDES%DATAKE(M,NK,ISPINOR+1)+SLOCAL-EVALUE_INI(NP),0._q)
!-MM- end of alterations
            PRECON(MM,NP)= REAL( 1._q/(X+ CMPLX( 0 , DE_ATT ,KIND=q) ) ,KIND=q) !new
          ENDDO
          ENDDO
        ELSE
          DO M=1,NPL
             PRECON(M,NP)=1
          ENDDO
        ENDIF
        DEIT(NP)=0
        IT(NP)  =0
!=======================================================================
      ENDIF
      ENDDO newband
!=======================================================================
! if the NB list is now empty end the bands DO loop
!=======================================================================
      LSTOP=.TRUE.
      LDO  =.FALSE.
      DO NP=1,NSIM
         IF ( NB(NP) /= 0 ) THEN
            LSTOP  =.FALSE.
            LDO(NP)=.TRUE.     ! band not finished yet
            IT(NP) =IT(NP)+1   ! increase iteration count
         ENDIF
      ENDDO
      IF (LSTOP) EXIT bands
!=======================================================================
! intra-band minimisation
!=======================================================================
      i1: DO NP=1,NSIM


      N=NB(NP); ITER=IT(NP); IF (.NOT. LDO(NP)) CYCLE i1
!-----------------------------------------------------------------------
! fill current wavefunctions into work arrays CF at position ITER
!-----------------------------------------------------------------------
      DO M=1,NPL
        CF(M,ITER,NP)=W1(NP)%CPTWFP(M)
      ENDDO

      DO NPRO=1,WDES%NPRO
        CPROF(NPRO,ITER,NP)= W1(NP)%CPROJ(NPRO)
      ENDDO
!-----------------------------------------------------------------------
! calculate the search-directions  (H-epsilon S) |phi_opt>
!-----------------------------------------------------------------------
      EVALUE(NP)=W%CELEN(N,NK,ISP)
      ENDDO i1

!      DO NP=1,NSIM
!      IF (.NOT. LDO(NP)) CYCLE
!      CALL HAMILT(WDES1,W1(NP),NONLR_S,NONL_S,GRID,  INFO%LREAL,EVALUE_INI(NP), &
!     &    LMDIM,CDIJ(1,1,1,ISP),CQIJ(1,1,1,ISP), SV(1,ISP), CF(1,2*NRES,NP))
!      ENDDO
      !  store H | psi > temporarily in upmost storage position (2*NRES)
      !  to have uniform stride for result array
      CALL HAMILTMU(WDES1,W1,NONLR_S,NONL_S,GRID,  INFO%LREAL,EVALUE_INI, &
     &     LMDIM,CDIJ(1,1,1,ISP),CQIJ(1,1,1,ISP), SV(1,ISP), CF(1,2*NRES,1),LD, NSIM, LDO)

      i2: DO NP=1,NSIM


      N=NB(NP); ITER=IT(NP); IF (.NOT. LDO(NP)) CYCLE i2

      FNORM=0
      FPRE_ =0
      DO ISPINOR=0,WDES%NRSPINORS-1
      DO M=1,NGVECTOR
        MM=M+ISPINOR*NGVECTOR
        CF(MM,NRES+ITER,NP)=CF(MM,2*NRES,NP)-EVALUE_INI(NP)*W1(NP)%CPTWFP(MM)
!-MM- changes to accommodate spin spirals
! original statement
!       IF (LDELAY .AND. WDES%DATAKE(M,NK)>INFO%ENINI) CF(MM,NRES+ITER,NP)=0
        IF (LDELAY .AND. WDES%DATAKE(M,NK,ISPINOR+1)>INFO%ENINI) CF(MM,NRES+ITER,NP)=0
!-MM- end of alterations
        FNORM =FNORM+CF(MM,NRES+ITER,NP)*CONJG(CF(MM,NRES+ITER,NP))
        FPRE_ =FPRE_+CF(MM,NRES+ITER,NP)*CONJG(CF(MM,NRES+ITER,NP)) &
     &              *PRECON(MM,NP)
      ENDDO
      ENDDO
      
!     norm of total error vector before start
!     norm smaller than EBREAK stop |e -e(app)| < | Residuum |
      IF (ABS(FNORM)<INFO%EBREAK) THEN
         LDO(NP)=.FALSE.
         CYCLE i2
      ENDIF

      IF (INFO%IALGO==6)  FPRE_=FNORM
      FPRE(NP)=FPRE_
      IF (IDUMP==2) WRITE(*,'(E9.2,"R")',ADVANCE='NO') SQRT(ABS(FNORM))
      IF (ITER==1) THEN
        RMS=RMS+WDES%RSPIN*WDES%WTKPT(NK)*W%FERWE(N,NK,ISP)* &
     &      SQRT(ABS(FNORM))/WDES%NB_TOT
      ENDIF

!----------------------------------------------------------------------

      optsubspace: IF (.NOT. LDELAY .AND. ITER > 1) THEN

! better conditioning for search
!DIR$ IVDEP
!OCL NOVREL
      DO M=1,NPL
         CF(M,ITER-1,NP)=CF(M,ITER,NP)-CF(M,ITER-1,NP)
         CF(M,NRES+ITER-1,NP)=CF(M,NRES+ITER,NP)-CF(M,NRES+ITER-1,NP)
      ENDDO

!DIR$ IVDEP
!OCL NOVREL
      DO NPRO=1,WDES%NPRO
         CPROF(NPRO,ITER-1,NP)=CPROF(NPRO,ITER,NP)-CPROF(NPRO,ITER-1,NP)
      ENDDO
!***********************************************************************
! RMM-DIIS step
!
! minimize norm of residual vector in the subspace spanned by
! the set of wavefunctions stored in CF (with the corresponding
! projections stored in CPROF)
!
!***********************************************************************
      optsub: IF (INFO%IALGO /=0 ) THEN
!----------------------------------------------------------------------
! calculate  matrix INFO%IALGO=7 or 8
! CHAM(n2,n1)= /R(n2)/ R(n1)/
! CTMP(n2,n1)= /a(n2)/ S /a(n1)/
!----------------------------------------------------------------------
    CTMP=0
    CHAM=0

    buildh: DO N1=1,ITER
      IF (INFO%IALGO==8 .OR. INFO%IALGO==9) THEN
      DO M=1,NPL
        WTMP(NP)%CPTWFP(M)=CF(M,NRES+N1,NP)*PRECON(M,NP)
      ENDDO
      ELSE
      DO M=1,NPL
        WTMP(NP)%CPTWFP(M)=CF(M,NRES+N1,NP)
      ENDDO
      ENDIF

      CALL ZGEMV( 'C', NPL, ITER+1-N1,(1._q,0._q) ,CF(1,NRES+N1,NP), &
     &            WDES%NRPLWV, WTMP(NP)%CPTWFP(1) , 1 , (0._q,0._q),  CWORK1(1), 1)

      DO N2=N1,ITER
        CHAM(N2,N1)=      (CWORK1(N2-N1+1))
        CHAM(N1,N2)=(CONJG(CWORK1(N2-N1+1)))
      ENDDO

      CALL ZGEMV( 'C', NPL, ITER+1-N1,(1._q,0._q) ,CF(1,N1,NP), &
     &            WDES%NRPLWV, CF(1,N1,NP) , 1 , (0._q,0._q),  CWORK1(1), 1)

      IF (INFO%LOVERL .AND. WDES%NPROD>0 ) THEN
        WDES1%NBANDS=1    ! is used only here not quite clean
        CALL OVERL(WDES1, INFO%LOVERL,LMDIM,CQIJ, CPROF(1,N1,NP),WTMP(NP)%CPROJ(1))
        CALL ZGEMV( 'C', WDES%NPRO, ITER+1-N1, (1._q,0._q) ,CPROF(1,N1,NP), &
     &           WDES%NPROD, WTMP(NP)%CPROJ(1) , 1 ,(1._q,0._q) ,  CWORK1(1), 1)

      ENDIF

      DO N2=N1,ITER
         CTMP(N2,N1)=      (CWORK1(N2-N1+1))
         CTMP(N1,N2)=(CONJG(CWORK1(N2-N1+1)))
      ENDDO
    ENDDO buildh

    
    

      DO N1=1,ITER
         IF (ABS(AIMAG(CHAM(N1,N1)))>1E-2_q) THEN
            WRITE(*,*)'WARNING: Sub-Space-Matrix is not hermitian in rmm', &
                       AIMAG(CHAM(N1,N1))
         ENDIF
         CHAM(N1,N1)= REAL( CHAM(N1,N1) ,KIND=q)
      ENDDO

! solve eigenvalue-problem and calculate lowest eigenvector
! this eigenvector corresponds to a minimal residuum
! CHAM(n1,n2) U(n2,1) = E(1) S(n1,n2)  U(n2,1)

      IF (.FALSE.) THEN
         
         NPL2=MIN(10,ITER)
         WRITE(6,*)
         DO N1=1,NPL2
            WRITE(6,1)N1,(REAL( CHAM(N1,N2) ,KIND=q) ,N2=1,NPL2)
         ENDDO
         WRITE(6,*)
         DO N1=1,NPL2
            WRITE(6,3)N1,(AIMAG(CHAM(N1,N2)),N2=1,NPL2)
         ENDDO
         WRITE(6,*)
         DO N1=1,NPL2
            WRITE(6,1)N1,(REAL( CTMP(N1,N2) ,KIND=q) ,N2=1,NPL2)
         ENDDO
         WRITE(6,*)
         DO N1=1,NPL2
            WRITE(6,3)N1,(AIMAG(CTMP(N1,N2)),N2=1,NPL2)
         ENDDO
         WRITE(6,*)

 1       FORMAT(1I2,3X,20F9.5)
 3       FORMAT(1I2,3X,20E9.1)
         
      ENDIF
!
! onion award of the year for IBM,
! who use a completely different DSYGV calling sequence
!
      CALL ZHEGV &
     &  (1,'V','U',ITER,CHAM,NRES,CTMP,NRES,R, &
     &           CWORK(1),LWORK*INFO%NDAV,RWORK(1),IFAIL)
      

      IF (IFAIL/=0) THEN
         IF (IU6>=0) &
         WRITE(IU6,219) IFAIL,ITER,N
         IF (IU0>=0) &
         WRITE(IU0,219) IFAIL,ITER,N
!  try to save things somehow, goto next band
         LDO(NP)=.FALSE.
         CYCLE i2
      ENDIF
  219 FORMAT('WARNING in EDDRMM: call to ZHEGV failed, returncode =',I4,I2,I2)
      FPRE(NP)=R(1)
      IF (IDUMP==2)  WRITE(*,'(E9.2,"R")',ADVANCE='NO') SQRT(ABS(FPRE(NP)))

!     write out 'optimal trial step' i.e step which would have minimized
!     the residuum
      IF (ITER==2 .AND. IDUMP==2) THEN
         OTRIAL= REAL( 1+CHAM(1,1)/CHAM(2,1) ,KIND=q) *TRIAL(NP)
         WRITE(*,'(1X,F7.4,"o")',ADVANCE='NO') OTRIAL
      ENDIF

!     some heuristic for numerical accuracy problems
!     small residuum and negative step -> stop immedately
      IF (ITER==2) THEN
         OTRIAL= REAL( 1+CHAM(1,1)/CHAM(2,1) ,KIND=q) *TRIAL(NP)
         IF (OTRIAL <0 .AND.  ABS(FPRE(NP))< 1E-9_q) THEN
            IF (IU0>=0) WRITE(IU0,'(" num prob ")',ADVANCE='NO')
            LDO(NP)=.FALSE.
            CYCLE i2
         ENDIF
      ENDIF

      IF (.FALSE.) THEN
      
         NPL2=MIN(10,ITER)

         WRITE(77,*)
         DO N1=1,NPL2
            WRITE(77,1)N1,R(N1),(REAL( CHAM(N2,N1) ,KIND=q) ,N2=1,NPL2)
         ENDDO
         WRITE(77,*)
         DO N1=1,NPL2
            WRITE(77,3)N1,R(N1),(AIMAG(CHAM(N2,N1)),N2=1,NPL2)
         ENDDO
         WRITE(77,*)
      
      ENDIF

      ELSE optsub
!***********************************************************************
! inverse interation step
! minimize
!    || ( H - e S)  | phi > - | phi_initial > ||
!                                                m
! in the subspace spanned by the set of wavefunctions stored in
! CF
!
!***********************************************************************
    AM=0
    B =0

!        < phi_j |  ( H - e S) ( H - e S)  | phi_i >
    builda: DO N1=1,ITER
      DO M=1,NPL
         WTMP(NP)%CPTWFP(M)=CF(M,NRES+N1,NP)
      ENDDO
      CALL ZGEMV( 'C', NPL, ITER+1-N1, (1._q,0._q), CF(1,NRES+N1,NP), &
     &            WDES%NRPLWV, WTMP(NP)%CPTWFP(1), 1, (0._q,0._q),  CWORK1(1), 1)

      DO N2=N1,ITER
        AM(N2,N1)=      (CWORK1(N2-N1+1))
        AM(N1,N2)=(CONJG(CWORK1(N2-N1+1)))
      ENDDO
     ENDDO builda

      DO M=1,NPL
         WTMP(NP)%CPTWFP(M)=CF_INI(M,NP)
      ENDDO
!        < phi_j |  ( H - e S) | phi_ini >
      CALL ZGEMV( 'C', NPL, ITER, (1._q,0._q), CF(1,NRES+1,NP), &
     &            WDES%NRPLWV, WTMP(NP)%CPTWFP(1), 1, (0._q,0._q),  CWORK1(1), 1)

      DO N1=1,ITER
         B(N1)=      (CWORK1(N1))
      ENDDO


    
    

      AM_=AM
      B_=B

      IF (.FALSE.) THEN
         
         NPL2=MIN(10,ITER)
         WRITE(6,*)
         DO N1=1,NPL2
            WRITE(6,'(I3,8E14.7)')N1, (AM(N1,N2) ,N2=1,NPL2)
         ENDDO

         WRITE(6,'(A3,8E14.7)') 'v', (B(N1) ,N1=1,NPL2)

         
      ENDIF
!      CALL ZGETRF( ITER, ITER, AM, NRES, IPIV, IFAIL )
      IF (IFAIL ==0) &
!      CALL ZGETRS('T', ITER, 1, AM, NRES, IPIV, B, NRES, IFAIL)
      T=0
      DO K=1,ITER
      DO KP=1,ITER
         T=T+B(K)*AM_(K,KP)*B(KP)
      ENDDO
      T=T-2*B_(K)*B(K)
      ENDDO
      IF (IDUMP>=2) WRITE(*,'(1X,F7.4,"a")',ADVANCE='NO') T


      IF (ITER == 1) B(1) = 1

      IF (IFAIL/=0) THEN
         IF (IU6>=0) &
         WRITE(IU6,219) IFAIL,ITER,N
         IF (IU0>=0) &
         WRITE(IU0,219) IFAIL,ITER,N
!  try to save things somehow, goto next band
         LDO(NP)=.FALSE.
         CYCLE i2
      ENDIF
!     write out 'optimal trial step' i.e step which would have minimized
!     the residuum
      IF (ITER==2 .AND. IDUMP==2) THEN
         OTRIAL= REAL( 1+B(1)/B(2) ,KIND=q) *TRIAL(NP)
         WRITE(*,'(1X,F7.4,"o")',ADVANCE='NO') OTRIAL
      ENDIF

      IF (.FALSE.) THEN
      
         NPL2=MIN(10,ITER)
         WRITE(6,*)
         WRITE(6,'(A3,8E14.7)') 'e',(B(N1),N1=1,NPL2)
      
      ENDIF

      DO I=1,ITER
        CHAM(I,1)=B(I)
      ENDDO

      IF (IDUMP == 2) THEN

         DO M=1,NPL
            WTMP(NP)%CPTWFP(M)=-CF_INI(M,NP)
            IF (ITER > 1) THEN
               DO I=1,ITER
                  WTMP(NP)%CPTWFP(M)=WTMP(NP)%CPTWFP(M)+CF(M,I+NRES,NP)*B(I)
               ENDDO
            ENDIF
            WTMP(NP)%CPTWFP(M)=WTMP(NP)%CPTWFP(M)
         ENDDO

         CALL ZGEMV( 'C', NPL, 1, (1._q,0._q),  WTMP(NP)%CPTWFP(1), &
               WDES%NRPLWV, WTMP(NP)%CPTWFP(1), 1, (0._q,0._q),  CWORK1(1), 1)
         IF (IDUMP>=2) WRITE(*,'(1X,F7.4,"A")',ADVANCE='NO') REAL(CWORK1(1))
      ENDIF

      ENDIF optsub
!=======================================================================
! now performe the trial step (ITER > 1 use previous trial step)
! but restrict trial step to 1.0
!=======================================================================
! transform wavefunction
      IF (TRIAL(NP)<0) TRIAL(NP)=ABS(TRIAL(NP))
         DO M=1,NPL
            WTMP(NP)%CPTWFP(M)=0
         ENDDO

         DO I=1,ITER
!DIR$ IVDEP
!OCL NOVREC
            DO M=1,NPL
               WTMP(NP)%CPTWFP(M)=WTMP(NP)%CPTWFP(M)+CF(M,I,NP)*CHAM(I,1)
            ENDDO
         ENDDO

! trial step on wavefunction
!DIR$ IVDEP
!OCL NOVREC
         DO I=1,ITER
            DO M=1,NPL
               WTMP(NP)%CPTWFP(M)=WTMP(NP)%CPTWFP(M)-CF(M,NRES+I,NP)*CHAM(I,1)*PRECON(M,NP)*TRIAL(NP)
            ENDDO
         ENDDO
! transform the wave-function to real space

         DO ISPINOR=0,WDES%NRSPINORS-1
            CALL FFTWAV(NGVECTOR,WDES%NINDPW(1,NK),WTMP(NP)%CR(1+ISPINOR*WDES1%MPLWV),WTMP(NP)%CPTWFP(1+ISPINOR*NGVECTOR),GRID)
         ENDDO
      ELSE optsubspace
!***********************************************************************
!
! minimize energy starting from current wavefunctions stored in
! W1 along the current searchdirection stored in CF(..,NRES+ITER)
!
!***********************************************************************
! trial vector in line minimization

      DO M=1,NPL
          WTMP(NP)%CPTWFP(M)=CF(M,NRES+ITER,NP)*PRECON(M,NP)
      ENDDO

      DO ISPINOR=0,WDES1%NRSPINORS-1
         CALL FFTWAV(NGVECTOR,WDES%NINDPW(1,NK),WTMP(NP)%CR(1+ISPINOR*WDES1%MPLWV),WTMP(NP)%CPTWFP(1+ISPINOR*NGVECTOR),GRID)
      ENDDO
      ENDIF optsubspace
      ENDDO i2

!-----------------------------------------------------------------------
! calculate results of projection operatores
! and new eigenenergy
!-----------------------------------------------------------------------
      IF ( INFO%LREAL ) THEN
         IF (NSIM >1 ) THEN
         CALL RPROMU(NONLR_S, WDES1, WTMP, NSIM, LDO)
         ELSE
         DO NP=1,NSIM
            IF (.NOT. LDO(NP)) CYCLE
            CALL RPRO1(NONLR_S, WDES1, WTMP(NP))
         ENDDO
         ENDIF
      ELSE
         DO NP=1,NSIM
            IF (.NOT. LDO(NP)) CYCLE
            CALL PROJ1(NONL_S,WDES1,WTMP(NP))
         ENDDO
      ENDIF


      i3: DO NP=1,NSIM
      N=NB(NP); ITER=IT(NP); IF (.NOT. LDO(NP)) CYCLE i3
!***********************************************************************
!
! finish trial step of RMM-DIIS/ inverse iteration by copying WTMP to W1
!
!***********************************************************************
      mine: IF (.NOT. LDELAY .AND. ITER > 1) THEN
!DIR$ IVDEP
!OCL NOVREC
      DO M=1,NPL
         W1(NP)%CPTWFP(M)=WTMP(NP)%CPTWFP(M)
      ENDDO
!DIR$ IVDEP
!OCL NOVREC
      DO NPRO=1,WDES%NPRO
        W1(NP)%CPROJ(NPRO)=WTMP(NP)%CPROJ(NPRO)
      ENDDO
!DIR$ IVDEP
!OCL NOVREC
      DO ISPINOR=0,WDES%NRSPINORS-1
      DO K=1,GRID%RL%NP
         KK=K+ISPINOR*WDES1%MPLWV
         W1(NP)%CR(KK)=WTMP(NP)%CR(KK)
      ENDDO
      ENDDO
!***********************************************************************
!
! finish line minimization of <phi| H | phi> along trial direction
!
!***********************************************************************
      ELSE mine
    ! < g | phi > 1. order energy change
      A2=0
      DO M=1,NPL
         A2 =A2+WTMP(NP)%CPTWFP(M)*CONJG(CF(M,NRES+ITER,NP))
      ENDDO
      
      A2=-A2*2
    ! get exact energy change along trial step (A1)

!DIR$ IVDEP
!OCL NOVREC
      DO M=1,NPL
         W1(NP)%CPTWFP(M)=W1(NP)%CPTWFP(M)-WTMP(NP)%CPTWFP(M)
      ENDDO
!DIR$ IVDEP
!OCL NOVREC
      DO NPRO=1,WDES%NPRO
        W1(NP)%CPROJ(NPRO)= W1(NP)%CPROJ(NPRO)-WTMP(NP)%CPROJ(NPRO)
      ENDDO
!DIR$ IVDEP
!OCL NOVREC
      DO ISPINOR=0,WDES%NRSPINORS-1
      DO K=1,GRID%RL%NP
         KK=K+ISPINOR*WDES1%MPLWV
         W1(NP)%CR(KK)=W1(NP)%CR(KK)-WTMP(NP)%CR(KK)
      ENDDO
      ENDDO

      CALL ECCP(WDES1,W1(NP),W1(NP),LMDIM,CDIJ(1,1,1,ISP),GRID,SV(1,ISP), W%CELEN(N,NK,ISP))
      CALL CNORMA(WDES1,W1(NP), INFO%LOVERL,LMDIM,CQIJ,WSCAL)

      W%CELEN(N,NK,ISP) =W%CELEN(N,NK,ISP)*WSCAL**2
      A1= W%CELEN(N,NK,ISP)-EVALUE(NP)
   ! quadratic interpolation to find the minimum

      TRIAL(NP)= -A2/(A1-A2)/2
      DE       = (A2+(A1-A2)*TRIAL(NP))*TRIAL(NP)
   ! avoid too large trial steps
      IF (IDUMP>=2) WRITE(*,'(1X,F7.4,"T")',ADVANCE='NO') TRIAL(NP)
      IF (.NOT. LDELAY) THEN
         IF (TRIAL(NP)>0 .AND.(TRIAL(NP) > 1)) TRIAL(NP)= 1
         IF (TRIAL(NP)>0 .AND.(TRIAL(NP)<0.1)) TRIAL(NP)=0.1
         IF (TRIAL(NP)<0 .AND.(TRIAL(NP) <-1)) TRIAL(NP)=-1
         IF (TRIAL(NP)<0 .AND.(TRIAL(NP)>-.1)) TRIAL(NP)=-0.1
      ENDIF

      IF (IDUMP>=2) WRITE(*,'(1X,F7.4,"T")',ADVANCE='NO') TRIAL(NP)

   ! set W1 finally

!DIR$ IVDEP
!OCL NOVREC
      DO M=1,NPL
         W1(NP)%CPTWFP(M)=W1(NP)%CPTWFP(M)-WTMP(NP)%CPTWFP(M)*(TRIAL(NP)-1)
      ENDDO
!DIR$ IVDEP
!OCL NOVREC
      DO NPRO=1,WDES%NPRO
         W1(NP)%CPROJ(NPRO)= W1(NP)%CPROJ(NPRO)-WTMP(NP)%CPROJ(NPRO)*(TRIAL(NP)-1)
      ENDDO
!DIR$ IVDEP
!OCL NOVREC
      DO ISPINOR=0,WDES%NRSPINORS-1
      DO K=1,GRID%RL%NP
         KK=K+ISPINOR*WDES1%MPLWV
         W1(NP)%CR(KK)=W1(NP)%CR(KK)-WTMP(NP)%CR(KK)*(TRIAL(NP)-1)
      ENDDO
      ENDDO

      ENDIF mine
!=======================================================================
! common code
!=======================================================================

      CALL CNORMN(WDES1,W1(NP), INFO%LOVERL,LMDIM,CQIJ,WSCAL)
    ! rescale also the stored real space wavefuntion

      DO ISPINOR=0,WDES%NRSPINORS-1
      DO K=1,GRID%RL%NP
         KK=K+ISPINOR*WDES1%MPLWV
         W1(NP)%CR(KK)=W1(NP)%CR(KK)*WSCAL
      ENDDO
      ENDDO

      CALL ECCP(WDES1,W1(NP),W1(NP),LMDIM,CDIJ(1,1,1,ISP),GRID,SV(1,ISP), W%CELEN(N,NK,ISP))
      DECEL =W%CELEN(N,NK,ISP)-EVALUE(NP)
      DE    =DECEL

      IF (IDUMP==2) WRITE(*,'(E10.2,2H |)',ADVANCE='NO') DECEL

      DESUM =DESUM +WDES%RSPIN*WDES%WTKPT(NK)*W%FERWE(N,NK,ISP)*DECEL
      ICOUEV=ICOUEV+1
!=======================================================================
! break of intra-band-minimisation
! at the moment we performe a break of the intra-band minimization if
! ) DE is less then INFO%DEPER % of the change in the first minimization
!     of this band (relative breakcondition)
! ) DE less then INFO%EBREAK (absolut breakcondition)
! ) if unoccupied band break after 2. iteration
!=======================================================================
      DE=ABS(DE)
      IF (DE<INFO%EBREAK) LDO(NP)=.FALSE.
      IF (ABS(W%FERWE(N,NK,ISP))<INFO%WEIMIN .AND. &
     &    ABS(W%FERWE(N,NK,ISP)*DE)<INFO%WEIMIN .AND. ITER >= 2) LDO(NP)=.FALSE.
      IF (ABS(FPRE(NP))<DEIT(NP)) LDO(NP)=.FALSE.
      IF (ITER==1) THEN
        DEIT(NP)=ABS(FPRE(NP))*INFO%DEPER
      ENDIF
      IF (ITER == NITER) LDO(NP)=.FALSE.

      ENDDO i3

! (1._q,0._q) band just finished ?, set NB(NP) also to 0 and finish everything
      DO NP=1,NSIM
         N=NB(NP)
         IF (.NOT. LDO(NP) .AND. N /=0 ) THEN
            NB(NP)=0
            IF (IDUMP==2)  WRITE(*,'(F9.4,2H q)')REAL( W%CELEN(N,NK,ISP) ,KIND=q)
            IF (IDUMP==10) WRITE(*,*)
         ENDIF
      ENDDO
!=======================================================================
! move onto the next Band
!=======================================================================
      ENDDO bands
!=======================================================================
      END DO kpoints
      ENDDO spin
!=======================================================================

      
      
      

      IF (W%OVER_BAND) THEN
         W%OVER_BAND=.FALSE.
         CALL REDIS_PW_DEALLOC(H_PW)
      ENDIF

      DO NP=1,NSIM
         DEALLOCATE(W1(NP)%CR)
         CALL DELWAV(WTMP(NP) ,.TRUE.)
      ENDDO
      DEALLOCATE(CF,CF_INI,PRECON,CPROF,CHAM,CTMP,CWORK1, &
                 AM,B,IPIV,AM_,B_)

      RETURN
      END SUBROUTINE
      END MODULE
