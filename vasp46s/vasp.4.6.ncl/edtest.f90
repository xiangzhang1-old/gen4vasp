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





!************************ SUBROUTINE EDSTEP*****************************
! RCS:  $Id: edtest.F,v 1.3 2002/08/14 13:59:38 kresse Exp $
!
! simple steepest descent relaxation of the wavefunctions
! allows duing several steepest descent steps on each
! wavefunctions to mimic Alavis algorithm
! it is curcial not to overrelax if several steepest decent steps
! are done
!
!***********************************************************************


      SUBROUTINE EDEXP(GRID,INFO,LATT_CUR,NONLR_S,NONL_S,W,WDES, &
        LMDIM,CDIJ,CQIJ, RMS,DESUM,ICOUEV, SV)
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
      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRID
      TYPE (info_struct) INFO
      TYPE (latt)        LATT_CUR
      TYPE (nonlr_struct) NONLR_S
      TYPE (nonl_struct) NONL_S
      TYPE (wavefun)     W
      TYPE (wavedes)     WDES

      COMPLEX(q)   SV(GRID%RL%NP)          ! local potential
      COMPLEX(q) CDIJ(LMDIM,LMDIM,WDES%NIONS),CQIJ(LMDIM,LMDIM,WDES%NIONS)

!----- local work arrays
      TYPE (wavedes1)    WDES1         ! descriptor for (1._q,0._q) k-point
      TYPE (wavefun1)    W1            ! current wavefunction
      TYPE (wavefun1)    WSEARCH       ! current search direction
      REAL(q),ALLOCATABLE::    PRECON(:)

      ALLOCATE(W1%CR(GRID%MPLWV),PRECON(WDES%NRPLWV))
      CALL NEWWAV(WSEARCH ,WDES,GRID%MPLWV,.FALSE.)

      DELTA=INFO%TIME

      RMS=0
      ICOUEV=0
      DESUM=0
!=======================================================================
      kpoints: DO NK=1,WDES%NKPTS
!=======================================================================
      CALL SETWDES(WDES,WDES1,NK); CALL SETWGRID(WDES1,GRID)
      IF (INFO%LREAL) THEN
        CALL PHASER(GRID,LATT_CUR,NONLR_S,NK,WDES,0.0_q,0.0_q,0.0_q)
      ELSE
        CALL PHASE(WDES,NONL_S,NK)
      ENDIF
!=======================================================================
      bands: DO N=WDES%NBANDS,1,-1
!=======================================================================

      NPL=WDES1%NPL
      IDUMP=0
      IF (MOD(N,5)==0) IDUMP=0
      IF (IDUMP==2) WRITE(*,'(I3,1X,$)') N

      CALL SETWAV(W,W1,N,NK)  ! allocation for W1%CR done above
!-----------------------------------------------------------------------
! start with the exact evaluation of the eigenenergy
!-----------------------------------------------------------------------
      IF (IDUMP==2) WRITE(*,'(F9.4,$)') REAL( W%CELEN(N,NK) ,KIND=q)

      CALL FFTWAV(NPL,WDES%NINDPW(1,NK),W1%CR(1),W1%CPTWFP(1),GRID)
      CALL ECCP(WDES1,W1,W1,LMDIM,CDIJ,GRID,SV, W%CELEN(N,NK))

      IF (IDUMP==2) WRITE(*,'(F9.4,$)') REAL( W%CELEN(N,NK) ,KIND=q)
!-----------------------------------------------------------------------
! calculate the preconditioning matrix
!-----------------------------------------------------------------------
      IF (INFO%IALGO==7 .OR.INFO%IALGO==8) THEN
        EKIN=0
        DO M=1,WDES1%NPL
          CPT=W1%CPTWFP(M)
!-MM- changes to accommodate spin spirals
! original statement
!         EKIN =EKIN+ REAL( CPT*CONJG(CPT) ,KIND=q) * WDES1%DATAKE(M)
          EKIN =EKIN+ REAL( CPT*CONJG(CPT) ,KIND=q) * WDES1%DATAKE(M,1)
!-MM- end of alterations
        ENDDO
        IF (EKIN<2.0_q) EKIN=2.0_q
        EKIN=EKIN*1.5_q
        IF (IDUMP==2) WRITE(*,'(F9.4,$)') EKIN

        FAKT=2._q/EKIN
        DO M=1,WDES1%NPL
!-MM- changes to accommodate spin spirals
! original statement
!         X=WDES1%DATAKE(M)/EKIN
          X=WDES1%DATAKE(M,1)/EKIN
!-MM- end of alterations
          X2= 27+X*(18+X*(12+8*X))
          PRECON(M)=X2/(X2+16*X*X*X*X)*FAKT
        ENDDO
      ELSE
        PRECON=1
      ENDIF

      NITER=4
      EVALUE=W%CELEN(N,NK)
!=======================================================================
      iter_band: DO ITER=1,NITER
!=======================================================================
      ICOUEV=ICOUEV+1
!
!  calculate the search-directions  H-epsilon S |phi>
!
      CALL HAMILT(WDES1,W1,NONLR_S,NONL_S,GRID,  INFO%LREAL,EVALUE, &
     &    LMDIM,CDIJ,CQIJ, SV,WSEARCH%CPTWFP(1))

      DETEST=0
      CFNORM=0
      DO M=1,WDES1%NPL
        CPT     =W1%CPTWFP(M)
        WSEARCH%CPTWFP(M)   =WSEARCH%CPTWFP(M)-EVALUE*CPT
        CFNORM =CFNORM+WSEARCH%CPTWFP(M)*CONJG(WSEARCH%CPTWFP(M))
      ENDDO
!     norm of total error vector before start
      IF (ITER==1) THEN
      IF (IDUMP==2) WRITE(*,'(F9.4,$)') SQRT(ABS(CFNORM))
        RMS=RMS+2*WDES%WTKPT(NK)*W%FERWE(N,NK)* &
     &       SQRT(ABS(CFNORM))/WDES%NBANDS
      ENDIF
!
! add preconditioned search direction to wavefunction
!
       DO M=1,NPL
          W1%CPTWFP(M)=W1%CPTWFP(M)-PRECON(M)*WSEARCH%CPTWFP(M)*DELTA
       ENDDO

       CALL FFTWAV(NPL,WDES%NINDPW(1,NK),W1%CR(1),W1%CPTWFP(1),GRID)
       IF (INFO%LREAL) THEN
         CALL RPRO1(NONLR_S,WDES1,W1)
       ELSE
         CALL PROJ1(NONL_S,WDES1,W1)
       ENDIF

       CALL CNORMN(WDES1,W1, INFO%LOVERL,LMDIM,CQIJ,WSCAL)
       DO K=1,WDES1%NPLWVL
         W1%CR(K)=W1%CR(K)*WSCAL
       ENDDO

!=======================================================================
      ENDDO iter_band
!=======================================================================

      CALL ECCP(WDES1,W1,W1,LMDIM,CDIJ,GRID,SV, W%CELEN(N,NK))

      DECEL=W%CELEN(N,NK)-EVALUE
      DESUM =DESUM+ 2*WDES%WTKPT(NK)*W%FERWE(N,NK)*DECEL
      IF (IDUMP==2) WRITE(*,'(F9.4)') REAL( W%CELEN(N,NK) ,KIND=q)

!=======================================================================
      ENDDO bands
      ENDDO kpoints
!=======================================================================

      DEALLOCATE(W1%CR,PRECON)
      CALL DELWAV(WSEARCH ,.FALSE.)

      RETURN
      END SUBROUTINE
