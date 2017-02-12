#include "symbol.inc"
      MODULE rot
      USE prec
      CONTAINS
!************************ SUBROUTINE EDWAV *****************************
! RCS:  $Id: rot.F,v 1.2 2001/02/20 14:45:00 kresse Exp $
!
! subroutine for optimizing all bands simultaneous
! the subroutine caclulates a search direction into which the
! total energy is minimized, than a trial step into this direction
! is done.
! On the next call the new exact energy must be supplied by the calling
! routine, EDWAV then tries to dertermine an optimum steplenght
! and steps into the minimum.
! IFLAG determines the status of the routine:
!  0   trial - steepest descent step
!  1   trial - conjugate gradient step(IALGO can overwrite this mode to 0)
!  2   go to the  to minimum
!      TOTEN must contain the old energy TOTEN2 the energy after the step
!  10  increase trial step (mainly for checking accuracy of routine)
! on exit
! >0   main routine performes anouther energy-evaluation
!  0   main routine ends the loop
!
! IALGO determines the specific algorithm
!  5   steepest descent without optimizing the stepsize
!  7   steepest descent with optimization of stepsize
!  8   conjugated gradient
!  9   conjugate gradient with gradient information only
!
! IROT determines if suspace rotation and/or wavefunction are optimized
!  1   wavefunctions
!  2   subspacerotation
!  3   both
! ICEL determines how partial occupancies are optimized
!  0   no optimization
!  1   standart optimization (optimize epsilon)
!  2   optimized partial occupancies directly
! DESUM1 expected energy change energy change
!
!
!***********************************************************************

      SUBROUTINE EDWAV(IALGO,IFLAG,IROT,ICEL,IU6,IU0,BTRIAL, &
             EFERMI,LOVERL,LREAL,NBLK,ORT,RMS,TOTEN,TOTEN2,DESUM1, &
             GRID,KPOINTS,LATT_CUR,NONLR_S,NONL_S,W,WDES,W_F,W_G, &
             LMDIM,CQIJ,CDIJ,CHAM,CHF,SV )
      USE prec

      USE wave
      USE dfast
      USE lattice
      USE mpimy
      USE mgrid
      USE nonl
      USE nonlr
      USE hamil
      USE constant
      USE mkpoints
      USE choleski

      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRID
      TYPE (latt)        LATT_CUR
      TYPE (nonlr_struct) NONLR_S
      TYPE (nonl_struct) NONL_S
      TYPE (wavefun)     W     ! current wavefunctions
      TYPE (wavefun)     W_F   ! last search direction
      TYPE (wavefun)     W_G   ! current gradient
          ! W_F%CELEN is special, and stores variational eigenvalues
          ! from which partial occupancies are set
      TYPE (wavefun1)    W1,WPRE
      TYPE (wavedes)     WDES
      TYPE (wavedes1)    WDES1
      TYPE (kpoints_struct) KPOINTS
      LOGICAL LREAL,LOVERL


      RGRID   SV(GRID%RL%NP)          ! local potential
!-----Arrays for Non-local Contribution
      OVERLAP CDIJ(LMDIM,LMDIM,WDES%NIONS),CQIJ(LMDIM,LMDIM,WDES%NIONS)
!-----Arrays for gradient
      GDEF CHAM(WDES%NBANDS,WDES%NBANDS,WDES%NKPTS),CHF(WDES%NBANDS,WDES%NBANDS,WDES%NKPTS)

      PARAMETER (NSTEP=10)
      REAL(q)    STEP(NSTEP),ENERGY(NSTEP)
! work array
      REAL(q)    HFEIG(WDES%NBANDS)
      COMPLEX(q) CGRAD(WDES%NRPLWV)
      GDEF COVL(WDES%NBANDS,WDES%NBANDS),CUNI(WDES%NBANDS,WDES%NBANDS), &
     &     CEIDB(WDES%NBANDS,WDES%NBANDS)

      SAVE GNORM,GNORML,BSTEP,DESUM0,ICOUNT,GAMMA,NCOUNT,GSL
      SAVE STRIAL1,STRIAL2

      DATA GNORML /0/
      DATA BSTEP  /0/
      DATA ILOEW  /1/
      DATA ICOUNT /-1/
      DATA AVERAG /0/

#ifndef noPAW
! we cant use this code for spin polarized mode yet
      IF (WDES%ISPIN==2) THEN
         WRITE(*,*) 'EDWAV: ISPIN=2 NOT YET ALLOWED (NOT IMPLEMETED)'
         STOP
      ENDIF

!    initialized average step size to trial step size
      IF (AVERAG==0) AVERAG=BTRIAL
      BTRIAL=AVERAG
      WEIMIN=0.000_q
      IDUMP=1

      CALL NEWWAV(WPRE ,WDES,GRID%MPLWV,.TRUE.)
      W1%CR => WPRE%CR  ! need storage only once

!*********************************************************************
!   TRIAL STEP (or PREDICTOR STEP)
! ) calculate the gradient vector and store
!   preconditioned gradient in W_G
! ) calculate change in energy if step into direction W_G is done DESUM
! ) calculate norm of gradient-vector GNORM
! ) test orthonormality with old direction ORT
! ) do trial step
!***********************************************************************
      main: IF (IFLAG< 2) THEN
      NCOUNT=0

      GNORM =0
      ORT   =0
      ORT1  =0
      ORT2  =0
      ORTCEL=0
      DE1   =0
      DE2   =0
      DECEL =0
!=======================================================================
! first loop
! calculate gradient (for wavefunctions)
!=======================================================================
      kpoint: DO NK=1,WDES%NKPTS
!=======================================================================
      NPL=WDES%NPLWKP(NK)

      CALL SETWDES(WDES,WDES1,NK); CALL SETWGRID(WDES1,GRID)
!-----------------------------------------------------------------------
! calculate the phase-factor on the compressed grid (only non local)
!-----------------------------------------------------------------------
      IF (LREAL) THEN
        CALL PHASER(GRID,LATT_CUR,NONLR_S,NK,WDES,0.0_q,0.0_q,0.0_q)
      ELSE
        CALL PHASE(WDES,NONL_S,NK)
      ENDIF

!=======================================================================
      bands: DO N=1,WDES%NBANDS
!=======================================================================
      CALL SETWAV(W,W1,N,NK)  ! allocation for W1%CR done above

      WEIGHT =2*WDES%WTKPT(NK)*W%FERWE(N,NK)
!-----------------------------------------------------------------------
! transform the wave-function to real space
! and calculate the result  of the Ham. acting onto the wavefunction
!-----------------------------------------------------------------------
      CALL FFTWAV(WDES1%NPL,WDES1%NINDPW(1),W1%CR(1),W1%CW(1),GRID)

      CALL HAMNNP(GRID,NONLR_S,NONL_S,W,WDES,W1%CR(1), LREAL,SV,N,NK, &
         LMDIM,CDIJ,CHAM)
!-----------------------------------------------------------------------
! start with the exact evaluation of the eigenenergy
!-----------------------------------------------------------------------
      IF (IDUMP>=3) WRITE(*,'(F9.4)',ADVANCE='NO') REAL( W%CELEN(N,NK) ,KIND=q)
      W%CELEN(N,NK)=CHAM(N,N,NK)
      IF (IDUMP>=3) WRITE(*,'(F9.4)',ADVANCE='NO') REAL( W%CELEN(N,NK) ,KIND=q)
      EVALUE=W%CELEN(N,NK)
!-----------------------------------------------------------------------
! calculate the gradient-vector  H |phi_i> -S H_ii' |phi_i'>
!-----------------------------------------------------------------------
      CALL HAMIL0(GRID,NONLR_S,NONL_S,W,WDES,W1%CR(1), LREAL,SV,N,NK, &
         LMDIM,CQIJ,CDIJ,CHAM(1,1,NK),  CGRAD)
!-----------------------------------------------------------------------
! test orthonormality to last search direction
!-----------------------------------------------------------------------
! avoid conjugation of those elements which do not contribute to energy
      IF (ABS(W%FERWE(N,NK))<WEIMIN) THEN
      DO M=1,NPL
        W_F%CW(M,N,NK)=0
      ENDDO
      DO NPRO=1,WDES%NPRO
        W_F%CPROJ(NPRO,N,NK)=0
      ENDDO
      ENDIF

      DORT=0
      DO M=1,NPL
        DORT  =DORT+ CGRAD(M)*CONJG(W_F%CW(M,N,NK))
      ENDDO
      ORT1  =ORT1+ DORT *WEIGHT *2
!-----------------------------------------------------------------------
! calculate kinetic energy for preconditioning
! and set WPRE%CW
!-----------------------------------------------------------------------
      EKIN=0
      ENORM=0
      DO M=1,NPL
        CPT=W%CPTWFP(M,N,NK)
        ENORM=ENORM+ REAL( CPT*CONJG(CPT) ,KIND=q)
        EKIN =EKIN+ REAL( CPT*CONJG(CPT) ,KIND=q) * WDES%DATAKE(M,NK)
        WPRE%CW(M)=CGRAD(M)
      ENDDO
      EKIN=EKIN/ENORM

      IF (EKIN<0.5_q) EKIN=0.5_q
!-----------------------------------------------------------------------
! performe preconditioning of gradient
!-----------------------------------------------------------------------
      FAKT=2._q/EKIN
      DO M=1,NPL
         X =WDES%DATAKE(M,NK)/EKIN
         X2=27+X*(18+X*(12+8*X))
         WPRE%CW(M)=WPRE%CW(M)*X2/(X2+16*X*X*X*X)*FAKT
      ENDDO
!----------------------------------------------------------------------
!  calculate result of projection operatores
!----------------------------------------------------------------------
      IF (LREAL) THEN
        CALL FFTWAV(NPL,WDES1%NINDPW(1),WPRE%CR(1),WPRE%CW(1),GRID)
        CALL RPRO1(NONLR_S,WDES1,WPRE)
      ELSE
        CALL PROJ1(NONL_S,WDES1,WPRE)
      ENDIF

!-----------------------------------------------------------------------
! calculate expected energy-change into search direction
! and store this (preconditioned) gradient in W_G%CW
!-----------------------------------------------------------------------
      IF (IROT==1 .OR. IROT==3) THEN
      DE=0
      DO M=1,NPL
         DE =DE+ CGRAD(M)*CONJG(WPRE%CW(M))
      ENDDO
      DE1   =DE1   +DE*WEIGHT *2
      GNORM =GNORM +DE*WEIGHT *2

      DO M=1,NPL
        W_G%CW(M,N,NK)=WPRE%CW(M)
      ENDDO

      DO NPRO=1,WDES%NPRO
        W_G%CPROJ(NPRO,N,NK)=WPRE%CPROJ(NPRO)
      ENDDO

      ENDIF
!-----------------------------------------------------------------------
      ENDDO bands
      ENDDO kpoint
!-----------------------------------------------------------------------
! dump CHAM
!-----------------------------------------------------------------------
      IF (.TRUE.) THEN
        NK=1
        NPL2=MIN(10,WDES%NBANDS)

        DO N1=1,NPL2
          WRITE(6,681)N1,(REAL( CHAM(N1,N2,NK) ,KIND=q) ,N2=1,NPL2)
        ENDDO
        WRITE(6,*)
#ifndef gammareal
        DO N1=1,NPL2
          WRITE(6,682)N1,(AIMAG(CHAM(N1,N2,NK)),N2=1,NPL2)
        ENDDO
        WRITE(6,*)
#endif
  681   FORMAT(1I2,3X,20F9.5)
  682   FORMAT(1I2,3X,20E9.1)
      ENDIF
!=======================================================================
! calculate energy-change due to subspacerotation
! and the orthonormality to the last rotation
! the gradient is  simply  i [H,F]
! ORT   = i [H,F] * CHF (CHF== last rotation matrix)
! GNORM = i [H,F] * new rotation matrix
!=======================================================================
      IF (IROT==2 .OR. IROT==3) THEN
         DO NK=1,WDES%NKPTS
         WEIGHT =WDES%WTKPT(NK)
!-----------------------------------------------------------------------
         DE_ATT=ABS(W%CELEN(WDES%NBANDS,NK)-W%CELEN(1,NK))/8

         IF (ILOEW==1 .OR.ILOEW==2) THEN
           DIFMAX=0.001_q
           DO N2=1,WDES%NBANDS
           DO N1=1,WDES%NBANDS
              DIFFER=    (W%FERWE(N1,NK)-W%FERWE(N2,NK))
              DIFCEL= REAL( W%CELEN(N2,NK)-W%CELEN(N1,NK) ,KIND=q)
              CEXPRE=WDES%RSPIN*CHAM(N2,N1,NK)*DIFFER
              ENHANCE=1
              ! slow rotation for elements with different occupancies
              ! and very close energy
              ! this is the point where sloshing occurs
              ! (i.e. reason for  divergence of dielectric matrix
              !  in semiconductors there are no such elements !!)
              IF (ABS(DIFCEL) < ABS(DE_ATT*DIFFER)) THEN
                   ENHANCE=ENHANCE*ABS(DIFCEL/(DE_ATT*DIFFER))
              ENDIF
             ! avoid rotation if ordering of weights is wrong
              IF (ABS(DIFCEL)<DIFMAX.OR. &
                DIFCEL* REAL( W_F%CELEN(N2,NK)-W_F%CELEN(N1,NK) ,KIND=q) <-1E-2_q)  &
              THEN
                CROT  =0
              ELSE
                CROT  =CHAM(N1,N2,NK)/DIFCEL*ENHANCE
                IF (ABS(CROT)>0.1_q) THEN
                   FAKT= 0.1_q/ABS(CROT)
                   CROT  = CROT*FAKT
                ENDIF
              ENDIF
             ! avoid conjugation of those elements which do not contribute to energy
              IF (ABS(DIFFER)<WEIMIN) THEN
                CHF(N1,N2,NK)=0
              ENDIF
              CEIDB(N1,N2)=CROT
              GNORM =GNORM +CEXPRE*CROT*WEIGHT
              ORT2  =ORT2  +CEXPRE*CHF(N1,N2,NK)*WEIGHT
              DE2   =DE2   +CEXPRE*CROT*WEIGHT
           ENDDO
           ENDDO
           DO N2=1,WDES%NBANDS
           DO N1=1,WDES%NBANDS
             CHAM(N1,N2,NK)=CEIDB(N1,N2)
           ENDDO; ENDDO
!-----------------------------------------------------------------------
         ENDIF
         ENDDO
      ENDIF
!=======================================================================
! calculate gradient for fermi-weights
!=======================================================================
      IF (KPOINTS%SIGMA>0 .AND. ICEL == 1) THEN
      SHIFT=0
      EF=EFERMI

      SHIFT=0
      WSUM =0
      DO NK=1,WDES%NKPTS
      WEIGHT =WDES%WTKPT(NK)*2
      DO N =1,WDES%NBANDS
        X=  (EF- REAL( W_F%CELEN(N,NK) ,KIND=q) )/KPOINTS%SIGMA
        CALL DELSTP(KPOINTS%ISMEAR,X,DFUN,SFUN)
        WSUM =WSUM +DFUN* WEIGHT
        SHIFT=SHIFT+DFUN*(W_F%CELEN(N,NK)-W%CELEN(N,NK))/KPOINTS%SIGMA*WEIGHT
      ENDDO
      ENDDO
      SHIFT=SHIFT/WSUM
      WMIN=WSUM/WDES%NBANDS*WEIMIN

      DO NK=1,WDES%NKPTS
      WEIGHT =WDES%WTKPT(NK)*2
      DO N =1,WDES%NBANDS
          X=  (EF- REAL( W_F%CELEN(N,NK) ,KIND=q) )/KPOINTS%SIGMA
          CALL DELSTP(KPOINTS%ISMEAR,X,DFUN,SFUN)
! avoid conjugation of those elements which do not contribute to energy
          IF (ABS(DFUN)<WMIN) THEN
             W_F%FERWE(N,NK)=0
          ENDIF
          DFUN= DFUN*((W_F%CELEN(N,NK)-W%CELEN(N,NK))/KPOINTS%SIGMA-SHIFT)
          DECEL =DECEL +DFUN* REAL( W_F%CELEN(N,NK)-W%CELEN(N,NK) ,KIND=q) *WEIGHT
          GNORM =GNORM +DFUN* REAL( W_F%CELEN(N,NK)-W%CELEN(N,NK) ,KIND=q) *WEIGHT
          ORTCEL=ORTCEL+DFUN*W_F%FERWE(N,NK)*WEIGHT
          W_G%FERWE(N,NK)= W_F%CELEN(N,NK)-W%CELEN(N,NK)
      ENDDO
      ENDDO

      ENDIF
      ORT=ORT+ORTCEL+ORT1+ORT2
!=======================================================================
! performe CORRECTOR step taking into account the new gradient only
!=======================================================================
      cor: IF (IALGO==9 .AND.IFLAG/=0 .AND. MOD(ICOUNT,2)==0) THEN

      BCSTEP= BTRIAL*ORT/(GSL-ORT)
      BSTEP = BTRIAL+BCSTEP

      AVERAG=AVERAG*0.8_q+BSTEP*0.2_q

      IF (IU0>=0) &
      WRITE(IU0,12,ADVANCE='NO') BCSTEP,GSL,ORT
   12 FORMAT('cor-step=',F10.5,' , ',F13.6,F13.6)
      IF (IU0>=0) WRITE(IU0,*)
      ICOUNT=ICOUNT+1
      DESUM1= -GSL*(BCSTEP+BTRIAL)/2

! step for wavefunction

      IF (IROT==1 .OR. IROT==3) THEN
        DO NK=1,WDES%NKPTS
        NPL=WDES%NPLWKP(NK)
        DO N=1,WDES%NBANDS
          DO M=1,NPL
            W%CPTWFP(M,N,NK)=W%CPTWFP(M,N,NK)   -W_F%CW(M,N,NK)      *BCSTEP
          ENDDO

          DO NPRO=1,WDES%NPRO
            W%CPROJ(NPRO,N,NK)=W%CPROJ(NPRO,N,NK)-W_F%CPROJ(NPRO,N,NK)*BCSTEP
          ENDDO
        ENDDO
        ENDDO
      ENDIF

! step for suspace-matrix

      IF (IROT==2 .OR. IROT==3) THEN
      DO NK=1,WDES%NKPTS
        NPL=WDES%NPLWKP(NK)

        IF (ILOEW==1) THEN
        CALL ROT1(WDES%NBANDS,NK,CHF,BCSTEP,HFEIG,CUNI,CEIDB,COVL)
        ELSE IF (ILOEW==2) THEN
        CALL ROT2(WDES%NBANDS,NK,CHF,BCSTEP,HFEIG,CUNI,CEIDB,COVL)
        ENDIF

        CALL LINCOM('F',W%CPTWFP(1,1,NK),W%CPROJ(1,1,NK),CEIDB, &
          WDES%NBANDS,WDES%NBANDS,NPL,WDES%NPRO,WDES%NRPLWV,WDES%NPROD,WDES%NBANDS, &
          NBLK,W%CPTWFP(1,1,NK),W%CPROJ(1,1,NK))

      ENDDO
      ENDIF

      CALL ORTHCH(WDES,W, LOVERL, LMDIM,CQIJ,NBLK)

! step for fermi-weights

      IF (KPOINTS%SIGMA>0 .AND. ICEL == 1) THEN
        DO NK=1,WDES%NKPTS
        DO N =1,WDES%NBANDS
          W_F%CELEN(N,NK)=W_F%CELEN(N,NK)-BCSTEP*W_F%FERWE(N,NK)
          W_G%FERWE(N,NK)=W%FERWE(N,NK)
        ENDDO
        ENDDO
      ENDIF

      IFLAG=0
      ORT=0
      GOTO 1000

      ENDIF cor
!=======================================================================
! performe trial step (possibly with conjugation of search directions)
! and calculate change in energy up to first order
! the conjugate gradient optimization is restarted
! ) after 6 steps
! ) based on a 10 %  criterion
!=======================================================================
      IF ( (IALGO <8) .OR.  IFLAG/=1 &
     &    .OR. ABS(ORT) > GNORM/5 .OR.ICOUNT==-1 ) THEN
        GAMMA=0
        ICOUNT=0
      ELSE
        GAMMA=GNORM/GNORML
        ICOUNT=ICOUNT+1
        IF (GAMMA*ABS(ORT) > GNORM/5) THEN
          GAMMA=0
          ICOUNT=0
        ENDIF
      ENDIF

      IF (IDUMP>=1.AND. IU0>=0) &
     &  WRITE(IU0,20) GAMMA,DE1,DE2,DECEL,ORT1,ORT2,ORTCEL
      GNORML=GNORM
      RMS   =GNORM

  20  FORMAT(/' gam=',F6.3, &
     &  ' g(H,U,f)= ',3E10.3,' ort(H,U,f) =',3E10.3)
!-----------------------------------------------------------------------
!  DESUM0 is expected first order energy-change into the trial direction
!  the trial direction is  f(N) = p(N) + GAMMA f(N-1)
!  (g(N) gradient, p(N) preconditioned gradient, f(N) search direction)
!  because ORT= g(N) f(N-1) and GNORM = g(N) p(N)
!  DESUM is given by DESUM = g(N) * f(N) = GNORM+GAMMA* ORT
!-----------------------------------------------------------------------
      DESUM0=GNORM+GAMMA* ORT
      GSL   =GNORM+GAMMA*ORT
!-----------------------------------------------------------------------
! conjugate directions
!-----------------------------------------------------------------------
      conjug: DO NK=1,WDES%NKPTS

      IF (IROT==1 .OR. IROT==3) THEN
      NPL=WDES%NPLWKP(NK)
      DO N=1,WDES%NBANDS
        DO M=1,NPL
           W_F%CW(M,N,NK)    =W_G%CW(M,N,NK)    +W_F%CW(M,N,NK)*GAMMA
        ENDDO
        DO NPRO=1,WDES%NPRO
          W_F%CPROJ(NPRO,N,NK)=W_G%CPROJ(NPRO,N,NK)+W_F%CPROJ(NPRO,N,NK)*GAMMA
        ENDDO
      ENDDO
      ENDIF

      IF (IROT==2 .OR. IROT==3) THEN
        DO N2=1,WDES%NBANDS
        DO N1=1,WDES%NBANDS
          CHF(N1,N2,NK)=CHAM(N1,N2,NK)+GAMMA*CHF(N1,N2,NK)
        ENDDO
        ENDDO
      ENDIF

      IF (KPOINTS%SIGMA>0 .AND. (ICEL == 1 .OR. ICEL==2)) THEN
      DO N =1,WDES%NBANDS
          W_F%FERWE(N,NK)=W_G%FERWE(N,NK)+W_F%FERWE(N,NK)*GAMMA
      ENDDO
      ENDIF

      ENDDO conjug

!-----------------------------------------------------------------------
! trial step
!-----------------------------------------------------------------------
      BSTEP=BTRIAL
      IF (IROT==1 .OR. IROT==3) THEN
        DO NK=1,WDES%NKPTS
        NPL=WDES%NPLWKP(NK)
        DO N=1,WDES%NBANDS
          DO M=1,NPL
             ! W_G%CW(M,N,NK)    =W%CPTWFP(M,N,NK)
             W%CPTWFP(M,N,NK)=W%CPTWFP(M,N,NK)-W_F%CW(M,N,NK)*BSTEP
          ENDDO
          DO NPRO=1,WDES%NPRO
            ! W_G%CPROJ(NPRO,N,NK)=W%CPROJ(NPRO,N,NK)
            W%CPROJ(NPRO,N,NK)=W%CPROJ(NPRO,N,NK)-W_F%CPROJ(NPRO,N,NK)*BSTEP
          ENDDO
        ENDDO
        ENDDO
      ENDIF

! trial step for suspace-matrix

      IF (IROT==2 .OR. IROT==3) THEN
      DO NK=1,WDES%NKPTS
        IF (ILOEW==1) THEN
        CALL ROT1(WDES%NBANDS,NK,CHF,BSTEP,HFEIG,CUNI,CEIDB,COVL)
        ELSE
        CALL ROT2(WDES%NBANDS,NK,CHF,BSTEP,HFEIG,CUNI,CEIDB,COVL)
        ENDIF

       NPL=WDES%NPLWKP(NK)
       CALL LINCOM('F',W%CPTWFP(1,1,NK),W%CPROJ(1,1,NK),CEIDB, &
           WDES%NBANDS,WDES%NBANDS,NPL,WDES%NPRO,WDES%NRPLWV,WDES%NPROD,WDES%NBANDS, &
           NBLK,W%CPTWFP(1,1,NK),W%CPROJ(1,1,NK))

      ENDDO
      ENDIF

      CALL ORTHCH(WDES,W, LOVERL, LMDIM,CQIJ,NBLK)

! fermi-weights
! update fermi-weights W%FERWE, store old fermi-weights in  W_G%FERWE

      IF (KPOINTS%SIGMA>0 .AND. ICEL == 1) THEN
      DO NK=1,WDES%NKPTS
      DO N =1,WDES%NBANDS
          W_F%CELEN(N,NK)=W_F%CELEN(N,NK)-BSTEP*W_F%FERWE(N,NK)
          W_G%FERWE(N,NK)=W%FERWE(N,NK)
      ENDDO
      ENDDO

      ELSE

      DO NK=1,WDES%NKPTS
      DO N =1,WDES%NBANDS
          W_G%FERWE(N,NK)=W%FERWE(N,NK)
      ENDDO
      ENDDO

      ENDIF
!-----------------------------------------------------------------------
! set IFLAG and return to main routine
! DESUM1 is the first order energy change
!-----------------------------------------------------------------------
      IF (IALGO<=6 .OR. IALGO==9) THEN
        IFLAG=0
      ELSE
        IFLAG=1
      ENDIF
      IF (IFLAG==0) DESUM1 =-DESUM0*BSTEP
      GOTO 1000
!***********************************************************************
! CORRECTOR STEP
! IFLAG=2  goto the minimum
! ) calculate the remaining distance we have to go
! ) and step to the minimum
! IFLAG>2  increase trial step
!***********************************************************************
      ELSE IF (IFLAG>=2) THEN main
! store step and energy for later use
      NCOUNT=NCOUNT+1
      STEP(NCOUNT)  =BSTEP
      ENERGY(NCOUNT)=TOTEN2

! calculate position of minium using a quadratic interpolation
      A1= TOTEN2-TOTEN
      A2=-DESUM0*BSTEP

      IF (IDUMP>=2 .AND. IU0>=0) WRITE(IU0,11,ADVANCE='NO') BSTEP,A1,A2
   11 FORMAT(/'step=',F10.5,' , ',F13.6,F13.6)

      BOPT  = -A2/(A1-A2)/2
      DESUM1= (A2+(A1-A2)*BOPT)*BOPT
!-----------------------------------------------------------------------
! IFLAG==2
!-----------------------------------------------------------------------
      IF (IFLAG==2) THEN

!---- first trial step, search near minimum

        IF (NCOUNT==1) THEN
! abort loop (if the next line is commented out the position of
!             the minimum is determined with higher precission
           IFLAG=0
           IF (IDUMP>=1 .AND.IU0>=0 ) &
           WRITE(IU0,10) GAMMA,BTRIAL,BOPT*BSTEP,AVERAG
   10 FORMAT(/' gam=',F6.3,' trial=',F6.3,'  step=',F6.3,' mean=',F6.3)

          STRIAL1=BSTEP*BOPT
          STRIAL2=BSTEP*BOPT*0.05_q

          IF (IFLAG==0) STRIAL2=0
          IF (BOPT<0)  THEN
             IF (IU0>=0) &
             WRITE(IU0,*)'BOPT < 0, no corrections made'
             GOTO 1000
          ELSE
             BCSTEP=(STRIAL1-STRIAL2)-BSTEP
             BSTEP =(STRIAL1-STRIAL2)
          ENDIF

!---- second trial step, increase BSTEP by factor 0.01

        ELSE IF (NCOUNT==2 .OR.NCOUNT==3) THEN
           BCSTEP=STRIAL2
           BSTEP =BSTEP+STRIAL2

!---- linear approximation between two search steps near minimum

        ELSE IF (NCOUNT==4) THEN
! abort loop
           IFLAG=0
           BOPT= STRIAL1+STRIAL2*(ENERGY(2)- ENERGY(4)) &
     &         / (2*ENERGY(2)-4*ENERGY(3)+2*ENERGY(4))
           A1  = (-ENERGY(2)+ENERGY(4))/2/STRIAL2
           A2  = ( ENERGY(2)-2*ENERGY(3)+ENERGY(4))/2/STRIAL2**2
           DESUM1=A1*(BOPT-STRIAL1)+A2*(BOPT-STRIAL1)**2
           IF (IU0>=0) &
           WRITE(IU0,10) GAMMA,BSTEP,BOPT,AVERAG
           BCSTEP=BOPT-BSTEP
           BSTEP =BOPT
        ENDIF
!-----------------------------------------------------------------------
! IFLAG>2
!-----------------------------------------------------------------------
      ELSE
!---- increase trial step by BTRIAL
        BCSTEP=BTRIAL
        BSTEP =BSTEP+BTRIAL
      ENDIF
      AVERAG=AVERAG*0.8_q+BSTEP*0.2_q

      NK=1
!-----------------------------------------------------------------------
! trial step
!-----------------------------------------------------------------------
      IF (IROT==1 .OR. IROT==3) THEN

      DO NK=1,WDES%NKPTS
      NPL=WDES%NPLWKP(NK)
      DO N=1,WDES%NBANDS

        DO M=1,NPL
          ! W%CPTWFP(M,N,NK)=W_G%CW(M,N,NK)   -W_F%CW(M,N,NK)      *BSTEP
          W%CPTWFP(M,N,NK)=W%CPTWFP(M,N,NK)   -W_F%CW(M,N,NK)      *BCSTEP
        ENDDO
        DO NPRO=1,WDES%NPRO
          ! W%CPROJ(NPRO,N,NK)=W_G%CPROJ(NPRO,N,NK)-W_F%CPROJ(NPRO,N,NK)*BSTEP
          W%CPROJ(NPRO,N,NK)=W%CPROJ(NPRO,N,NK)-W_F%CPROJ(NPRO,N,NK)*BCSTEP
        ENDDO

      ENDDO
      ENDDO

      ENDIF

! step for suspace-matrix

      IF (IROT==2 .OR. IROT==3) THEN
      DO NK=1,WDES%NKPTS

        NPL=WDES%NPLWKP(NK)

        IF (ILOEW==1) THEN
        ! CALL ROT1(WDES%NBANDS,NK,CHF,BSTEP,HFEIG,CUNI,CEIDB,COVL)
          CALL ROT1(WDES%NBANDS,NK,CHF,BCSTEP,HFEIG,CUNI,CEIDB,COVL)
        ELSE IF (ILOEW==2) THEN
        ! CALL ROT2(WDES%NBANDS,NK,CHF,BSTEP,HFEIG,CUNI,CEIDB,COVL)
          CALL ROT2(WDES%NBANDS,NK,CHF,BCSTEP,HFEIG,CUNI,CEIDB,COVL)
        ENDIF

        CALL LINCOM('F',W%CPTWFP(1,1,NK),W%CPROJ(1,1,NK),CEIDB, &
          WDES%NBANDS,WDES%NBANDS,NPL,WDES%NPRO,WDES%NRPLWV,WDES%NPROD,WDES%NBANDS, &
          NBLK,W%CPTWFP(1,1,NK),W%CPROJ(1,1,NK))

      ENDDO
      ENDIF

      CALL ORTHCH(WDES,W, LOVERL, LMDIM,CQIJ,NBLK)

! step for fermi-weights

      IF (KPOINTS%SIGMA>0 .AND. ICEL == 1) THEN
        DO NK=1,WDES%NKPTS
        DO N =1,WDES%NBANDS
            W_F%CELEN(N,NK)=W_F%CELEN(N,NK)-BCSTEP*W_F%FERWE(N,NK)
        ENDDO
        ENDDO
      ENDIF

!-----------------------------------------------------------------------
      ENDIF main

 1000 CONTINUE
      CALL DELWAV(WPRE,.TRUE.)
#endif
      RETURN
      END SUBROUTINE

!************************ SUBROUTINE ROTDIA ****************************
!
! this subroutine calculates the rotation matrix  for one k-point
! the commutator CHF=[F,H] must be supplied by the calling routine
! on return the rotation matrix is stored in   CEIDB
! the rotation matrix is
!    -i d [F,H]       -i d E +                     +
!   e            = U e      U    with [F,H] = U E U
!
! i.e. U is the matrix which diagonolizes [F,H] and E is the diagonal
!      matrix containing the eigenvalues
!  routine is no longer used
!***********************************************************************

      SUBROUTINE ROTDIA(NBANDS,NK,CHF,BSTEP,HFEIG,W,NW,CUNI,CEIDB,CTMP)
      USE prec

      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)
      DIMENSION CHF(NBANDS,NBANDS,NK)
      DIMENSION CTMP(NBANDS,NBANDS),CUNI(NBANDS,NBANDS), &
     &          CEIDB(NBANDS,NBANDS)

      DIMENSION HFEIG(NBANDS),W(NW)

!=======================================================================
! Diagononalize the resulting Matrix (search direction)
!=======================================================================
      DO N1=1,NBANDS
      DO N2=1,NBANDS
        CUNI(N1,N2)=CHF(N1,N2,NK)
      ENDDO
      ENDDO

      CALL ZHEEV('V','U',NBANDS,CUNI,NBANDS,HFEIG,CTMP,NBANDS*NBANDS, &
     &         W,  IFAIL)
      IF (IFAIL/=0) THEN
         WRITE(*,*) 'ERROR ROTDIA: Call to routine ZHEEV failed! '// &
     &              'Error code was ',IFAIL
         STOP
      ENDIF
!=======================================================================
! set up the unitary transformation CEIDB for wavefunctions
!=======================================================================
      DO N1=1,NBANDS
      DO N2=1,NBANDS
         CTMP(N1,N2)=CUNI(N1,N2)*EXP(-BSTEP*(0._q,1._q)*HFEIG(N2))
      ENDDO
      ENDDO
      CALL ZGEMM( 'N', 'C', NBANDS, NBANDS, NBANDS, (1.,0.), CTMP, &
     &             NBANDS, CUNI, NBANDS, (0._q,0._q), CEIDB, NBANDS)
!=======================================================================
! rotate the wavefunction (with CEIDB): warning about LINCOM! We must
! take the transposed of CEIDB because MATMUL uses the transposed ...
!=======================================================================
      DO N1=1,NBANDS
      DO N2=1,NBANDS
         CTMP(N1,N2)=CEIDB(N2,N1)
      ENDDO
      ENDDO

      DO N1=1,NBANDS
      DO N2=1,NBANDS
         CEIDB(N2,N1)=CTMP(N2,N1)
      ENDDO
      ENDDO
      RETURN
      END SUBROUTINE


!************************ SUBROUTINE ROT1  ****************************
!
! this subroutine calculates the rotation matrix
!  CEIDB = 1 - CHF* BSTEP
! for one k-point, where CHF is the Loewdin-perturbation matrix
! which must be supplied by the calling routine
!
!***********************************************************************

      SUBROUTINE ROT1   (NBANDS,NK,CHF,BSTEP,HFEIG,CUNI,CEIDB,CTMP)
      USE prec

      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)
      GDEF CHF(NBANDS,NBANDS,NK)
      GDEF CTMP(NBANDS,NBANDS),CUNI(NBANDS,NBANDS), &
     &     CEIDB(NBANDS,NBANDS)

      DIMENSION HFEIG(NBANDS)

!=======================================================================
! CDEIB contains the searchdirection
! take the transposed of CEIDB because MATMUL uses the transposed ...
!=======================================================================
      DO N1=1,NBANDS
      DO N2=1,NBANDS
         CEIDB(N2,N1)=CHF(N1,N2,NK)*(-BSTEP)
      ENDDO
      ENDDO

      DO N1=1,NBANDS
         CEIDB(N1,N1)=1
      ENDDO

      RETURN
      END SUBROUTINE

!************************ SUBROUTINE ROT2   ****************************
!
! this subroutine calculates the rotation matrix
!    - CHF d          - E d  +                   +
!   e            = U e      U    with CHF = U E U
!
! for one k-point, the Loewdin-perturbation CHF matrix must be supplied
! the resulting matrix is unitar
!
! i.e. U is the matrix which diagonolizes CHF   and E is the diagonal
!      matrix containing the eigenvalues
! up to first order the result from ROT1 and ROT2 is the same !
!
!***********************************************************************

      SUBROUTINE ROT2(NBANDS,NK,CHF,BSTEP,HFEIG,CUNI,CEIDB,CTMP)
      USE prec

      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)
      GDEF CHF(NBANDS,NBANDS,NK)
      GDEF CTMP(NBANDS,NBANDS),CUNI(NBANDS,NBANDS), &
     &     CEIDB(NBANDS,NBANDS)

      DIMENSION HFEIG(NBANDS),W(3*NBANDS)

!=======================================================================
! Diagononalize the resulting Matrix (search direction)
!=======================================================================
      DO N1=1,NBANDS
      DO N2=1,NBANDS
        CUNI(N1,N2)=CHF(N1,N2,NK)
      ENDDO
      ENDDO

#ifdef  gammareal
      CALL DSYEV &
     &         ('V','U',NBANDS,CUNI,NBANDS,HFEIG,CTMP,NBANDS*NBANDS, &
     &           IFAIL)
#else
      CALL ZHEEV &
     &         ('V','U',NBANDS,CUNI,NBANDS,HFEIG,CTMP,NBANDS*NBANDS, &
     &           W,  IFAIL)
#endif
      IF (IFAIL/=0) THEN
         WRITE(*,*) 'ERROR ROTDIA: Call to routine ZHEEV failed! '// &
     &              'Error code was ',IFAIL
         STOP
      ENDIF
!=======================================================================
! set up the unitary transformation CEIDB for wavefunctions
!=======================================================================
      DO N1=1,NBANDS
      DO N2=1,NBANDS
         CTMP(N1,N2)=CUNI(N1,N2)*EXP(-BSTEP*HFEIG(N2))
      ENDDO
      ENDDO
      CALL GGEMM( 'N',trans, NBANDS, NBANDS, NBANDS, one, CTMP, &
     &             NBANDS, CUNI, NBANDS, zero, CEIDB, NBANDS)
!=======================================================================
! rotate the wavefunction (with CEIDB): warning about MATMUL! We must
! take the transposed of CEIDB because MATMUL uses the transposed ...
!=======================================================================
      DO N1=1,NBANDS
      DO N2=1,NBANDS
         CTMP(N1,N2)=CEIDB(N2,N1)
      ENDDO
      ENDDO
      DO N1=1,NBANDS
      DO N2=1,NBANDS
         CEIDB(N2,N1)=CTMP(N2,N1)
      ENDDO
      ENDDO
      RETURN
      END SUBROUTINE
      END MODULE
