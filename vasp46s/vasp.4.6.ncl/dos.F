
#include "symbol.inc"
!********************DENINI**********************************************
! RCS:  $Id: dos.F,v 1.3 2002/08/14 13:59:37 kresse Exp $
!
! initialise the fermi-weights
! if the total number of bands is larger than the number of occupied
! bands:
! initialize the fermi occupation function
! FERWE  1  for the lowest fully occupied nelect/2 bands
!       .5  for the (nelect/2+1)th band if nelect odd
!        0  for the empty bands
!
! non-metallic systems set all components of FERWE equal 1
!
!************************************************************************

      SUBROUTINE  DENINI(FERWE,NBANDS,NKPTS,ELECT)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      DIMENSION FERWE(NBANDS,NKPTS)

       IF (ABS(NBANDS-ELECT/2._q)>0.000001_q) THEN
          NBE=INT(ELECT/2._q)
          DO NB=1,NBE
          DO KP=1,NKPTS
            FERWE(NB,KP)=1._q
          ENDDO
          ENDDO
          IF (ABS(ELECT/2._q-NBE)>0.000001_q) THEN
            DO KP=1,NKPTS
              FERWE(NBE+1,KP)=0.5_q*(ELECT/2._q-NBE)
            ENDDO
          ELSE
            DO KP=1,NKPTS
              FERWE(NBE+1,KP)=0
            ENDDO
          ENDIF
          DO NB=NBE+2,NBANDS
          DO KP=1,NKPTS
             FERWE(NB,KP)=0
          ENDDO
          ENDDO
        ELSE
          DO NB=1,NBANDS
          DO KP=1,NKPTS
            FERWE(NB,KP)=1
          ENDDO
          ENDDO
        ENDIF
      RETURN
      END SUBROUTINE  DENINI



!***********************SUBROUTINE DENSTA*******************************
! DENSTA is the dispatcher to the routines calculating the
! partial occupancies
! depending on ISMEAR the required algorithm is chosen
! ISMEAR=0     gaussian smearing (DENMP)
! ISMEAR>0     the generalized form of Methfessel and Paxton of order
!              N=ISMEAR will be used instead of Gaussians to get the dos ...
! ISMEAR=-1    fermi-smearing DENSTA calls DENFER
! ISMEAR=-4,-5 the tetrahedron method will be used
!
! The width of the smearing is given by SIGMA
! corrections to the total energy are calculated and returned in ENTROPY
! if these corrections are added to the total energy the free
! variational energy is obtained
!
! in addition the routine also calculates the DOS and the integrated
! and if required also the partial DOS (NIOND and LDIMP must be set)
! 
! this is the worst spagethi code in VASP, but still it works
!
!***********************************************************************


      SUBROUTINE DENSTA( IU0, IU, WDES, W, KPOINTS, NELECT, &
              NUP_DOWN, ENTROPY, EFERMI, SIGMA, &
              NEDOS, LDIMP, NIOND, DOS, DOSI, PAR, DOSPAR)

      USE prec
      USE wave
      USE mkpoints

      TYPE (wavedes)        WDES
      TYPE (wavespin)       W
      TYPE (kpoints_struct) KPOINTS
      INTEGER NEDOS, LDIMP, NIOND, IU, IU0, JOBPAR
      REAL(q) NELECT,NUP_DOWN,ENTROPY,SIGMA,EFERMI,ENTROPY_
      REAL(q) DOS(NEDOS,WDES%ISPIN),DOSI(NEDOS,WDES%ISPIN)
      REAL(q) PAR(WDES%NB_TOT, WDES%NKDIM, LDIMP, NIOND, WDES%NCDIJ)
      REAL(q) DOSPAR(NEDOS,LDIMP,NIOND,WDES%NCDIJ)
! local variables
      REAL(q) EADD,EPS,NEL,NEL_SHIFT
      REAL(q) RSPIN
      LOGICAL  :: LAUTO=.FALSE.
      INTEGER  :: ICALLS=0,ISP,ISP2,K,N,ISPIN,ISPIN_MAX

      ENTROPY=0
      ENTROPY_=0

      IF (NIOND == 0 .OR. LDIMP==0) THEN
        JOBPAR=0
      ELSE
        JOBPAR=1
      ENDIF

! set KPOINTS%EMAX and KPOINTS%EMIN if required
      ! first call initialize LAUTO
      IF (ICALLS==0) THEN
         LAUTO=(KPOINTS%EMAX<=KPOINTS%EMIN)
         ICALLS=1
      ENDIF

      IF (LAUTO) THEN
        KPOINTS%EMIN=1.E30_q
        KPOINTS%EMAX=-KPOINTS%EMIN
        DO ISP=1,WDES%ISPIN
        DO K=1,WDES%NKPTS
        DO N=1,WDES%NB_TOT
          EPS=W%CELTOT(N,K,ISP)
          KPOINTS%EMAX=MAX(KPOINTS%EMAX,EPS)
          KPOINTS%EMIN=MIN(KPOINTS%EMIN,EPS)
        ENDDO
        ENDDO
        ENDDO

        EADD=(KPOINTS%EMAX-KPOINTS%EMIN)*0.05_q
        EADD=MAX(EADD,10 *ABS(SIGMA))
        KPOINTS%EMIN=KPOINTS%EMIN-EADD
        KPOINTS%EMAX=KPOINTS%EMAX+EADD
      ENDIF

      ISPIN_MAX=1

      ISPIN=WDES%ISPIN
      RSPIN=WDES%RSPIN

      NEL_SHIFT=0
    ! calculation for a specific number of electrons in up and down component
    ! this actually requires me to do the worst fiddling ;(
      IF (NUP_DOWN >= 0 .AND. WDES%ISPIN>1 ) THEN
         ISPIN_MAX=2
         ISPIN    =1
         RSPIN    =2
         NEL_SHIFT=NUP_DOWN
    !
    ! non collinear calculation and partial DOS required
    ! also fiddle a little bit (call the DOS routines 4 time each time
    ! with a different PAR and DOSPAR)
      ELSE IF (WDES%LNONCOLLINEAR .AND. JOBPAR==1) THEN
         ISPIN_MAX=4
      ENDIF
    !
    ! now set the occupancies and calculate the DOS
    !
    spin: DO ISP=1,ISPIN_MAX
      ISP1=ISP
      IF (WDES%LNONCOLLINEAR) ISP1=1
      ISP2=ISP1+ISPIN-1

      IF (ISP==1) THEN
         NEL=NELECT+NEL_SHIFT
      ELSE
         NEL=NELECT-NEL_SHIFT
      ENDIF

      IF (KPOINTS%ISMEAR==-1) THEN
        CALL DENFER(IU0, WDES, ISPIN, RSPIN, KPOINTS%EMIN,KPOINTS%EMAX, NEL, ENTROPY_, &
           EFERMI, SIGMA, W%FERWE(:,:,ISP1:ISP2), W%CELEN(:,:,ISP1:ISP2),&
           NEDOS, LDIMP, NIOND, DOS(1,ISP1), DOSI(1,ISP1), &
           PAR(1,1,1,1,ISP), DOSPAR(1,1,1,ISP), JOBPAR)
        CALL MRG_FERWE(WDES,W)
      ELSE IF (KPOINTS%ISMEAR==-4) THEN
        CALL DENTET(IU0, W%CELTOT(1,1,ISP1),WDES%WTKPT(1),WDES%NB_TOT,WDES%NKPTS, &
           DOS(1,ISP1),DOSI(1,ISP1), &
           NEDOS,ISPIN,RSPIN,KPOINTS%EMIN,KPOINTS%EMAX,KPOINTS%IDTET(0,1),KPOINTS%NTET, &
           KPOINTS%VOLWGT,NEL,EFERMI,W%FERTOT(1,1,ISP1), &
           ENTROPY_,-2,IU,PAR(1,1,1,1,ISP), DOSPAR(1,1,1,ISP),WDES%NKDIM,LDIMP,NIOND,JOBPAR)
      ELSE IF (KPOINTS%ISMEAR==-5) THEN
        CALL DENTET(IU0, W%CELTOT(1,1,ISP1),WDES%WTKPT(1),WDES%NB_TOT,WDES%NKPTS, &
           DOS(1,ISP1),DOSI(1,ISP1), &
           NEDOS,ISPIN,RSPIN,KPOINTS%EMIN,KPOINTS%EMAX,KPOINTS%IDTET(0,1),KPOINTS%NTET,  &
           KPOINTS%VOLWGT,NEL,EFERMI,W%FERTOT(1,1,ISP1), &
           ENTROPY_,2,IU,PAR(1,1,1,1,ISP), DOSPAR(1,1,1,ISP),WDES%NKDIM,LDIMP,NIOND,JOBPAR)
      ELSE
        CALL DENMP(IU0, WDES, ISPIN, RSPIN, KPOINTS%EMIN,KPOINTS%EMAX, NEL, ENTROPY_, &
           EFERMI, KPOINTS%ISMEAR, SIGMA, &
           W%FERWE(:,:,ISP1:ISP2) , W%CELEN(:,:,ISP1:ISP2),&
           NEDOS, LDIMP, NIOND,  DOS(1,ISP1), DOSI(1,ISP1), &
           PAR(1,1,1,1,ISP), DOSPAR(1,1,1,ISP), JOBPAR)

        CALL MRG_FERWE(WDES,W)
      ENDIF
    ! for ISMEAR>=30, Methfessel-Paxton smearing of order ISMEAR-30 is used to
    ! calculate the  Fermi weights, but the tetrahedron method is used to get
    ! the density of states (one of the secret flag settings in VASP ...)
      ENTROPY=ENTROPY+ENTROPY_

      IF (KPOINTS%LTET .AND. KPOINTS%ISMEAR /=-4 .AND. KPOINTS%ISMEAR /= -5 ) THEN
        CALL DENTET(IU0,W%CELTOT(1,1,ISP1),WDES%WTKPT(1),WDES%NB_TOT,WDES%NKPTS, &
           DOS(1,ISP1),DOSI(1,ISP1), &
           NEDOS,ISPIN,RSPIN,KPOINTS%EMIN,KPOINTS%EMAX,KPOINTS%IDTET(0,1),KPOINTS%NTET, &
           KPOINTS%VOLWGT,NEL,EFERMI,W%FERTOT(1,1,ISP1:ISP2), &
           ENTROPY_,0,IU,PAR(1,1,1,1,ISP), DOSPAR(1,1,1,ISP),WDES%NKDIM,LDIMP,NIOND,JOBPAR)
      ENDIF

      ENDDO spin

    ! due to all this fiddling the entropy term is not right
    ! we need to correct this now
      IF (NUP_DOWN >= 0 .AND. WDES%ISPIN>1 ) THEN
         ENTROPY=ENTROPY/2
      ELSE IF (WDES%LNONCOLLINEAR .AND. JOBPAR==1) THEN
         ENTROPY=ENTROPY/4
      ENDIF

      CONTAINS

!***********************SUBROUTINE DENNP *******************************
!
! if ISMEAR=0
! subroutine DENSTA calculates a continuous density of states in the
! interval (EMIN,EMAX) by applying a gaussian broadening to the discrete
! eigenvalue spectrum contained in CELEN(NBANDS,NKPTS). The width of the
! gaussians is SIGMA. The fermi energy EFERMI is calculated from the
! integrated dos
! correction to the variational energy is calculated (-SIGMA S)
! according to A.de Vita
! if ISMEAR>0 the generalized form of Methfessel and Paxton of order
!        N=ISMEAR will be used instead of Gaussians to get the dos ...
! routine is parallelized to get full speed ...
! initially it took 2 seconds to calculate occupancies (gK)
!
!***********************************************************************

     SUBROUTINE DENMP(IU0, WDES, ISPIN, RSPIN, EMIN, EMAX, NELECT, ECORR, EFERMI, &
              ISMEAR, SIGMA, FERWE, CELEN , &
              NEDOS, LDIMP, NIOND, DOS, DOSI, PAR, DOSPAR, JOBPAR)

      USE prec
      USE constant
      USE wave
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (wavedes)        WDES
      REAL(q) DOS(NEDOS,ISPIN),DOSI(NEDOS,ISPIN),NELECT
      REAL(q) PAR(WDES%NB_TOT,WDES%NKDIM,LDIMP,NIOND,ISPIN)
      REAL(q) DOSPAR(NEDOS,LDIMP,NIOND,ISPIN)
      REAL(q) FERWE(:,:,:)
      COMPLEX(q) CELEN(:,:,:)
! local variables
      LOGICAL LOWB,HIGHB
      INTEGER ISP

      SIGMA_=ABS(SIGMA)
      IF (SIGMA_==0) RETURN

      DELTAE=(EMAX-EMIN)/(NEDOS-1)
!=======================================================================
! initialize arrays for dos and integr. dos
!=======================================================================
      IF (JOBPAR==1) THEN
         DOSPAR=0
      ENDIF

      DOS =0
      DOSI=0
!=======================================================================
! accumulate dos and integrated dos
!=======================================================================
      NCALL=0
      DO ISP=1,ISPIN
      DO K=1,WDES%NKPTS
      DO N=1,WDES%NBANDS
        EPS=CELEN(N,K,ISP)
        WEIGHT= RSPIN*WDES%WTKPT(K)

        NELOW=(EPS-8._q*SIGMA_-EMIN)/DELTAE+1
        NEHIG=(EPS+8._q*SIGMA_-EMIN)/DELTAE+1
        IF (NELOW<1)     NELOW=1
        IF (NELOW>NEDOS) NELOW=NEDOS
        IF (NEHIG<1)     NEHIG=1
        IF (NEHIG>NEDOS) NEHIG=NEDOS

        SFUN_DONE=0
        DO I=NELOW,NEHIG
          E=EMIN+DELTAE*(I-1)-EPS
          CALL DELSTP(ISMEAR,(E/SIGMA_),DFUN,SFUN)
          EPSDOS=DFUN/SIGMA_
!gK fix the DOS so that the integrated DOS yields accurate results
          EPSDOS=(SFUN-SFUN_DONE)/DELTAE
          SFUN_DONE=SFUN

          DOS(I,ISP) =DOS(I,ISP) +(WEIGHT*EPSDOS)
          DOSI(I,ISP)=DOSI(I,ISP)+WEIGHT*SFUN
          IF (JOBPAR==1) THEN
             DO NI=1,NIOND
             DO LP=1,LDIMP
#ifdef MPI
                NB_GLOBAL=(N-1)*WDES%COMM_INTER%NCPU+WDES%COMM_INTER%NODE_ME
#else
                NB_GLOBAL=N
#endif
                DOSPAR(I,LP,NI,ISP)=DOSPAR(I,LP,NI,ISP)+ &
                         (WEIGHT*EPSDOS)*PAR(NB_GLOBAL,K,LP,NI,ISP)
             ENDDO; ENDDO
          ENDIF
        ENDDO
        DO I=NEHIG+1,NEDOS
          DOSI(I,ISP)=DOSI(I,ISP)+WEIGHT
        ENDDO
      ENDDO
      ENDDO
      ENDDO
      CALLMPI (M_sum_d(WDES%COMM_INTER, DOS(1,1),  NEDOS*ISPIN))
      CALLMPI (M_sum_d(WDES%COMM_INTER, DOSI(1,1), NEDOS*ISPIN))
      IF (JOBPAR==1) THEN
        CALLMPI (M_sum_d(WDES%COMM_INTER, DOSPAR(1,1,1,1), NEDOS*LDIMP*NIOND*ISPIN))
      ENDIF
!=======================================================================
! calculate approximated fermi energy
!=======================================================================
      DO I=1,NEDOS
        DOSTOT=DOSI(I,1)
        IF (ISPIN==2) DOSTOT=DOSTOT+DOSI(I,2)
        IF (ABS(DOSTOT-NELECT)<0.01_q .OR.DOSTOT>NELECT) EXIT
      ENDDO
      EFERMI= EMIN+(I-1)*DELTAE

      IF (SIGMA<1E-5_q) RETURN
!=======================================================================
! search now for exact Fermi-level using bisectioning
!=======================================================================
      EF1= EMIN+(I-2)*DELTAE
      LOWB =.FALSE.
      EF2= EMIN+(I-1)*DELTAE
      HIGHB=.FALSE.
      NITER=0

      setfermi: DO

      EFERMI=(EF1+EF2)/2
      NITER=NITER+1

      ELECT=0
      DO ISP=1,ISPIN
      DO K=1,WDES%NKPTS
      DO N=1,WDES%NBANDS
        EPS=CELEN(N,K,ISP)
        X1=(EFERMI-EPS)/SIGMA_
        CALL DELSTP(ISMEAR,X1,DFUN,SFUN)
        FERWE(N,K,ISP)=SFUN
        ELECT=ELECT+FERWE(N,K,ISP)*WDES%WTKPT(K)
      ENDDO
      ENDDO
      ENDDO
      ELECT=ELECT*RSPIN
      CALLMPI( M_sum_d( WDES%COMM_INTER, ELECT, 1))

      ! compare with number of electrons

      IF ( ABS(ELECT-NELECT)<1E-10_q) GOTO 110
      IF ( (ABS(EF1-EF2)/(ABS(EFERMI)+1.E-10_q))<1E-14_q) GOTO 120
      IF ( ELECT>NELECT) THEN
        IF (.NOT.LOWB)  EF1=EF1-DELTAE
        HIGHB=.TRUE.
        EF2  =EFERMI
      ELSE
        IF (.NOT.HIGHB) EF2=EF2+DELTAE
        LOWB =.TRUE.
        EF1  =EFERMI
      ENDIF
      ENDDO setfermi

  120 CONTINUE
      IF (IU0>=0) THEN
      WRITE(IU0,*)' WARNING: DENMP: can''t reach specified precision'
      WRITE(IU0,*)' Number of Electrons is NELECT =',ELECT
      ENDIF
  110 CONTINUE

!=======================================================================
! calculate entropy -SIGMA_ * S
!=======================================================================
      ECORR=0
      DO ISP=1,ISPIN
      DO K=1,WDES%NKPTS
      DO N=1,WDES%NBANDS
        X1=(EFERMI-CELEN(N,K,ISP))/SIGMA_
        CALL DELSTP(ISMEAR,X1,DFUN,SFUN)
        ECORR=ECORR+0.5_q*DFUN*WDES%WTKPT(K)
        IF (ISMEAR>0) THEN
          CALL DELSTP((ISMEAR-1),X1,DFUN,SFUN)
          ECORR=ECORR-0.5_q*DFUN*WDES%WTKPT(K)
        ENDIF
      ENDDO
      ENDDO
      ENDDO

      ECORR=-ECORR*SIGMA_*RSPIN
      CALLMPI( M_sum_d( WDES%COMM_INTER, ECORR, 1))

      RETURN
      END SUBROUTINE

!***********************SUBROUTINE DENFER******************************
!
! subroutine DENFER calculates a continuous density of states in the
! interval (EMIN,EMAX) by applying a fermi broadening to the discrete
! eigenvalue spectrum contained in CELEN(NBANDS,NKPTS). The width of the
! broadening is SIGMA_ (=1/BETA). EFERMI is calculated so that
! sum over occupation-numbers is equal to number of electrons
! in addition the correction to the variational energy is calculated
! as proposed by M. Weinert J.W. Davenport, Phys Rev B 45, 13709 (1992)
!
!**********************************************************************

     SUBROUTINE DENFER(IU0, WDES, ISPIN, RSPIN, EMIN, EMAX, NELECT,  &
              ECORR, EFERMI, SIGMA, FERWE, CELEN , &
              NEDOS, LDIMP, NIOND, DOS, DOSI, PAR, DOSPAR, JOBPAR)

      USE prec
      USE wave
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)
      INTEGER ISPIN

      TYPE (wavedes)        WDES
      REAL(q) DOS(NEDOS,ISPIN),DOSI(NEDOS,ISPIN),NELECT
      REAL(q) PAR(WDES%NB_TOT,WDES%NKDIM,LDIMP,NIOND,ISPIN)
      REAL(q) DOSPAR(NEDOS,LDIMP,NIOND,ISPIN)
      REAL(q) FERWE(:,:,:)
      COMPLEX(q) CELEN(:,:,:)
! local variables
      LOGICAL LOWB,HIGHB
      INTEGER ISP

! fermi function and its derivative
      F(E,SIG)=  1/(1 + EXP(E /SIG))
      G(E,SIG)=  1/EXP(E/SIG)/(1 + EXP(-E/SIG))**2/SIG

      SIGMA_=ABS(SIGMA)
      IF (SIGMA_==0) RETURN

      DELTAE=(EMAX-EMIN)/(NEDOS-1)
!=======================================================================
! initialize arrays for dos and integr. dos
!=======================================================================
      IF (JOBPAR==1) THEN
         DOSPAR=0
      ENDIF

      DOS =0
      DOSI=0

!=======================================================================
! accumulate dos and integrated dos
!=======================================================================

      DO ISP=1,ISPIN
      DO K=1,WDES%NKPTS
      DO N=1,WDES%NBANDS
        EPS=CELEN(N,K,ISP)
        WEIGHT= RSPIN*WDES%WTKPT(K)

        NELOW=(EPS-8._q*SIGMA_-EMIN)/DELTAE+1
        NEHIG=(EPS+8._q*SIGMA_-EMIN)/DELTAE+1
        IF (NELOW<1)     NELOW=1
        IF (NELOW>NEDOS) NELOW=NEDOS
        IF (NEHIG<1)     NEHIG=1
        IF (NEHIG>NEDOS) NEHIG=NEDOS

        SFUN_DONE=0
        DO I=NELOW,NEHIG
          E=EMIN+DELTAE*(I-1)-EPS
          SFUN=F(-E,SIGMA_)
          EPSDOS=G(E,SIGMA_)
!gK fix the DOS so that the integrated DOS yields accurate results
          EPSDOS=(SFUN-SFUN_DONE)/DELTAE
          SFUN_DONE=SFUN

          DOS(I,ISP) =DOS(I,ISP) +(WEIGHT*EPSDOS)
          DOSI(I,ISP)=DOSI(I,ISP)+WEIGHT*SFUN
          IF (JOBPAR==1) THEN
             DO NI=1,NIOND
             DO LP=1,LDIMP
#ifdef MPI
                NB_GLOBAL=(N-1)*WDES%COMM_INTER%NCPU+WDES%COMM_INTER%NODE_ME
#else
                NB_GLOBAL=N
#endif
                DOSPAR(I,LP,NI,ISP)=DOSPAR(I,LP,NI,ISP)+ &
                          (WEIGHT*EPSDOS)*PAR(NB_GLOBAL,K,LP,NI,ISP)
             ENDDO; ENDDO
          ENDIF
        ENDDO
        DO I=NEHIG+1,NEDOS
          DOSI(I,ISP)=DOSI(I,ISP)+WEIGHT
        ENDDO
      ENDDO
      ENDDO
      ENDDO

      CALLMPI (M_sum_d(WDES%COMM_INTER, DOS(1,1), NEDOS*ISPIN))
      CALLMPI (M_sum_d(WDES%COMM_INTER, DOSI(1,1), NEDOS*ISPIN))
      IF (JOBPAR==1) THEN
        CALLMPI (M_sum_d(WDES%COMM_INTER, DOSPAR(1,1,1,1), NEDOS*LDIMP*NIOND*ISPIN))
      ENDIF
!=======================================================================
! calculate approximated fermi energy
!=======================================================================
      DO I=1,NEDOS
        DOSTOT=DOSI(I,1)
        IF (ISPIN==2) DOSTOT=DOSTOT+DOSI(I,2)
        IF (ABS(DOSTOT-NELECT)<0.01_q .OR.DOSTOT>NELECT) EXIT
      ENDDO

      EFERMI= EMIN+(I-1)*DELTAE
      IF (SIGMA<1E-5_q) RETURN

!=======================================================================
! search now for exact Fermi-level using bisectioning
!=======================================================================
      EF1= EMIN+(I-2)*DELTAE
      LOWB =.FALSE.
      EF2= EMIN+(I-1)*DELTAE
      HIGHB=.FALSE.

      setfermi: DO
      EFERMI=(EF1+EF2)/2

      ELECT=0
      DO ISP=1,ISPIN
      DO K=1,WDES%NKPTS
      DO N=1,WDES%NBANDS
        EPS=CELEN(N,K,ISP)
!HKim modification begin
        XX = (EPS-EFERMI)/SIGMA_
        IF ( XX < -40.0_q ) THEN
           FERWE(N,K,ISP) = 1.0_q
        ELSE IF ( XX > 40.0_q ) THEN
           FERWE(N,K,ISP) = 0.0_q
        ELSE
           FERWE(N,K,ISP) = F(EPS-EFERMI,SIGMA_)
        ENDIF
!HKim modification end
        ELECT=ELECT+FERWE(N,K,ISP)*WDES%WTKPT(K)
      ENDDO
      ENDDO
      ENDDO

      ELECT=ELECT*RSPIN
      CALLMPI( M_sum_d( WDES%COMM_INTER, ELECT, 1))

      ! compare with number of electrons

      IF ( ABS(ELECT-NELECT)<1E-10_q) GOTO 110
      IF ( (ABS(EF1-EF2)/(ABS(EFERMI)+1.E-10_q))<1E-14_q) GOTO 120
      IF ( ELECT>NELECT) THEN
        IF (.NOT.LOWB)  EF1=EF1-DELTAE
        HIGHB=.TRUE.
        EF2  =EFERMI
      ELSE
        IF (.NOT.HIGHB) EF2=EF2+DELTAE
        LOWB =.TRUE.
        EF1  =EFERMI
      ENDIF
      ENDDO setfermi

  120 CONTINUE
      IF (IU0>=0) THEN
      WRITE(IU0,*)' WARNING: DENFER: can''t reach specified precision'
      WRITE(IU0,*)' Number of Electrons is NELECT =',ELECT
      ENDIF

  110 CONTINUE
!=======================================================================
! calculate entropy
!=======================================================================
      ECORR= 0

      DO ISP=1,ISPIN
      DO K=1,WDES%NKPTS
      DO N=1,WDES%NBANDS
        IF (FERWE(N,K,ISP) /=0 .AND. FERWE(N,K,ISP) /=1) &
     &  ECORR=ECORR+FERWE(N,K,ISP)    *LOG(FERWE(N,K,ISP))*WDES%WTKPT(K)+ &
     &            (1-FERWE(N,K,ISP))*LOG(1-FERWE(N,K,ISP))*WDES%WTKPT(K)
      ENDDO
      ENDDO
      ENDDO

      ECORR=ECORR*SIGMA_*RSPIN
      CALLMPI( M_sum_d( WDES%COMM_INTER, ECORR, 1))

      RETURN
      END SUBROUTINE DENFER 

      END SUBROUTINE DENSTA
!***********************SUBROUTINE DENTET******************************
!
! subroutine DENTET calculates a continuous density of states in the
! interval (EMIN,EMAX) applying the tetrahedron method to the discrete
! eigenvalue spectrum in CELEN(NBANDS,NKPTS). EFERMI is calculated so
! that sum over occupation-numbers is equal to number of electrons
! executed on all nodes
!**********************************************************************

      SUBROUTINE DENTET(IU0,CELEN,WTKPT,NBANDS,NKPTS,DOS,DOSI,NEDOS, &
     &           ISPIN,RSPIN,EMIN,EMAX,IDTET,NTET,VOLWGT,NELECT,EFERMI,FERWE, &
     &         ECORR,JOB,IU6,PAR,DOSPAR,NKDIM,LDIMP,NIOND,JOBPAR)
      USE prec
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)
      LOGICAL LOWB,HIGHB
      REAL(q)   NELECT
      DIMENSION CELEN(NBANDS,NKDIM,ISPIN),FERWE(NBANDS,NKDIM,ISPIN)
      DIMENSION DOS(NEDOS,ISPIN),DOSI(NEDOS,ISPIN),IDTET(0:4,NTET)
      DIMENSION DOSPAR(NEDOS,LDIMP,NIOND,ISPIN),WTKPT(NKPTS)
      DIMENSION PAR(NBANDS,NKDIM,LDIMP,NIOND,ISPIN)

      DELTAE=(EMAX-EMIN)/(NEDOS-1)
!=======================================================================
! initialize arrays for dos and integr. dos
!=======================================================================

      IF (JOBPAR==1) THEN
       DOSPAR=0
      ENDIF
      DOS =0
      DOSI=0

!=======================================================================
! calculate dos and integrated dos
!=======================================================================

      CALL BZINTS(0,FERWE,CELEN,WTKPT,NBANDS,NBANDS,NKPTS,IDTET,NTET, &
     &     ISPIN,RSPIN,VOLWGT,EMIN,EMAX,DOS,DOSI,NEDOS,EFERMI,SUMWEI, &
     &     SUME,100,PAR,DOSPAR,NKDIM,LDIMP,NIOND,NIOND,JOBPAR)

!=======================================================================
! calculate approximated fermi energy
!=======================================================================

      DO I=1,NEDOS
        DOSTOT=DOSI(I,1)
        IF (ISPIN==2) DOSTOT=DOSTOT+DOSI(I,2)
        IF (ABS(DOSTOT-NELECT)<0.01_q .OR.DOSTOT>NELECT) EXIT
      ENDDO

      EFERMI= EMIN+(I-1)*DELTAE
      IF (JOB==0) RETURN

!=======================================================================
! search now for exact Fermi-level
!=======================================================================
      EF1= EMIN+(I-2)*DELTAE
      LOWB =.FALSE.
      EF2= EMIN+(I-1)*DELTAE
      HIGHB=.FALSE.

!=======================================================================
! calculate fermi-weighting function, and their sum
!=======================================================================
   calcfermi: DO
      EFERMI=(EF1+EF2)/2

      CALL BZINTS(JOB,FERWE,CELEN,WTKPT,NBANDS,NBANDS,NKPTS,IDTET,NTET, &
     &       ISPIN,RSPIN,VOLWGT,EMIN,EMAX,DOS,DOSI,NEDOS,EFERMI,SUMWEI, &
     &       SUME,100,PAR,DOSPAR,NKDIM,LDIMP,NIOND,NIOND,JOBPAR)
      ELECT=SUMWEI*RSPIN

!=======================================================================
! compare now with Number of Electrons
!=======================================================================
      IF ( ABS(ELECT-NELECT)<1E-10_q) GOTO 110
      IF ( (ABS(EF1-EF2)/(ABS(EFERMI)+1.E-10_q))<1E-14_q) GOTO 120
      IF ( ELECT>NELECT) THEN
        IF (.NOT.LOWB)  EF1=EF1-DELTAE
        HIGHB=.TRUE.
        EF2  =EFERMI
      ELSE
        IF (.NOT.HIGHB) EF2=EF2+DELTAE
        LOWB =.TRUE.
        EF1  =EFERMI
      ENDIF
      ENDDO calcfermi

  120 CONTINUE
      IF (IU0>=0) THEN
         WRITE(*,*)' WARNING: DENTET: can''t reach specified precision'
         WRITE(*,*)' Number of Electrons is NELECT =',ELECT
      ENDIF

  110 CONTINUE
! Final call to BZINTS (not really absolutely necessary ...) - can be
! commented out if no informational output/debugging desired ... !
      CALL BZINTS(JOB,FERWE,CELEN,WTKPT,NBANDS,NBANDS,NKPTS,IDTET,NTET, &
     &      ISPIN,RSPIN,VOLWGT,EMIN,EMAX,DOS,DOSI,NEDOS,EFERMI,SUMWEI, &
     &      SUME,IU6,PAR,DOSPAR,NKDIM,LDIMP,NIOND,NIOND,JOBPAR)

!=======================================================================
! How to calculate the correction term to the total Energy???? Is there
! a correction term at all ('no smearing', analytical interpolation and
! then integration with 'delta-function sampling'!!!) ?????????????????
! People say generally: NO! THERE IS NO ENTROPY! --> believe it or not!
!=======================================================================
      ECORR= 0

      RETURN
      END SUBROUTINE DENTET


!******************** DELSTP    ****************************************
!
! Returns generalised delta and step functions (Methfessel & Paxton)
!
!  Input:
!      n > -1 : order of approximant; x : argument
!  Output:
!      D_n (x) ,  S_n (x)
!  Remarks:
!      D_n (x) = exp -x^2 * sum_i=0^n A_i H_2i(x)
!      S_n (x) = (1 - erf x)/2 + exp -x^2 * sum_i=1^n A_i H_{2i-1}(x)
!      where H is a Hermite polynomial and
!      A_i = (-1)^i / ( i! 4^i sqrt(pi) )
!
!***********************************************************************

      SUBROUTINE DELSTP(N,X,D,S)
      USE prec
      USE constant
      IMPLICIT REAL(q) (A-H,O-Z)

      IF (X<-1.E5_q) THEN
         D=0._q
         S=0._q
         RETURN
      END IF
      IF (X>1.E5_q) THEN
         D=0._q
         S=1._q
         RETURN
      END IF
!=======================================================================
!  If n < 0 : assume Gaussian type smearing
!  (must return  same as N=0 or ... )
!=======================================================================
      IF (N<0) THEN
         D=EXP(-(X*X))/SQRT(PI)
         S=0.5_q+0.5_q*ERRF(X)
         RETURN
      END IF
!=======================================================================
! Methfessel & Paxton
!=======================================================================
      EX2=EXP(-(X*X))
      S0=0.5_q*ERRF(X)
      A=1._q/SQRT(PI)
      K=0
      H1=1._q
      H2=2._q*X
      S=0._q
      D=A
      DO I=1,N
         A=A/((-4._q)*I)
         K=K+1
         H3=H1
         H1=H2
         H2=2._q*X*H2-2*K*H3
         S=S+A*H1
         K=K+1
         H3=H1
         H1=H2
         H2=2._q*X*H2-2*K*H3
         D=D+A*H1
      ENDDO
      D=D*EX2
      S=0.5_q+S0-S*EX2
      RETURN
      END SUBROUTINE DELSTP
