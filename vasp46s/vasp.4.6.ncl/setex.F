#include "symbol.inc"
      MODULE SETEXM
      USE prec
      INCLUDE "setexm.inc"
      ! interpolation of correlation from paramagnetic to 
      ! ferromagnetic case according to 
      ! Vosko, Wilk and Nusair, CAN. J. PHYS. 58, 1200 (1980)
      INTEGER :: LFCI=1

      CONTAINS
!************************ PROGRAM RD_EXCH*******************************
! RCS:  $Id: setex.F,v 1.3 2001/01/19 11:50:40 kresse Exp $
!
!  VASP interpolate the XC-potential from a table (at least
!    the plane wave part)
!  the required table is generated here
!  up to vasp.4.4 it was possible to read a file EXHCAR,
!  this facility has been removed however in this version
!
!***********************************************************************
      SUBROUTINE RD_EX(E,ISPIN,LEXCH,LEXCHF,IU6,IU0,IDIOT)
      USE prec
      USE ini
!      USE pseudo
      IMPLICIT NONE

      INTEGER  LEXCH            ! xc-type
      LOGICAL  LEXCHF
      INTEGER  ISPIN            ! spin

      TYPE (EXCTABLE) E

! arrays for tutor call
      INTEGER IU6,IU0,IDIOT,ITUT(3),RDUM,CDUM,LDUM
! temporary
      INTEGER N,NDUMMY,LOXCH
      REAL(q) AMARG
      CHARACTER*1 CSEL
      CHARACTER*2 CEXCH

      LEXCHF=.FALSE.
      IF (IU6>=0) WRITE(IU6,*) 'EXHCAR: internal setup'

      CALL SET_EX(LEXCH,NEXCH,E%EXCTAB,E%NEXCHF,E%RHOEXC,IU0)

      AMARG=1E30_q  ! natural boundary conditions required
      CALL SPLCOF(E%EXCTAB(1,1,1),E%NEXCHF(2),NEXCH,AMARG)

      IF (ISPIN==2) THEN
         CALL SPLCOF(E%EXCTAB(1,1,2),E%NEXCHF(2),NEXCH,AMARG)
         CALL SPLCOF(E%EXCTAB(1,1,3),E%NEXCHF(2),NEXCH,AMARG)
         CALL SPLCOF(E%EXCTAB(1,1,4),E%NEXCHF(2),NEXCH,AMARG)
         CALL SPLCOF(E%EXCTAB(1,1,5),E%NEXCHF(2),NEXCH,AMARG)
         CALL SPLCOF(E%EXCTAB(1,1,6),E%NEXCHF(2),NEXCH,AMARG)
      ENDIF

 7002 CONTINUE
      IF (IU6>=0) WRITE(IU6,7004)LEXCH,E%RHOEXC(1),E%NEXCHF(1), &
     &                     E%RHOEXC(2),E%NEXCHF(2)
 7004 FORMAT(' exchange correlation table for  LEXCH = ',I8/ &
     &       '   RHO(1)= ',F8.3,5X,'  N(1)  = ',I8/ &
     &       '   RHO(2)= ',F8.3,5X,'  N(2)  = ',I8)

      RETURN
      END SUBROUTINE


!************************ PROGRAM SETEX  *******************************
!
! set up the default xc-table
! i.e.
!  Ceperly Alder
!  with standard interpolation to spin
!  with relativistic correction
!
!***********************************************************************

      SUBROUTINE SET_EX(LEXCHV,N,EXCTAB,NEXCHF,RHOEXC,IU0)
      USE prec

      IMPLICIT REAL(q) (A-H,O-Z)
      DIMENSION EXCTAB(N,5,6),NEXCHF(2),RHOEXC(2)
      LOGICAL TREL
      CHARACTER*13 CEXCH

      TREL=.TRUE.
! convert LEXCHV to local format
      IF (LEXCHV==0) THEN
        LEXCH=0
      ELSE IF (LEXCHV==1) THEN
    ! Hedin Lundquist
        LEXCH=4
      ELSE IF (LEXCHV==2) THEN
    ! Ceperly-Alder
        LEXCH=1
      ELSE IF (LEXCHV==3) THEN
    ! Wigner
        LEXCH=6
      ELSE IF (LEXCHV==8 .OR. LEXCHV==9) THEN
    ! Pade approximation of Perdew
        LEXCH=7
      ELSE
        LEXCH=1
      ENDIF

      IF (LEXCH==0) THEN
         CEXCH='  '
      ELSE IF (LEXCH==1) THEN
         CEXCH='Ceperly-Alder'
      ELSE IF (LEXCH==2) THEN
         CEXCH='VW'
      ELSE IF (LEXCH==3) THEN
         CEXCH='GL'
      ELSE IF (LEXCH==4) THEN
         CEXCH='HL'
      ELSE IF (LEXCH==5) THEN
         CEXCH='BH'
      ELSE IF (LEXCH==6) THEN
         CEXCH='WI'
      ELSE IF (LEXCH==7) THEN
         CEXCH='PB'
      ELSE
         IF (IU0>=0) &
         WRITE(IU0,*) 'Wrong exchange correlation type!'
         STOP
      ENDIF
      

      IF (IU0>=0) THEN
         IF (LEXCH==7) THEN
            WRITE(IU0,*)'LDA part: xc-table for Pade appr. of Perdew'
         ELSE
            IF (LFCI==1) THEN 
               WRITE(IU0,*)'LDA part: xc-table for ',CEXCH, ', Vosko type interpolation para-ferro'
            ELSE
               WRITE(IU0,*)'LDA part: xc-table for ',CEXCH, ', standard interpolation'
            ENDIF
         ENDIF
      ENDIF

! Slater parameter
      SLATER=1._q
! standard interpolation   from para- to ferromagnetic corr
      RHOSMA=INT(0.5_q*1000)/1000._q
      NSMA  =2000
#ifdef hugeXCtable
      RHOMAX=INT(2000.5_q*1000)/1000._q
      WRITE(0,*) '******* RHOMAX set to 2000 instead of 100 *****'
#else
      RHOMAX=INT(100.5_q*1000)/1000._q
#endif


      RHOEXC(1)=RHOSMA
      RHOEXC(2)=RHOMAX
      NEXCHF(1)=NSMA
      NEXCHF(2)=N

      IF (NSMA/=0) THEN
         RH=RHOSMA/NSMA/2
      ELSE
         RH=RHOMAX/N/100
      ENDIF
      J=1
      CALL EXCHG(LEXCH,RH,EXCP, &
     &           DEXF,DECF,ALPHA,LFCI,SLATER,TREL)
        EXCTAB(J,1,1)=RH
        EXCTAB(J,1,2)=RH
        EXCTAB(J,1,3)=RH
        EXCTAB(J,1,4)=RH

        EXCTAB(J,2,1)=EXCP
        EXCTAB(J,2,2)=DEXF
        EXCTAB(J,2,3)=DECF
        EXCTAB(J,2,4)=ALPHA


      IF (NSMA/=0) THEN
         DRHO=RHOSMA/NSMA
         DO 100 I=1,NSMA-1
            J=I+1
            RH=DRHO*I
            CALL EXCHG(LEXCH,RH,EXCP, &
     &                 DEXF,DECF,ALPHA,LFCI,SLATER,TREL)
          EXCTAB(J,1,1)=RH
          EXCTAB(J,1,2)=RH
          EXCTAB(J,1,3)=RH
          EXCTAB(J,1,4)=RH

          EXCTAB(J,2,1)=EXCP
          EXCTAB(J,2,2)=DEXF
          EXCTAB(J,2,3)=DECF
          EXCTAB(J,2,4)=ALPHA

  100    CONTINUE

         J=NSMA+1
         RH=RHOSMA
         CALL EXCHG(LEXCH,RH,EXCP, &
     &              DEXF,DECF,ALPHA,LFCI,SLATER,TREL)
         EXCTAB(J,1,1)=RH
         EXCTAB(J,1,2)=RH
         EXCTAB(J,1,3)=RH
         EXCTAB(J,1,4)=RH

         EXCTAB(J,2,1)=EXCP
         EXCTAB(J,2,2)=DEXF
         EXCTAB(J,2,3)=DECF
         EXCTAB(J,2,4)=ALPHA
      ENDIF

      DRHO=(RHOMAX-RHOSMA)/(N-NSMA)
      DO 200 I=1,N-NSMA-1
         J=I+NSMA+1
         RH=DRHO*I+RHOSMA
         CALL EXCHG(LEXCH,RH,EXCP, &
     &              DEXF,DECF,ALPHA,LFCI,SLATER,TREL)
         EXCTAB(J,1,1)=RH
         EXCTAB(J,1,2)=RH
         EXCTAB(J,1,3)=RH
         EXCTAB(J,1,4)=RH

         EXCTAB(J,2,1)=EXCP
         EXCTAB(J,2,2)=DEXF
         EXCTAB(J,2,3)=DECF
         EXCTAB(J,2,4)=ALPHA
  200 CONTINUE

      J=1
      ZETA=0
      FZA =.854960467080682810_q
      FZB =.854960467080682810_q
      EXCTAB(J,1,5)=ZETA
      EXCTAB(J,1,6)=ZETA
      EXCTAB(J,2,5)=FZA
      EXCTAB(J,2,6)=FZB

      DO 800 I=1,N-1
         ZETA=FLOAT(I)/FLOAT(N-1)
         FZA=FZ0(ZETA)/ZETA/ZETA
         FZB=FZ0(ZETA)/ZETA/ZETA
         J=I+1
         EXCTAB(J,1,5)=ZETA
         EXCTAB(J,1,6)=ZETA
         EXCTAB(J,2,5)=FZA
         EXCTAB(J,2,6)=FZB
  800 CONTINUE

      IF (.FALSE.) THEN
      DO 300 I=1,N
         WRITE(97,20) EXCTAB(I,1,1),EXCTAB(I,2,1)
  300 CONTINUE
      DO 400 I=1,N
         WRITE(97,20)  EXCTAB(I,1,2),EXCTAB(I,2,2)
  400 CONTINUE
      DO 500 I=1,N
         WRITE(97,20)  EXCTAB(I,1,3),EXCTAB(I,2,3)
  500 CONTINUE
      DO 600 I=1,N
         WRITE(97,20)  EXCTAB(I,1,4),EXCTAB(I,2,4)
  600 CONTINUE
      DO 700 I=1,N
         WRITE(97,20)  EXCTAB(I,1,5),EXCTAB(I,2,5)
  700 CONTINUE
      DO 710 I=1,N
         WRITE(97,20)  EXCTAB(I,1,6),EXCTAB(I,2,6)
  710 CONTINUE
   20 FORMAT((3(E24.16,2X)))
      ENDIF

      RETURN
      END SUBROUTINE


!**************** SUBROUTINE EXCHG *************************************
!   EXCHG calculated xc-energy
!   uses xclib
!***********************************************************************

      SUBROUTINE EXCHG(LEXCH,RHO,EXCP,DEXF,DECF,ALPH,LFCI,SLATER,TREL)
      USE prec
      USE constant

      IMPLICIT REAL(q) (A-H,O-Z)
      LOGICAL TREL
      INTEGER LFCI

      IF (RHO==0) THEN
         EXCP=0._q
         DEXF=0._q
         DECF=0._q
         ALPH=0._q
         RETURN
      ENDIF

      RHOTHD = RHO**(1/3._q)
      RS = (3._q/(4._q*PI)/RHO)**(1/3._q) /AUTOA

      IF (LEXCH==0) THEN
         EXCP=SLATER*EX(RS,1,TREL)
         DEXF=SLATER*EX(RS,2,TREL)-EXCP
         DECF=0._q
         ALPH=0._q
      ELSE IF (LEXCH==1) THEN
         EXCP=EX(RS,1,TREL)+ECCA(RS,1)
         DEXF=EX(RS,2,TREL)-EX(RS,1,TREL)
         DECF=ECCA(RS,2)-ECCA(RS,1)
         ALPH=0._q
      ELSE IF (LEXCH==2) THEN
         EXCP=EX(RS,1,TREL)+ECVO(RS,1)
         DEXF=EX(RS,2,TREL)-EX(RS,1,TREL)
         DECF=ECVO(RS,2)-ECVO(RS,1)
         ALPH=0._q
      ELSE IF (LEXCH==3) THEN
         EXCP=EX(RS,1,TREL)+ECGL(RS,1)
         DEXF=EX(RS,2,TREL)-EX(RS,1,TREL)
         DECF=ECGL(RS,2)-ECGL(RS,1)
         ALPH=0._q
      ELSE IF (LEXCH==4) THEN
         EXCP=EX(RS,1,TREL)+ECHL(RS,1)
         DEXF=EX(RS,2,TREL)-EX(RS,1,TREL)
         DECF=ECHL(RS,2)-ECHL(RS,1)
         ALPH=0._q
      ELSE IF (LEXCH==5) THEN
         EXCP=EX(RS,1,TREL)+ECBH(RS,1)
         DEXF=EX(RS,2,TREL)-EX(RS,1,TREL)
         DECF=ECBH(RS,2)-ECBH(RS,1)
         ALPH=0._q
      ELSE IF (LEXCH==6) THEN
         EXCP=EX(RS,1,TREL)+ECWI(RS,1)
         DEXF=EX(RS,2,TREL)-EX(RS,1,TREL)
         DECF=ECWI(RS,2)-ECWI(RS,1)
         ALPH=0._q
      ELSE IF (LEXCH==7) THEN
         ZETA=0  ! paramagnetic result
         CALL CORPBE_LDA(RS,ZETA,ECLDA,ECD1LDA,ECD2LDA)
         ZETA=1  ! ferromagnetic result
         CALL CORPBE_LDA(RS,ZETA,ECLDA_MAG,ECD1LDA,ECD2LDA)

         EXCP=EX(RS,1,.FALSE.)+ECLDA
         DEXF=EX(RS,2,.FALSE.)-EX(RS,1,.FALSE.)
         DECF=ECLDA_MAG-ECLDA
!        WRITE(77,'(5F14.7)') RS, ECCA(RS,1),ECLDA,ECCA(RS,2),ECLDA_MAG
         ALPH=0._q
      ELSE
         EXCP=0._q
         DEXF=0._q
         DECF=0._q
         ALPH=0._q
      ENDIF
      !
      ! for Perdew, Burke Ernzerhof we allways use the 
      ! recommended  interpolation from nm to magnetic
      ! 
      IF (LFCI==1 .OR. LEXCH==7 ) THEN
         IF (LEXCH==7) THEN
            A0=PBE_ALPHA(RS)
!           WRITE(78,'(5F14.7)') RS,A0,ALPHA0(RS)
         ELSE
            A0=ALPHA0(RS)
         ENDIF
         ALPH=DECF-A0
         DECF=A0
      ENDIF

      EXCP=EXCP*RYTOEV/RHOTHD
      DEXF=DEXF*RYTOEV/RHOTHD
      DECF=DECF*RYTOEV/RHOTHD
      ALPH=ALPH*RYTOEV/RHOTHD

      RETURN
      END SUBROUTINE
      END MODULE
