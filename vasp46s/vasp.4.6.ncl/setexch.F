#include "symbol.inc"
!************************ PROGRAM FEXCSP *******************************
! RCS:  $Id: setexch.F,v 1.1 2000/11/15 08:13:54 kresse Exp $
!
! this subroutine set up file EXHCAR, which are necessary to
! calculate the exchange-correlation energy in VAMP ...
! Warning: a factor rho^(1/3) is not contained in the tables
!
!***********************************************************************

      PROGRAM FEXCSP

      PARAMETER (N=4000)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      PARAMETER (AUTOA=0.529177249_q,RYTOEV=13.605826_q, &
     &           PI=3.141_q 592 653 589 793 238)
      CHARACTER*2 CEXCH
      LOGICAL TREL
      DIMENSION RH(N),EXCP(N),DEXF(N),DECF(N),FINTX(N),FINTC(N),ALPHA(N)

!=======================================================================
! ask the user for the type of exchange correlation he wants to use
!=======================================================================
      WRITE(*,10)
   10 FORMAT(' Exchange-Correlation Type:'/ &
     & ' 0   Slater exchange correlation (usual setting is SLATER=2/3)'/ &
     & ' 1   Perdew and Zunger, PHYS. REV. B23, 5048 (1982)'/ &
     & ' 2   Vosko, Wilk and Nusair, CAN. J. PHYS. 58, 1200 (1980)'/ &
     & ' 3   Gunnarson and Lundqvist'/ &
     & ' 4   Hedin and Lundqvist, J. PHYS. C4, 2064 (1971)'/ &
     & ' 5   Barth and Hedin'/ &
     & ' 6   Wigner-interpolation')
      READ(*,*) LEXCH

      WRITE(*,'(A)') ' Relativistic corrections? (.T. of .F.)'
      READ(*,*) TREL

      WRITE(*,11)
   11 FORMAT(' Interpolation type from para- to ferromagnetic corr.'/ &
     & ' 0   exchange-like ''standard interpolation'''/ &
     & ' 1   Vosko-type function (CAN. J. PHYS. 58, 1200 (1980)')
      READ(*,*) LFCI

      SLATER=1._q
      IF (LEXCH==0) THEN
         WRITE(*,*)'Slater-parameter SLATER?'
         READ(*,*) SLATER
         SLATER=SLATER*3._q/2._q
      ENDIF

      WRITE(*,*)'maximal small electron density RHO(1) ?'
      READ(*,*) RHOSMA
      RHOSMA=INT(RHOSMA*1000)/1000._q
      WRITE(*,*)'number of points N(1) between 0 and RHO(1) ?'
      READ(*,*) NSMA
      WRITE(*,*)'maximal electron density RHO(2) ?'
      READ(*,*) RHOMAX
      RHOMAX=INT(RHOMAX*1000)/1000._q

!=======================================================================
! write out exchange-correlation table
!=======================================================================

      OPEN(10,FILE='EXHCAR')

      IF (LEXCH==0) THEN
         CEXCH='  '
      ELSE IF (LEXCH==1) THEN
         CEXCH='CA'
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
      ELSE
         WRITE(*,*) 'Wrong exchange correlation type!'
         STOP
      ENDIF

      WRITE(10,30) CEXCH,RHOSMA,NSMA,RHOMAX,N
   30 FORMAT('exchange1 correlation table'/ &
     &       '  LEXCH = ',A2/ &
     &       '  RHO(1)= ',F8.3/ &
     &       '  N(1)  = ',I8/ &
     &       '  RHO(2)= ',F8.3/ &
     &       '  N(2)  = ',I8)

      IF (NSMA/=0) THEN
         RH(1)=RHOSMA/NSMA/2
      ELSE
         RH(1)=RHOMAX/N/100
      ENDIF
      CALL EXCHG(LEXCH,RH(1),EXCP(1), &
     &           DEXF(1),DECF(1),ALPHA(1),LFCI,SLATER,TREL)
   20 FORMAT((3(E24.16,2X)))

      IF (NSMA/=0) THEN
         DRHO=RHOSMA/NSMA
         DO 100 I=1,NSMA-1
            J=I+1
            RH(J)=DRHO*I
            CALL EXCHG(LEXCH,RH(J),EXCP(J), &
     &                 DEXF(J),DECF(J),ALPHA(J),LFCI,SLATER,TREL)
  100    CONTINUE

         J=NSMA+1
         RH(J)=RHOSMA
         CALL EXCHG(LEXCH,RH(J),EXCP(J), &
     &              DEXF(J),DECF(J),ALPHA(J),LFCI,SLATER,TREL)
      ENDIF

      DRHO=(RHOMAX-RHOSMA)/(N-NSMA)
      DO 200 I=1,N-NSMA-1
         J=I+NSMA+1
         RH(J)=DRHO*I+RHOSMA
         CALL EXCHG(LEXCH,RH(J),EXCP(J), &
     &              DEXF(J),DECF(J),ALPHA(J),LFCI,SLATER,TREL)
  200 CONTINUE

      DO 300 I=1,N
         WRITE(10,20) RH(I),EXCP(I)
  300 CONTINUE
      DO 400 I=1,N
         WRITE(10,20) RH(I),DEXF(I)
  400 CONTINUE
      DO 500 I=1,N
         WRITE(10,20) RH(I),DECF(I)
  500 CONTINUE
      DO 600 I=1,N
         WRITE(10,20) RH(I),ALPHA(I)
  600 CONTINUE
      WRITE(10,20) 0._q,.854960467080682810_q
      DO 700 I=1,N-1
         ZETA=FLOAT(I)/FLOAT(N-1)
         WRITE(10,20) ZETA,FZ0(ZETA)/ZETA/ZETA
  700 CONTINUE
      WRITE(10,20) 0._q,.854960467080682810_q
      DO 800 I=1,N-1
         ZETA=FLOAT(I)/FLOAT(N-1)
         WRITE(10,20) ZETA,FZ0(ZETA)/ZETA/ZETA
  800 CONTINUE

      CLOSE(10)

      END


!**************** SUBROUTINE EXCHG *************************************
!
!     EXCHG calculated xc-energy
!     uses xclib
!***********************************************************************

      SUBROUTINE EXCHG(LEXCH,RHO,EXCP,DEXF,DECF,ALPH,LFCI,SLATER,TREL)
      USE prec

      IMPLICIT REAL(q) (A-H,O-Z)
      PARAMETER (AUTOA=0.529177249_q,RYTOEV=13.605826_q, &
     &           PI=3.141_q 592 653 589 793 238)

      LOGICAL TREL

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
      ELSE
         EXCP=0._q
         DEXF=0._q
         DECF=0._q
         ALPH=0._q
      ENDIF

      IF (LFCI==1) THEN
         A0=ALPHA0(RS)
!elio Moroni
!         ALPH=DECF+DEXF-A0
         ALPH=DECF-A0
         DECF=A0
      ENDIF

      EXCP=EXCP*RYTOEV/RHOTHD
      DEXF=DEXF*RYTOEV/RHOTHD
      DECF=DECF*RYTOEV/RHOTHD
      ALPH=ALPH*RYTOEV/RHOTHD

      RETURN
      END
