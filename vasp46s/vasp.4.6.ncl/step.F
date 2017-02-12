#include "symbol.inc"
!*********************************************************************
! RCS:  $Id: step.F,v 1.1 2000/11/15 08:13:54 kresse Exp $
!
!  subroutine STEP performs one predictor-corrector step  on the
!  position of teh ions
!
!  This part of the program uses slightly different reduced units:
!
!  particel coordinates are in direct lattice   range [0,1]
!  unit time     UT  =  POTIM * 1E-15
!    in reduced units one timestep for the ion movents takes 1 UT
!  unit energy   UE  = 1 eV
!
!  these parameters must be supplied by the main porgram
!  NIONS    number of IONS
!  A        direct lattice
!  TS       soll-temperature
!  POTIM    length of timestep
!  X        positions after the predictor step
!           (on first call initial positions)
!  V        velocities of ions in r.u.
!           these are the velocities after the corrector step
!           (on first call initial velocities)
!  D2C      D2C must contain half the accelations of the ions in r.u.
!           evaluated at the predictor coordinates X,
!           these accelerations serve as
!           corrector-derivatives of  second order
!           (on first call accelerations for initial positions)
!  SMASS    mass Parameter for nose dynamic, it has the dimension
!           of a Energy*time**2 = mass*length**2 and is supplied in
!           atomic mass * ULX**2 (this makes it  easy to define
!           the parameter)
!         0,-1   microcanonical ensemble
!                if ISCALE is 1 scale velocities to TS
!         -2     dont change velocities at all
!
!  execpt for initialisation only D2C must be recalculated by
!  the main program
!
!  Following varibles are used by step, and should not be changed
!  S        nose-parameter
!  D2       predictor-derivatives of  second order
!  D3       predictor-derivatives of   third order
!
!  following Quantities are return
!  XC       'exact positions' after the corrector step (C-step)
!           position are equal to those at the previous timestep
!  X        updated predictor position of ions (see above)
!  EKIN     kinetic energy of system evaluated after C-step
!  EPS      kinetic energy of fictivious Nose-parameter after C-step
!  ES       potential energy of fictivious Nose-parameter after C-step
!  DISMAX   maximum distance the ions have moved in r.u.
!
!  NOTE:
!  if Energy is calculated as EKIN+EPS+ES+Pot. ENERGY at predictor
!  coordinates, then there is allways a small error in the total
!  Energy, but there should be no overall drift !
!
!*********************************************************************

      SUBROUTINE STEP(INIT,ISCALE,NIONS,A,ANORM,D2C,SMASS,X,XC, &
     &      POTIM,POMASS,NTYP,ITYP,TS,V,D2,D3,S, &
     &      EKIN,EPS,ES,DISMAX,NDEGREES_OF_FREEDOM,IU6)
      USE prec

      IMPLICIT REAL(q) (A-H,O-Z)
      DIMENSION  X(3,NIONS),XC(3,NIONS), &
     &           V(3,NIONS),D2(3,NIONS),D3(3,NIONS)
      DIMENSION  ITYP(NIONS),POMASS(NTYP)
      DIMENSION  S(4)
      DIMENSION  A(3,3),ANORM(3)
!=======================================================================
!   Some important Parameters, to convert to UL,UT,UE
!   BOLK is the Bolzmannk.,BOLKEV is the Bolzmannk. in EV
!=======================================================================
      PARAMETER (EVTOJ=1.60217733E-19_q,AMTOKG=1.6605402E-27_q, &
     &           BOLK=1.380658E-23_q,BOLKEV=8.6173857E-5_q)

!=======================================================================
!   D2C ..are the calculated vectors for the accelerations
!=======================================================================
      DIMENSION  D2C(3,NIONS)

!=======================================================================
!   set unit length,unit time ...
!=======================================================================
      UL =1E-10_q*ANORM(1)
      UT =POTIM*1E-15_q
!=======================================================================
!   convert SMASS to more usefull units
!   SMASS is in         (atomic mass* UL**2)
!   what we want now is (eV*UT**2)
!
!=======================================================================
      FACT=(AMTOKG/EVTOJ)*(UL/UT)**2

      IF (SMASS>0) THEN
        SQQ=SMASS *FACT
      ELSE
        SQQ=0
        IF (SMASS==-2) THEN
! SMASS=-2 set D2C and D2 and D3 to 0 (simulate no accelerations at all)
        DO 50 N=1,3
        DO 50 I=1,NIONS
          D2C(N,I)=0
          D2(N,I) =0
          D3(N,I)=0
 50    CONTINUE
      ENDIF
      ENDIF

!=======================================================================
!   INIT=0 initialize variables and perform first predictor step
!   copy D2C to     D2
!        X   to     XC
!   set  D3  to     0
!   also initialize S
!=======================================================================
      IF (INIT==0) THEN
        INIT=1
        DO 100 N=1,3
        DO 100 I=1,NIONS
          D2(N,I)=D2C(N,I)
          D3(N,I)=0
          XC(N,I)=X(N,I)
 100    CONTINUE

!-------use kinetic Energy to calculate initial values for S
        CALL EKINC(EKIN,NIONS,NTYP,ITYP,POMASS,POTIM,A,V)

        S(1)=1._q
        S(2)=0._q
        S(3)=0._q
        IF (SMASS>0) S(3)=(EKIN-NDEGREES_OF_FREEDOM*BOLKEV*TS/2)*S(1)/SQQ
        S(4)=0._q

        GOTO 1000
      ENDIF

!*********************************************************************
!   c o r r e c t o r  -  s t e p
!   calculation of: corrector-derivative of third order (D3X,D3Y,D3Z)
!                   predictor & corrector for velocities (VX,VY,VZ)
!                   corrector for particles coordinates
!*********************************************************************

      IF (SMASS>0) THEN
        SC=(S(2)+2._q*S(3)+3._q*S(4))/(S(1)+S(2)+S(3)+S(4))
      ELSE
        SC=0
      ENDIF

      DO 200 N=1,3
      DO 200 I=1,NIONS

         D2C(N,I)=D2C(N,I) &
     &     -SC*(0.5_q*(V(N,I)+D3(N,I))+D2(N,I)+D3(N,I))
         D3X1=(D2C(N,I)-D2(N,I))/3._q
         V(N,I)=V(N,I)+D2C(N,I)+D2(N,I)-(D3X1-D3(N,I))*0.5_q
         X(N,I)=X(N,I)+(D3X1-D3(N,I))*0.5_q
         D3(N,I)=D3X1

         IF(X(N,I)>=1._q) X(N,I)=X(N,I)-1._q
         IF(X(N,I)<0._q) X(N,I)=X(N,I)+1._q
         D2(N,I)=D2C(N,I)
         XC(N,I)=X(N,I)

  200 CONTINUE

!=======================================================================
!   calculate kintic energy and temperature
!=======================================================================

      CALL EKINC(EKIN,NIONS,NTYP,ITYP,POMASS,POTIM,A,V)

!=======================================================================
!   c o r r e c t o r  -  i t e r a t i o n   o f  S
!   calculation of: corrector of S,S(1)
!                   1.,2.,3.-derivatives of S,S(2),S(3),S(4)
!=======================================================================

      IF (SMASS>0) THEN
      SC=(1._q/6._q)*((S(1)+S(2)+S(3)+S(4)) &
     &  *((EKIN-NDEGREES_OF_FREEDOM*BOLKEV*TS/2)/SQQ+0.5_q*SC**2)-S(3))-0.5_q*S(4)
      S(1)=S(1)+S(2)+S(3)+S(4)+SC
      S(2)=S(2)+2._q*S(3)+3._q*S(4)+5._q*SC
      S(3)=S(3)+3._q*S(4)+6._q*SC
      S(4)=S(4)+2._q*SC
      ENDIF

!=======================================================================
!   scale velocities if requested
!   the 3. derivatives must be scaled for consistency of the
!   system of differential equations.
!=======================================================================

 1000 CONTINUE

      IF (ISCALE==1) THEN
      SCALE=SQRT(TS/(2*EKIN/(BOLKEV*NDEGREES_OF_FREEDOM)))
      IF (IU6>=0) &
      WRITE(IU6,11) SCALE
   11 FORMAT('scaling velocities: factor=',F10.4)

      DO 500 I =1,3
      DO 500 NI=1,NIONS
           V(I,NI) =V(I,NI)*SCALE
           D3(I,NI)=D3(I,NI)*SCALE
  500 CONTINUE

      CALL EKINC(EKIN,NIONS,NTYP,ITYP,POMASS,POTIM,A,V)

      ENDIF
!*********************************************************************
!   p r e d i c t o r  -  s t e p
!   calculation of predictor-coordinates (x,y,z) for all particles
!   calculation of predictor-derivatives of first up to third order
!   is included in the corrector step. this data is not required
!   for calculating the accelerations, and hence the correction terms.
!   in the same way S and its derivatives are calculated.
!*********************************************************************

      EPS=0.5_q*SQQ*(S(2)/S(1))**2
      ES =NDEGREES_OF_FREEDOM*BOLKEV*TS*LOG(S(1))

      DO 110 N=1,3
      DO 110 I=1,NIONS
         X(N,I)=X(N,I)+V(N,I)+D2(N,I)+D3(N,I)
         IF(X(N,I)>=1._q) X(N,I)=X(N,I)-1._q
         IF(X(N,I)<0._q) X(N,I)=X(N,I)+1._q
  110 CONTINUE

!=======================================================================
!   calculate the maximum distance moved by the ions
!=======================================================================

      DISMAX=0
      DO 120 I=1,NIONS
         DISMAX=MAX(DISMAX, &
     &        (V(1,I)+D2(1,I)+D3(1,I))**2+ &
     &        (V(2,I)+D2(2,I)+D3(2,I))**2+ &
     &        (V(3,I)+D2(3,I)+D3(3,I))**2)
  120 CONTINUE

      DISMAX=SQRT(DISMAX)

!=======================================================================
!   now the routine drops back to MAIN
!   MAIN evaluates the accelerations D2C at the predictor
!   coordinates X
!   (these accelerations will serve as corrector-derivatives of second
!   order).
!=======================================================================

      RETURN
      END
