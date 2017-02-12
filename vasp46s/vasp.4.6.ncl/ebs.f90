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





      MODULE ebs
      USE prec
      CONTAINS
!************************ subroutine FEWALD ****************************
! RCS:  $Id: ebs.F,v 1.2 2002/04/30 15:36:32 kresse Exp $
!
! rewritten by Georg  Kresse (original version Mike Payne?)
!
! This subroutine calculates the ewald energy due to coulomb energy
! between the ions and the neutralising background.
! Arbitrary cellshape and arbitrary atomic species are allowed.
! In addition the subroutine calculates forces and stress corresponding
! to the ewald energy.
!
! rather than computing the values of the ewald integrals this program
! interpolates the values of the integrals from a spline fit of
! the follwing two functions:
!   EWRSPL= f1(x) = Sqrt(Pi)/2 Erfc(x)/x
!   EWQSPL= f2(x) = Exp (-x^2) / 2 x^2)
!
! important parameters for these tables:
!   ARGMR  maximum value of real space table  EWRSPL
!   ARGMQ  maximum value for the reciprocal table EWQSPL
!   NEWPTS number of values in both arrays
!
! some temporary work arrays are used:
!
! RFORCE(:,A,B) = real space part of
!        force on ion A from periodic array of ions B
! GFORCE(:,A,B) =  reciprocal space parts same as above
! RSIF1,2,3 = force on unit cell from the real space parts
! GSIF1,2,3 = force on unit cell from the reciprocal space parts
!
! MIND: that this routine assumes that all positions are within
! [0,1], if this is not the case the real space summation might not
! be done for sufficient cells
!
!***********************************************************************
      SUBROUTINE FEWALD &
     &    (POSION,EWIFOR,A,B,ANORM,BNORM, &
     &     OMEGA,EWSIF,TEWEN,NTYP,ZVAL,NIONS,NIOND,ITYP,NITYP,IU)
      USE prec
      USE ini
      USE constant
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      PARAMETER(NEWPTS=4000,ARGMR=4.0_q,ARGMQ=4.0_q)

      DIMENSION POSION(3,NIONS)
      DIMENSION EWIFOR(3,NIONS)
! old tables supplied by main program (no longer used)
      DATA INIT  /0/
      SAVE INIT,EWRSPL,EWQSPL,ARGMIR,ARGMIQ
      DIMENSION EWRSPL(NEWPTS,5),EWQSPL(NEWPTS,5)

      DIMENSION EWSIF(3,3)
      DIMENSION NITYP(NTYP),ZVAL(NTYP),ITYP(NIONS)

! work arrys
      DIMENSION RFORCE(3,NIONS,NIONS),GFORCE(3,NIONS,NIONS)
      DIMENSION A(3,3),B(3,3),ANORM(3),BNORM(3)
      DIMENSION CPHASE(NIONS)
!=======================================================================
! set up splinetable
! first position in array corresponds to  ARGMR/NEWPTS
!=======================================================================
      IF (INIT==0) THEN
        DO N=1,NEWPTS
          EWRSPL(N,1) =  ARGMR*N/NEWPTS
          EWRSPL(N,2) =  ER(EWRSPL(N,1))

          EWQSPL(N,1) =  ARGMQ*N/NEWPTS
          EWQSPL(N,2) =  EX(EWQSPL(N,1))
       ENDDO

        X =   ARGMR/NEWPTS
        CALL SPLCOF(EWRSPL,NEWPTS,NEWPTS, &
     &   -1/(EXP(X*X)*X) - SQRT(PI)*ERRFC(X)/(2*X*X))
        ARGMIR=ARGMR*40/NEWPTS

        X =   ARGMQ/NEWPTS
        CALL SPLCOF(EWQSPL,NEWPTS,NEWPTS, &
     &   -EXP(-X*X)/X*(1+1/(X*X)) )
        ARGMIQ=ARGMQ*40/NEWPTS

      ENDIF

!=======================================================================
! initialise the values of the work arrays
!=======================================================================
      GFORCE=0.0_q
      RFORCE=0.0_q
      EWIFOR=0.0_q

      RSIF1 =0.0_q
      RSIF2 =0.0_q
      RSIF3 =0.0_q
      RSIF12=0.0_q
      RSIF23=0.0_q
      RSIF31=0.0_q

      GSIF1 =0.0_q
      GSIF2 =0.0_q
      GSIF3 =0.0_q
      GSIF12=0.0_q
      GSIF23=0.0_q
      GSIF31=0.0_q

      TEWEN=0.0_q
      REWEN=0.0_q
      GEWEN=0.0_q
!=======================================================================
! Choose the variable SCALE that determines the rate of decay of the
! real and reciprocal parts of the ewald sums. The choice of SCALE used
! below is good for systems where the lengths of the lattice vectors
! are similar, for odder shape systems a different choice for the value
! of SCALE may reduce the amount of work involved in calculating the
! coulomb lattice sums.
!=======================================================================
      SIZMIN=OMEGA**(1/3._q)
      SCALE=SQRT(PI)/SIZMIN
      PIDSCA=PI/SCALE
!=======================================================================
! MAXC1,2,3 defines the maximum number of cells over which the
!  real space sum is done
! MAXC1,2,3 is chose so that all ions within the sphere
!  ARGMR/SCALE are summed
!=======================================================================
      RCUT=ARGMR/SCALE
      MAXC1=RCUT*BNORM(1)+.99_q
      MAXC2=RCUT*BNORM(2)+.99_q
      MAXC3=RCUT*BNORM(3)+.99_q
!-----------------------------------------------------------------------
! MAXGP1,2,3 defines the maximum number of cells over which the
!  reciprocal space sum is done
! MAXGP1,2,3 is chose so that all contributions within the sphere
!  ARGMQ*PIDSCA are summed
!-----------------------------------------------------------------------
      QCUT=ARGMQ/PIDSCA
      MAXGP1=QCUT*ANORM(1)+.99_q
      MAXGP2=QCUT*ANORM(2)+.99_q
      MAXGP3=QCUT*ANORM(3)+.99_q
      IF (INIT==0) THEN
        IF (IU>=0) WRITE(IU,5) SCALE,MAXC1,MAXC2,MAXC3,MAXGP1,MAXGP2,MAXGP3
   5    FORMAT(' First call to EWALD:  gamma=',F8.3/ &
     &         ' Maximum number of real-space cells',I2,'x',I2,'x',I2/ &
     &         ' Maximum number of reciprocal cells',I2,'x',I2,'x',I2/)

        INIT=1
      ENDIF
!=======================================================================
!
!                      ***********
! caclulate the real space part
!                      ***********
!
!=======================================================================
      CALL VTIME(TV0,TC0)
!-----------------------------------------------------------------------
! scaling factors
!-----------------------------------------------------------------------
      RSCALE=(SCALE**3)*4*FELECT/SQRT(PI)
      RENSCA=SCALE*2*FELECT/SQRT(PI)
!=======================================================================
! calculate the real space contribution to the total energy etc by
! summing the contributions due to ion NI interacting with all the
! ions NNI in the block of unit cells retained in the real space sum
!=======================================================================
dircell: DO N1=-MAXC1,MAXC1
      DO N2=-MAXC2,MAXC2
      DO N3=-MAXC3,MAXC3

      ion1: DO NI=1,NIONS
      ion2: DO NNI=NI,NIONS
      ZZ=ZVAL(ITYP(NI))*ZVAL(ITYP(NNI))
!=======================================================================
! the energy of interaction between each pair of ions must be divided
! equally between the two ions. Since the loops over the ions give the
! energy of interation of ion NI with the periodic array of ions NNI
! where the index NNI is greater than or equal to NNI we include half
! the energy if NI is equal to NNI
!=======================================================================
      IF(NNI==NI) THEN
         DCOUNT=0.5_q
      ELSE
         DCOUNT=1.0_q
      ENDIF
!-----------------------------------------------------------------------
! initialise the distances in the X,Y,Z directions between ion NI and
! ion NNI
! X1,2,3 direct lattice X,Y,Z DIFF cartesian coordinates
!-----------------------------------------------------------------------
      X1=MOD(POSION(1,NI)-POSION(1,NNI)+100.5_q,1._q)-0.5_q+N1
      X2=MOD(POSION(2,NI)-POSION(2,NNI)+100.5_q,1._q)-0.5_q+N2
      X3=MOD(POSION(3,NI)-POSION(3,NNI)+100.5_q,1._q)-0.5_q+N3
      XDIFF= X1*A(1,1)+X2*A(1,2)+X3*A(1,3)
      YDIFF= X1*A(2,1)+X2*A(2,2)+X3*A(2,3)
      ZDIFF= X1*A(3,1)+X2*A(3,2)+X3*A(3,3)

      DIS=SQRT(XDIFF**2+YDIFF**2+ZDIFF**2)
      ARG=DIS*SCALE
!-----------------------------------------------------------------------
! convert distance to an address in the table
! and evaluate f1(x) and g1(x) using spline table
!-----------------------------------------------------------------------
      IF( ARG<ARGMR .AND. ABS(ARG)>1E-10_q ) THEN
      IF (ARG>ARGMIR) THEN
         I  =INT(ARG*NEWPTS/ARGMR)
         REM=ARG-EWRSPL(I,1)

         EW    = EWRSPL(I,2)+REM*(EWRSPL(I,3)+ &
              REM*(EWRSPL(I,4)+REM*EWRSPL(I,5)))
         EW    = EW*ZZ
         EWD   = EWRSPL(I,3)+REM*(EWRSPL(I,4)*2+REM*EWRSPL(I,5)*3)
         EWD   = -EWD/2/ARG*ZZ
      ELSE
         EW    =  ZZ*ER(ARG)
         EWD   = -ERD(ARG)/2/ARG*ZZ
      ENDIF
!-----------------------------------------------------------------------
! add the contribution to the total force on ion ni due to the periodic
! array of ions nni from ion ni in the present real space cell
!-----------------------------------------------------------------------
      RFORCE(1,NNI,NI)=RFORCE(1,NNI,NI)+(XDIFF*EWD)
      RFORCE(2,NNI,NI)=RFORCE(2,NNI,NI)+(YDIFF*EWD)
      RFORCE(3,NNI,NI)=RFORCE(3,NNI,NI)+(ZDIFF*EWD)
!-----------------------------------------------------------------------
! add contribution to  total force on the unit cell
!-----------------------------------------------------------------------
      EWD= EWD*DCOUNT
      RSIF1 =RSIF1 +(EWD*XDIFF)*XDIFF
      RSIF2 =RSIF2 +(EWD*YDIFF)*YDIFF
      RSIF3 =RSIF3 +(EWD*ZDIFF)*ZDIFF
      RSIF12=RSIF12+(EWD*XDIFF)*YDIFF
      RSIF23=RSIF23+(EWD*YDIFF)*ZDIFF
      RSIF31=RSIF31+(EWD*ZDIFF)*XDIFF
!-----------------------------------------------------------------------
! add contribution to total energy
!-----------------------------------------------------------------------
      REWEN=REWEN+EW*DCOUNT
      ENDIF
!=======================================================================
! + end of the loop over ion NNI and NI move onto the next ion
! + move onto ion NNI in the next real space cell in the x,y,z direction
!=======================================================================
      ENDDO ion2
      ENDDO ion1

      ENDDO; ENDDO; ENDDO dircell

!      CALL VTIME(TV,TC)
!      WRITE(0,6)'R-EWALD',TV-TV0,TC-TC0
!6     FORMAT(2X,A8,':  VPU time',F7.2,': CPU time',F7.2)
!      CALL VTIME(TV0,TC0)
!=======================================================================
!
!                      ***********
! now move onto the reciprocal space part of the ewald summations
!                      ***********
!
!=======================================================================
!-----------------------------------------------------------------------
! calculate the scaling factors
!-----------------------------------------------------------------------
      GSISC1=-4*(PI**3)*FELECT/(OMEGA*SCALE**4)
      GSISC2= 2*PI*FELECT/(OMEGA*SCALE**2)

      GSCALE=((2*PI)**2)*FELECT/(OMEGA*SCALE**2)
      GENSCA=2*PI*FELECT/(OMEGA*SCALE**2)
!-----------------------------------------------------------------------
! Loop over all reciprocal cells
!-----------------------------------------------------------------------
reccell: DO N1=-MAXGP1,MAXGP1
      DO N2=-MAXGP2,MAXGP2
      DO N3=-MAXGP3,MAXGP3
!-----------------------------------------------------------------------
! calculate the dimensionless length of the reciprocal lattice vector
!-----------------------------------------------------------------------
      GX= N1*B(1,1)+N2*B(1,2)+N3*B(1,3)
      GY= N1*B(2,1)+N2*B(2,2)+N3*B(2,3)
      GZ= N1*B(3,1)+N2*B(3,2)+N3*B(3,3)

      GARG=SQRT(GX**2+GY**2+GZ**2)*PIDSCA
!-----------------------------------------------------------------------
! convert distance to an address in the ewald integral
! and evaluate f2(x) and g2(x) using spline table
!-----------------------------------------------------------------------
      IF(GARG<ARGMQ.AND.ABS(GARG)>1E-10_q) THEN

      IF (GARG>ARGMIQ) THEN
         I  =INT(GARG*NEWPTS/ARGMQ)
         REM=GARG-EWQSPL(I,1)

         EW    = EWQSPL(I,2)+REM*(EWQSPL(I,3)+ &
              REM*(EWQSPL(I,4)+REM*EWQSPL(I,5)))
         EWD   = EWQSPL(I,3)+REM*(EWQSPL(I,4)*2+REM*EWQSPL(I,5)*3)
         EWD   = -EWD/2/GARG
      ELSE
         EW    = EX(GARG)
         EWD   = -EXD(GARG)/2/GARG
      ENDIF
!-----------------------------------------------------------------------
! using the values of the ewald integrals for this reciprocal lattice
! vector calculate the functions of the reciprocal lattice vector and
! the ewald integrals used to calculate the forces etc
!-----------------------------------------------------------------------
      FORCG1=GX*EW
      FORCG2=GY*EW
      FORCG3=GZ*EW
      GPRO1= GX*B(1,1)+GY*B(2,1)+GZ*B(3,1)
      GPRO2= GX*B(1,2)+GY*B(2,2)+GZ*B(3,2)
      GPRO3= GX*B(1,3)+GY*B(2,3)+GZ*B(3,3)

      TSI1 =EWD*GX*GX*GSISC1+EW*GSISC2
      TSI2 =EWD*GY*GY*GSISC1+EW*GSISC2
      TSI3 =EWD*GZ*GZ*GSISC1+EW*GSISC2
      TSI12=EWD*GX*GY*GSISC1
      TSI23=EWD*GY*GZ*GSISC1
      TSI31=EWD*GZ*GX*GSISC1

!-----------------------------------------------------------------------
! calculate phase factor for corrent G vector
!-----------------------------------------------------------------------
      DO NI=1,NIONS
         CPHASE(NI)=EXP(-CITPI*( &
     &            POSION(1,NI)*N1+POSION(2,NI)*N2+POSION(3,NI)*N3))
      END DO

!-----------------------------------------------------------------------
! loop over all pair of ions to get structure factor and forces
!-----------------------------------------------------------------------
      SFACT=0
      ion3: DO NI=1,NIONS
      ion4: DO NNI=NI,NIONS
      ZZ=ZVAL(ITYP(NI))*ZVAL(ITYP(NNI))

      ! as before, dcount accounts for the energy of interaction between pairs
      ! of ions for the restricted sum nni greater or equal to ni

      IF(NI==NNI) THEN
        DCOUNT=0.5_q
      ELSE
        DCOUNT=1.0_q
      ENDIF
      ! phase factor
      CEXPHF=CPHASE(NI)*CONJG(CPHASE(NNI))

      ! calculate the contribution to the force between ion NI and ions NNI
      GFORCE(1,NNI,NI)=GFORCE(1,NNI,NI)+FORCG1*(AIMAG(CEXPHF)*ZZ)
      GFORCE(2,NNI,NI)=GFORCE(2,NNI,NI)+FORCG2*(AIMAG(CEXPHF)*ZZ)
      GFORCE(3,NNI,NI)=GFORCE(3,NNI,NI)+FORCG3*(AIMAG(CEXPHF)*ZZ)

      ! sum structure factor
      SFACT=SFACT+CEXPHF*DCOUNT*ZZ
      ENDDO ion4
      ENDDO ion3
!-----------------------------------------------------------------------
! add to total energy and stress
!-----------------------------------------------------------------------
      GEWEN =GEWEN +EW*SFACT
      GSIF1 =GSIF1 +TSI1 *SFACT
      GSIF2 =GSIF2 +TSI2 *SFACT
      GSIF3 =GSIF3 +TSI3 *SFACT
      GSIF12=GSIF12+TSI12*SFACT
      GSIF23=GSIF23+TSI23*SFACT
      GSIF31=GSIF31+TSI31*SFACT

      ENDIF
      ENDDO; ENDDO; ENDDO reccell

!      CALL VTIME(TV,TC)
!      WRITE(0,6)'G-EWALD',TV-TV0,TC-TC0
!=======================================================================
! use inversion through a point midway between ion NI and NNI to obtain
! the forces on ion NNI due to interaction with the periodic array of
! ions NNI
!=======================================================================
      DO NI=1,NIONS
      DO NNI=1,NI-1
         DO M=1,3
            GFORCE(M,NNI,NI)=-GFORCE(M,NI,NNI)
            RFORCE(M,NNI,NI)=-RFORCE(M,NI,NNI)
         ENDDO
      ENDDO
      ENDDO
!=======================================================================
! calculate the total forces on the ions
!=======================================================================
      DO NI=1,NIONS
      DO NNI=1,NIONS
         DO M=1,3
            EWIFOR(M,NI)=EWIFOR(M,NI)+RFORCE(M,NNI,NI)*RSCALE- &
                 &    GFORCE(M,NNI,NI)*GSCALE
         ENDDO
      ENDDO
      ENDDO
!=======================================================================
! mean Z and Z^2
!=======================================================================
      ZION=0
      ZZION=0
      DO N=1,NTYP
        ZION =ZION +ZVAL(N)        *NITYP(N)
        ZZION=ZZION+ZVAL(N)*ZVAL(N)*NITYP(N)
     ENDDO
!=======================================================================
! total forces on the unit cell
!=======================================================================
      EWSIF(1,1)=RSIF1 *RSCALE+ GSIF1  -ZION**2*GENSCA/4
      EWSIF(2,2)=RSIF2 *RSCALE+ GSIF2  -ZION**2*GENSCA/4
      EWSIF(3,3)=RSIF3 *RSCALE+ GSIF3  -ZION**2*GENSCA/4
      EWSIF(1,2)=RSIF12*RSCALE+ GSIF12
      EWSIF(2,3)=RSIF23*RSCALE+ GSIF23
      EWSIF(3,1)=RSIF31*RSCALE+ GSIF31
      EWSIF(2,1)=EWSIF(1,2)
      EWSIF(3,2)=EWSIF(2,3)
      EWSIF(1,3)=EWSIF(3,1)
!=======================================================================
! calculate the total coulomb energy per cell including the energy due
! the interaction between the ions and the uniform neutralising
! background (given by the final two terms in the sum below)
!=======================================================================
      TEWEN=RENSCA*REWEN+ GEWEN *GENSCA-(ZION**2*GENSCA/4) &
     &   -(SCALE*FELECT*ZZION/SQRT(PI))

      RETURN
      CONTAINS


      FUNCTION ER(X)
        REAL(q) X,ER
        ER =  (SQRT(PI)/2.)*ERRFC(X)/X
      END FUNCTION ER
      FUNCTION ERD(X)
        REAL(q) X,ERD
        ERD=-1/(EXP(X*X)*X) - SQRT(PI)*ERRFC(X)/(2*(X*X))
      END FUNCTION ERD

      FUNCTION EX(X)
        REAL(q) X,EX
        EX =  EXP(-(X*X))/(2*(X*X))
      END FUNCTION EX

      FUNCTION EXD(X)
        REAL(q) X,EXD
        EXD=-EXP(-(X*X))/X*(1+1/(X*X))
      END FUNCTION EXD

      END SUBROUTINE

      END MODULE
