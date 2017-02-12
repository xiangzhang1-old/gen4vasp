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





!************************ SUBROUTINE FEXCGC *****************************
! RCS:  $Id: xcgrad.F,v 1.4 2001/01/31 11:52:00 kresse Exp $
!
!  the latest version of the routine requires as input
!  the charge density in real space (CHTOT) and returns
!  the potential in real space (CWORK)
!
!  get GGA potential   (mind not LDA contribution is calculated)
!  this version supports vectorization if
!#define vector
!  is used.
!  In this case only PW-91 is supported, vectorization should be
!  possible by inlining up to 200 lines.
!  If you want to use other GGAs use
!#undef vector
!  Routine was written by jF, and rewritten  by aE and gK
!  to get the potential the algorithm proposed by
!  White and Bird Phys.Rev.B 50,7 (1994) 4954) is used
!  stress is also calculated according to this algorithm
!
! ATTENTION ATTENTION ATTENTION ATTENTION ATTENTION ATTENTION
! We use a quite dangerous construction
! to support  REAL(q) <-> COMPLEX(q)   fft s
! several arrays are passed twice to the routine FEXCG_
! on some compilers this makes troubles,
! we call an external subroutine OPSYNC to avoid that compilers
! move DO Loops around violating our assumption that
! DWORK and CWORK point ot the same location
! (the OPSYNC subroutine actually does nothing at all)
! ATTENTION ATTENTION ATTENTION ATTENTION ATTENTION ATTENTION
!***********************************************************************
      MODULE xcgrad
      USE prec
      CONTAINS

      SUBROUTINE FEXCG(LEXCHG,GRIDC,LATT_CUR,XCENC,EXC,CVZERO,XCSIF, &
                  CHTOT,CWORK,DENCOR)
      USE prec
      USE lattice
      USE mpimy
      USE mgrid

      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRIDC
      TYPE (latt)        LATT_CUR

      COMPLEX(q)  CHTOT(GRIDC%MPLWV),CWORK(GRIDC%MPLWV)
      COMPLEX(q)       DENCOR(GRIDC%RL%NP)
      REAL(q)     XCSIF(3,3)
! work arrays
      COMPLEX(q),ALLOCATABLE:: CWGRAD(:)
      REAL(q),ALLOCATABLE   :: DWORKG(:),DWORK1(:),DWORK2(:),DWORK3(:),DCHARG(:)

      NP1=GRIDC%RL%NP
      ALLOCATE(CWGRAD(GRIDC%MPLWV), &
               DWORKG(NP1),DWORK1(NP1),DWORK2(NP1),DWORK3(NP1),DCHARG(NP1))

      CALL FEXCG_(LEXCHG,GRIDC,LATT_CUR,XCENC,EXC,CVZERO,XCSIF, &
     &            CWGRAD,CHTOT,CWORK, &
     &            CWGRAD,CHTOT,CWORK, &
     &            DENCOR,DWORKG,DWORK1,DWORK2,DWORK3,DCHARG)

      DEALLOCATE(CWGRAD,DWORKG,DWORK1,DWORK2,DWORK3,DCHARG)

      RETURN
      END SUBROUTINE
      END MODULE

      SUBROUTINE FEXCG_(LEXCHG,GRIDC,LATT_CUR,XCENC,EXC,CVZERO,XCSIF, &
     &            CWGRAD,CHTOT,CWORK, &
     &            DWGRAD,DHTOT,DWORK, &
     &            DENCOR,DWORKG,DWORK1,DWORK2,DWORK3,DCHARG)
      USE prec

      USE lattice
      USE mpimy
      USE mgrid
      USE constant

      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRIDC
      TYPE (latt)        LATT_CUR

!  ATTENTION ATTENTION ATTENTION ATTENTION ATTENTION ATTENTION

! Mind CWORK and DWORK point actually to the same storagelocation
! similar to e EQUIVALENCE (CWORK(1),DWORK(1))
! same is true for (CWGRAD,DWGRAD) and   (CHTOT,DHTOT)
! so we can interchange both arrays arbitrarily
!
      COMPLEX(q) CHTOT(GRIDC%MPLWV),CWGRAD(GRIDC%MPLWV),CWORK(GRIDC%MPLWV)
      COMPLEX(q)      DHTOT(GRIDC%MPLWV),DWGRAD(GRIDC%MPLWV),DWORK(GRIDC%MPLWV)
      COMPLEX(q)      DENCOR(GRIDC%RL%NP)
      REAL(q) DWORKG(GRIDC%RL%NP),DWORK1(GRIDC%RL%NP),DWORK2(GRIDC%RL%NP), &
              DWORK3(GRIDC%RL%NP),DCHARG(GRIDC%RL%NP)
      REAL(q) XCSIF(3,3)

      NODE_ME=0
      IONODE =0
      IDUMP=0

! in vector mode only PW91 is supported

      RINPL=1._q/GRIDC%NPLWV
!=======================================================================
! First phase: Transform DENCOR (core charge) and
!  CHTOT (pseudo chargedensity) to real space
!=======================================================================
! get real charge density + core charge in reciprocal space
      DO I=1,GRIDC%RL%NP
         DWORK(I)=(DENCOR(I)+DHTOT(I))*RINPL/LATT_CUR%OMEGA
      ENDDO
      CALL OPSYNC(DWORK,CWORK,GRIDC%NPLWV)
      CALL FFT3RC(CWORK,GRIDC,-1)
      CALL TRUNC_HIGH_FREQU(LATT_CUR, GRIDC, CWORK)
!=======================================================================
! now calculate the gradient of the chargedensity
!=======================================================================
      DO  I=1,GRIDC%RC%NP
        CWGRAD(I)=CWORK(I)
      ENDDO
! x-component:
!DIR$ IVDEP
!$DIR FORCE_VECTOR
!OCL NOVREC
      DO I=1,GRIDC%RC%NP
         N1= MOD((I-1),GRIDC%RC%NROW) +1
         NC= (I-1)/GRIDC%RC%NROW+1
         N2= GRIDC%RC%I2(NC)
         N3= GRIDC%RC%I3(NC)
         GX=(GRIDC%LPCTX(N1)*LATT_CUR%B(1,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(1,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(1,3))
         CWORK(I)=CWORK(I)*GX*CITPI
      ENDDO
! grad_x in real space:
      CALL SETUNB(CWORK,GRIDC)
      CALL FFT3RC(CWORK,GRIDC,1)
      CALL OPSYNC(DWORK,CWORK,GRIDC%NPLWV)
      DO I=1,GRIDC%RL%NP
        DWORK1(I)= REAL( DWORK(I) ,KIND=q)
      ENDDO

! y-component:
      DO I=1,GRIDC%RC%NP
        CWORK(I)=CWGRAD(I)
      ENDDO
!DIR$ IVDEP
!$DIR FORCE_VECTOR
!OCL NOVREC
      DO I=1,GRIDC%RC%NP
         N1= MOD((I-1),GRIDC%RC%NROW) +1
         NC= (I-1)/GRIDC%RC%NROW+1
         N2= GRIDC%RC%I2(NC)
         N3= GRIDC%RC%I3(NC)
         GY=(GRIDC%LPCTX(N1)*LATT_CUR%B(2,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(2,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(2,3))
         CWORK(I)=CWORK(I)*GY*CITPI
      ENDDO
! grad_y in real space:
      CALL SETUNB(CWORK,GRIDC)
      CALL FFT3RC(CWORK,GRIDC,1)
      CALL OPSYNC(DWORK,CWORK,GRIDC%NPLWV)
      DO I=1,GRIDC%RL%NP
        DWORK2(I)= REAL( DWORK(I) ,KIND=q)
      ENDDO

! z-component:
      DO I=1,GRIDC%RC%NP
        CWORK(I)=CWGRAD(I)
      ENDDO
!DIR$ IVDEP
!$DIR FORCE_VECTOR
!OCL NOVREC
      DO I=1,GRIDC%RC%NP
         N1= MOD((I-1),GRIDC%RC%NROW) +1
         NC= (I-1)/GRIDC%RC%NROW+1
         N2= GRIDC%RC%I2(NC)
         N3= GRIDC%RC%I3(NC)
         GZ=(GRIDC%LPCTX(N1)*LATT_CUR%B(3,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(3,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(3,3))
         CWORK(I)=CWORK(I)*GZ*CITPI
      ENDDO
! grad_z in real space:
      CALL SETUNB(CWORK,GRIDC)
      CALL FFT3RC(CWORK,GRIDC,1)
      CALL OPSYNC(DWORK,CWORK,GRIDC%NPLWV)
      DO I=1,GRIDC%RL%NP
        DWORK3(I)= REAL( DWORK(I) ,KIND=q)
      ENDDO
! calculate total charge in real space, and abs nabla rho
!DIR$ IVDEP
!$DIR FORCE_VECTOR
!OCL NOVREC
      DO I=1,GRIDC%RL%NP
         DCHARG(I)=(DHTOT(I)+DENCOR(I))/LATT_CUR%OMEGA
         G2=DWORK1(I)*DWORK1(I)+DWORK2(I)*DWORK2(I)+DWORK3(I)*DWORK3(I)
         DWORKG(I)=SQRT(G2)
      ENDDO
!=======================================================================
!  grad rho    d    f_xc
! ---------- * ------------      (Phys.Rev.B 50,7 (1994) 4954)
! |grad rho|   d |grad rho|
!
!  MIND: the factor OMEGA is problematic
!   1/N sum_r energy_density * rho *OMEGA = Energy
!   1/N sum_r energy_density * \bar rho   = Energy (\bar rho=rho*OMEGA)
!  the array DCHARG(I) is the real charge density (incl. part. core)
!=======================================================================

      EXC=0
      DO I=1,GRIDC%RL%NP
         RHO= DCHARG(I)
         CALL GGAALL(LEXCHG,RHO*AUTOA3,DWORKG(I)*AUTOA4,EXCL,DEXC,DVXC,.FALSE.)
         EXC=EXC+EXCL*RHO*RYTOEV*LATT_CUR%OMEGA
         DVXC=DVXC*RYTOEV*AUTOA
        !  store d f/ d (|d rho| ) / |d rho|  in DWORK
         DWORK(I)  = DVXC / MAX(DWORKG(I),1.E-10_q)
        !  store d f/ d rho  in DWORKG
         DWORKG(I) = DEXC*RYTOEV
      ENDDO
!=======================================================================
! gradient terms in stress tensor
!          d    f_xc     grad rho  x grad rho
! sum_r   ------------   --------------------- * LATT_CUR%OMEGA
!         d |grad rho|        |grad rho|
!=======================================================================
      SIF11=0
      SIF22=0
      SIF33=0
      SIF12=0
      SIF23=0
      SIF31=0
      DO I=1,GRIDC%RL%NP
        SIF11=SIF11+DWORK1(I)*DWORK1(I)*DWORK(I)
        SIF22=SIF22+DWORK2(I)*DWORK2(I)*DWORK(I)
        SIF33=SIF33+DWORK3(I)*DWORK3(I)*DWORK(I)
        SIF12=SIF12+DWORK1(I)*DWORK2(I)*DWORK(I)
        SIF23=SIF23+DWORK2(I)*DWORK3(I)*DWORK(I)
        SIF31=SIF31+DWORK3(I)*DWORK1(I)*DWORK(I)
      ENDDO
      SIF11=SIF11*RINPL*LATT_CUR%OMEGA
      SIF22=SIF22*RINPL*LATT_CUR%OMEGA
      SIF33=SIF33*RINPL*LATT_CUR%OMEGA
      SIF12=SIF12*RINPL*LATT_CUR%OMEGA
      SIF23=SIF23*RINPL*LATT_CUR%OMEGA
      SIF31=SIF31*RINPL*LATT_CUR%OMEGA

      DO I=1,GRIDC%RL%NP
         DWORK1(I) = DWORK1(I)* REAL( DWORK(I) ,KIND=q)
         DWORK2(I) = DWORK2(I)* REAL( DWORK(I) ,KIND=q)
         DWORK3(I) = DWORK3(I)* REAL( DWORK(I) ,KIND=q)
      ENDDO
!=======================================================================
! times i G_k in reciprocal space...
!=======================================================================
! x-component:
      DO I=1,GRIDC%RL%NP
        DWORK(I) = DWORK1(I)
      ENDDO
      CALL OPSYNC(DWORK,CWORK,GRIDC%NPLWV)
      CALL FFT3RC(CWORK,GRIDC,-1)
      CALL TRUNC_HIGH_FREQU(LATT_CUR, GRIDC, CWORK)

      DO I=1,GRIDC%RC%NP
         N1= MOD((I-1),GRIDC%RC%NROW) +1
         NC= (I-1)/GRIDC%RC%NROW+1
         N2= GRIDC%RC%I2(NC)
         N3= GRIDC%RC%I3(NC)
         GX=(GRIDC%LPCTX(N1)*LATT_CUR%B(1,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(1,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(1,3))
         CWGRAD(I)=CWORK(I)*GX*CITPI
      ENDDO

! y-component:
      DO I=1,GRIDC%RL%NP
         DWORK(I) = DWORK2(I)
      ENDDO
      CALL OPSYNC(DWORK,CWORK,GRIDC%NPLWV)
      CALL FFT3RC(CWORK,GRIDC,-1)
      CALL TRUNC_HIGH_FREQU(LATT_CUR, GRIDC, CWORK)

      DO I=1,GRIDC%RC%NP
         N1= MOD((I-1),GRIDC%RC%NROW) +1
         NC= (I-1)/GRIDC%RC%NROW+1
         N2= GRIDC%RC%I2(NC)
         N3= GRIDC%RC%I3(NC)
         GY=(GRIDC%LPCTX(N1)*LATT_CUR%B(2,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(2,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(2,3))
         CWGRAD(I)=CWGRAD(I)+CWORK(I)*GY*CITPI
      ENDDO

! z-component:
      DO I=1,GRIDC%RL%NP
         DWORK(I) = DWORK3(I)
      ENDDO
      CALL OPSYNC(DWORK,CWORK,GRIDC%NPLWV)
      CALL FFT3RC(CWORK,GRIDC,-1)
      CALL TRUNC_HIGH_FREQU(LATT_CUR, GRIDC, CWORK)

      DO I=1,GRIDC%RC%NP
         N1= MOD((I-1),GRIDC%RC%NROW) +1
         NC= (I-1)/GRIDC%RC%NROW+1
         N2= GRIDC%RC%I2(NC)
         N3= GRIDC%RC%I3(NC)
         GZ=(GRIDC%LPCTX(N1)*LATT_CUR%B(3,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(3,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(3,3))
         CWGRAD(I)=CWGRAD(I)+CWORK(I)*GZ*CITPI
      ENDDO

      CALL SETUNB(CWGRAD,GRIDC)
      CALL FFT3RC(CWGRAD,GRIDC,1)
      CALL OPSYNC(CWGRAD,DWGRAD,GRIDC%NPLWV)
!=======================================================================
! Now prepare the rest:
! (store rho in DWORK3 and quantity of above in DWORK1)
!=======================================================================
      XCENC=0._q
      CVZERO=0._q
      XCENCC=0._q

      DO I=1,GRIDC%RL%NP
         RHO= REAL( DCHARG(I) ,KIND=q)
         VXC=DWORKG(I)- REAL( DWGRAD(I) ,KIND=q) *RINPL
         DWORK(I)=VXC
         CVZERO=CVZERO+VXC
         XCENCC=XCENCC-VXC*RHO*LATT_CUR%OMEGA
         XCENC=XCENC  -VXC* REAL( DHTOT(I) ,KIND=q)
      ENDDO
! array reduction

      EXC   =EXC
      CVZERO=CVZERO*RINPL
      XCENC =(XCENC+EXC)*RINPL
      XCENCC=(XCENCC+EXC)*RINPL
      EXC=EXC*RINPL

      SIF11=SIF11-XCENCC
      SIF22=SIF22-XCENCC
      SIF33=SIF33-XCENCC
      XCSIF(1,1)=SIF11
      XCSIF(2,2)=SIF22
      XCSIF(3,3)=SIF33
      XCSIF(1,2)=SIF12
      XCSIF(2,1)=SIF12
      XCSIF(2,3)=SIF23
      XCSIF(3,2)=SIF23
      XCSIF(3,1)=SIF31
      XCSIF(1,3)=SIF31

      
      
      

! Test dumps:
      IF (IDUMP/=0) THEN
         WRITE(*,'(A,F24.14)') '<rho*excgc> =',EXC
         WRITE(*,'(A,F24.14)') '<rho*vxcgc> =',EXC-XCENC
         WRITE(*,'(A,F24.14)') '    xcencgc =',XCENC
      ENDIF

      RETURN
      END


!************************ SUBROUTINE GGAALL *****************************
!
!  switch between different GGAs
!  presently only PW91, PBE and RPBE are implemented
!  (i.e. d exc / d rho and  d exc / d | grad rho | are calculated
!  directly 
!  for other GGA functional finite differences are used to calculate
!  the required derivatives
!  LLDA allows to include the LDA contribution directly in this 
!  routine.
!  This only works for PBE and RPBE and fails in all other cases
!
!***********************************************************************

      SUBROUTINE GGAALL(LEXCHG,D,DD, EXC,EXCD,EXCDD,LLDA)
      USE prec
      USE constant
      IMPLICIT REAL(q) (A-H,O-Z)
      PARAMETER (DDELTA=1E-4_q)
      PARAMETER (THRD=1._q/3._q)
      LOGICAL LLDA

      IF (LEXCHG==4) THEN

        ! PW91 using the routines of Bird and White

        CALL GGA91_WB(D,DD,EXC,EXCD,EXCDD)
        EXC = 2*EXC / D
        EXCD =2*EXCD
        EXCDD=2*EXCDD
!        WRITE(*,'(10F14.7)') D, EXC,EXCD,EXCDD

      ELSE IF (LEXCHG==5 .OR. LEXCHG==6) THEN

        ! Perdew Burke Ernzerhof and revised functional

        IF (LEXCHG==5) THEN
           ukfactor=1.0_q
        ELSE
           ukfactor=0.0_q
        ENDIF

        IF (D<=0) THEN
           EXC   = 0._q
           EXCD  = 0._q
           EXCDD = 0._q
           RETURN
        ENDIF

        DTHRD=exp(log(D)*THRD)
        RS=(0.75_q/PI)**THRD/DTHRD
        FK=(3._q*PI*PI)**THRD*DTHRD
        SK = SQRT(4.0_q*FK/PI)
        IF(D>1.E-10_q)THEN
           S=DD/(D*FK*2._q)
           T=DD/(D*SK*2._q)
        ELSE
           S=0.0_q
           T=0.0_q
        ENDIF

        CALL EXCHPBE(D,DTHRD,S,EXLDA,EXC,EXDLDA,EXCD,EXCDD, &
             ukfactor)

        CALL CORunspPBE(RS,ECLDA,ECDLDA,SK, &
             T,EC,ECD,ECDD,.TRUE.)

!        WRITE(*,'(10F14.7)') D,S,(EXC-EXLDA)/D,EXCD-EXDLDA,EXCDD
!        WRITE(*,'(10F14.7)') D,RS,SK,T,EC,ECD,ECDD

        IF (LLDA) THEN
          EXC  =EC  +EXC/D+ECLDA
          EXCD =ECD +EXCD +ECDLDA
          EXCDD=ECDD+EXCDD
        ELSE
          EXC  =EC  +(EXC-EXLDA)/D
          EXCD =ECD +EXCD-EXDLDA
          EXCDD=ECDD+EXCDD
        ENDIF

        ! Hartree -> Rydberg conversion
        EXC = 2*EXC
        EXCD =2*EXCD
        EXCDD=2*EXCDD

!        WRITE(*,'(10F14.7)') D, EXC,EXCD,EXCDD
      ELSE

        ! for all other functionals
        ! we use finite differences to calculate the required 
        ! quantities 
        ! presently no other functional are supported
        ! but in case (1._q,0._q) needs these routines
         
        DELTA=MIN(DDELTA,ABS(D)/100)
        D1=D-DELTA
        D2=D+DELTA
        CALL GGAEALL(LEXCHG,D1,DD,1._q,1._q,VXC,EXC1)
        CALL GGAEALL(LEXCHG,D2,DD,1._q,1._q,VXC,EXC2)
        EXCD=(EXC2*D2-EXC1*D1)/MAX((D2-D1),1E-10_q)

        DELTA=MIN(DDELTA,ABS(DD)/100)
        DD1=DD-DELTA
        DD2=DD+DELTA
        CALL GGAEALL(LEXCHG,D,DD1,1._q,1._q,VXC,EXC1)
        CALL GGAEALL(LEXCHG,D,DD2,1._q,1._q,VXC,EXC2)
        EXCDD=(EXC2*D-EXC1*D)/MAX((DD2-DD1),1E-10_q)
        CALL GGAEALL(LEXCHG,D,DD,1._q,1._q,VXC,EXC)
      ENDIF
      RETURN
      END

!
!  Common interface to all gradient correction routines:
!  returns energy + potential (if defined)
!
      SUBROUTINE GGAEALL(LEXCHG,RHO,DRHO,DABGRH,DDRHO,VXC,EXC)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      
      WRITE(*,*) 'internal ERROR GGAEALL: Wrong LEXCHG, scheme not implemented!'
      STOP

      RETURN
      END



!************************ SUBROUTINE GGA91_WB **************************
!
!  calculates the PW 91 xc-functional using the modifications of
!  J.A. White and D.M. Bird
!  written by  Perdew
!  modified by J.A. White and D.M. Bird
!    (Phys.Rev.B 50,7 (1994) 4954)
!
!***********************************************************************

!
!  interface for PW-91 only (used in vector mode only)
!  this subroutine converts to Rydberg units
!
      SUBROUTINE GGA91_WB_RY(D,DD,EXC,EXCD,EXCDD)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      CALL GGA91_WB(D,DD,EXC,EXCD,EXCDD)
      EXC = 2._q *EXC / D
      EXCD =2._q *EXCD
      EXCDD=2._q *EXCDD
      RETURN
      END

      SUBROUTINE GGA91_WB(D,DD,EXC,EXCD,EXCDD)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      REAL(q),PARAMETER :: PI=3.1415926536_q
      REAL(q),PARAMETER :: THRD=0.333333333333333_q
      IF (D<0) THEN
       EXC   = 0._q
       EXCD  = 0._q
       EXCDD = 0._q
       RETURN
      ENDIF

      G = 1.0_q
      RS=(0.75_q/(PI*D))**THRD
      FK=(3._q*PI*PI*D)**THRD
      SK = SQRT(4.0_q*FK/PI)
      IF(D>1.E-10_q)THEN
        S=DD/(D*FK*2._q)
        T=DD/(D*SK*2._q)
      ELSE
        S=0.0_q
        T=0.0_q
      ENDIF
      CALL EXCH2(D,S,EX,EXD,EXDD)
      CALL GCOR1(0.0310907_q,0.21370_q,7.5957_q,3.5876_q,1.6382_q, &
     &    0.49294_q,1.00_q,RS,EU,EURS)
      EC = 0._q
      ECD = 0._q
      CALL CORGGA1(D,RS,T,EC1,EC1D,EC1DD, &
     &             FK,SK,G,EU,EURS,ECZET)
      EXC = EX+EC+EC1
!      WRITE(*,'(10F14.6)') D,S,EX/D,EXD,EXDD
!      WRITE(*,'(10F14.6)') D,RS,SK,T,(EC+EC1)/D,ECD+EC1D,EC1DD
      EXCD = EXD+ECD+EC1D
      EXCDD = EXDD+EC1DD
      RETURN
      END
!=======================================================================
! differecne betwenn EXCH1 and EXCH2 is that EXCH1 includes
! LDA part wheres EXCH2 excludes this part
!
      SUBROUTINE EXCH2(D,S,EX,EXD,EXDD)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
!  gga91 exchange for a spin-unpolarized electronic system
!  input d : density
!  input s:  abs(grad d)/(2*kf*d)
!  output:  exchange energy per electron (ex) and ites derivatives
!           w.r.t. d (exd) and dd (exdd)
      REAL(q),PARAMETER :: A1=0.19645_q,A2=0.27430_q,A3=0.15084_q,A4=100._q
      REAL(q),PARAMETER :: AX=-0.7385588_q,A=7.7956_q,B1=0.004_q
      REAL(q),PARAMETER :: THPITH=3.0936677262801_q
      THRD = 1._q/3._q
      THRD4 = 4._q/3._q
      FAC = AX*D**THRD
      S2 = S*S
      S3 = S2*S
      S4 = S2*S2
      P0 = 1.0_q/SQRT(1.0_q+A*A*S2)
      P1 = LOG(A*S+1.0_q/P0)
      P2 = EXP(-A4*S2)
      P3 = 1.0_q/(1.0_q+A1*S*P1+B1*S4)
      P4 = 1.0_q+A1*S*P1+(A2-A3*P2)*S2
      F = P3*P4
      EX = FAC*(F-1)*D
!     IF (S==0)
!    &         PRINT *,P1,P2,P3,P4,FAC,F
      P5 = 2.0_q*(S*(A2-A3*P2)+A3*A4*S3*P2-2.0_q*B1*S3)
      P6 = (A1*(P1+A*S*P0)+4.0_q*B1*S3)*((A2-A3*P2)*S2-B1*S4)
      FS = (P5*P3-P6*P3*P3)
      EXD = THRD4*FAC*(F-S*FS-1)
      EXDD = AX*FS*0.5_q/THPITH
      RETURN
      END
!=======================================================================
      SUBROUTINE GCOR1(A,A1,B1,B2,B3,B4,P,RS,GG,GGRS)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
!  called by subroutine corlsd
      P1 = P + 1.0_q
      Q0 = -2.0_q*A*(1.0_q+A1*RS)
      RS12 = SQRT(RS)
      RS32 = RS12**3
      RSP = RS**P
      Q1 = 2.0_q*A*(B1*RS12+B2*RS+B3*RS32+B4*RS*RSP)
      Q2 = LOG(1.0_q+1.0_q/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/RS12+2.0_q*B2+3.0_q*B3*RS12+2.0_q*B4*P1*RSP)
      GGRS = -2.0_q*A*A1*Q2-Q0*Q3/(Q1**2+Q1)
      RETURN
      END
!================================================================
      SUBROUTINE CORGGA1(D,RS,T,EC1,EC1D,EC1DD, &
     &                  FK,SK,G,EC,ECRS,ECZET)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
!  gga91 correlation
!  input rs: seitz radius
!  input t: abs(grad d)/(d*2.*ks)
!  output hn: nonlocal part of correlation energy
!  hnd,hndd : derivatives of hn w.r.t. d and dd
      REAL(q),PARAMETER :: XNU=15.75592_q,CC0=0.004235_q,CX=-0.001667212_q,ALF=0.09_q
      REAL(q),PARAMETER :: C1=0.002568_q,C2=0.023266_q,C3=7.389E-6_q,C4=8.723_q
      REAL(q),PARAMETER :: C5=0.472_q,C6=7.389E-2_q,A4=100.0_q
      REAL(q),PARAMETER :: THRD=0.33333333333333_q,SIXTH7=1.16666666666667_q
      REAL(q),PARAMETER :: PI=3.1415926536_q
      BET = XNU*CC0
      DELT = 2.0_q*ALF/BET
      G3 = G**3
      G4 = G3*G
      PON = -DELT*EC/(G3*BET)
      B = DELT/(EXP(PON)-1.0_q)
      B2 = B*B
      T2 = T*T
      T4 = T2*T2
      T6 = T4*T2
      RS2 = RS*RS
      RS3 = RS2*RS
      Q4 = 1.0_q+B*T2
      Q5 = 1.0_q+B*T2+B2*T4
      Q6 = C1+C2*RS+C3*RS2
      Q7 = 1.0_q+C4*RS+C5*RS2+C6*RS3
      CC = -CX + Q6/Q7
      R0 = (SK/FK)**2
      R1 = A4*R0*G4
      COEFF = CC-CC0-3.0_q*CX/7.0_q
      R2 = XNU*COEFF*G3
      R3 = EXP(-R1*T2)
      H0 = G3*(BET/DELT)*LOG(1.0_q+DELT*Q4*T2/Q5)
      H1 = R3*R2*T2
      H = H0+H1
!============================================================
      Q8 = Q5*Q5+DELT*Q4*Q5*T2
      H0T = 2.0_q*BET*T*(1.0_q+2.0_q*B*T2)/Q8
      H0B = -BET*T6*(2.0_q*B+B2*T2)/Q8
      H0RS = H0B*B*ECRS*(B+DELT)/BET
      H1T = 2.0_q*R3*R2*T*(1.0_q-R1*T2)
      CCRS = (C2+2._q*C3*RS)/Q7 - Q6*(C4+2._q*C5*RS+3._q*C6*RS2)/Q7**2
      R1RS = 100.0_q*R0/RS
      H1RS = XNU*T2*R3*(CCRS - COEFF*T2*R1RS)
! = = = = = = = = = = = =
      HT = H0T + H1T
      HRS = H0RS + H1RS
      EC1 = D*H
      EC1D = H-THRD*RS*HRS-SIXTH7*T*HT
      EC1DD = 0.5_q*HT/SK
      RETURN
      END

!---------------------------------------------------------------------
! The P B E
! this routine was kindly supplied by Bjork Hammer
! it is essentially identical to Burkes routine but includes
! the revised P B E routines
!---------------------------------------------------------------------
!######################################################################
!----------------------------------------------------------------------
      SUBROUTINE EXCHPBE(rho,rhothrd,s,exlda,expbe,exdlda,exd,exdd, &
     &                   ukfactor)
!----------------------------------------------------------------------
!  PBE EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
!  K Burkes modification of PW91 codes, May 14, 1996
!  Modified again by K. Burke, June 29, 1996, with simpler Fx(s)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!  INPUT rho : DENSITY
!  INPUT rhothrd : DENSITY^(1/3)
!  INPUT S:  ABS(GRAD rho)/(2*KF*rho), where kf=(3 pi^2 rho)^(1/3)
!  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX) 
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! References:
! [a]J.P.~Perdew, K.~Burke, and M.~Ernzerhof, submiited to PRL, May96
! [b]J.P. Perdew and Y. Wang, Phys. Rev.  B {\bf 33},  8800  (1986);
!     {\bf 40},  3399  (1989) (E).
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Formulas:
!       e_x[unif]=ax*rho^(4/3)  [LDA]
! ax = -0.75*(3/pi)^(1/3)
!       e_x[PBE]=e_x[unif]*FxPBE(s)
!       FxPBE(s)=1+uk-uk/(1+ul*s*s)                 [a](13)
! uk, ul defined after [a](13) 
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      parameter(thrd=1._q/3._q,thrd4=4._q/3._q)
      parameter(pi=3.14159265358979323846264338327950_q)
      parameter(ax=-0.738558766382022405884230032680836_q)
      parameter(um=0.2195149727645171_q,uk1=0.8040_q,ul1=um/uk1)
!----------------------------------------------------------------------
! construct LDA exchange energy density
      exunif=AX*rhothrd
      exlda=exunif*rho
      exdlda=exunif*thrd4
!----------------------------------------------------------------------
! construct PBE enhancement factor
      S2 = S*S
!----------------------------------------------------------------------
      if (ukfactor.ne.0.0_q) then
! These are the PBE96 and revPBE98 functionals
! scale uk with a factor
         uk = uk1*ukfactor
         ul = ul1/ukfactor 
         P0=1._q+ul*S2
         FxPBE = 1._q+uk-uk/P0
         expbe = exlda*FxPBE
!----------------------------------------------------------------------
!  ENERGY DONE. NOW THE POTENTIAL:
!  find first derivatives of Fx w.r.t s.
!  Fs=(1/s)*d FxPBE/ ds
         Fs=2._q*um/(P0*P0)
      else 
! This is the RPBE functional [Hammer et al, PRB 59, 7413 (1999)]
         P0=exp(-ul1*S2)
         FxPBE = 1._q+uk1*(1.0_q-P0)
         expbe = exlda*FxPBE
!----------------------------------------------------------------------
!  ENERGY DONE. NOW THE POTENTIAL:
!  find first derivatives of Fx w.r.t s.
!  Fs=(1/s)*d FxPBE/ ds
         Fs=2._q*um*P0
      endif
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! calculate the partial derivatives of ex wrt n and |grad(n)|
!  0.3232409194=(3*pi^2)^(-1/3)
      exd =exunif*THRD4*(FxPBE-S2*Fs)
      exdd=0.5_q*ax*0.3232409194_q*S*Fs
      RETURN
      END
!----------------------------------------------------------------------
!----------------------------------------------------------------------

!######################################################################
!----------------------------------------------------------------------
      SUBROUTINE CORunspPBE(RS,EC,VC,sk, &
     &                  T,H,DVC,ecdd,lgga)
!----------------------------------------------------------------------
!  Official PBE correlation code. K. Burke, May 14, 1996.
!  INPUT: RS=SEITZ RADIUS=(3/4pi rho)^(1/3)
!       : ZET=RELATIVE SPIN POLARIZATION = (rhoup-rhodn)/rho
!       : t=ABS(GRAD rho)/(rho*2.*KS*G)  -- only needed for PBE
!       : lgga=flag to do gga (0=>LSD only)
!       : lpot=flag to do potential (0=>energy only)
!  output: ec=lsd correlation energy from [a]
!        : vcup=lsd up correlation potential
!        : vcdn=lsd dn correlation potential
!        : h=NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
!        : dvcup=nonlocal correction to vcup
!        : dvcdn=nonlocal correction to vcdn
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! References:
! [a] J.P.~Perdew, K.~Burke, and M.~Ernzerhof, 
!     {\sl Generalized gradient approximation made simple}, sub.
!     to Phys. Rev.Lett. May 1996.
! [b] J. P. Perdew, K. Burke, and Y. Wang, {\sl Real-space cutoff
!     construction of a generalized gradient approximation:  The PW91
!     density functional}, submitted to Phys. Rev. B, Feb. 1996.
! [c] J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
!----------------------------------------------------------------------
!----------------------------------------------------------------------
       USE prec
       IMPLICIT REAL(q) (A-H,O-Z)
! thrd*=various multiples of 1/3
! numbers for use in LSD energy spin-interpolation formula, [c](9).
!      GAM= 2^(4/3)-2
!      FZZ=f''(0)= 8/(9*GAM)
! numbers for construction of PBE
!      gamma=(1-log(2))/pi^2
!      bet=coefficient in gradient expansion for correlation, [a](4).
!      eta=small number to stop d phi/ dzeta from blowing up at 
!          |zeta|=1.
      logical*4 lgga
      parameter(thrd=1._q/3._q,thrdm=-thrd,thrd2=2._q*thrd)
      parameter(sixthm=thrdm/2._q)
      parameter(thrd4=4._q*thrd)
      parameter(GAM=0.5198420997897463295344212145565_q)
      parameter(fzz=8._q/(9._q*GAM))
      parameter(gamma=0.03109069086965489503494086371273_q)
      parameter(bet=0.06672455060314922_q,delt=bet/gamma)
      parameter(eta=1.e-12_q)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! find LSD energy contributions, using [c](10) and Table I[c].
! EU=unpolarized LSD correlation energy
! EURS=dEU/drs
! EP=fully polarized LSD correlation energy
! EPRS=dEP/drs
! ALFM=-spin stiffness, [c](3).
! ALFRSM=-dalpha/drs
! F=spin-scaling factor from [c](9).
! construct ec, using [c](8)
      rtrs=dsqrt(rs)
      CALL gcor2(0.0310907_q,0.21370_q,7.5957_q,3.5876_q,1.6382_q, &
     &    0.49294_q,rtrs,EU,EURS)
      EC = EU
! check for (0._q,0._q) energy, immediate return if true
      IF (EC==0.) THEN
        H=0; DVC=0; ecdd=0
        RETURN
      ENDIF
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! LSD potential from [c](A1)
! ECRS = dEc/drs [c](A2)
! ECZET=dEc/dzeta [c](A3)
! FZ = dF/dzeta [c](A4)
      ECRS = EURS
      VC = EC -RS*ECRS/3._q
      if (.not.lgga) return
!----------------------------------------------------------------------
! PBE correlation energy
! G=phi(zeta), given after [a](3)
! DELT=bet/gamma
! B=A of [a](8)
      PON=-EC/(gamma)
      B = DELT/(DEXP(PON)-1._q)
      B2 = B*B
      T2 = T*T
      T4 = T2*T2
      RS2 = RS*RS
      RS3 = RS2*RS
      Q4 = 1._q+B*T2
      Q5 = 1._q+B*T2+B2*T4
      H = (BET/DELT)*DLOG(1._q+DELT*Q4*T2/Q5)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! ENERGY DONE. NOW THE POTENTIAL, using appendix E of [b].
      T6 = T4*T2
      RSTHRD = RS/3._q
      FAC = DELT/B+1._q
      BEC = B2*FAC/(BET)
      Q8 = Q5*Q5+DELT*Q4*Q5*T2
      Q9 = 1._q+2._q*B*T2
      hB = -BET*B*T6*(2._q+B*T2)/Q8
      hRS = -RSTHRD*hB*BEC*ECRS
      hT = 2._q*BET*Q9/Q8
      DVC = H+HRS-7.0_q*T2*HT/6._q
      ecdd=0.5_q/sk*t*ht
      RETURN
      END
!----------------------------------------------------------------------


!######################################################################
!----------------------------------------------------------------------
      SUBROUTINE CORPBE(RS,ZET,EC,VCUP,VCDN,g,sk, &
     &                  T,H,DVCUP,DVCDN,ecdd,lgga)
!----------------------------------------------------------------------
!  Official PBE correlation code. K. Burke, May 14, 1996.
!  INPUT: RS=SEITZ RADIUS=(3/4pi rho)^(1/3)
!       : ZET=RELATIVE SPIN POLARIZATION = (rhoup-rhodn)/rho
!       : t=ABS(GRAD rho)/(rho*2.*KS*G)  -- only needed for PBE
!       : lgga=flag to do gga (0=>LSD only)
!       : lpot=flag to do potential (0=>energy only)
!  output: ec=lsd correlation energy from [a]
!        : vcup=lsd up correlation potential
!        : vcdn=lsd dn correlation potential
!        : h=NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
!        : dvcup=nonlocal correction to vcup
!        : dvcdn=nonlocal correction to vcdn
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! References:
! [a] J.P.~Perdew, K.~Burke, and M.~Ernzerhof, 
!     {\sl Generalized gradient approximation made simple}, sub.
!     to Phys. Rev.Lett. May 1996.
! [b] J. P. Perdew, K. Burke, and Y. Wang, {\sl Real-space cutoff
!     construction of a generalized gradient approximation:  The PW91
!     density functional}, submitted to Phys. Rev. B, Feb. 1996.
! [c] J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
! thrd*=various multiples of 1/3
! numbers for use in LSD energy spin-interpolation formula, [c](9).
!      GAM= 2^(4/3)-2
!      FZZ=f''(0)= 8/(9*GAM)
! numbers for construction of PBE
!      gamma=(1-log(2))/pi^2
!      bet=coefficient in gradient expansion for correlation, [a](4).
!      eta=small number to stop d phi/ dzeta from blowing up at 
!          |zeta|=1.
      logical*4 lgga
      parameter(thrd=1._q/3._q,thrdm=-thrd,thrd2=2._q*thrd)
      parameter(sixthm=thrdm/2._q)
      parameter(thrd4=4._q*thrd)
      parameter(GAM=0.5198420997897463295344212145565_q)
      parameter(fzz=8._q/(9._q*GAM))
      parameter(gamma=0.03109069086965489503494086371273_q)
      parameter(bet=0.06672455060314922_q,delt=bet/gamma)
      parameter(eta=1.e-12_q)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! find LSD energy contributions, using [c](10) and Table I[c].
! EU=unpolarized LSD correlation energy
! EURS=dEU/drs
! EP=fully polarized LSD correlation energy
! EPRS=dEP/drs
! ALFM=-spin stiffness, [c](3).
! ALFRSM=-dalpha/drs
! F=spin-scaling factor from [c](9).
! construct ec, using [c](8)
      rtrs=dsqrt(rs)
      CALL gcor2(0.0310907_q,0.21370_q,7.5957_q,3.5876_q,1.6382_q, &
     &    0.49294_q,rtrs,EU,EURS)
      CALL gcor2(0.01554535_q,0.20548_q,14.1189_q,6.1977_q,3.3662_q, &
     &    0.62517_q,rtRS,EP,EPRS)
      CALL gcor2(0.0168869_q,0.11125_q,10.357_q,3.6231_q,0.88026_q, &
     &    0.49671_q,rtRS,ALFM,ALFRSM)
      ALFC = -ALFM
      Z4 = ZET**4
      F=((1._q+ZET)**THRD4+(1._q-ZET)**THRD4-2._q)/GAM
      EC = EU*(1._q-F*Z4)+EP*F*Z4-ALFM*F*(1._q-Z4)/FZZ
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! LSD potential from [c](A1)
! ECRS = dEc/drs [c](A2)
! ECZET=dEc/dzeta [c](A3)
! FZ = dF/dzeta [c](A4)
      ECRS = EURS*(1._q-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1._q-Z4)/FZZ
      FZ = THRD4*((1._q+ZET)**THRD-(1._q-ZET)**THRD)/GAM
      ECZET = 4._q*(ZET**3)*F*(EP-EU+ALFM/FZZ)+FZ*(Z4*EP-Z4*EU &
     &        -(1._q-Z4)*ALFM/FZZ)
      COMM = EC -RS*ECRS/3._q-ZET*ECZET
      VCUP = COMM + ECZET
      VCDN = COMM - ECZET
      if (.not.lgga) return
!----------------------------------------------------------------------
! PBE correlation energy
! G=phi(zeta), given after [a](3)
! DELT=bet/gamma
! B=A of [a](8)
      G3 = G**3
      PON=-EC/(G3*gamma)
      B = DELT/(DEXP(PON)-1._q)
      B2 = B*B
      T2 = T*T
      T4 = T2*T2
      RS2 = RS*RS
      RS3 = RS2*RS
      Q4 = 1._q+B*T2
      Q5 = 1._q+B*T2+B2*T4
      H = G3*(BET/DELT)*DLOG(1._q+DELT*Q4*T2/Q5)
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! ENERGY DONE. NOW THE POTENTIAL, using appendix E of [b].
      G4 = G3*G
      T6 = T4*T2
      RSTHRD = RS/3._q
      GZ=(((1._q+zet)**2+eta)**sixthm- &
     &((1._q-zet)**2+eta)**sixthm)/3._q
      FAC = DELT/B+1._q
      BG = -3._q*B2*EC*FAC/(BET*G4)
      BEC = B2*FAC/(BET*G3)
      Q8 = Q5*Q5+DELT*Q4*Q5*T2
      Q9 = 1._q+2._q*B*T2
      hB = -BET*G3*B*T6*(2._q+B*T2)/Q8
      hRS = -RSTHRD*hB*BEC*ECRS
      hZ = 3._q*GZ*h/G + hB*(BG*GZ+BEC*ECZET)
      hT = 2._q*BET*G3*Q9/Q8
      COMM = H+HRS-7.0_q*T2*HT/6._q
      PREF = HZ-GZ*T2*HT/G
      COMM = COMM-PREF*ZET
      DVCUP = COMM + PREF
      DVCDN = COMM - PREF
      ecdd=0.5_q/(sk*g)*t*ht
      RETURN
      END
!----------------------------------------------------------------------
!######################################################################
!----------------------------------------------------------------------
      SUBROUTINE GCOR2(A,A1,B1,B2,B3,B4,rtrs,GG,GGRS)
! slimmed down version of GCOR used in PW91 routines, to interpolate
! LSD correlation energy, as given by (10) of
! J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
! K. Burke, May 11, 1996.
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      Q0 = -2._q*A*(1._q+A1*rtrs*rtrs)
      Q1 = 2._q*A*rtrs*(B1+rtrs*(B2+rtrs*(B3+B4*rtrs)))
      Q2 = DLOG(1._q+1._q/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/rtrs+2._q*B2+rtrs*(3._q*B3+4._q*B4*rtrs))
      GGRS = -2._q*A*A1*Q2-Q0*Q3/(Q1*(1._q+Q1))
      RETURN
      END
!----------------------------------------------------------------------
!######################################################################
!----------------------------------------------------------------------



