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





!************************ SUBROUTINE FEXGCS *****************************
! RCS:  $Id: xcspin.F,v 1.6 2003/06/27 13:22:24 kresse Exp kresse $
!
!  the latest version of the routine requires as input
!  the charge density in real space (CHTOT) and returns
!  the potential in real space (CWORK)
!
!  Routine was written by Elio Moroni, and rewritten to f90 by gK
!  to get the potentials the algorithm proposed by
!  White and Bird Phys.Rev.B 50,7 (1994) 4954) is used
!  stress is also calculated according to this algorithm
!
! ATTENTION ATTENTION ATTENTION ATTENTION ATTENTION ATTENTION
! We use a quite dangerous construction
! to support  REAL(q) <-> COMPLEX(q)   ffts
! several arrays are passed twice to the routine FEXCG_
! on some compilers this makes troubles,
! we call an external subroutine OPSYNC to avoid that compilers
! move DO Loops around violating our assumption that
! DWORK and CWORK point to the same location in storage
! (the OPSYNC subroutine actually does nothing at all)
! ATTENTION ATTENTION ATTENTION ATTENTION ATTENTION ATTENTION
!*********** ************************************************************

      SUBROUTINE FEXCGS(NCDIJ,GRIDC,LATT_CUR,XCENC,EXC,CVZERO,XCSIF, &
                  CHTOT,CWORK,DENCOR, LEXCHG)
      USE prec
      USE lattice
      USE mpimy
      USE mgrid

      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRIDC
      TYPE (latt)        LATT_CUR

      COMPLEX(q) CHTOT(GRIDC%MPLWV,NCDIJ),CWORK(GRIDC%MPLWV,NCDIJ)
      COMPLEX(q)      DENCOR(GRIDC%RL%NP)
      DIMENSION XCSIF(3,3)
! work arrays
      COMPLEX(q),ALLOCATABLE:: CWGRAD(:,:)
      REAL(q),ALLOCATABLE   :: DWORKG(:,:),DWORK1(:,:),DWORK2(:,:),DWORK3(:,:), &
     &                      DVC(:)

      NP1=GRIDC%RL%NP
      IF (NCDIJ==2) THEN
         ALLOCATE(CWGRAD(GRIDC%MPLWV,NCDIJ), DWORKG(NP1,NCDIJ), &
           DWORK1(NP1,NCDIJ),DWORK2(NP1,NCDIJ),DWORK3(NP1,NCDIJ),DVC(NP1))

         CALL FEXCGS_(NCDIJ,GRIDC,LATT_CUR,XCENC,EXC,CVZERO,XCSIF, &
        &            CWGRAD,CHTOT,CWORK, &
        &            CWGRAD,CHTOT,CWORK, &
        &            DENCOR,DWORKG,DWORK1,DWORK2,DWORK3,DVC, LEXCHG)

         DEALLOCATE(CWGRAD,DWORKG,DWORK1,DWORK2,DWORK3,DVC)               
      ELSEIF (NCDIJ==4) THEN
!-MM- gradient corrections in the noncollinear case are calculated
!     a bit differently than in the collinear case
         ALLOCATE(CWGRAD(GRIDC%MPLWV,NCDIJ), DWORKG(NP1,NCDIJ/2), &
           DWORK1(NP1,NCDIJ),DWORK2(NP1,NCDIJ),DWORK3(NP1,NCDIJ),DVC(NP1))

         CALL FEXCGS_NONCOL_(NCDIJ,GRIDC,LATT_CUR,XCENC,EXC,CVZERO,XCSIF, &
        &            CWGRAD,CHTOT,CWORK, &
        &            CWGRAD,CHTOT,CWORK, &
        &            DENCOR,DWORKG,DWORK1,DWORK2,DWORK3,DVC, LEXCHG)

         DEALLOCATE(CWGRAD,DWORKG,DWORK1,DWORK2,DWORK3,DVC)         
!-MM- end of changes to calculation of gga in noncollinear case
      ENDIF
      
      RETURN
      END


!  ATTENTION ATTENTION ATTENTION ATTENTION ATTENTION ATTENTION
! Mind CWORK and DWORK point actually to the same storagelocation
! similar to an EQUIVALENCE (CWORK(1),DWORK(1))
! the same for  (CWGRAD,DWGRAD) and   (CHTOT,DHTOT)
! so we can interchange both arrays arbitrarily

      SUBROUTINE FEXCGS_(ISPIN,GRIDC,LATT_CUR,XCENC,EXC,CVZERO,XCSIF, &
     &            CWGRAD,CHTOT,CWORK, &
     &            DWGRAD,DHTOT,DWORK, &
     &            DENCOR,DWORKG,DWORK1,DWORK2,DWORK3,DVC, LEXCHG)
      USE prec
      USE lattice
      USE mpimy
      USE mgrid
      USE constant

      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRIDC
      TYPE (latt)        LATT_CUR

      COMPLEX(q) CHTOT(GRIDC%MPLWV,ISPIN),CWORK(GRIDC%MPLWV,ISPIN), &
              CWGRAD(GRIDC%MPLWV,ISPIN)
      COMPLEX(q)   DHTOT(GRIDC%MPLWV,ISPIN),DWORK(GRIDC%MPLWV,ISPIN), &
              DWGRAD(GRIDC%MPLWV,ISPIN)
      COMPLEX(q)   DENCOR(GRIDC%RL%NP)
      REAL(q) XCSIF(3,3)
      REAL(q) DWORKG(GRIDC%RL%NP,ISPIN),DWORK1(GRIDC%RL%NP,ISPIN), &
              DWORK2(GRIDC%RL%NP,ISPIN),DWORK3(GRIDC%RL%NP,ISPIN), &
              DVC(GRIDC%RL%NP)

      REAL(q) :: CHGMIN=1E-10
! set to (1._q,0._q) for error-dumps
      IDUMP=0
! 
! only for the analytical PBE functional we can use a low charge density threshold
! since numerical differentiation leads to errors for the other functionals
!
      IF (LEXCHG==5 .OR. LEXCHG==6) THEN
         CHGMIN=1E-10
      ELSE
         CHGMIN=1E-7
      ENDIF

      RINPL=1._q/GRIDC%NPLWV
!=======================================================================
! First phase: Transform DENCOR (core charge) and
!  CHTOT (pseudo chargedensity) to real space
!=======================================================================
      spin: DO ISP=1,ISPIN
! set CWORK to total real charge in reciprocal space
      DO I=1,GRIDC%RL%NP
         DWORK(I,ISP)=(DENCOR(I)/2+DHTOT(I,ISP))*RINPL/LATT_CUR%OMEGA
      ENDDO
      CALL OPSYNC(DWORK(1,ISP),CWORK(1,ISP),GRIDC%NPLWV)
      CALL FFT3RC(CWORK(1,ISP),GRIDC,-1)
      CALL TRUNC_HIGH_FREQU(LATT_CUR, GRIDC, CWORK(1,ISP))

!=======================================================================
! now calculate the gradient of the chargedensity
!=======================================================================
      DO  I=1,GRIDC%RC%NP
        CWGRAD(I,ISP)=CWORK(I,ISP)
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
         CWORK(I,ISP)=CWORK(I,ISP)*GX*CITPI
       ENDDO
! grad_x in real space:
      CALL SETUNB(CWORK(1,ISP),GRIDC)
      CALL FFT3RC(CWORK(1,ISP),GRIDC,1)
      CALL OPSYNC(DWORK(1,ISP),CWORK(1,ISP),GRIDC%NPLWV)
      DO I=1,GRIDC%RL%NP
        DWORK1(I,ISP)= REAL( DWORK(I,ISP) ,KIND=q)
      ENDDO
!
! y-component:
      DO  I=1,GRIDC%RC%NP
        CWORK(I,ISP)=CWGRAD(I,ISP)
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
         CWORK(I,ISP)=CWORK(I,ISP)*GY*CITPI
      ENDDO
! grad_y in real space:
      CALL SETUNB(CWORK(1,ISP),GRIDC)
      CALL FFT3RC(CWORK(1,ISP),GRIDC,1)
      CALL OPSYNC(DWORK(1,ISP),CWORK(1,ISP),GRIDC%NPLWV)
      DO I=1,GRIDC%RL%NP
        DWORK2(I,ISP)= REAL( DWORK(I,ISP) ,KIND=q)
      ENDDO


! z-component:
      DO  I=1,GRIDC%RC%NP
        CWORK(I,ISP)=CWGRAD(I,ISP)
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
         CWORK(I,ISP)=CWORK(I,ISP)*GZ*CITPI
      ENDDO
! grad_z in real space:
      CALL SETUNB(CWORK(1,ISP),GRIDC)
      CALL FFT3RC(CWORK(1,ISP),GRIDC,1)
      CALL OPSYNC(DWORK(1,ISP),CWORK(1,ISP),GRIDC%NPLWV)
      DO I=1,GRIDC%RL%NP
        DWORK3(I,ISP)= REAL( DWORK(I,ISP) ,KIND=q)
      ENDDO
      ENDDO spin
!=======================================================================
!  grad rho    d    f_xc
! ---------- * ------------      (Phys.Rev.B 50,7 (1994) 4954)
! |grad rho|   d |grad rho|
!
!  MIND: the factor OMEGA is difficult to understand:
!   1/N sum_r energy_density * rho *OMEGA = Energy
!   1/N sum_r energy_density * \bar rho   = Energy (\bar rho=rho*LATT_CUR%OMEGA)
!=======================================================================
      EXC=0
      DO I=1,GRIDC%RL%NP
         RHO1= MAX(REAL((DHTOT(I,1)+DENCOR(I)/2)/LATT_CUR%OMEGA ,KIND=q), CHGMIN)
         RHO2= MAX(REAL((DHTOT(I,2)+DENCOR(I)/2)/LATT_CUR%OMEGA ,KIND=q), CHGMIN)

         ABSNABUP=SQRT(DWORK1(I,1)*DWORK1(I,1)+DWORK2(I,1)*DWORK2(I,1) &
                   +DWORK3(I,1)*DWORK3(I,1))

         ABSNABDW=SQRT(DWORK1(I,2)*DWORK1(I,2)+DWORK2(I,2)*DWORK2(I,2) &
                   +DWORK3(I,2)*DWORK3(I,2))

         ABSNAB= (DWORK1(I,1)+DWORK1(I,2))*(DWORK1(I,1)+DWORK1(I,2))+ &
                 (DWORK2(I,1)+DWORK2(I,2))*(DWORK2(I,1)+DWORK2(I,2))+ &
                 (DWORK3(I,1)+DWORK3(I,2))*(DWORK3(I,1)+DWORK3(I,2))
         ABSNAB=SQRT(ABSNAB)
!
! presently VASP coverges very badly if correlation_ABS_DRHOUP_ABS_DRHOD
! is not set
! it is probably the exchange potential, which is partly canceled by
! the correlation energy if ABSNAB=ABSNABUP+ABSNABDW
! so let us stick to ABSNAB=ABSNABUP+ABSNABDW for the time beeing

!#define  correlation_ABS_DRHOUP_ABS_DRHOD
!
         CALL GGASPIN(RHO1*AUTOA3,RHO2*AUTOA3, &
     &               ABSNABUP*AUTOA4, ABSNABDW*AUTOA4,ABSNAB*AUTOA4, &
     &           EXCL,DEXC1,DEXC2,DVXC1,DVXC2,DVC_,LEXCHG,.FALSE.)
!
         RHO=RHO1+RHO2

         EXC=EXC+EXCL*RHO*RYTOEV*LATT_CUR%OMEGA
         DVXC1=DVXC1*RYTOEV*AUTOA
         DVXC2=DVXC2*RYTOEV*AUTOA
         DVC_ =DVC_ *RYTOEV*AUTOA
!
!   store d f/ d (|d rho| ) / |d rho|  in DWORK
!
         DWORK(I,1)  = DVXC1 / MAX(ABSNABUP,1.E-10_q)
         DWORK(I,2)  = DVXC2 / MAX(ABSNABDW,1.E-10_q)
         DVC(I)      = DVC_  / MAX(ABSNAB,1.E-10_q)
!
!   store d f/ d rho  in DWORKG
!
         DWORKG(I,1) = DEXC1*RYTOEV
         DWORKG(I,2) = DEXC2*RYTOEV
      ENDDO
!=======================================================================
! gradient terms in stress tensor
!          d    f_xc     grad rho  x grad rho
! sum_r   ------------   --------------------- * LATT_CUR%OMEGA
!         d |grad rho|        |grad rho|
!=======================================================================
!
!    verify that this sum is correct also when ISPIN=2
!
      SIF11=0
      SIF22=0
      SIF33=0
      SIF12=0
      SIF23=0
      SIF31=0
      DO ISP=1,ISPIN
      DO I=1,GRIDC%RL%NP
        SIF11=SIF11+DWORK1(I,ISP)*DWORK1(I,ISP)*DWORK(I,ISP)
        SIF22=SIF22+DWORK2(I,ISP)*DWORK2(I,ISP)*DWORK(I,ISP)
        SIF33=SIF33+DWORK3(I,ISP)*DWORK3(I,ISP)*DWORK(I,ISP)
        SIF12=SIF12+DWORK1(I,ISP)*DWORK2(I,ISP)*DWORK(I,ISP)
        SIF23=SIF23+DWORK2(I,ISP)*DWORK3(I,ISP)*DWORK(I,ISP)
        SIF31=SIF31+DWORK3(I,ISP)*DWORK1(I,ISP)*DWORK(I,ISP)

        SIF11=SIF11+DWORK1(I,ISP)*(DWORK1(I,1)+DWORK1(I,2))*DVC(I)
        SIF22=SIF22+DWORK2(I,ISP)*(DWORK2(I,1)+DWORK2(I,2))*DVC(I)
        SIF33=SIF33+DWORK3(I,ISP)*(DWORK3(I,1)+DWORK3(I,2))*DVC(I)
        SIF12=SIF12+DWORK1(I,ISP)*(DWORK2(I,1)+DWORK2(I,2))*DVC(I)
        SIF23=SIF23+DWORK2(I,ISP)*(DWORK3(I,1)+DWORK3(I,2))*DVC(I)
        SIF31=SIF31+DWORK3(I,ISP)*(DWORK1(I,1)+DWORK1(I,2))*DVC(I)
      ENDDO
      ENDDO
      SIF11=SIF11*RINPL*LATT_CUR%OMEGA
      SIF22=SIF22*RINPL*LATT_CUR%OMEGA
      SIF33=SIF33*RINPL*LATT_CUR%OMEGA
      SIF12=SIF12*RINPL*LATT_CUR%OMEGA
      SIF23=SIF23*RINPL*LATT_CUR%OMEGA
      SIF31=SIF31*RINPL*LATT_CUR%OMEGA

!=======================================================================
! calculate 
!              d    f_xc     grad rho
!        div  (------------  --------  ) 
!              d |grad rho| |grad rho|
!
! in reciprocal space
!=======================================================================

      DO I=1,GRIDC%RL%NP
         ANAB1U= DWORK1(I,1)
         ANAB2U= DWORK2(I,1)
         ANAB3U= DWORK3(I,1)
         ANAB1D= DWORK1(I,2)
         ANAB2D= DWORK2(I,2)
         ANAB3D= DWORK3(I,2)

         DWORK1(I,1) = ANAB1U* DWORK(I,1) + (ANAB1U+ANAB1D) * DVC(I)
         DWORK2(I,1) = ANAB2U* DWORK(I,1) + (ANAB2U+ANAB2D) * DVC(I)
         DWORK3(I,1) = ANAB3U* DWORK(I,1) + (ANAB3U+ANAB3D) * DVC(I)

         DWORK1(I,2) = ANAB1D* DWORK(I,2) + (ANAB1U+ANAB1D) * DVC(I)
         DWORK2(I,2) = ANAB2D* DWORK(I,2) + (ANAB2U+ANAB2D) * DVC(I)
         DWORK3(I,2) = ANAB3D* DWORK(I,2) + (ANAB3U+ANAB3D) * DVC(I)
      ENDDO


      spin2: DO ISP=1,ISPIN
! x-component:
!DIR$ IVDEP
!$DIR FORCE_VECTOR
!OCL NOVREC
      DO I=1,GRIDC%RL%NP
        DWORK(I,ISP) = DWORK1(I,ISP)
      ENDDO
      CALL OPSYNC(DWORK(1,ISP),CWORK(1,ISP),GRIDC%NPLWV)
      CALL FFT3RC(CWORK(1,ISP),GRIDC,-1)
      CALL TRUNC_HIGH_FREQU(LATT_CUR, GRIDC, CWORK(1,ISP))

      DO I=1,GRIDC%RC%NP
         N1= MOD((I-1),GRIDC%RC%NROW) +1
         NC= (I-1)/GRIDC%RC%NROW+1
         N2= GRIDC%RC%I2(NC)
         N3= GRIDC%RC%I3(NC)

         GX=(GRIDC%LPCTX(N1)*LATT_CUR%B(1,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(1,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(1,3))
         CWGRAD(I,ISP)=CWORK(I,ISP)*GX*CITPI
      ENDDO

! y-component:
      DO I=1,GRIDC%RL%NP
        DWORK(I,ISP) = DWORK2(I,ISP)
      ENDDO
      CALL OPSYNC(DWORK(1,ISP),CWORK(1,ISP),GRIDC%NPLWV)
      CALL FFT3RC(CWORK(1,ISP),GRIDC,-1)
      CALL TRUNC_HIGH_FREQU(LATT_CUR, GRIDC, CWORK(1,ISP))

      DO I=1,GRIDC%RC%NP
         N1= MOD((I-1),GRIDC%RC%NROW) +1
         NC= (I-1)/GRIDC%RC%NROW+1
         N2= GRIDC%RC%I2(NC)
         N3= GRIDC%RC%I3(NC)

         GY=(GRIDC%LPCTX(N1)*LATT_CUR%B(2,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(2,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(2,3))
         CWGRAD(I,ISP)=CWGRAD(I,ISP)+CWORK(I,ISP)*GY*CITPI
      ENDDO

! z-component:
      DO I=1,GRIDC%RL%NP
        DWORK(I,ISP) = DWORK3(I,ISP)
      ENDDO
      CALL OPSYNC(DWORK(1,ISP),CWORK(1,ISP),GRIDC%NPLWV)
      CALL FFT3RC(CWORK(1,ISP),GRIDC,-1)
      CALL TRUNC_HIGH_FREQU(LATT_CUR, GRIDC, CWORK(1,ISP))

      DO I=1,GRIDC%RC%NP
         N1= MOD((I-1),GRIDC%RC%NROW) +1
         NC= (I-1)/GRIDC%RC%NROW+1
         N2= GRIDC%RC%I2(NC)
         N3= GRIDC%RC%I3(NC)

         GZ=(GRIDC%LPCTX(N1)*LATT_CUR%B(3,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(3,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(3,3))
         CWGRAD(I,ISP)=CWGRAD(I,ISP)+CWORK(I,ISP)*GZ*CITPI
      ENDDO

      CALL SETUNB(CWGRAD(1,ISP),GRIDC)
      CALL FFT3RC(CWGRAD(1,ISP),GRIDC,1)
      CALL OPSYNC(CWGRAD(1,ISP),DWGRAD(1,ISP),GRIDC%NPLWV)
!
      ENDDO spin2


!
!=======================================================================
! Now prepare the rest:
! (store rho in DWORK3 and quantity of above in DWORK1)
!=======================================================================
      XCENC=0._q
      CVZERO=0._q
      XCENCC=0._q

      DO I=1,GRIDC%RL%NP
         RHO1= MAX(REAL((DHTOT(I,1)+DENCOR(I)/2)/LATT_CUR%OMEGA ,KIND=q), CHGMIN)
         RHO2= MAX(REAL((DHTOT(I,2)+DENCOR(I)/2)/LATT_CUR%OMEGA ,KIND=q), CHGMIN)

         VXC1=DWORKG(I,1)- REAL( DWGRAD(I,1) ,KIND=q) *RINPL
         VXC2=DWORKG(I,2)- REAL( DWGRAD(I,2) ,KIND=q) *RINPL
         DWORK(I,1)=VXC1
         DWORK(I,2)=VXC2
         VXC = 0.5_q*(VXC1+VXC2)
         CVZERO=CVZERO+VXC
         XCENCC=XCENCC-VXC1*RHO1*LATT_CUR%OMEGA-VXC2*RHO2*LATT_CUR%OMEGA
         XCENC=XCENC  -VXC1* REAL( DHTOT(I,1) ,KIND=q) -VXC2* REAL( DHTOT(I,2) ,KIND=q)
      ENDDO

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
      END SUBROUTINE


!  ATTENTION ATTENTION ATTENTION ATTENTION ATTENTION ATTENTION
! Mind CWORK and DWORK point actually to the same storagelocation
! similar to an EQUIVALENCE (CWORK(1),DWORK(1))
! the same for  (CWGRAD,DWGRAD) and   (CHTOT,DHTOT)
! so we can interchange both arrays arbitrarily

      SUBROUTINE FEXCGS_NONCOL_(NCDIJ,GRIDC,LATT_CUR,XCENC,EXC,CVZERO,XCSIF, &
     &            CWGRAD,CHTOT,CWORK, &
     &            DWGRAD,DHTOT,DWORK, &
     &            DENCOR,DWORKG,DWORK1,DWORK2,DWORK3,DVC, LEXCHG)
      USE prec
      USE lattice
      USE mpimy
      USE mgrid
      USE constant

      IMPLICIT COMPLEX(q) (C)

      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRIDC
      TYPE (latt)        LATT_CUR

      COMPLEX(q) CHTOT(GRIDC%MPLWV,NCDIJ),CWORK(GRIDC%MPLWV,NCDIJ), &
              CWGRAD(GRIDC%MPLWV,NCDIJ)
      COMPLEX(q)   DHTOT(GRIDC%MPLWV,NCDIJ),DWORK(GRIDC%MPLWV,NCDIJ), &
              DWGRAD(GRIDC%MPLWV,NCDIJ)
      COMPLEX(q)   DENCOR(GRIDC%RL%NP)
      REAL(q) XCSIF(3,3)
      REAL(q) DWORKG(GRIDC%RL%NP,NCDIJ/2),DWORK1(GRIDC%RL%NP,NCDIJ), &
              DWORK2(GRIDC%RL%NP,NCDIJ),DWORK3(GRIDC%RL%NP,NCDIJ), &
              DVC(GRIDC%RL%NP)
      REAL(q) MAG_NORM,NABMAG(3)
      REAL(q) :: CHGMIN=1E-10

! set to (1._q,0._q) for error-dumps
      IDUMP=0
! 
! only for the analytical PBE functional we can use a low charge density threshold
! since numerical differentiation leads to errors for the other functionals
!
      IF (LEXCHG==5 .OR. LEXCHG==6) THEN
         CHGMIN=1E-10
      ELSE
         CHGMIN=1E-7
      ENDIF

      RINPL=1._q/GRIDC%NPLWV
!=======================================================================
      spin: DO ISP=1,NCDIJ
      IF (ISP==1) THEN
! Add core charge density to psuedo charge density
         DO I=1,GRIDC%RL%NP
            DWORK(I,ISP)=(DENCOR(I)+DHTOT(I,ISP))*RINPL/LATT_CUR%OMEGA
         ENDDO
      ELSE
         DO I=1,GRIDC%RL%NP
            DWORK(I,ISP)=DHTOT(I,ISP)*RINPL/LATT_CUR%OMEGA
         ENDDO
      ENDIF
      CALL OPSYNC(DWORK(1,ISP),CWORK(1,ISP),GRIDC%NPLWV)
      CALL FFT3RC(CWORK(1,ISP),GRIDC,-1)
      CALL TRUNC_HIGH_FREQU(LATT_CUR, GRIDC, CWORK(1,ISP))
!=======================================================================
! now calculate the gradients
!=======================================================================
! x-component:
      DO I=1,GRIDC%RC%NP
         CWGRAD(I,ISP)=CWORK(I,ISP)
      ENDDO
!DIR$ IVDEP
!$DIR FORCE_VECTOR
!OCL NOVREC
      DO I=1,GRIDC%RC%NP
         N1= MOD((I-1),GRIDC%RC%NROW) +1
         NC= (I-1)/GRIDC%RC%NROW+1
         N2= GRIDC%RC%I2(NC)
         N3= GRIDC%RC%I3(NC)
         GX=(GRIDC%LPCTX(N1)*LATT_CUR%B(1,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(1,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(1,3))
         CWORK(I,ISP)=CWORK(I,ISP)*GX*CITPI
      ENDDO
! grad_x in real space:
      CALL SETUNB(CWORK(1,ISP),GRIDC)
      CALL FFT3RC(CWORK(1,ISP),GRIDC,1)
      CALL OPSYNC(DWORK(1,ISP),CWORK(1,ISP),GRIDC%NPLWV)
      DO I=1,GRIDC%RL%NP
         DWORK1(I,ISP)= REAL( DWORK(I,ISP) ,KIND=q)
      ENDDO
!
! y-component:
      DO I=1,GRIDC%RC%NP
         CWORK(I,ISP)=CWGRAD(I,ISP)
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
         CWORK(I,ISP)=CWORK(I,ISP)*GY*CITPI
      ENDDO
! grad_y in real space:
      CALL SETUNB(CWORK(1,ISP),GRIDC)
      CALL FFT3RC(CWORK(1,ISP),GRIDC,1)
      CALL OPSYNC(DWORK(1,ISP),CWORK(1,ISP),GRIDC%NPLWV)
      DO I=1,GRIDC%RL%NP
         DWORK2(I,ISP)= REAL( DWORK(I,ISP) ,KIND=q)
      ENDDO
!
! z-component:
      DO I=1,GRIDC%RC%NP
         CWORK(I,ISP)=CWGRAD(I,ISP)
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
         CWORK(I,ISP)=CWORK(I,ISP)*GZ*CITPI
      ENDDO
! grad_z in real space:
      CALL SETUNB(CWORK(1,ISP),GRIDC)
      CALL FFT3RC(CWORK(1,ISP),GRIDC,1)
      CALL OPSYNC(DWORK(1,ISP),CWORK(1,ISP),GRIDC%NPLWV)
      DO I=1,GRIDC%RL%NP
         DWORK3(I,ISP)= REAL( DWORK(I,ISP) ,KIND=q)
      ENDDO
      ENDDO spin

!=======================================================================
!  grad rho    d    f_xc
! ---------- * ------------      (Phys.Rev.B 50,7 (1994) 4954)
! |grad rho|   d |grad rho|
!
!  MIND: the factor OMEGA is difficult to understand:
!   1/N sum_r energy_density * rho *OMEGA = Energy
!   1/N sum_r energy_density * \bar rho   = Energy (\bar rho=rho*LATT_CUR%OMEGA)
!=======================================================================
      EXC=0
      DO I=1,GRIDC%RL%NP
!
! | m |
!
         MAG_NORM=MAX(SQRT(ABS( & 
        &    DHTOT(I,2)*DHTOT(I,2)+DHTOT(I,3)*DHTOT(I,3)+DHTOT(I,4)*DHTOT(I,4))), 1E-10_q)        
!
! RHO1: \rho_up   = ( \rho_tot + | m | )/2
! RHO2: \rho_down = ( \rho_tot - | m | )/2
!
         RHO1= MAX(REAL((DHTOT(I,1)+DENCOR(I)+MAG_NORM)/2/LATT_CUR%OMEGA ,KIND=q), CHGMIN)
         RHO2= MAX(REAL((DHTOT(I,1)+DENCOR(I)-MAG_NORM)/2/LATT_CUR%OMEGA ,KIND=q), CHGMIN)
!
! \nabla | m |
!
         NABMAG(1)=(DWORK1(I,2)*DHTOT(I,2)+DWORK1(I,3)*DHTOT(I,3)+DWORK1(I,4)*DHTOT(I,4))/MAG_NORM
         NABMAG(2)=(DWORK2(I,2)*DHTOT(I,2)+DWORK2(I,3)*DHTOT(I,3)+DWORK2(I,4)*DHTOT(I,4))/MAG_NORM
         NABMAG(3)=(DWORK3(I,2)*DHTOT(I,2)+DWORK3(I,3)*DHTOT(I,3)+DWORK3(I,4)*DHTOT(I,4))/MAG_NORM
!
! | ( \nabla \rho + \nabla | m | )/2 |
!
         ABSNABUP=SQRT( &
        &            (DWORK1(I,1)+NABMAG(1))*(DWORK1(I,1)+NABMAG(1)) + &
        &             (DWORK2(I,1)+NABMAG(2))*(DWORK2(I,1)+NABMAG(2)) + &
        &              (DWORK3(I,1)+NABMAG(3))*(DWORK3(I,1)+NABMAG(3)) ) / 2
!
! | ( \nabla \rho - \nabla | m | )/2 |
!
         ABSNABDW=SQRT( &
        &            (DWORK1(I,1)-NABMAG(1))*(DWORK1(I,1)-NABMAG(1)) + &
        &             (DWORK2(I,1)-NABMAG(2))*(DWORK2(I,1)-NABMAG(2)) + &
        &              (DWORK3(I,1)-NABMAG(3))*(DWORK3(I,1)-NABMAG(3)) ) / 2
!
! | \nabla \rho |
!
         ABSNAB=SQRT(DWORK1(I,1)*DWORK1(I,1)+DWORK2(I,1)*DWORK2(I,1)+DWORK3(I,1)*DWORK3(I,1))
!
! Refill DWORK[1..3](:,2) with ( \nabla \rho - \nabla | m | )/2
!
         DWORK1(I,2)=(DWORK1(I,1)-NABMAG(1))/2
         DWORK2(I,2)=(DWORK2(I,1)-NABMAG(2))/2
         DWORK3(I,2)=(DWORK3(I,1)-NABMAG(3))/2
!
! Refill DWORK[1..3](:,1) with ( \nabla \rho + \nabla | m | )/2
!
         DWORK1(I,1)=(DWORK1(I,1)+NABMAG(1))/2
         DWORK2(I,1)=(DWORK2(I,1)+NABMAG(2))/2
         DWORK3(I,1)=(DWORK3(I,1)+NABMAG(3))/2
!
! presently VASP converges very badly if correlation_ABS_DRHOUP_ABS_DRHOD
! is not set
! it is probably the exchange potential, which is partly canceled by
! the correlation energy if ABSNAB=ABSNABUP+ABSNABDW
! so let us stick to ABSNAB=ABSNABUP+ABSNABDW for the time beeing

!#define  correlation_ABS_DRHOUP_ABS_DRHOD

         CALL GGASPIN(RHO1*AUTOA3,RHO2*AUTOA3, &
     &               ABSNABUP*AUTOA4, ABSNABDW*AUTOA4,ABSNAB*AUTOA4, &
     &           EXCL,DEXC1,DEXC2,DVXC1,DVXC2,DVC_,LEXCHG,.FALSE.)

         RHO=RHO1+RHO2

         EXC=EXC+EXCL*RHO*RYTOEV*LATT_CUR%OMEGA
         DVXC1=DVXC1*RYTOEV*AUTOA
         DVXC2=DVXC2*RYTOEV*AUTOA
         DVC_ =DVC_ *RYTOEV*AUTOA
!
!   store d f/ d (|d rho| ) / |d rho|  in DWORK
!
         DWORK(I,1)  = DVXC1 / MAX(ABSNABUP,1.E-10_q)
         DWORK(I,2)  = DVXC2 / MAX(ABSNABDW,1.E-10_q)
         DVC(I)      = DVC_  / MAX(ABSNAB,1.E-10_q)
!
!   store d f/ d rho  in DWORKG
!
         DWORKG(I,1) = DEXC1*RYTOEV
         DWORKG(I,2) = DEXC2*RYTOEV
      ENDDO
!      stop
!=======================================================================
! gradient terms in stress tensor
!          d    f_xc     grad rho  x grad rho
! sum_r   ------------   --------------------- * LATT_CUR%OMEGA
!         d |grad rho|        |grad rho|
!=======================================================================
!
!    verify that this sum is correct also when ISPIN=2
!
      SIF11=0
      SIF22=0
      SIF33=0
      SIF12=0
      SIF23=0
      SIF31=0
      DO ISP=1,2
      DO I=1,GRIDC%RL%NP
        SIF11=SIF11+DWORK1(I,ISP)*DWORK1(I,ISP)*DWORK(I,ISP)
        SIF22=SIF22+DWORK2(I,ISP)*DWORK2(I,ISP)*DWORK(I,ISP)
        SIF33=SIF33+DWORK3(I,ISP)*DWORK3(I,ISP)*DWORK(I,ISP)
        SIF12=SIF12+DWORK1(I,ISP)*DWORK2(I,ISP)*DWORK(I,ISP)
        SIF23=SIF23+DWORK2(I,ISP)*DWORK3(I,ISP)*DWORK(I,ISP)
        SIF31=SIF31+DWORK3(I,ISP)*DWORK1(I,ISP)*DWORK(I,ISP)

        SIF11=SIF11+DWORK1(I,ISP)*(DWORK1(I,1)+DWORK1(I,2))*DVC(I)
        SIF22=SIF22+DWORK2(I,ISP)*(DWORK2(I,1)+DWORK2(I,2))*DVC(I)
        SIF33=SIF33+DWORK3(I,ISP)*(DWORK3(I,1)+DWORK3(I,2))*DVC(I)
        SIF12=SIF12+DWORK1(I,ISP)*(DWORK2(I,1)+DWORK2(I,2))*DVC(I)
        SIF23=SIF23+DWORK2(I,ISP)*(DWORK3(I,1)+DWORK3(I,2))*DVC(I)
        SIF31=SIF31+DWORK3(I,ISP)*(DWORK1(I,1)+DWORK1(I,2))*DVC(I)
      ENDDO
      ENDDO
      SIF11=SIF11*RINPL*LATT_CUR%OMEGA
      SIF22=SIF22*RINPL*LATT_CUR%OMEGA
      SIF33=SIF33*RINPL*LATT_CUR%OMEGA
      SIF12=SIF12*RINPL*LATT_CUR%OMEGA
      SIF23=SIF23*RINPL*LATT_CUR%OMEGA
      SIF31=SIF31*RINPL*LATT_CUR%OMEGA

!=======================================================================
! calculate 
!              d    f_xc     grad rho
!        div  (------------  --------  ) 
!              d |grad rho| |grad rho|
!
! in reciprocal space
!=======================================================================

      DO I=1,GRIDC%RL%NP
         ANAB1U= DWORK1(I,1)
         ANAB2U= DWORK2(I,1)
         ANAB3U= DWORK3(I,1)
         ANAB1D= DWORK1(I,2)
         ANAB2D= DWORK2(I,2)
         ANAB3D= DWORK3(I,2)

         DWORK1(I,1) = ANAB1U* DWORK(I,1) + (ANAB1U+ANAB1D) * DVC(I)
         DWORK2(I,1) = ANAB2U* DWORK(I,1) + (ANAB2U+ANAB2D) * DVC(I)
         DWORK3(I,1) = ANAB3U* DWORK(I,1) + (ANAB3U+ANAB3D) * DVC(I)

         DWORK1(I,2) = ANAB1D* DWORK(I,2) + (ANAB1U+ANAB1D) * DVC(I)
         DWORK2(I,2) = ANAB2D* DWORK(I,2) + (ANAB2U+ANAB2D) * DVC(I)
         DWORK3(I,2) = ANAB3D* DWORK(I,2) + (ANAB3U+ANAB3D) * DVC(I)
      ENDDO

      spin2: DO ISP=1,2
! x-component:
!DIR$ IVDEP
!$DIR FORCE_VECTOR
!OCL NOVREC
      DO I=1,GRIDC%RL%NP
        DWORK(I,ISP) = DWORK1(I,ISP)
      ENDDO
      CALL OPSYNC(DWORK(1,ISP),CWORK(1,ISP),GRIDC%NPLWV)
      CALL FFT3RC(CWORK(1,ISP),GRIDC,-1)
      CALL TRUNC_HIGH_FREQU(LATT_CUR, GRIDC, CWORK(1,ISP))

      DO I=1,GRIDC%RC%NP
         N1= MOD((I-1),GRIDC%RC%NROW) +1
         NC= (I-1)/GRIDC%RC%NROW+1
         N2= GRIDC%RC%I2(NC)
         N3= GRIDC%RC%I3(NC)

         GX=(GRIDC%LPCTX(N1)*LATT_CUR%B(1,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(1,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(1,3))
         CWGRAD(I,ISP)=CWORK(I,ISP)*GX*CITPI
      ENDDO

! y-component:
      DO I=1,GRIDC%RL%NP
        DWORK(I,ISP) = DWORK2(I,ISP)
      ENDDO
      CALL OPSYNC(DWORK(1,ISP),CWORK(1,ISP),GRIDC%NPLWV)
      CALL FFT3RC(CWORK(1,ISP),GRIDC,-1)
      CALL TRUNC_HIGH_FREQU(LATT_CUR, GRIDC, CWORK(1,ISP))

      DO I=1,GRIDC%RC%NP
         N1= MOD((I-1),GRIDC%RC%NROW) +1
         NC= (I-1)/GRIDC%RC%NROW+1
         N2= GRIDC%RC%I2(NC)
         N3= GRIDC%RC%I3(NC)

         GY=(GRIDC%LPCTX(N1)*LATT_CUR%B(2,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(2,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(2,3))
         CWGRAD(I,ISP)=CWGRAD(I,ISP)+CWORK(I,ISP)*GY*CITPI
      ENDDO

! z-component:
      DO I=1,GRIDC%RL%NP
        DWORK(I,ISP) = DWORK3(I,ISP)
      ENDDO
      CALL OPSYNC(DWORK(1,ISP),CWORK(1,ISP),GRIDC%NPLWV)
      CALL FFT3RC(CWORK(1,ISP),GRIDC,-1)
      CALL TRUNC_HIGH_FREQU(LATT_CUR, GRIDC, CWORK(1,ISP))

      DO I=1,GRIDC%RC%NP
         N1= MOD((I-1),GRIDC%RC%NROW) +1
         NC= (I-1)/GRIDC%RC%NROW+1
         N2= GRIDC%RC%I2(NC)
         N3= GRIDC%RC%I3(NC)

         GZ=(GRIDC%LPCTX(N1)*LATT_CUR%B(3,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(3,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(3,3))
         CWGRAD(I,ISP)=CWGRAD(I,ISP)+CWORK(I,ISP)*GZ*CITPI
      ENDDO

      CALL SETUNB(CWGRAD(1,ISP),GRIDC)
      CALL FFT3RC(CWGRAD(1,ISP),GRIDC,1)
      CALL OPSYNC(CWGRAD(1,ISP),DWGRAD(1,ISP),GRIDC%NPLWV)
!
      ENDDO spin2


!
!=======================================================================
! Now prepare the rest:
! (store rho in DWORK3 and quantity of above in DWORK1)
!=======================================================================
      XCENC=0._q
      CVZERO=0._q
      XCENCC=0._q

      DO I=1,GRIDC%RL%NP
!
! | m |
!
         MAG_NORM=SQRT(ABS(DHTOT(I,2)*DHTOT(I,2)+ DHTOT(I,3)*DHTOT(I,3) + DHTOT(I,4)*DHTOT(I,4)))
!
! RHO1: \rho_up   = ( \rho_tot + | m | )/2
! RHO2: \rho_down = ( \rho_tot - | m | )/2
!
         RHO1= MAX(REAL((DHTOT(I,1)+DENCOR(I)+MAG_NORM)/2/LATT_CUR%OMEGA ,KIND=q), CHGMIN)
         RHO2= MAX(REAL((DHTOT(I,1)+DENCOR(I)-MAG_NORM)/2/LATT_CUR%OMEGA ,KIND=q), CHGMIN)

         VXC1=DWORKG(I,1)- REAL( DWGRAD(I,1) ,KIND=q) *RINPL
         VXC2=DWORKG(I,2)- REAL( DWGRAD(I,2) ,KIND=q) *RINPL
         DWORK(I,1)=VXC1
         DWORK(I,2)=VXC2
         VXC = 0.5_q*(VXC1+VXC2)
         CVZERO=CVZERO+VXC
         XCENCC=XCENCC-VXC1*RHO1*LATT_CUR%OMEGA-VXC2*RHO2*LATT_CUR%OMEGA
         XCENC=XCENC  -VXC1* REAL( (DHTOT(I,1)+MAG_NORM)/2 ,KIND=q) &
        &               -VXC2* REAL( (DHTOT(I,1)-MAG_NORM)/2 ,KIND=q)
      ENDDO

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
      END SUBROUTINE


!************************ SUBROUTINE GGASPIN *****************************
!
!  switch between different GGAs
!  presently only PBE and RPBE are implemented efficiently
!  (i.e. d exc / d rho and  d exc / d | grad rho | are calculated
!  directly 
!  for other GGA functional finite differences are used to calculate
!  the required derivatives
!  LLDA allows to include the LDA contribution directly in this 
!  routine.
!  This only works for PBE and RPBE and fails in all other cases
!
!***********************************************************************

      SUBROUTINE GGASPIN(D1,D2,DD1,DD2,DDA,EXC, &
         excd1,excd2,excq1,excq2,ecq,LEXCHG,LLDA)

!     D1   density up
!     D2   density down
!     DD1  |gradient of density up|
!     DD2  |gradient of density down|
!     DDA  |gradient of the total density|
!     LLDA add lda contributions

      USE prec
      USE constant
      IMPLICIT REAL(q) (A-H,O-Z)
!
      LOGICAL LLDA
      PARAMETER (THRD=1._q/3._q)

      ecq=0

      IF (LEXCHG==4) THEN
         CALL GGAXCS(D1,D2,DD1,DD2,EXC,EXCD1,EXCD2,EXCQ1,EXCQ2)
!
         D = D1 + D2
         EXC    = 2._Q* EXC / D
         EXCD1  = 2._Q* EXCD1
         EXCD2  = 2._Q* EXCD2
         EXCQ1  = 2._Q* EXCQ1
         EXCQ2  = 2._Q* EXCQ2
      ELSE IF (LEXCHG==5 .OR. LEXCHG==6) THEN
        IF (LEXCHG==5) THEN
           ukfactor=1.0_q
        ELSE
           ukfactor=0.0_q
        ENDIF

        D=2*D1


        DTHRD=exp(log(D)*THRD)
        FK=(3._q*PI*PI)**THRD*DTHRD
        IF(D>1.E-10_q)THEN

           S=DD1/(D*FK)

        ELSE
           S=0.0_q
        ENDIF

        CALL EXCHPBE(D,DTHRD,S,EXLDA1,EXC1,EXDLDA1,EXCD1,EXCQ1, &
             ukfactor)
!        WRITE(*,'(10F14.7)') D,S,(EXC1-EXLDA1)/D,EXCD1-EXDLDA1,EXCQ1
!        EXLDA1=0; EXC1=0; EXDLDA1=0; EXCD1=0; EXCQ1=0

        D=2*D2

        DTHRD=exp(log(D)*THRD)
        FK=(3._q*PI*PI)**THRD*DTHRD
        IF(D>1.E-10_q)THEN
           S=DD2/(D*FK)
        ELSE
           S=0.0_q
        ENDIF

        CALL EXCHPBE(D,DTHRD,S,EXLDA2,EXC2,EXDLDA2,EXCD2,EXCQ2, &
             ukfactor)
!        EXLDA2=0; EXC2=0; EXDLDA2=0; EXCD2=0; EXCQ2=0
!        WRITE(*,'(10F14.7)') D,S,(EXC2-EXLDA2)/D,EXCD2-EXDLDA2,EXCQ2

        D=D1+D2
        DTHRD=exp(log(D)*THRD)
        RS=(0.75_q/PI)**THRD/DTHRD
        ZETA=(D1-D2)/D
        ZETA=MIN(MAX(ZETA,-0.9999999999999_q),0.9999999999999_q)

        FK=(3._q*PI*PI)**THRD*DTHRD
        SK = SQRT(4.0_q*FK/PI)

        G = (exp((2*THRD)*log(1._q+ZETA)) &
                +exp((2*THRD)*log(1._q-ZETA)))/2._q
        T = DDA/(D*2._q*SK*G)

        CALL corpbe(RS,ZETA,ECLDA,ECD1LDA,ECD2LDA,G,SK, &
                     T,EC,ECD1,ECD2,ECQ,.TRUE.)
!        ECLDA=0 ; ECD1LDA=0 ; ECD2LDA=0; EC=0 ; ECD1=0 ; ECD2=0 ; ECQ=0 

!        WRITE(*,'(10F14.7)') D,RS,SK,T,EC,ECD1,ECD2,ECQ

        IF (LLDA) THEN
          EXC  =EC  +(EXC1+EXC2)/(2*D)+ECLDA
          EXCD1=ECD1 +EXCD1+ECD1LDA
          EXCD2=ECD2 +EXCD2+ECD2LDA
        ELSE
          EXC  =EC  +(EXC1-EXLDA1+EXC2-EXLDA2)/(2*D)
          EXCD1=ECD1 +EXCD1-EXDLDA1
          EXCD2=ECD2 +EXCD2-EXDLDA2
        ENDIF
        ! Hartree -> Rydberg conversion
        EXC    = 2._q*EXC
        EXCD1  = 2._q* EXCD1
        EXCD2  = 2._q* EXCD2
        EXCQ1  = 2._q* EXCQ1
        EXCQ2  = 2._q* EXCQ2
        ECQ    = 2._q* ECQ
      ELSE
         WRITE(*,*) 'internal ERROR GGASPIN: Wrong LEXCHG, scheme not implemented!'
         STOP
      ENDIF

      RETURN
      END

!************************ SUBROUTINE GGAXCS *****************************
!
!  calculates the PW 91 xc-functional for spin-dependent potentials
!  using the algorithm of J.A. White and D.M. Bird  
!    (Phys.Rev.B 50,7 (1994) 4954)
!  presently finite differences are used to calculate
!  the required derivatives
!
!***********************************************************************

      SUBROUTINE GGAXCS(ro1,ro2,q1,q2,exc, &
     &   excd1,excd2,excq1,excq2)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      REAL(q),PARAMETER :: PI=3.1415926536_q
      REAL(q),PARAMETER :: THRD=0.333333333333333_q

      ro   = ro1 + ro2
      exc  = fun(ro1,ro2,q1,q2)

!*****************************************
!  The derivatives of f_xc necessary to
!  construct spin up and down potentials are
!  computed numerically
!*****************************************
      excd1 = dfro1(ro1,ro2,q1,q2)
      excd2 = dfro2(ro1,ro2,q1,q2)
      excq1 = dfq1 (ro1,ro2,q1,q2)
      excq2 = dfq2 (ro1,ro2,q1,q2)

      RETURN
      END


      function fun(ro1,ro2,q1,q2)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)

!**************************************
!     f_xc PW91 construction
!**************************************
      ro   = ro1 + ro2
      qq   = q1  + q2
      xi   = (ro1 - ro2)/ro
!
      call expw(ro1,q1,ex1)
      call expw(ro2,q2,ex2)
      call cpw(xi,ro,qq,ec)
      despin=0
!     call cpwsp(ro1,ro2,q1,q2,despin)
!
      fx1  = ro1*ex1
      fx2  = ro2*ex2
      fcg  = ro*ec
      fsg  = ro*despin
      fun  = fx1 + fx2 + fcg + fsg
!
      return
      end

      function dfq1(ro1,ro2,q1,q2)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)

      parameter (delta=1E-2_q)
      eps = min(delta,abs(0.001_q*q1))+1E-15_q
      x1  = q1 - eps
      x3  = q1 + eps
      f1  = fun(ro1,ro2,x1,q2)
      f3  = fun(ro1,ro2,x3,q2)
      dfq1 = (f3-f1) / (2*eps)
      return
      end

      function dfq2(ro1,ro2,q1,q2)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)

      parameter (delta=1E-2_q)
      eps = min(delta,abs(0.001_q*q2))+1E-15_q
      x1  = q2 - eps
      x3  = q2 + eps
      f1  = fun(ro1,ro2,q1,x1)
      f3  = fun(ro1,ro2,q1,x3)
      dfq2 = (f3-f1) / (2*eps)
      return
      end

      function dfro1(ro1,ro2,q1,q2)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)

      parameter (delta=1E-4_q)
      eps = min(delta,abs(0.001_q*ro1))+1E-15_q
      x1  = ro1 - eps
      x3  = ro1 + eps
      f1  = fun(x1,ro2,q1,q2)
      f3  = fun(x3,ro2,q1,q2)
      dfro1 = (f3-f1) / (2*eps)
      return
      end

      function dfro2(ro1,ro2,q1,q2)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)

!     eps = 0.001*ro2 + 1.0d-15
      parameter (delta=1E-4_q)
      eps = min(delta,abs(0.001_q*ro2))+1E-15_q
      x1  = ro2 - eps
      x3  = ro2 + eps
      f1  = fun(ro1,x1,q1,q2)
      f3  = fun(ro1,x3,q1,q2)
      dfro2 = (f3-f1) / (2*eps)
      return
      end

      subroutine expw(rhos,rhops,ex)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)

!*****************************************
!     GGA91 EXCHANGE
!  gga91 exchange for a spin-unpolarized electronic system
!  input d : density
!  input s:  abs(grad d)/(2*kf*d)
!  output:  exchange energy per electron (ex) and its derivatives
!           w.r.t. d (exd) and dd (exdd)
!*****************************************
!
      REAL(q),PARAMETER :: A1=0.19645_q,A2=0.27430_q,A3=0.15084_q,A4=100.0_q
      REAL(q),PARAMETER :: AX=-0.7385588_q,A=7.7956_q,B1=0.004_q
      REAL(q),PARAMETER :: THRD=0.33333333333333_q,THRD4=1.333333333333333_q
      REAL(q),PARAMETER :: THPITH=3.0936677262801_q
!
      rho=2._q*rhos
      rhop=2._q*rhops
      FAC = AX*rho**THRD
!
      w=0.16162046_q
      s=w*abs(rhop)/rho**(4._q/3._q)
      S2 = S*S
      S3 = S2*S
      S4 = S2*S2
!
      S2 = S*S
      S3 = S2*S
      S4 = S2*S2
      P0 = 1.0_q/SQRT(1.0_q+A*A*S2)
      P1 = LOG(A*S+1.0_q/P0)
      P2 = EXP(-A4*S2)
      P3 = 1.0_q/(1.0_q+A1*S*P1+B1*S4)
      P4 = 1.0_q+A1*S*P1+(A2-A3*P2)*S2
      F = P3*P4
      EX = FAC*(F-1)
      RETURN
      END
!
      subroutine cpw(xi,ro,qq,ecc)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)

!*****************************************
!     GGA91 CORRELATION
!*****************************************
      REAL(q),PARAMETER :: C1=0.002568_q,C2=0.023266_q,C3=7.389E-6_q,C4=8.723_q
      REAL(q),PARAMETER :: C5=0.472_q,C6=7.389E-2_q,A4=100.0_q
      REAL(q),PARAMETER :: XNU=15.75592_q,CC0=0.004235_q,CX=-0.001667212_q,ALF=0.09_q
      REAL(q),PARAMETER :: THRDM=-0.333333333333_q,THRD2=0.666666666667_q
      REAL(q),PARAMETER :: THRD=0.33333333333333_q,SIXTH7=1.16666666666667_q
      REAL(q),PARAMETER :: PI=3.1415926536_q
!
!     LOCAL CORRELATION is needed to compute EC:
!
      rs   = (4._q*pi*ro/3._q)**(-1._q/3._q)
      zet   = xi
      CALL CORLSD(RS,ZET,EC,VCUP,VCDN,ECRS,ECZET,ALFC)
!
!     NONLOCAL CORRELATION:
!
      FK = 1.91915829_q/RS
      SK = SQRT(4.0_q*FK/PI)
      G = ((1.0_q+ZET)**THRD2+(1.0_q-ZET)**THRD2)/2.0_q
      t =  abs(qq)/(2.0_q*G*sk*ro)
!
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
!
      ECC=H
!
      RETURN
      END


      subroutine cpwsp(ro1,ro2,q1,q2,despin)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)

      REAL(q),PARAMETER :: cc0=0.004235e0_q
      REAL(q),PARAMETER :: fryd=2.e0_q
      REAL(q),PARAMETER :: thrd=0.333333333333e0_q
!
      ro   = ro1 + ro2
      qq   = q1  + q2
      xi   = (ro1 - ro2)/ro
!
!     local polarization gradient
!
      qxi  = (q1-q2)/ro - xi*qq/ro
!
      deno1 = (ro*(1.e0_q-xi**2))**thrd
      term1 = -0.458e0_q*xi*qxi*(qq/ro)/deno1
      deno2 = (1.e0_q-xi**2)*ro**thrd
      term2 = (-0.037e0_q + 0.10e0_q*xi**2)*(qxi**2)/deno2
!
      despin = cc0*(term1 + term2)
!
      return
      end


      SUBROUTINE CORLSD(RS,ZET,EC,VCUP,VCDN,ECRS,ECZET,ALFC)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
!  UNIFORM-GAS CORRELATION OF PERDEW AND WANG 1991
!  INPUT: SEITZ RADIUS (RS), RELATIVE SPIN POLARIZATION (ZET)
!  OUTPUT: CORRELATION ENERGY PER ELECTRON (EC), UP- AND DOWN-SPIN
!     POTENTIALS (VCUP,VCDN), DERIVATIVES OF EC WRT RS (ECRS) & ZET (ECZET)
!  OUTPUT: CORRELATION CONTRIBUTION (ALFC) TO THE SPIN STIFFNESS
!     IMPLICIT REAL*8 (A-H,O-Z)
      REAL(q),PARAMETER :: GAM=0.5198421E0_q,FZZ=1.709921E0_q
      REAL(q),PARAMETER :: THRD=0.333333333333E0_q,THRD4=1.333333333333E0_q
      F = ((1.E0_q+ZET)**THRD4+(1.E0_q-ZET)**THRD4-2.E0_q)/GAM
      CALL GCOR(0.0310907E0_q,0.21370E0_q,7.5957E0_q,3.5876E0_q,1.6382E0_q, &
     &    0.49294E0_q,1.00E0_q,RS,EU,EURS)
      CALL GCOR(0.01554535E0_q,0.20548E0_q,14.1189E0_q,6.1977E0_q,3.3662E0_q, &
     &    0.62517E0_q,1.00E0_q,RS,EP,EPRS)
      CALL GCOR(0.0168869E0_q,0.11125E0_q,10.357E0_q,3.6231E0_q,0.88026E0_q, &
     &    0.49671E0_q,1.00E0_q,RS,ALFM,ALFRSM)
!  ALFM IS MINUS THE SPIN STIFFNESS ALFC
      ALFC = -ALFM
      Z4 = ZET**4
      EC = EU*(1.E0_q-F*Z4)+EP*F*Z4-ALFM*F*(1.E0_q-Z4)/FZZ
!  ENERGY DONE. NOW THE POTENTIAL:
      ECRS = EURS*(1.E0_q-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1.E0_q-Z4)/FZZ
      FZ = THRD4*((1.E0_q+ZET)**THRD-(1.E0_q-ZET)**THRD)/GAM
      ECZET = 4.E0_q*(ZET**3)*F*(EP-EU+ALFM/FZZ)+FZ*(Z4*EP-Z4*EU &
     &        -(1.E0_q-Z4)*ALFM/FZZ)
      COMM = EC -RS*ECRS/3.E0_q-ZET*ECZET
      VCUP = COMM + ECZET
      VCDN = COMM - ECZET
      RETURN
      END


      SUBROUTINE GCOR(A,A1,B1,B2,B3,B4,P,RS,GG,GGRS)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
!  CALLED BY SUBROUTINE CORLSD
!     IMPLICIT REAL*8 (A-H,O-Z)
      P1 = P + 1.E0_q
      Q0 = -2.E0_q*A*(1.E0_q+A1*RS)
      RS12 = SQRT(RS)
      RS32 = RS12**3
      RSP = RS**P
      Q1 = 2.E0_q*A*(B1*RS12+B2*RS+B3*RS32+B4*RS*RSP)
      Q2 = LOG(1.E0_q+1.E0_q/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/RS12+2.E0_q*B2+3.E0_q*B3*RS12+2.E0_q*B4*P1*RSP)
      GGRS = -2.E0_q*A*A1*Q2-Q0*Q3/(Q1**2+Q1)
      RETURN
      END

