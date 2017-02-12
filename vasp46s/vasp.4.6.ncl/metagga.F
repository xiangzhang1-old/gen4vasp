#include "symbol.inc"
!************************** M E T A G G A . F ******************************
! Implementation of the METAGGA according to 
! Perdew, Kurth, Zupan and Blaha (PRL 82, 2544)
!
! All subroutines in this File were written by Robin Hirschl in Dec. 2000
! using templates supplied by Georg Kresse
! Thanks to Georg for his encouragement and support
!
!***************************************************************************


!************************ SUBROUTINE METAEXC_PW ****************************
!
! this subroutine calculates  the total local potential CVTOT
! which is the sum of the hartree potential, the exchange-correlation
! potential and the ionic local potential
! the routine also calculates the total local potential SV on the small
! grid
! on entry: 
!  CHTOT(:,1) density
!  CHTOT(:,2) respectively CHTOT(:,2:4) contain the magnetization
! on return (LNONCOLLINEAR=.FALSE.):
!  CVTOT(:,1) potential for up
!  CVTOT(:,2) potential for down
! on return (LNONCOLLINEAR=.TRUE.):
!  CVTOT(:,1) 
!
! Robin Hirschl 20001220 (template: POTLOK from pot.F)
!***************************************************************************

      SUBROUTINE METAEXC_PW(GRIDC,WDES,INFO,E,LATT_CUR,CHTOT,TAU,TAUW,DENCOR)
      USE prec
      USE mpimy
      USE mgrid
      USE pseudo
      USE lattice
      USE poscar
      USE setexm
      USE base
      USE xcgrad
      USE wave
!#define robdeb
      IMPLICIT NONE

      TYPE (grid_3d)     GRIDC
      TYPE (wavedes)     WDES
      TYPE (info_struct) INFO
      TYPE (energy)      E
      TYPE (latt)        LATT_CUR

      COMPLEX(q) CHTOT(GRIDC%MPLWV, WDES%NCDIJ)
      RGRID ::  TAU(DIMREAL(GRIDC%MPLWV),WDES%NCDIJ)   ! kinetic energy density
      RGRID ::  TAUW(DIMREAL(GRIDC%MPLWV),WDES%NCDIJ)  ! Weizsaecker kinetic energy density
      RGRID      DENCOR(GRIDC%RL%NP)
      REAL(q)    TMPSIF(3,3)
! work arrays 
      COMPLEX(q), ALLOCATABLE::  CWORK(:,:),TMPWORK(:,:)
      REAL(q) EXC,TMP1,TMP2
      INTEGER ISP,I
#ifdef robdeb
      REAL(q) TOTKINE,TOTKINESP
#endif  
      
      ALLOCATE(CWORK(GRIDC%MPLWV,WDES%NCDIJ),TMPWORK(GRIDC%MPLWV,WDES%ISPIN))
!-----------------------------------------------------------------------
!
!  calculate the exchange correlation energy
!
!-----------------------------------------------------------------------
      EXC     =0._q
      E%EXCM =0._q

      xc: IF (INFO%LEXCHG >= 0) THEN
     ! transform the charge density to real space
         TMPSIF=0

         DO ISP=1,WDES%NCDIJ
            CALL FFT3RC(CHTOT(1,ISP),GRIDC,1)
         ENDDO
         IF (WDES%ISPIN==2 .OR. WDES%LNONCOLLINEAR) THEN

            ! get the charge and the total magnetization
            CALL MAG_DENSITY(CHTOT, CWORK, GRIDC, WDES%NCDIJ)

            IF (INFO%LEXCHG >0) THEN
               ! unfortunately METAGGA_PW requires (up,down) density
               ! instead of (rho,mag)
               CALL RL_FLIP(CWORK, GRIDC, 2, .TRUE.)
               ! GGA potential
               CALL METAGGA_PW(2, GRIDC, LATT_CUR, CWORK, TMPWORK, TAU, TAUW, &
                    DENCOR, TMP1, EXC, TMP2, TMPSIF)
               CALL RL_FLIP(CWORK, GRIDC, 2, .FALSE.)
            ENDIF

         ELSE
            IF (INFO%LEXCHG >0) THEN
               CALL METAGGA_PW(1, GRIDC, LATT_CUR, CHTOT, TMPWORK, TAU, TAUW, &
                    DENCOR, TMP1, EXC, TMP2, TMPSIF)
            ENDIF
                
         ENDIF

         E%EXCM=EXC
         
#ifdef robdeb         
         TOTKINE=0;
         DO ISP=1,WDES%NCDIJ
            TOTKINESP=0.
            DO I=1,GRIDC%RL%NP
               TOTKINESP=TOTKINESP+TAU(I,ISP)
            ENDDO
            TOTKINESP=TOTKINESP*LATT_CUR%OMEGA/GRIDC%RL%NP
            TOTKINE=TOTKINE+TOTKINESP
            WRITE(*,*) 'ISP: ',ISP,' TotKinE PW grid: ',TOTKINESP
         ENDDO
         WRITE(*,*) 'TotKinEW PW grid: ',TOTKINE
         WRITE(*,*) '--------------------------------------------------------'
         WRITE(*,*) 'Exc PW grid:     ',E%EXCM
         WRITE(*,*) '--------------------------------------------------------'
#endif

      ELSE xc
         DO ISP=1,WDES%NCDIJ
            CALL FFT3RC(CHTOT(1,ISP),GRIDC,1)
         ENDDO
      ENDIF xc

      ! CHTOT back to reciprocal space
      DO ISP=1,WDES%NCDIJ
         CALL FFT_RC_SCALE(CHTOT(1,ISP),CHTOT(1,ISP),GRIDC)
         CALL SETUNB_COMPAT(CHTOT(1,ISP),GRIDC)
      ENDDO 
      
      DEALLOCATE(CWORK,TMPWORK)
#undef robdeb
      RETURN
    END SUBROUTINE METAEXC_PW


!************************ SUBROUTINE METAGGA_PW *****************************
! RCS:  $Id: metagga.F,v 1.9 2003/06/27 13:22:19 kresse Exp kresse $
!
!  This routine calculates the meta GGA according to 
!  Perdew, Kurth, Zupan and Blaha (PRL 82, 2544)
!
! the charge density must be passed in real 
!  
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

!************************************************************************
!
! calculate the meta exchange correlation on a plane wave grid
!
!************************************************************************



 SUBROUTINE METAGGA_PW(ISPIN, GRIDC, LATT_CUR, CHTOT, CPOT, TAU, TAUW, DENCOR, &
               XCENC, EXC, CVZERO, XCSIF)
   USE prec
   USE lattice
   USE mpimy
   USE mgrid

   IMPLICIT COMPLEX(q) (C)

   IMPLICIT REAL(q) (A-B,D-H,O-Z)


   INTEGER    ISPIN                     ! ISPIN 1 or 2
   TYPE (grid_3d)     GRIDC             ! descriptor for grid
   TYPE (latt)        LATT_CUR          ! lattice descriptor
   COMPLEX(q) CHTOT(GRIDC%MPLWV,ISPIN)  ! charge density in real space
   COMPLEX(q) TAU(GRIDC%MPLWV,ISPIN)    ! kinetic energy density
   COMPLEX(q) TAUW(GRIDC%MPLWV,ISPIN)   ! Weizsaecker kinetic energy density
   COMPLEX(q) CPOT(GRIDC%MPLWV,ISPIN)   ! exhcange correlation potential
   COMPLEX(q) CMU (GRIDC%MPLWV,ISPIN)   ! exhcange correlation kinetic energy density
   RGRID      DENCOR(GRIDC%RL%NP)       ! pseudo core charge density in real sp

   REAL(q) :: XCENC                     ! double counting correction (unsupported)
   REAL(q) :: EX,EC,EXC                 ! exchange correlation energy
   REAL(q) :: CVZERO                    ! average xc potential
   REAL(q) :: XCSIF(3,3)                ! stress tensor (unsupported)
      
! work arrays
   COMPLEX(q),ALLOCATABLE:: CWGRAD(:,:)
   REAL(q),ALLOCATABLE   :: DWORKG(:,:),DWORK1(:,:),DWORK2(:,:),DWORK3(:,:),DVC(:)

   NP1=GRIDC%RL%NP
   ALLOCATE(CWGRAD(GRIDC%MPLWV,ISPIN), DWORKG(NP1,ISPIN), &
        DWORK1(NP1,ISPIN),DWORK2(NP1,ISPIN),DWORK3(NP1,ISPIN),DVC(NP1))
   
   CALL METAGGA_PW_(ISPIN,GRIDC,LATT_CUR, XCENC,EXC,CVZERO,XCSIF, &
        CWGRAD,CHTOT,CPOT,TAU,TAUW, &
        CWGRAD,CHTOT,CPOT,TAU,TAUW, &
        DENCOR,DWORKG,DWORK1,DWORK2,DWORK3,DVC)
   DEALLOCATE(CWGRAD,DWORKG,DWORK1,DWORK2,DWORK3,DVC)
 RETURN
 END SUBROUTINE METAGGA_PW


!  ATTENTION ATTENTION ATTENTION ATTENTION ATTENTION ATTENTION
! Mind CWORK and DWORK point actually to the same storagelocation
! similar to an EQUIVALENCE (CWORK(1),DWORK(1))
! the same for  (CWGRAD,DWGRAD) and   (CHTOT,DHTOT)
! so we can interchange both arrays arbitrarily

 SUBROUTINE METAGGA_PW_(ISPIN,GRIDC,LATT_CUR,XCENC,EXC,CVZERO,XCSIF, &
      CWGRAD,CHTOT,CWORK,CKE,CWKE, &
      DWGRAD,DHTOT,DWORK,DKE,DWKE, &
      DENCOR,DWORKG,DWORK1,DWORK2,DWORK3,DVC)

   USE prec
   USE lattice
   USE mpimy
   USE mgrid
   USE constant

   IMPLICIT COMPLEX(q) (C)
   IMPLICIT REAL(q) (A-B,D-H,O-Z)

   INTEGER    ISPIN                     ! ISPIN 1 or 2
   TYPE (grid_3d)     GRIDC             ! descriptor for grid
   TYPE (latt)        LATT_CUR          ! lattice descriptor

   
   RGRID      DENCOR(GRIDC%RL%NP)       ! pseudo core charge density in real sp
   REAL(q) :: XCENC                     ! double counting correction (unsupported)
   REAL(q) :: EXC                       ! exchange correlation energy
   REAL(q) :: CVZERO                    ! average xc potential
   REAL(q) :: XCSIF(3,3)                ! stress tensor (unsupported)

   COMPLEX(q) CHTOT (GRIDC%MPLWV,ISPIN),CWORK(GRIDC%MPLWV,ISPIN), &
              CWGRAD(GRIDC%MPLWV,ISPIN),CKE(GRIDC%MPLWV,ISPIN),CWKE(GRIDC%MPLWV,ISPIN)
   RGRID      DHTOT (DIMREAL(GRIDC%MPLWV),ISPIN),DWORK(DIMREAL(GRIDC%MPLWV),ISPIN), &
              DWGRAD(DIMREAL(GRIDC%MPLWV),ISPIN),DKE(DIMREAL(GRIDC%MPLWV),ISPIN), &
              DWKE(DIMREAL(GRIDC%MPLWV),ISPIN)
   REAL(q)    DWORKG(GRIDC%RL%NP,ISPIN),DWORK1(GRIDC%RL%NP,ISPIN), &
              DWORK2(GRIDC%RL%NP,ISPIN),DWORK3(GRIDC%RL%NP,ISPIN), &
              DVC(GRIDC%RL%NP)
   ! set to one for error-dumps
 #ifdef MPI
   NODE_ME=GRIDC%COMM%NODE_ME
   IONODE =GRIDC%COMM%IONODE
   IDUMP=0
   #ifdef debug
       IF (NODE_ME==IONODE) IDUMP=1
   #endif
 #else
   IDUMP=0
   #ifdef debug
     IDUMP=1
   #endif
 #endif

! important constants
   RINPL=1._q/GRIDC%NPLWV                       ! Scaling of Energy
   EVTOH=1._q/(2.*HSQDTM)*AUTOA5                ! KinEDens eV to Hartree
!=======================================================================
! First phase: Transform DENCOR (core charge) and
!  CHTOT (pseudo chargedensity) to reciprocal space space
!=======================================================================
   spin: DO ISP=1,ISPIN
      ! set CWORK to total real charge in reciprocal space
      DO I=1,GRIDC%RL%NP
         DWORK(I,ISP)=(DENCOR(I)/ISPIN+DHTOT(I,ISP))*RINPL/LATT_CUR%OMEGA
      ENDDO
      CALL OPSYNC(DWORK(1,ISP),CWORK(1,ISP),GRIDC%NPLWV)
      CALL FFT3RC(CWORK(1,ISP),GRIDC,-1)
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
   EX=0; EC=0
   DO I=1,GRIDC%RL%NP
      IF (ISPIN==2) THEN
! spin polarized calculation
         RHO1= MAX(REAL((DHTOT(I,1)+DENCOR(I)/ISPIN)/LATT_CUR%OMEGA ,KIND=q),1.E-10_q)
         RHO2= MAX(REAL((DHTOT(I,2)+DENCOR(I)/ISPIN)/LATT_CUR%OMEGA ,KIND=q),1.E-10_q)
         
         ABSNABUP=SQRT(DWORK1(I,1)*DWORK1(I,1)+DWORK2(I,1)*DWORK2(I,1) &
              +DWORK3(I,1)*DWORK3(I,1))
         
         ABSNABDW=SQRT(DWORK1(I,2)*DWORK1(I,2)+DWORK2(I,2)*DWORK2(I,2) &
              +DWORK3(I,2)*DWORK3(I,2))
         
         ABSNAB= SQRT((DWORK1(I,1)+DWORK1(I,2))*(DWORK1(I,1)+DWORK1(I,2))+ &
              (DWORK2(I,1)+DWORK2(I,2))*(DWORK2(I,1)+DWORK2(I,2))+ &
              (DWORK3(I,1)+DWORK3(I,2))*(DWORK3(I,1)+DWORK3(I,2)))
#ifdef correlation_ABS_DRHOUP_ABS_DRHOD
         ABSNAB=ABSNABUP+ABSNABDW
#endif
! kinetic energy density
         TAUU=MAX(REAL(DKE(I,1),KIND=q),1.E-10_q)
         TAUD=MAX(REAL(DKE(I,2),KIND=q),1.E-10_q)
! Weizsaecker kinetic energy density 
         TAUWU=MAX(REAL(DWKE(I,1),KIND=q),1.E-10_q)
         TAUWD=MAX(REAL(DWKE(I,2),KIND=q),1.E-10_q)
! correct kinetic energy densities
! charge density is not the same as the one for which kinetic energy was
! calculated (e.g. augmentation charge)
! use difference in Weizsaecker kinetic energy density for correction
!!$         IF (RHO1>0) THEN
!!$            TAUWTOT=MAX(0.25*HSQDTM*ABSNABUP**2/RHO1,1E-10_q)
!!$            TAUDIFF=TAUWTOT-TAUWU
!!$            TAUWU=TAUWU+TAUDIFF
!!$            TAUU=TAUU+TAUDIFF
!!$         ENDIF
!!$         IF (RHO2>0) THEN
!!$            TAUWTOT=MAX(0.25*HSQDTM*ABSNABDW**2/RHO2,1E-10_q)
!!$            TAUDIFF=TAUWTOT-TAUWD
!!$            TAUWD=TAUWD+TAUDIFF
!!$            TAUD=TAUD+TAUDIFF
!!$         ENDIF
      ELSE
! non spin polarized calculation
! up and down-spin variables get half of the total values
         RHO1= MAX(REAL((DHTOT(I,1)+DENCOR(I))/2._q/LATT_CUR%OMEGA ,KIND=q),1.E-10_q)
         RHO2= RHO1
         ABSNAB=SQRT(DWORK1(I,1)*DWORK1(I,1)+DWORK2(I,1)*DWORK2(I,1) &
              +DWORK3(I,1)*DWORK3(I,1))
         ABSNABUP= 0.5_q*ABSNAB
         ABSNABDW= 0.5_q*ABSNAB
! kinetic energy density
         TAUU= MAX(0.5_q*REAL(DKE(I,1),KIND=q),1.E-10_q)
! Weizsaecker kinetic energy density 
         TAUWU= MAX(0.5_q*REAL(DWKE(I,1),KIND=q),1.E-10_q)
! correct kinetic energy densities
! charge density is not the same as the one for which kinetic energy was
! calculated (e.g. augmentation charge)
! use difference in Weizsaecker kinetic energy density for correction
!!$         IF (RHO1>0) THEN
!!$            TAUWTOT=MAX(0.25*HSQDTM*ABSNABUP**2/RHO1,1E-10_q)
!!$            TAUDIFF=TAUWTOT-TAUWU
!!$            TAUWU=TAUWU+TAUDIFF
!!$            TAUU=TAUU+TAUDIFF
!!$         ENDIF
         TAUD= TAUU
         TAUWD=TAUWU        
      ENDIF


! All parameters for subroutine Metagga must be passed in Hartree
!      WRITE(77,'(I6,9E14.4)') I,RHO1*AUTOA3, &
!           TAUU*EVTOH,TAUD*EVTOH,TAUWU,TAUWD,TAUWU/(TAUU*EVTOH), &
!           TAUWD/(TAUD*EVTOH),ECL*LATT_CUR%OMEGA/AUTOA3,DENCOR(I)/(2*LATT_CUR%OMEGA)    

      CALL METAGGA(RHO1*AUTOA3,RHO2*AUTOA3, &
     &               ABSNABUP*AUTOA4, ABSNABDW*AUTOA4,ABSNAB*AUTOA4, &
     &           TAUU*EVTOH,TAUD*EVTOH, &
     &           TAUWU*EVTOH,TAUWD*EVTOH,EXL,ECL,I)

      DEXC1=0; DEXC2=0; DVXC1=0; DVXC2=0; DVC_=0
!     RHO=RHO1+RHO2

! Conversion back to eV
      EX=EX+2*EXL*RYTOEV*LATT_CUR%OMEGA/AUTOA3
      EC=EC+2*ECL*RYTOEV*LATT_CUR%OMEGA/AUTOA3

! ATTENTION ATTENTION ATTENTION
! DWORKG is   N O T   properly defined in this function
! should be |nabla rho|

      DVXC1=DVXC1*RYTOEV*AUTOA
      DVXC2=DVXC2*RYTOEV*AUTOA
      DVC_ =DVC_ *RYTOEV*AUTOA
#ifdef correlation_ABS_DRHOUP_ABS_DRHOD
      DVXC1=DVXC1+DVC_ 
      DVXC2=DVXC2+DVC_ 
      DVC_=0
#endif
!
!   store d f/ d (|d rho| ) / |d rho|  in DWORK
!
!      DWORK(I,1)  = DVXC1 / MAX(DWORKG(I,1),1.E-10_q)
!      DWORK(I,2)  = DVXC2 / MAX(DWORKG(I,2),1.E-10_q)
!      DVC(I)      = DVC_  / MAX(ABSNAB,1.E-10_q)
!
!   store d f/ d rho  in DWORKG
!
!      DWORKG(I,1) = DEXC1*RYTOEV
!      DWORKG(I,2) = DEXC2*RYTOEV
   ENDDO
   
   ! OUTPUT OF MEGGA EXCHANGE AND CORRELATION
!   WRITE(*,*)
!   WRITE(*,'(2(A,F14.6))') 'Exchange energy    eV:',EX*RINPL,'  Hartree:',EX*RINPL/(2*RYTOEV)
!   WRITE(*,'(2(A,F14.6))') 'Correlation energy eV:',EC*RINPL,'  Hartree:',EC*RINPL/(2*RYTOEV)
!   WRITE(*,*)

! collect results from all nodes. Can be deleted as soon as the
! return command is shifted further down.   
   CALLMPI( M_sum_s(GRIDC%COMM, 2, EX, EC, 0._q , 0._q ))

   EXC=(EX+EC)*RINPL
   RETURN
!
! FOR LATER USE WE HAVE INCLUDED THE RELEVANT ROUTINES FOR THE
! EVALUATION OF THE POTENTIALS AND OF THE STRESS-TENSOR
!

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
      RHO1= MAX(REAL((DHTOT(I,1)+DENCOR(I)/2)/LATT_CUR%OMEGA ,KIND=q), 1E-10_q)
      RHO2= MAX(REAL((DHTOT(I,2)+DENCOR(I)/2)/LATT_CUR%OMEGA ,KIND=q), 1E-10_q)

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

   CALLMPI( M_sum_s(GRIDC%COMM, 2, EXC, XCENC, 0._q , 0._q ))
   CALLMPI( M_sum_d(GRIDC%COMM, XCSIF, 9))
   CALLMPI( M_sum_z(GRIDC%COMM,CVZERO,1))

! Test dumps:
   IF (IDUMP/=0) THEN
      WRITE(*,'(A,F24.14)') '<rho*excgc> =',EXC
      WRITE(*,'(A,F24.14)') '<rho*vxcgc> =',EXC-XCENC
      WRITE(*,'(A,F24.14)') '    xcencgc =',XCENC
   ENDIF
   RETURN
 END SUBROUTINE METAGGA_PW_



!************************ SUBROUTINE GGASPINCOR ************************
!
!  calculate the correlation energy density according to the
!  Perdew, Burke and Ernzerhof functional
!
!***********************************************************************

    SUBROUTINE GGASPINCOR(D1,D2,DDA, EC)

!     D1   density up
!     D2   density down
!     DDA  |gradient of the total density|

      USE prec
      USE constant
      IMPLICIT REAL(q) (A-H,O-Z)
      PARAMETER (THRD=1._q/3._q)

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

      EC  =(EC  +ECLDA)

      RETURN
    END SUBROUTINE GGASPINCOR

!************************ SUBROUTINE GGACOR *****************************
!
!  calculate the correlation energy density according to the
!  Perdew, Burke and Ernzerhof functional
!
!***********************************************************************

    SUBROUTINE GGACOR(D, DD, EC)
      USE prec
      USE constant
      IMPLICIT REAL(q) (A-H,O-Z)
      PARAMETER (THRD=1._q/3._q)

      IF (D<0) THEN
         EC   = 0._q
         RETURN
      ENDIF

      DTHRD=exp(log(D)*THRD)
      RS=(0.75_q/PI)**THRD/DTHRD
      FK=(3._q*PI*PI)**THRD*DTHRD
      SK = SQRT(4.0_q*FK/PI)

      IF(D>1.E-10_q)THEN
         T=DD/(D*SK*2._q)
      ELSE
         T=0.0_q
      ENDIF

      CALL CORunspPBE(RS,ECLDA,ECDLDA,SK, &
           T,EC,ECD,ECDD,.TRUE.)

      EC = (EC+ECLDA)

      RETURN
    END SUBROUTINE GGACOR


!=======================================================================
! 
! SUBROUTINE TAU_PW
!
! This subroutine calculates the kinetic energy of the PW part of the
! wavefunctions (0.5*|grad psi|**2), output on GRIDC
! INPUT: GRID, GRID_SOFT, GRIDC, SOFT_TO_C, LATT_CUR, SYMM, NIOND, W, WDES
! OUTPUT: (GRIDC%MPLWC, WDES%NCDIJ)
!
! |grad psi|**2=0.5*grad**2 rho+psi* grad**2 psi
!
!
! Robin Hirschl 20001119
!=======================================================================

SUBROUTINE TAU_PW(GRID,GRID_SOFT,GRIDC,SOFT_TO_C,LATT_CUR,SYMM,NIOND, &
        W,WDES,TAU)
      USE prec
      USE lattice
      USE mgrid
      USE msymmetry
      USE base
      USE wave
      USE mpimy
      USE constant
      
      IMPLICIT NONE

      TYPE (grid_3d)     GRID,GRID_SOFT,GRIDC
      TYPE (transit)     SOFT_TO_C
      TYPE (latt)        LATT_CUR
      TYPE (symmetry)    SYMM
      TYPE (wavespin)    W
      TYPE (wavedes)     WDES

      REAL(q)   WEIGHT,GX,GY,GZ,GSQU
      INTEGER   NIOND,MPLWV,ISP,NK,N,NPL,I,N1,NC,N2,N3
! result
      RGRID :: TAU(DIMREAL(GRIDC%MPLWV),WDES%NCDIJ)    ! kinetic energy density
! dynamic work array
      COMPLEX(q),ALLOCATABLE :: CW(:), CW2(:)
      RGRID, ALLOCATABLE :: CW3(:), CW4(:) 
      COMPLEX(q) :: CF(WDES%NRPLWV)

      MPLWV=MAX(GRID%MPLWV, GRID_SOFT%MPLWV)
      ALLOCATE(CW(MPLWV), CW2(MPLWV), CW3(DIMREAL(MPLWV)), CW4(DIMREAL(MPLWV)))
      
      TAU=0

      IF (WDES%NCDIJ==4) THEN
         WRITE(*,*) 'WARNING: kinetic energy density not implemented for non collinear case.'
         WRITE(*,*) 'exiting TAU_PW; sorry for the inconveniences.'
         RETURN
      ENDIF


spin: DO ISP=1,WDES%NCDIJ

!=======================================================================
! calculate psi* grad**2 psi -> result in CW3 on GRID_SOFT
! and recalculate chargedensity (psi*psi*) in CW (complex!) on GRID_SOFT
!=======================================================================
      CW3=0
      CW4=0
      
      band: DO N=1,WDES%NBANDS
         kpoints: DO NK=1,WDES%NKPTS
            WEIGHT=WDES%RSPIN*WDES%WTKPT(NK)*W%FERWE(N,NK,ISP)/LATT_CUR%OMEGA
            NPL=WDES%NPLWKP(NK)
            DO I=1,NPL
!-MM- changes to accommodate spin spirals
! original statement
!              CF(I)=WDES%DATAKE(I,NK)*W%CPTWFP(I,N,NK,ISP)
               CF(I)=WDES%DATAKE(I,NK,ISP)*W%CPTWFP(I,N,NK,ISP)
!-MM- end of alterations
            ENDDO
           
! fourier trafo of wave-function
            CALL FFTWAV(NPL,WDES%NINDPW(1,NK),CW(1),W%CPTWFP(1,N,NK,ISP),GRID)
! fourier trafo of k**2 psi
            CALL FFTWAV(NPL,WDES%NINDPW(1,NK),CW2(1),CF(1),GRID)

            DO I=1,GRID%RL%NP
               CW3(I)=CW3(I)+REAL(CW(I)*CONJG(CW2(I)), KIND=q)*WEIGHT
               CW4(I)=CW4(I)+REAL(CW(I)*CONJG(CW(I)), KIND=q)*WEIGHT
            ENDDO
            
         ENDDO kpoints
      ENDDO band
! merge results from nodes
#ifdef realmode
      CALLMPI( M_sum_d(WDES%COMM_INTER, CW3(1), GRID%RL%NP))
      CALLMPI( M_sum_d(WDES%COMM_INTER, CW4(1), GRID%RL%NP))
#else
      CALLMPI( M_sum_z(WDES%COMM_INTER, CW3(1), GRID%RL%NP))
      CALLMPI( M_sum_z(WDES%COMM_INTER, CW4(1), GRID%RL%NP))
#endif

! reset arrays CW and CW2 (will be reused)
      CW=0
      CW2=0

! charge density in complex array in rec. space
      CALL FFT_RC_SCALE(CW4,CW,GRID_SOFT)

!=======================================================================
! calculate grad**2 rho -> result in CW2 on GRID_SOFT
!=======================================================================
      
      DO I=1,GRID_SOFT%RC%NP
         N1= MOD((I-1),GRID_SOFT%RC%NROW) +1
         NC= (I-1)/GRID_SOFT%RC%NROW+1
         N2= GRID_SOFT%RC%I2(NC)
         N3= GRID_SOFT%RC%I3(NC)
         GX=(GRID_SOFT%LPCTX(N1)*LATT_CUR%B(1,1)+GRID_SOFT%LPCTY(N2)*LATT_CUR%B(1,2)+GRID_SOFT%LPCTZ(N3)*LATT_CUR%B(1,3))
         GY=(GRID_SOFT%LPCTX(N1)*LATT_CUR%B(2,1)+GRID_SOFT%LPCTY(N2)*LATT_CUR%B(2,2)+GRID_SOFT%LPCTZ(N3)*LATT_CUR%B(2,3))
         GZ=(GRID_SOFT%LPCTX(N1)*LATT_CUR%B(3,1)+GRID_SOFT%LPCTY(N2)*LATT_CUR%B(3,2)+GRID_SOFT%LPCTZ(N3)*LATT_CUR%B(3,3))
         GSQU=GX*GX+GY*GY+GZ*GZ
         CW2(I)=-CW(I)*GSQU*TPI*TPI
      ENDDO
      
! result in real space
      CALL FFT3RC(CW2,GRID_SOFT,1)


!==============================================================================
! sum up kinetic energy density, symmetrize and transform to fine grid (GRIDC)
!==============================================================================
      CALL RL_ADD(CW3,1.0_q,CW2,0.5_q*HSQDTM,CW3,GRID_SOFT)

! rescaling      
      DO I=1,GRID_SOFT%RL%NP
         CW3(I)=CW3(I)/GRID_SOFT%NPLWV
      ENDDO

! transition to finer grid only in reciprocal space     
      CALL FFT3RC(CW3,GRID_SOFT,-1)

! symmetrization of result CW2 
      IF (SYMM%ISYM>0) THEN
         ! symmetrization
         CALL RHOSYM(CW3,GRID_SOFT,SYMM%PTRANS,NIOND,SYMM%MAGROT,ISP)
      ENDIF
      CALL CPB_GRID(GRIDC,GRID_SOFT,SOFT_TO_C,CW3(1),TAU(1,ISP))
      CALL FFT3RC(TAU(1,ISP),GRIDC,1)
    
      ENDDO spin
      DEALLOCATE(CW,CW2,CW3,CW4)
      RETURN
END SUBROUTINE TAU_PW


!=======================================================================
! 
! SUBROUTINE TAU_PW_DIRECT
!
! This subroutine calculates the kinetic energy of the PW part of the
! wavefunctions (0.5*|grad psi|**2) 
! and the Weizsaecker kinetic energy density, output on GRIDC
! INPUT: GRIDC,  LATT_CUR, SYMM, NIOND, W, WDES
! OUTPUT: TAU(GRIDC%MPLWC, WDES%NCDIJ)
!         TAUW(GRIDC%MPLWC, WDES%NCDIJ)
! evaluation should be more accurate than TAU_PW due to avoidance of 
! wrap arounf errors
!
! Robin Hirschl 20001221 
!=======================================================================

SUBROUTINE TAU_PW_DIRECT(GRID,GRID_SOFT,GRIDC,SOFT_TO_C,LATT_CUR,SYMM,NIOND, &
        W,WDES,TAU,TAUW)
      USE prec
      USE lattice
      USE mgrid
      USE msymmetry
      USE base
      USE wave
      USE mpimy
      USE constant
      
      IMPLICIT NONE

      TYPE (grid_3d)     GRIDC,GRID,GRID_SOFT
      TYPE (latt)        LATT_CUR
      TYPE (transit)     SOFT_TO_C
      TYPE (symmetry)    SYMM
      TYPE (wavespin)    W
      TYPE (wavedes)     WDES
! result
      RGRID ::  TAU(DIMREAL(GRIDC%MPLWV),WDES%NCDIJ)  ! kinetic energy density
      RGRID ::  TAUW(DIMREAL(GRIDC%MPLWV),WDES%NCDIJ) ! weiz kin edens
! dynamic work array
      COMPLEX(q),ALLOCATABLE :: CW(:),CW3(:)
      RGRID,ALLOCATABLE :: CW2(:),CW4(:),CW5(:)
      COMPLEX(q) :: CF(WDES%NRPLWV,3)
      
      INTEGER :: MPLWV,ISP,N,NK,NPL,IDIR,I,NIOND
      REAL(q) :: WEIGHT,G1,G2,G3
      COMPLEX(q) :: GC

      MPLWV=MAX(GRID%MPLWV, GRID_SOFT%MPLWV)
      ALLOCATE(CW(MPLWV),CW2(DIMREAL(MPLWV)),CW3(MPLWV),CW4(DIMREAL(MPLWV)), &
           CW5(DIMREAL(MPLWV)))
      
      TAU=0; TAUW=0
      
      IF (WDES%NCDIJ==4) THEN
         WRITE(*,*) 'WARNING: kinetic energy density not implemented for non collinear case.'
         WRITE(*,*) 'exiting TAU_PW_DIRECT; sorry for the inconveniences.'
         RETURN
      ENDIF

spin: DO ISP=1,WDES%NCDIJ

      CW2=0; CW3=0; CW4=0; CW5=0
      band: DO N=1,WDES%NBANDS
         kpoints: DO NK=1,WDES%NKPTS
            WEIGHT=WDES%RSPIN*WDES%WTKPT(NK)*W%FERWE(N,NK,ISP)/LATT_CUR%OMEGA
            NPL=WDES%NPLWKP(NK)
            CF=0
! loop over plane wave coefficients            
            DO I=1,NPL
! get k-vector of respective k-point and coefficient
               G1=WDES%IGX(I,NK)+WDES%VKPT(1,NK)
               G2=WDES%IGY(I,NK)+WDES%VKPT(2,NK)
               G3=WDES%IGZ(I,NK)+WDES%VKPT(3,NK)
! loop over cartesian directions
               DO IDIR=1,3
                  GC=(G1*LATT_CUR%B(IDIR,1)+G2*LATT_CUR%B(IDIR,2)+G3*LATT_CUR%B(IDIR,3))
                  CF(I,IDIR)=GC*W%CPTWFP(I,N,NK,ISP)*CITPI
               ENDDO
            ENDDO

! fourier trafo of wave-function (result in CW3)
            CALL FFTWAV(NPL,WDES%NINDPW(1,NK),CW3(1),W%CPTWFP(1,N,NK,ISP),GRID)


! fourier trafo of gradient of wave-function 
! loop over cartesian directions
            DO IDIR=1,3
               CALL FFTWAV(NPL,WDES%NINDPW(1,NK),CW(1),CF(1,IDIR),GRID)
! update kinetic energy density (in CW2) and grad rho^2 (in CW4)
               DO I=1,GRID%RL%NP
                  CW2(I)=CW2(I)+HSQDTM*REAL(CW(I)*CONJG(CW(I)),KIND=q)*WEIGHT
                  CW4(I)=CW4(I)+HSQDTM*REAL(CW(I)*CW3(I)*CONJG(CW(I)*CW3(I)),KIND=q)* &
                       WEIGHT*WEIGHT/WDES%WTKPT(NK)
               ENDDO
            ENDDO

! update Weizsaecker KinEDens Denominator (charge density) (in CW5)           
           DO I=1,GRID%RL%NP
              CW5(I)=CW5(I)+REAL(CW3(I)*CONJG(CW3(I)),KIND=q)*WEIGHT
           ENDDO
        ENDDO kpoints
      ENDDO band


#ifdef realmode
      CALLMPI( M_sum_d(WDES%COMM_INTER, CW5(1), GRID%RL%NP))
#else
      CALLMPI( M_sum_z(WDES%COMM_INTER, CW5(1), GRID%RL%NP))
#endif

! calculate Weizsaecker KinEdens
      DO I=1,GRID%RL%NP
         IF (CW5(I)==0._q) THEN
            CW4(I)=0._q
         ELSE
            CW4(I)=CW4(I)/CW5(I)
         ENDIF
      ENDDO

! merge results from nodes
#ifdef realmode
      CALLMPI( M_sum_d(WDES%COMM_INTER, CW2(1), GRID%RL%NP))
      CALLMPI( M_sum_d(WDES%COMM_INTER, CW4(1), GRID%RL%NP))
#else
      CALLMPI( M_sum_z(WDES%COMM_INTER, CW2(1), GRID%RL%NP))
      CALLMPI( M_sum_z(WDES%COMM_INTER, CW4(1), GRID%RL%NP))
#endif

! rescaling
      DO I=1,GRID_SOFT%RL%NP
         CW2(I)=CW2(I)/GRID_SOFT%NPLWV
         CW4(I)=CW4(I)/GRID_SOFT%NPLWV         
      ENDDO
    
! to rec space
      CALL FFT3RC(CW2(1),GRID_SOFT,-1)
      CALL FFT3RC(CW4(1),GRID_SOFT,-1)
      
! transition to finer grid
      CALL CPB_GRID(GRIDC,GRID_SOFT,SOFT_TO_C,CW2(1),TAU(1,ISP))
      CALL CPB_GRID(GRIDC,GRID_SOFT,SOFT_TO_C,CW4(1),TAUW(1,ISP))

   ENDDO spin

! symmetrization of result TAU(:,ISP)
! needs (total,mag) instead of up,dw
   IF (SYMM%ISYM>0) THEN
      IF (WDES%NCDIJ==2) THEN
         CALL RC_FLIP(TAU,GRIDC,2,.FALSE.)
         CALL RC_FLIP(TAUW,GRIDC,2,.FALSE.)
      ENDIF
      ! symmetrization
      CALL RHOSYM(TAU(1,1),GRIDC,SYMM%PTRANS,NIOND,SYMM%MAGROT,1)
      CALL RHOSYM(TAUW(1,1),GRIDC,SYMM%PTRANS,NIOND,SYMM%MAGROT,1)
      IF (WDES%NCDIJ==2) THEN
         CALL RHOSYM(TAU(1,2),GRIDC,SYMM%PTRANS,NIOND,SYMM%MAGROT,2)
         CALL RHOSYM(TAUW(1,2),GRIDC,SYMM%PTRANS,NIOND,SYMM%MAGROT,2)
         CALL RC_FLIP(TAU,GRIDC,2,.TRUE.)
         CALL RC_FLIP(TAUW,GRIDC,2,.TRUE.)
      ENDIF
   ENDIF
   
   DO ISP=1,WDES%NCDIJ
! back to real space
      CALL FFT3RC(TAU(1,ISP),GRIDC,1)
      CALL FFT3RC(TAUW(1,ISP),GRIDC,1)

! ATTENTION:
! the transition to finer grid may cause instabilities with tauw being larger than tau
! this has to be corrected !
      TAUW(:,ISP)=MIN(REAL(TAUW(:,ISP),q),REAL(TAU(:,ISP),q))
   ENDDO
   DEALLOCATE(CW,CW2,CW3,CW4,CW5)
   
   RETURN
END SUBROUTINE TAU_PW_DIRECT



!************************ SUBROUTINE METAGGA *****************************
!
! calculates local contribution to metagga Exc according to 
! Perdew et. al. PRL 82, 12 (1999)
!
! RH 20001119
!
! everything in Hartree units
!
! ATTANTION: Every values are passed "as they are", i.e. including
! possibly unphysical numerical errors (e.g. negative charge densities)
! values need to be checked accordingly
!***********************************************************************

SUBROUTINE METAGGA(RU,RD,DRU,DRD,DRT,TAUU,TAUD,TAUWU,TAUWD,EX,EC,I)

! RU,RD      density up,down
! DRU, DRD   abs. val. gradient of density up/down
! DRT        abs. val. gradient of total density
! TAUU,TAUD  kinetic energy density up/down
! TAUWU,TAUWD Weizsaecker kinetic energy density up/down
! EXC        return value
  
  USE prec
  USE constant
  IMPLICIT REAL(q) (A-H,O-Z)
  
  INTEGER I
! the following parameters are given by Perdew et.al.
  PARAMETER (RKAPPA=0.804_q)
  PARAMETER (D=0.113_q)
  PARAMETER (C=0.53_q)
! other parameters
  PARAMETER (THRD=1._q/3._q)
  PARAMETER (TTHRD=2._q*THRD)
  PARAMETER (FTHRD=1._q+TTHRD)
  PARAMETER (ETHRD=1._q+FTHRD)
  PARAMETER (PISQ=PI*PI)


  EX=0._q;EC=0._q
! exchange energy
! spin up
     P=(2._q*DRU)**2._q/(4._q*(3._q*PISQ)**TTHRD*(2._q*RU)**ETHRD)
     QQS=6._q*TAUU/(2._q*(3._q*PISQ)**TTHRD*(2._q*RU)**FTHRD)-9._q/20._q-P/12._q
     X=10._q/81._q*P+146._q/2025._q*QQS*QQS-73._q/405._q*QQS*P
     X=X+(D+1._q/RKAPPA*(10._q/81._q)**2._q)*P*P
     FX=1._q+RKAPPA-RKAPPA/(1._q+(X/RKAPPA))
     EX=EX-RU*(3._q/(4._q*PI))*(3._q*PISQ*2._q*RU)**THRD*FX
! spin down
     P=(2._q*DRD)**2._q/(4._q*(3._q*PISQ)**TTHRD*(2._q*RD)**ETHRD)
     QQS=6._q*TAUD/(2._q*(3._q*PISQ)**TTHRD*(2._q*RD)**FTHRD)-9._q/20._q-P/12._q
     X=10._q/81._q*P+146._q/2025._q*QQS*QQS-73._q/405._q*QQS*P
     X=X+(D+1._q/RKAPPA*(10._q/81._q)**2._q)*P*P
     FX=1._q+RKAPPA-RKAPPA/(1._q+(X/RKAPPA))
     EX=EX-RD*(3._q/(4._q*PI))*(3._q*PISQ*2._q*RD)**THRD*FX

! correlation energy
     CALL GGASPINCOR(RU,RD,DRT,ECT)
     TAUK=(TAUWU+TAUWD)/(TAUU+TAUD)
     ECM1=(RU+RD)*ECT*(1._q+C*TAUK**2._q)
     
!     CALL GGACOR(RU,DRU,ECU)
     CALL GGASPINCOR(RU,0.0_q,DRU,ECU)
     TAUK=TAUWU/TAUU
     ECM2=TAUK**2._q*RU*ECU
     
!     CALL GGACOR(RD,DRD,ECD)
     CALL GGASPINCOR(RD,0.0_q,DRD,ECD)
     TAUK=TAUWD/TAUD
     ECM3=TAUK**2._q*RD*ECD
     EC=ECM1-(1._q+C)*(ECM2+ECM3)
     
  RETURN
END SUBROUTINE METAGGA
