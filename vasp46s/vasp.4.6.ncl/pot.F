#include "symbol.inc"
      MODULE pot
      USE prec
      USE charge
      CONTAINS
!************************ SUBROUTINE POTLOK ****************************
! RCS:  $Id: pot.F,v 1.5 2003/06/27 13:22:22 kresse Exp kresse $
!
! this subroutine calculates  the total local potential CVTOT
! which is the sum of the hartree potential, the exchange-correlation
! potential and the ionic local potential
! the routine also calculates the total local potential SV on the small
! grid
! on entry: 
!  CHTOT(:,1)    density
!  CHTOT(:,2)    respectively CHTOT(:,2:4) contain the magnetization
! on return (LNONCOLLINEAR=.FALSE.):
!  CVTOT(:,1)    potential for up
!  CVTOT(:,2)    potential for down
! on return (LNONCOLLINEAR=.TRUE.):
!  CVTOT(:,1:4)  spinor representation of potential
!
!***********************************************************************

      SUBROUTINE POTLOK(GRID,GRIDC,GRID_SOFT, COMM_INTER, WDES,  &
                  EXCTAB,INFO,P,T_INFO,E,LATT_CUR, DIP, &
                  CHTOT,CSTRF,CVTOT,DENCOR,SV, SOFT_TO_C,XCSIF )
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
!-MM- changes to accomodate constrained moments
      USE Constrained_M_modular
!-MM- end of addition

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRID,GRIDC,GRID_SOFT
      TYPE (wavedes)     WDES
      TYPE (transit)     SOFT_TO_C
      TYPE (exctable)    EXCTAB
      TYPE (info_struct) INFO
      TYPE (type_info)   T_INFO
      TYPE (potcar)      P (T_INFO%NTYP)
      TYPE (energy)      E
      TYPE (latt)        LATT_CUR
      TYPE (dipol)       DIP
      TYPE (communic)    COMM_INTER

      RGRID   SV(DIMREAL(GRID%MPLWV), WDES%NCDIJ)
      COMPLEX(q) CSTRF(GRIDC%MPLWV,T_INFO%NTYP), &
                 CHTOT(GRIDC%MPLWV, WDES%NCDIJ), CVTOT(GRIDC%MPLWV,WDES%NCDIJ)
      RGRID      DENCOR(GRIDC%RL%NP)
      REAL(q)    XCSIF(3,3),TMPSIF(3,3)
! work arrays (allocated after call to FEXCG)
      COMPLEX(q), ALLOCATABLE::  CWORK1(:),CWORK(:,:)
      REAL(q) ELECTROSTATIC
      LOGICAL, EXTERNAL :: L_NO_LSDA_GLOBAL
      
      MWORK1=MAX(GRIDC%MPLWV,GRID_SOFT%MPLWV)
      ALLOCATE(CWORK1(MWORK1),CWORK(GRIDC%MPLWV,WDES%NCDIJ))
!-----------------------------------------------------------------------
!
!  calculate the exchange correlation potential and the dc. correction
!
!-----------------------------------------------------------------------
      EXC     =0
      E%XCENC =0
      E%EXCG  =0
      E%CVZERO=0
      XCSIF   =0

      CVTOT=0
  xc: IF (INFO%LEXCHG >= 0) THEN
     ! transform the charge density to real space
        EXCG  =0
        XCENCG=0
        CVZERG=0
        TMPSIF=0

        DO ISP=1,WDES%NCDIJ
           CALL FFT3RC(CHTOT(1,ISP),GRIDC,1)
        ENDDO
        IF (WDES%ISPIN==2) THEN

          ! get the charge and the total magnetization
          CALL MAG_DENSITY(CHTOT, CWORK, GRIDC, WDES%NCDIJ)
! do LDA+U instead of LSDA+U
          IF (L_NO_LSDA_GLOBAL()) CWORK(:,2)=0
!
          IF (INFO%LEXCHG >0) THEN
             ! gradient corrections to LDA
             ! unfortunately FEXCGS requires (up,down) density
             ! instead of (rho,mag)
             CALL RL_FLIP(CWORK, GRIDC, 2, .TRUE.)
             ! GGA potential
             CALL FEXCGS(2, GRIDC, LATT_CUR, XCENCG, EXCG, CVZERG, TMPSIF, &
                  CWORK, CVTOT, DENCOR, INFO%LEXCHG)
             CALL RL_FLIP(CWORK, GRIDC, 2, .FALSE.)
          ENDIF

          ! add LDA part of potential
          CALL FEXCF(EXCTAB,GRIDC,LATT_CUR%OMEGA, &
             CWORK(1,1), CWORK(1,2), DENCOR, CVTOT(1,1), CVTOT(1,2), &
             E%CVZERO,EXC,E%XCENC,XCSIF, .TRUE.)
          ! we have now the potential for up and down stored in CVTOT(:,1) and CVTOT(:,2)

          ! get the proper direction vx = v0 + hat m delta v
          CALL MAG_DIRECTION(CHTOT(1,1), CVTOT(1,1), GRIDC, WDES%NCDIJ)
        ELSEIF (WDES%LNONCOLLINEAR) THEN
!-MM- gradient corrections in the noncollinear case are calculated
!     a bit differently than in the collinear case
          IF (INFO%LEXCHG >0) THEN
             ! GGA potential
             CALL FEXCGS(4, GRIDC, LATT_CUR, XCENCG, EXCG, CVZERG, TMPSIF, &
                  CHTOT, CVTOT, DENCOR, INFO%LEXCHG)
          ENDIF

          ! FEXCF requires (up,down) density instead of (rho,mag)
          CALL MAG_DENSITY(CHTOT, CWORK, GRIDC, WDES%NCDIJ)
! quick hack to do LDA+U instead of LSDA+U
          IF (L_NO_LSDA_GLOBAL()) CWORK(:,2)=0
! end of hack
          ! add LDA part of potential
          CALL FEXCF(EXCTAB,GRIDC,LATT_CUR%OMEGA, &
             CWORK(1,1), CWORK(1,2), DENCOR, CVTOT(1,1), CVTOT(1,2), &
             E%CVZERO,EXC,E%XCENC,XCSIF, .TRUE.)
          ! we have now the potential for up and down stored in CVTOT(:,1) and CVTOT(:,2)
          ! get the proper direction vx = v0 + hat m delta v
                    
          CALL MAG_DIRECTION(CHTOT(1,1), CVTOT(1,1), GRIDC, WDES%NCDIJ)
!-MM- end of changes to calculation of gga in noncollinear case
       ELSE
          IF (INFO%LEXCHG >0) THEN
             ! gradient corrections to LDA
             CALL FEXCG(INFO%LEXCHG,GRIDC,LATT_CUR,XCENCG,EXCG,CVZERG,TMPSIF, &
                  CHTOT,CVTOT,DENCOR)
          ENDIF
                
          ! LDA part of potential
          CALL FEXCP(EXCTAB,GRIDC,LATT_CUR%OMEGA, &
               CHTOT,DENCOR,CVTOT,CWORK,E%CVZERO,EXC,E%XCENC,XCSIF,.TRUE.)
       ENDIF

       XCSIF=XCSIF+TMPSIF
       E%EXCG=EXC+EXCG
       E%XCENC=E%XCENC+XCENCG
       E%CVZERO=E%CVZERO+CVZERG

      ELSE xc
         DO ISP=1,WDES%NCDIJ
            CALL FFT3RC(CHTOT(1,ISP),GRIDC,1)
         ENDDO
      ENDIF xc
!-MM- changes to accomodate constrained moments
!-----------------------------------------------------------------------
! add constraining potential
!-----------------------------------------------------------------------
#if ! defined(NGXhalf) && ! defined(NGZhalf)
      IF (M_CONSTRAINED()) THEN
      ! NB. at this point both CHTOT and CVTOT must be given
      ! in (charge,magnetization) convention in real space
         CALL M_INT(CHTOT,GRIDC,WDES)
         CALL ADD_CONSTRAINING_POT(CVTOT,GRIDC,WDES)
      ENDIF
#endif
!-MM- end of addition

!-----------------------------------------------------------------------
! calculate the total potential
!-----------------------------------------------------------------------
! add external electrostatic potential
      DIP%EDIPOL=0
      DIP%E_ION_EXTERN=0

      IF (DIP%LCOR_DIP) THEN
          ! get the total charge
          IF  ( WDES%NCDIJ > 1) THEN
             CALL MAG_DENSITY(CHTOT,CWORK, GRIDC, WDES%NCDIJ)
          ELSE
             CALL RL_ADD(CHTOT,1.0_q,CHTOT,0.0_q,CWORK,GRIDC)
          ENDIF

           CALL CDIPOL(GRIDC, LATT_CUR,P,T_INFO, DIP, &
             CWORK,CSTRF,CVTOT(1,1), WDES%NCDIJ, INFO%NELECT )

         CALL EXTERNAL_POT(GRIDC, LATT_CUR, CVTOT(1,1))
      ENDIF

      DO ISP=1,WDES%NCDIJ
         CALL FFT_RC_SCALE(CHTOT(1,ISP),CHTOT(1,ISP),GRIDC)
         CALL SETUNB_COMPAT(CHTOT(1,ISP),GRIDC)
      ENDDO
!-----------------------------------------------------------------------
! FFT of the exchange-correlation potential to reciprocal space
!-----------------------------------------------------------------------
      RINPL=1._q/GRIDC%NPLWV
      DO  ISP=1,WDES%NCDIJ 
         CALL RL_ADD(CVTOT(1,ISP),RINPL,CVTOT(1,ISP),0.0_q,CVTOT(1,ISP),GRIDC)
         CALL FFT3RC(CVTOT(1,ISP),GRIDC,-1)
      ENDDO
!-----------------------------------------------------------------------
! add the hartree potential and the double counting corrections
!-----------------------------------------------------------------------
      CALL POTHAR(GRIDC, LATT_CUR, CHTOT, CWORK,E%DENC)
      DO I=1,GRIDC%RC%NP
         CVTOT(I,1)=CVTOT(I,1)+CWORK(I,1)
      ENDDO
!-----------------------------------------------------------------------
!  add local pseudopotential potential
!-----------------------------------------------------------------------
      CALL POTION(GRIDC,P,LATT_CUR,T_INFO,CWORK,CWORK1,CSTRF,E%PSCENC)

      CALL CL_TEST_POT(CVTOT(1,1), CHTOT(1,1), CWORK, GRIDC, E_CL_ION)
      E%PSCENC=E%PSCENC+E_CL_ION

      ELECTROSTATIC=0
      NG=1
      col: DO NC=1,GRIDC%RC%NCOL
      N2= GRIDC%RC%I2(NC)
      N3= GRIDC%RC%I3(NC)
      row: DO N1=1,GRIDC%RC%NROW
        SETFACT1
        SETFACT

        ELECTROSTATIC=ELECTROSTATIC+ MULFACT CWORK(NG,1)*CONJG(CHTOT(NG,1))
        NG=NG+1
      ENDDO row
      ENDDO col
      ELECTROSTATIC=ELECTROSTATIC+E%PSCENC-E%DENC+E%TEWEN

      E%PSCENC=E%PSCENC + DIP%EDIPOL + DIP%E_ION_EXTERN

      DO I=1,GRIDC%RC%NP
         CVTOT(I,1)=CVTOT(I,1)+CWORK(I,1)
      ENDDO
      CALL POT_FLIP(CVTOT, GRIDC,WDES%NCDIJ )
!=======================================================================
! if overlap is used :
! copy CVTOT to SV and set contribution of unbalanced lattice-vectors
! to zero,  then  FFT of SV and CVTOT to real space
!=======================================================================

      DO ISP=1,WDES%NCDIJ
         CALL SETUNB_COMPAT(CVTOT(1,ISP),GRIDC)
         CALL CP_GRID(GRIDC,GRID_SOFT,SOFT_TO_C,CVTOT(1,ISP),CWORK1)
         CALL SETUNB(CWORK1,GRID_SOFT)
         CALL FFT3RC(CWORK1,GRID_SOFT, 1)
         CALL RL_ADD(CWORK1,1.0_q,CWORK1,0.0_q,SV(1,ISP),GRID_SOFT)

    !  final result is only correct for first in-band-group
    ! (i.e. proc with nodeid 1 in COMM_INTER)
    !  copy to other in-band-groups using COMM_INTER
    ! (see SET_RL_GRID() in mgrid.F, and M_divide() in mpi.F)
    #ifdef realmode
         CALLMPI( M_bcast_d(COMM_INTER, SV(1,ISP), GRID%RL%NP))
    #else
         CALLMPI( M_bcast_z(COMM_INTER, SV(1,ISP), GRID%RL%NP))
    #endif
         CALL FFT3RC(CVTOT(1,ISP),GRIDC,1)
      ENDDO

      DEALLOCATE(CWORK1,CWORK)

      RETURN
      END SUBROUTINE POTLOK


!************************ SUBROUTINE POTXC  ****************************
!
! this subroutine to calculate the XC-potential including gradient
! corrections
! for ITYPE
!     0   no gradient correction
!     1   gradient corrections on
! this routine is required to calculate the parital core
! corrections to the forces
!***********************************************************************

      SUBROUTINE POTXC(EXCTAB, GRIDC, INFO, WDES, LATT_CUR, CVTOT,CHTOT,DENCOR)
      USE prec

      USE xcgrad
      USE setexm
      USE mpimy
      USE mgrid
      USE lattice
      USE base
      USE wave

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRIDC
      TYPE (wavedes)     WDES
      TYPE (exctable)    EXCTAB
      TYPE (info_struct) INFO
      TYPE (latt)        LATT_CUR

      COMPLEX(q) CHTOT(GRIDC%MPLWV,WDES%NCDIJ),CVTOT(GRIDC%MPLWV,WDES%NCDIJ)
      RGRID      DENCOR(GRIDC%RL%NP)
! work arrays

      REAL(q)    XCSIF(3,3)
      COMPLEX(q), ALLOCATABLE:: CWORK(:,:)

      CVTOT = 0 
      ALLOCATE(CWORK(GRIDC%MPLWV,WDES%NCDIJ))

      DO ISP=1,WDES%NCDIJ
         CALL FFT3RC(CHTOT(1,ISP),GRIDC,1)
      ENDDO
        IF (WDES%ISPIN==2) THEN

          ! get the charge and the total magnetization
          CALL MAG_DENSITY(CHTOT, CWORK, GRIDC, WDES%NCDIJ)

          IF (INFO%LEXCHG >0) THEN
             ! gradient corrections to LDA
             ! unfortunately FEXCGS requires (up,down) density
             ! instead of (rho,mag)
             CALL RL_FLIP(CWORK, GRIDC, 2, .TRUE.)
             ! GGA potential
             CALL FEXCGS(2, GRIDC, LATT_CUR, XCENCG, EXCG, CVZERG, XCSIF, &
                  CWORK, CVTOT, DENCOR, INFO%LEXCHG)
             CALL RL_FLIP(CWORK, GRIDC, 2, .FALSE.)
          ENDIF

          ! add LDA part of potential
          CALL FEXCF(EXCTAB,GRIDC,LATT_CUR%OMEGA, &
             CWORK(1,1), CWORK(1,2), DENCOR, CVTOT(1,1), CVTOT(1,2), &
             CVZERO,EXC,XCENC,XCSIF, .TRUE.)
          ! we have now the potential for up and down stored in CVTOT(:,1) and CVTOT(:,2)

          ! get the proper direction vx = v0 + hat m delta v
          CALL MAG_DIRECTION(CHTOT(1,1), CVTOT(1,1), GRIDC, WDES%NCDIJ)
        ELSEIF (WDES%LNONCOLLINEAR) THEN
!-MM- gradient corrections in the noncollinear case are calculated
!     a bit differently than in the collinear case
          IF (INFO%LEXCHG >0) THEN
             ! GGA potential
             CALL FEXCGS(4, GRIDC, LATT_CUR, XCENCG, EXCG, CVZERG, XCSIF, &
                  CHTOT, CVTOT, DENCOR, INFO%LEXCHG)
          ENDIF

          ! FEXCF requires (up,down) density instead of (rho,mag)
          CALL MAG_DENSITY(CHTOT, CWORK, GRIDC, WDES%NCDIJ)
          ! add LDA part of potential
          CALL FEXCF(EXCTAB,GRIDC,LATT_CUR%OMEGA, &
             CWORK(1,1), CWORK(1,2), DENCOR, CVTOT(1,1), CVTOT(1,2), &
             CVZERO,EXC,XCENC,XCSIF, .TRUE.)
          ! we have now the potential for up and down stored in CVTOT(:,1) and CVTOT(:,2)
          ! get the proper direction vx = v0 + hat m delta v
          CALL MAG_DIRECTION(CHTOT(1,1), CVTOT(1,1), GRIDC, WDES%NCDIJ)
!-MM- end of changes to calculation of gga in noncollinear case

       ELSE
          IF (INFO%LEXCHG >0) THEN
             ! gradient corrections to LDA
             CALL FEXCG(INFO%LEXCHG,GRIDC,LATT_CUR,XCENCG,EXCG,CVZERG,XCSIF, &
                  CHTOT,CVTOT,DENCOR)
          ENDIF

          ! LDA part of potential
          CALL FEXCP(EXCTAB,GRIDC,LATT_CUR%OMEGA, &
               CHTOT,DENCOR,CVTOT,CWORK,CVZERO,EXC,XCENC,XCSIF,.TRUE.)
       ENDIF

       DO ISP=1,WDES%NCDIJ
          CALL FFT_RC_SCALE(CHTOT(1,ISP),CHTOT(1,ISP),GRIDC)
          CALL SETUNB_COMPAT(CHTOT(1,ISP),GRIDC)
       ENDDO
!-----------------------------------------------------------------------
! FFT of the exchange-correlation potential to reciprocal space
!-----------------------------------------------------------------------
      RINPL=1._q/GRIDC%NPLWV
      DO  ISP=1,WDES%NCDIJ 
         CALL RL_ADD(CVTOT(1,ISP),RINPL,CVTOT(1,ISP),0.0_q,CVTOT(1,ISP),GRIDC)
         CALL FFT3RC(CVTOT(1,ISP),GRIDC,-1)
         CALL SETUNB_COMPAT(CVTOT(1,ISP),GRIDC)
      ENDDO
      DEALLOCATE(CWORK)

      END SUBROUTINE


      END MODULE

!************************ SUBROUTINE POTHAR ****************************
!
! this subroutine calculates the hartree potential from the electronic
! charge density. The correction to the
! total energy due to overcounting the hartree energy on summing the
! electronic eigenvalues is also computed (hartree contribution to the
! total energy = 0.5*sum (vh(g)*rho(-g)). the sum of eigenvalues gives
! sum (vh(g)*rho(-g)) where vh(g) is the hartree potential at wavevector
! g and rho(g) is the charge density at wavevector g)
!
!***********************************************************************

      SUBROUTINE POTHAR(GRIDC,LATT_CUR, CHTOT,CVD,DENC)
      USE prec
      USE mpimy
      USE mgrid
      USE lattice
      USE constant
      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRIDC
      TYPE (latt)        LATT_CUR
      COMPLEX(q) CVD(GRIDC%RC%NP),CHTOT(GRIDC%RC%NP)

      DENC=0._q
!=======================================================================
! scale the hartree potential by edeps divided by the volume of the unit
! cell
!=======================================================================
      SCALE=EDEPS/LATT_CUR%OMEGA/TPI**2
!=======================================================================
! calculate the hartree potential on the grid of reciprocal lattice
! vectors and the correction to the total energy
!=======================================================================
      NI=0
      col: DO NC=1,GRIDC%RC%NCOL
      N2= GRIDC%RC%I2(NC)
      N3= GRIDC%RC%I3(NC)
      row: DO N1=1,GRIDC%RC%NROW

        NI=NI+1
        SETFACT1
        SETFACT

        GX= (GRIDC%LPCTX(N1)*LATT_CUR%B(1,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(1,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(1,3))
        GY= (GRIDC%LPCTX(N1)*LATT_CUR%B(2,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(2,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(2,3))
        GZ= (GRIDC%LPCTX(N1)*LATT_CUR%B(3,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(3,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(3,3))

        GSQU=GX**2+GY**2+GZ**2
!=======================================================================
! since the G=0 coulomb contributions to the hartree, ewald and
! electron-ion energies are individually divergent but together sum to
! zero, set the hartree potential at G=0 to zero.
!=======================================================================
        IF ((GRIDC%LPCTX(N1)==0).AND.(GRIDC%LPCTY(N2)==0).AND.(GRIDC%LPCTZ(N3)==0)) &
     & THEN
          CVD(NI)=(0.0_q,0.0_q)
        ELSE
          CVD(NI)=CHTOT(NI)/GSQU*SCALE
        ENDIF
      ENDDO row
      ENDDO col
      CALL SETUNB(CVD,GRIDC)
!=======================================================================
! calculate the correction to the total energy
!=======================================================================
      NI=0
      col2: DO NC=1,GRIDC%RC%NCOL
      N2= GRIDC%RC%I2(NC)
      N3= GRIDC%RC%I3(NC)
      row2: DO N1=1,GRIDC%RC%NROW

        NI=NI+1
        SETFACT1
        SETFACT

        DUM=MULFACT CVD(NI)*CONJG(CHTOT(NI))
        DENC=DENC+DUM
      ENDDO row2
      ENDDO col2
      DENC=-DENC/2

      CALLMPI( M_sum_d(GRIDC%COMM,DENC,1))

      RETURN
      END SUBROUTINE

!************************ SUBROUTINE POTION ****************************
!
! this subroutine calculates the pseudopotential and its derivatives
! multiplied by the partial structur-factors
! on the grid of  reciprocal lattice vectors
!
!***********************************************************************

      SUBROUTINE POTION(GRIDC,P,LATT_CUR,T_INFO,CVPS,CDVPS,CSTRF,PSCENC)
      USE prec

      USE mpimy
      USE mgrid
      USE pseudo
      USE lattice
      USE poscar
      USE constant

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRIDC
      TYPE (type_info)   T_INFO
      TYPE (potcar)      P (T_INFO%NTYP)
      TYPE (latt)        LATT_CUR

      COMPLEX(q) CSTRF(GRIDC%MPLWV,T_INFO%NTYP)
      COMPLEX(q) CVPS(GRIDC%RC%NP),CDVPS(GRIDC%RC%NP)

!=======================================================================
! calculate the contribution to the total energy from the non-coulomb
! part of the g=0 component of the pseudopotential and the force on the
! unit cell due to the change in this energy as the size of the cell
! changes
!=======================================================================
      ZVSUM=0
      DO NT=1,T_INFO%NTYP
         ZVSUM=ZVSUM+P(NT)%ZVALF*T_INFO%NITYP(NT)
      ENDDO

      PSCENC=0
      DO NT=1,T_INFO%NTYP
         PSCENC=PSCENC+P(NT)%PSCORE*T_INFO%NITYP(NT)*(ZVSUM/LATT_CUR%OMEGA)
      ENDDO

      CVPS =0
      CDVPS=0
!=======================================================================
! loop over all types of atoms
! multiply structur factor by local pseudopotential
!=======================================================================
      typ: DO NT=1,T_INFO%NTYP

      ARGSC=NPSPTS/P(NT)%PSGMAX
      PSGMA2=P(NT)%PSGMAX-P(NT)%PSGMAX/NPSPTS
      ZZ=  -4*PI*P(NT)%ZVALF*FELECT

      N=0
      col: DO NC=1,GRIDC%RC%NCOL
      N2= GRIDC%RC%I2(NC)
      N3= GRIDC%RC%I3(NC)
      row: DO N1=1,GRIDC%RC%NROW
        N=N+1
!=======================================================================
! calculate the magnitude of the reciprocal lattice vector
!=======================================================================
        GX= GRIDC%LPCTX(N1)*LATT_CUR%B(1,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(1,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(1,3)
        GY= GRIDC%LPCTX(N1)*LATT_CUR%B(2,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(2,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(2,3)
        GZ= GRIDC%LPCTX(N1)*LATT_CUR%B(3,1)+GRIDC%LPCTY(N2)*LATT_CUR%B(3,2)+GRIDC%LPCTZ(N3)*LATT_CUR%B(3,3)

        G=SQRT(GX**2+GY**2+GZ**2)*2*PI
        IF ( ( (GRIDC%LPCTX(N1)/=0) .OR. (GRIDC%LPCTY(N2)/=0) .OR. &
     &         (GRIDC%LPCTZ(N3)/=0) ) .AND. (G<PSGMA2) ) THEN
!=======================================================================
! convert the magnitude of the reciprocal lattice vector to a position
! in the pseudopotential arrays and interpolate the pseudopotential and
! its derivative
!=======================================================================
        I  =INT(G*ARGSC)+1
        REM=G-P(NT)%PSP(I,1)
        VPST =(P(NT)%PSP(I,2)+REM*(P(NT)%PSP(I,3)+ &
     &                     REM*(P(NT)%PSP(I,4)  +REM*P(NT)%PSP(I,5))))
        DVPST= P(NT)%PSP(I,3)+REM*(P(NT)%PSP(I,4)*2+REM*P(NT)%PSP(I,5)*3)

        CVPS (N)=CVPS (N)+( VPST+ ZZ / G**2)  /LATT_CUR%OMEGA*CSTRF(N,NT)
        CDVPS(N)=CDVPS(N)+( DVPST- 2*ZZ /G**3)/LATT_CUR%OMEGA*CSTRF(N,NT)
        ELSE
        CVPS (N)=0._q
        CDVPS(N)=0._q
        ENDIF

      ENDDO row
      ENDDO col
      ENDDO typ

      CALL SETUNB(CVPS,GRIDC)
      CALL SETUNB(CDVPS,GRIDC)

      RETURN
      END SUBROUTINE


!************************ SUBROUTINE  MAG_DIRECTION  *******************
!
! on entry CVTOT must contain the v_xc(up) and v_xc(down)
! on return CVTOT contains
! collinear case:
!  CVTOT(:,1) =  (v_xc(up) + v_xc(down))/2
!  CVTOT(:,2) =  (v_xc(up) - v_xc(down))/2
! non collinear case:
!  CVTOT(:,1) =  (v_xc(up) + v_xc(down))/2
!  CVTOT(:,2) =  hat m_x (v_xc(up) - v_xc(down))/2
!  CVTOT(:,3) =  hat m_y (v_xc(up) - v_xc(down))/2
!  CVTOT(:,4) =  hat m_z (v_xc(up) - v_xc(down))/2
! where hat m is the unit vector of the local magnetization density
!
!***********************************************************************

#undef strict_compatibility

      SUBROUTINE MAG_DIRECTION(CHTOT, CVTOT, GRID, NCDIJ)

      USE prec
      USE mgrid

      IMPLICIT NONE
      TYPE (grid_3d)     GRID
      INTEGER NCDIJ
      
      RGRID CHTOT(DIMREAL(GRID%MPLWV), NCDIJ), &
            CVTOT(DIMREAL(GRID%MPLWV), NCDIJ)
      ! local
      INTEGER K
      REAL(q) :: NORM2,DELTAV,V0
      
      IF (NCDIJ==2) THEN
         DO K=1,GRID%RL%NP
            V0    =(CVTOT(K,1)+CVTOT(K,2))/2
            DELTAV=(CVTOT(K,1)-CVTOT(K,2))/2
            CVTOT(K,1) = V0
#ifndef strict_compatibility
            CVTOT(K,2) = DELTAV
#else
            CVTOT(K,2) = DELTAV*REAL(CHTOT(K,2),KIND=q) / MAX(ABS(CHTOT(K,2)),1.E-20_q)
#endif
         ENDDO
      ELSE IF (NCDIJ==4) THEN
         DO K=1,GRID%RL%NP
            V0    =(CVTOT(K,1)+CVTOT(K,2))/2
            DELTAV=(CVTOT(K,1)-CVTOT(K,2))/2
#ifdef realmode
            NORM2 = MAX(SQRT(ABS(CHTOT(K,2)*CHTOT(K,2)+ CHTOT(K,3)*CHTOT(K,3)+ CHTOT(K,4)*CHTOT(K,4))),1.E-20_q)
#else
            NORM2 = MAX(SQRT(ABS(CHTOT(K,2)*CONJG(CHTOT(K,2))+ &
           &           CHTOT(K,3)*CONJG(CHTOT(K,3)) + CHTOT(K,4)*CONJG(CHTOT(K,4)))),1.E-20_q)
#endif

            CVTOT(K,1) = V0
            CVTOT(K,2) = DELTAV * REAL(CHTOT(K,2),KIND=q) / NORM2
            CVTOT(K,3) = DELTAV * REAL(CHTOT(K,3),KIND=q) / NORM2
            CVTOT(K,4) = DELTAV * REAL(CHTOT(K,4),KIND=q) / NORM2
         ENDDO
         
      ELSE
         WRITE(*,*) 'internal error: MAG_DIRECTION called with NCDIJ=',NCDIJ
         STOP
      ENDIF
         
      END SUBROUTINE MAG_DIRECTION


!************************ SUBROUTINE MAG_DENSITY ***********************
!
! this subroutine calculates the total charge density and the 
! absolute magnitude of the magnetization density
! on entry: 
!  CHTOT  rho, m_x, m_y, m_z
! on exit:
!  CWORK  rho, sqrt(m_x^2 + m_y^2 + m_z^2)
! in the collinear case, it this means a simple copy CHTOT to CWORK
!
!***********************************************************************

      SUBROUTINE MAG_DENSITY(CHTOT, CWORK, GRID, NCDIJ)

      USE prec
      USE mgrid

      IMPLICIT NONE
      TYPE (grid_3d)     GRID
      INTEGER NCDIJ
      
      RGRID CHTOT(DIMREAL(GRID%MPLWV), NCDIJ), &
            CWORK(DIMREAL(GRID%MPLWV), NCDIJ)
      ! local
      INTEGER K

      
      IF (NCDIJ==2) THEN
         DO K=1,GRID%RL%NP
            CWORK(K,1)=CHTOT(K,1)
#ifndef strict_compatibility
            CWORK(K,2)=CHTOT(K,2)
#else
            CWORK(K,2)=ABS(CHTOT(K,2))
#endif
         ENDDO
      ELSE IF (NCDIJ==4) THEN
         DO K=1,GRID%RL%NP
            CWORK(K,1)=CHTOT(K,1)
            CWORK(K,2)=SQRT(ABS(CHTOT(K,2)*CHTOT(K,2)+ CHTOT(K,3)*CHTOT(K,3) + CHTOT(K,4)*CHTOT(K,4)))
         ENDDO
      ELSE
         WRITE(*,*) 'internal error: MAG_DENSITY called with NCDIJ=',NCDIJ
         STOP
      ENDIF
      END SUBROUTINE MAG_DENSITY


!************************ SUBROUTINE POT_FLIP **************************
!
!
! rearranges the storage mode for spin components of potentials:
! for the collinear case calculate:
!  v0 1 + v_z
!  v0 1 - v_z
! for the non collinear case calculate
!  v  = v0 1 + sigma_x v_x + simga_y v_y + sigma_z v_z
! 
!***********************************************************************

      SUBROUTINE POT_FLIP(CVTOT, GRID, NCDIJ)
      USE prec
      USE mgrid
      IMPLICIT NONE
      INTEGER NCDIJ
      TYPE (grid_3d)     GRID
      COMPLEX(q) :: CVTOT(GRID%MPLWV, NCDIJ)

      ! local
      COMPLEX(q) :: C00,CX,CY,CZ
      REAL(q) :: FAC
      INTEGER K
      
      IF (NCDIJ==2) THEN
         FAC=1.0_q
         DO K=1,GRID%RC%NP
            C00=CVTOT(K,1)
            CZ =CVTOT(K,2)

            CVTOT(K,1)= (C00+CZ)*FAC           
            CVTOT(K,2)= (C00-CZ)*FAC           
         ENDDO
      ELSE IF (NCDIJ==4) THEN
         FAC=1.0_q
         DO K=1,GRID%RC%NP
            C00=CVTOT(K,1)
            CX =CVTOT(K,2)
            CY =CVTOT(K,3)
            CZ =CVTOT(K,4)

            CVTOT(K,1)= (C00+CZ)*FAC           
            CVTOT(K,2)= (CX-CY*(0._q,1._q))*FAC
            CVTOT(K,3)= (CX+CY*(0._q,1._q))*FAC
            CVTOT(K,4)= (C00-CZ)*FAC           
         ENDDO
      ELSE IF (NCDIJ==1) THEN
      ENDIF

    END SUBROUTINE POT_FLIP


!************************ SUBROUTINE EXTERNAL_POT **********************
!
! this subroutine can be used to add an external potential
! the units of the potential are eV
!
!***********************************************************************

      SUBROUTINE EXTERNAL_POT(GRIDC, LATT_CUR, CVTOT)

      USE prec
      USE base
      USE lattice
      USE mpimy
      USE mgrid
      USE poscar
      USE constant

      IMPLICIT COMPLEX(q) (C)
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      TYPE (grid_3d)     GRIDC
      TYPE (latt)        LATT_CUR
      RGRID      CVTOT(DIMREAL(GRIDC%MPLWV))

      RETURN

      NG=0

      IF (GRIDC%RL%NFAST==3) THEN
         ! mpi version: x-> N2, y-> N3, z-> N1
         N2MAX=GRIDC%NGX
         N3MAX=GRIDC%NGY
         N1MAX=GRIDC%NGZ

         DO NC=1,GRIDC%RL%NCOL
            N2= GRIDC%RL%I2(NC)
            N3= GRIDC%RL%I3(NC)
            DO N1=1,GRIDC%RL%NROW
               NG=NG+1
               CVTOT(NG)=CVTOT(NG)
            ENDDO
         ENDDO
      ELSE
         ! conventional version: x-> N1, y-> N2, z-> N3
         N1MAX=GRIDC%NGX
         N2MAX=GRIDC%NGY
         N3MAX=GRIDC%NGZ

         DO NC=1,GRIDC%RL%NCOL
            N2= GRIDC%RL%I2(NC)
            N3= GRIDC%RL%I3(NC)
            DO N1=1,GRIDC%RL%NROW
               NG=NG+1
               CVTOT(NG)=CVTOT(NG)
            ENDDO
         ENDDO
      ENDIF

      END SUBROUTINE
